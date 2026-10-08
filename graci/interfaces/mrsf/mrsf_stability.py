"""
Internal stability check of the ROKS/ROHF triplet reference of
MRSF-TDDFT (keyword stability = True in the $mrsftddft section).

The lowest eigenvalues of the orbital Hessian of the reference
(mrsf_hessian.RoksHessian: exact second derivative of the reference
energy with respect to the C->O, C->V and O->V rotations) are found by
a block Davidson started on the lowest diagonal elements. A Ritz value
is an upper bound of the lowest eigenvalue, so the check stops as soon
as one falls below -THRESH. An unstable reference is followed: the
orbitals are rotated along the eigenvector (line search on the
reference energy), the SCF is re-run from them, the integrals are
transformed again and the check is repeated. With point-group symmetry
only totally symmetric rotations are considered, so that the reference
keeps its symmetry.
"""

import sys
import time
import types
import numpy as np
import scipy.linalg

import graci.io.output as output
import graci.io.chkpt as chkpt
import graci.core.ao2mo as ao2mo
import graci.interfaces.mrsf.mrsf_init as mrsf_init
import graci.interfaces.mrsf.mrsf_hessian as mrsf_hessian

# lowest-eigenvalue threshold (Hartree, PySCF's convention for the
# Hessian of the same parametrisation)
THRESH = 1.e-5
# residual-norm threshold of the Davidson solver
TOL = 1.e-3
# maximum number of check / re-optimisation rounds
MAX_ROUNDS = 3
# trial step lengths along the unstable direction
STEPS = (0.25, 0.5, 1.0, 1.5)


def _print(msg):
    with output.output_file(output.file_names['out_file'], 'a+') as f:
        f.write(msg + '\n')


def davidson_lowest(apply, diag, mask, nroots=2, nguess=4, tol=TOL,
                    maxit=40, maxvec=48, thresh=THRESH):
    """
    lowest eigenpairs of the symmetric operator apply (vectors of length
    len(diag)) restricted to the components where mask is True.
    Returns (eigenvalues, eigenvectors, residual norms, iterations,
    status), status in 'unstable', 'stable', 'not converged'.
    """
    d = np.where(mask, diag, np.inf)
    nguess = min(nguess, int(mask.sum()))
    nroots = min(nroots, nguess)
    idx = np.argsort(d)[:nguess]
    V = np.zeros((len(diag), nguess))
    V[idx, np.arange(nguess)] = 1.0

    def op(X):
        return apply(X*mask[:, None])*mask[:, None]

    AV = op(V)
    for it in range(1, maxit + 1):
        G = V.T @ AV
        w, u = np.linalg.eigh(0.5*(G + G.T))
        X = V @ u[:, :nroots]
        AX = AV @ u[:, :nroots]
        R = AX - X*w[:nroots]
        rn = np.linalg.norm(R, axis=0)
        # a Ritz value bounds the lowest eigenvalue from above
        if w[0] < -thresh and rn[0] < max(tol, abs(w[0])/3.0):
            return w[:nroots], X, rn, it, 'unstable'
        if np.all(rn < tol):
            status = 'unstable' if w[0] < -thresh else 'stable'
            return w[:nroots], X, rn, it, status
        # preconditioned residuals of the unconverged roots
        new = []
        for j in range(nroots):
            if rn[j] < tol:
                continue
            den = diag - w[j]
            den[np.abs(den) < 1.e-4] = 1.e-4
            new.append(-R[:, j]/den*mask)
        N = np.array(new).T
        if V.shape[1] + N.shape[1] > maxvec:
            V, AV = X, AX
        for _ in range(2):
            N -= V @ (V.T @ N)
        keep = []
        for j in range(N.shape[1]):
            for i in keep:
                N[:, j] -= np.dot(N[:, i], N[:, j])*N[:, i]
            nrm = np.linalg.norm(N[:, j])
            if nrm > 1.e-8:
                N[:, j] /= nrm
                keep.append(j)
        if not keep:
            break
        N = N[:, keep]
        V = np.hstack([V, N])
        AV = np.hstack([AV, op(N)])
    return w[:nroots], X, rn, maxit, 'not converged'


def _check(ci, mf, fock_mo):
    """one stability check of the reference loaded in the library"""
    t0 = time.time()
    hess = mrsf_hessian.RoksHessian(mf, ci.mos, ci.occ_ref, fock_mo[0], fock_mo[1],
                                    float(ci.scf.hyb[2]), ci.eri_file(),
                                    (ci.nocca, ci.nvirb, ci.naux), ci.mem_budget)
    t_setup = time.time() - t0
    mosym = np.asarray(ci.mosym)
    mask = mosym[hess.Rp] == mosym[hess.Rq]
    t0 = time.time()
    w, X, rn, it, status = davidson_lowest(hess.apply, hess.diag(), mask)
    res = {'status': status, 'eig': w, 'vec': X[:, 0], 'Rp': hess.Rp, 'Rq': hess.Rq,
           'nprod': hess.nprod, 'niter': it, 'lz': int(mask.sum()),
           't_setup': t_setup, 't_check': time.time() - t0,
           'cached': None if hess.xcgrid is None else bool(hess.xcgrid.cache)}
    hess.free()
    return res


def _rotate(C, Rp, Rq, z, t):
    """C exp(t K), K_pq = z_pq = -K_qp for the rotation pairs (Rp, Rq)"""
    K = np.zeros((C.shape[1], C.shape[1]))
    K[Rp, Rq] = t*z
    K[Rq, Rp] = -t*z
    return C @ scipy.linalg.expm(K)


def _follow(ci, mf, res):
    """line search on the reference energy along the unstable direction;
    returns the rotated orbitals (all MOs of the Scf object) and the step"""
    scf = ci.scf
    C = np.array(scf.orbs)
    occ = np.asarray(scf.orb_occ)
    nmo = ci.nmo
    z = res['vec']/np.linalg.norm(res['vec'])
    best = (scf.energy, None, None)
    for t in STEPS:
        Ct = C.copy()
        Ct[:, :nmo] = _rotate(C[:, :nmo], res['Rp'], res['Rq'], z, t)
        e = mf.energy_tot(mf.make_rdm1(Ct, occ))
        if e < best[0]:
            best = (e, t, Ct)
    if best[1] is None:
        t = STEPS[0]
        Ct = C.copy()
        Ct[:, :nmo] = _rotate(C[:, :nmo], res['Rp'], res['Rq'], z, t)
        best = (None, t, Ct)
    return best[2], best[1], best[0]


def ensure_stable(ci, mo_ints):
    """
    check the internal stability of the reference of ci (library session
    without the extended space) and, if it is unstable, re-optimise the
    shared Scf object along the unstable direction, transform the
    integrals again and repeat (at most MAX_ROUNDS checks). Returns True
    if the final reference is stable.
    """
    scf = ci.scf
    if mo_ints is None:
        mo_ints = ao2mo.Ao2mo()
        mo_ints.emo_cut = ci.mo_cutoff
        mo_ints.run(scf, ci.precision)
    if ci.precision != 'double':
        output.print_message('stability check skipped: it requires precision = double')
        return False

    _print('\n Internal stability of the ROKS triplet reference')
    _print(' ------------------------------------------------')
    stable = False
    for rnd in range(1, MAX_ROUNDS + 1):
        fock_mo = np.array([ci.mos.T @ scf.fock_ao[0] @ ci.mos,
                            ci.mos.T @ scf.fock_ao[1] @ ci.mos])
        mrsf_init.set_extended(False)
        mrsf_init.init(ci, fock_mo)
        mrsf_init.init_ints(ci, ci.eri_file())
        ci.nocca, ci.nvirb, ci.xdim, ci.naux = [int(x) for x in mrsf_init.dims()]
        mf = scf.pyscf_obj()
        res = _check(ci, mf, fock_mo)
        _print('  check %d: lowest Hessian eigenvalues %s Eh -> %s' %
               (rnd, ' '.join('%.5f' % e for e in res['eig']), res['status']))
        _print('           %d rotations, %d Hessian products in %d iterations, '
               'set-up %.1f s, check %.1f s%s' %
               (res['lz'], res['nprod'], res['niter'], res['t_setup'], res['t_check'],
                '' if res['cached'] is None else
                (', AO values cached' if res['cached'] else ', AO values streamed '
                 '(increase mem_budget to cache them)')))
        if res['status'] != 'unstable':
            stable = res['status'] == 'stable'
            break
        if rnd == MAX_ROUNDS:
            break

        # follow the instability: line search, then a new SCF from there
        e0 = scf.energy
        Ct, t, et = _follow(ci, mf, res)
        _print('  following the unstable direction: step %.2f, E = %s Eh' %
               (t, 'n/a' if et is None else '%.8f' % et))
        guess = types.SimpleNamespace(mol=scf.mol, orbs=Ct, orb_occ=scf.orb_occ)
        if scf.run(scf.mol, guess) is None:
            sys.exit('stability check: the re-optimised SCF did not converge')
        _print('  re-optimised reference: E = %.8f Eh (%.4f eV below the previous one)' %
               (scf.energy, (e0 - scf.energy)*27.211386245988))
        chkpt.write(scf)
        mo_ints.run(scf, ci.precision)
        # the integral file keeps its name: drop the integrals held by the
        # library so that the next initialisation reads the new ones
        mrsf_init.finalize()
        ci.set_scf(scf, mo_ints=mo_ints)
        ci.check_input()
    if not stable:
        _print('  WARNING: the reference is not confirmed to be stable')
    if rnd > 1:
        _print('  NOTE: CI sections that ran before with the same $scf section '
               'used the previous reference')
    return stable
