"""
Trust-region augmented Hessian (TRAH) SCF solver

Follows B. Helmich-Paris, J. Chem. Phys. 154, 164104 (2021). The orbital
gradient, Hessian-vector products and orbital rotations are taken from
PySCF's second-order (CIAH) machinery, so any reference PySCF's
newton() supports (RHF/ROHF/RKS/ROKS, density fitting, X2C) works here.
"""

import numpy as np
from pyscf import scf, dft
from pyscf.soscf import newton_ah
from pyscf.lib import logger


def kernel(mf, dm0=None, conv_tol=1e-8, conv_tol_grad=None,
           max_cycle=50, trust0=0.5, trust_max=1.0, max_micro=30,
           hess_grid_level=0, verbose=None):
    """
    Converge the SCF of the PySCF mean-field object mf with TRAH.
    On return, mf carries the converged orbitals, energies and e_tot,
    as after mf.kernel()

    hess_grid_level: for KS references, a (coarser) DFT grid level for
    the XC kernel in the Hessian-vector products only. The energy and
    gradient keep the full grid, so the converged result is unchanged
    """
    log = logger.new_logger(mf, verbose)

    if conv_tol_grad is None:
        conv_tol_grad = np.sqrt(conv_tol)

    # second-order wrapper: supplies gen_g_hop, update_rotate_matrix
    # and rotate_mo for the reference type of mf
    so  = mf.newton()
    if hess_grid_level is not None and hasattr(mf, 'grids'):
        grids       = dft.gen_grid.Grids(mf.mol)
        grids.level = hess_grid_level
        grids.prune = mf.grids.prune
        so.grids    = grids.build()
    mol = mf.mol
    h1e = mf.get_hcore(mol)
    s1e = mf.get_ovlp(mol)

    # starting orbitals from a single diagonalisation of the guess Fock
    if dm0 is None:
        dm0 = mf.get_init_guess(mol, mf.init_guess)
    vhf      = mf.get_veff(mol, dm0)
    fock     = mf.get_fock(h1e, s1e, vhf, dm0)
    mo_e, mo = mf.eig(fock, s1e)
    mo_occ   = mf.get_occ(mo_e, mo)
    so.mo_coeff, so.mo_occ = mo, mo_occ

    def energy(mo):
        dm  = mf.make_rdm1(mo, mo_occ)
        vhf = mf.get_veff(mol, dm)
        return mf.energy_tot(dm, h1e, vhf), dm, vhf

    e_tot, dm, vhf = energy(mo)
    log.info('TRAH initial E= %.15g', e_tot)

    trust     = trust0
    converged = False
    de        = 0.
    new_point = True
    x_last    = None

    for imacro in range(max_cycle):

        # gradient and Hessian at the current orbitals (only
        # recomputed after an accepted step)
        if new_point:
            fock = mf.get_fock(h1e, s1e, vhf, dm)
            g, h_op, h_diag = so.gen_g_hop(mo, mo_occ, fock)
            if isinstance(mf, scf.rohf.ROHF) and _ROHF_HESS_INCOMPLETE:
                h_op = _rohf_hop(h_op, mo, mo_occ, fock)
            gnorm = np.linalg.norm(g)
            space = None

        if gnorm < conv_tol_grad and abs(de) < conv_tol:
            converged = True
            break

        # after a rejected step the subspace is reused: g and H are
        # unchanged, only the trust radius has shrunk
        x, pred, nmicro, space = _ah_step(g, h_op, h_diag, trust, gnorm,
                                          max_micro, space, x_last)

        u      = so.update_rotate_matrix(x, mo_occ, mo_coeff=mo)
        mo_new = so.rotate_mo(mo, u)
        e_new, dm_new, vhf_new = energy(mo_new)

        de_try = e_new - e_tot
        ratio  = de_try / pred if pred < 0. else -1.
        log.info('TRAH macro= %d  E= %.15g  dE= %.3e  |g|= %.3e  '
                 '|x|= %.3e  trust= %.3e  ratio= %.3f  micro= %d',
                 imacro, e_new, de_try, gnorm, np.linalg.norm(x), trust,
                 ratio, nmicro)

        # trust-radius update
        if ratio < 0.25:
            trust *= 0.7
        elif ratio > 0.75:
            trust  = min(1.2 * trust, trust_max)

        # step acceptance: reject anything that raises the energy
        if de_try < 0. or abs(de_try) < conv_tol:
            mo, dm, vhf, de = mo_new, dm_new, vhf_new, de_try
            e_tot     = e_new
            so.mo_coeff = mo
            x_last    = x
            new_point = True
        else:
            log.info('TRAH step rejected')
            new_point = False

    # canonicalise and store the solution on mf
    fock = mf.get_fock(h1e, s1e, vhf, dm)
    mo_e, mo_can = mf.canonicalize(mo, mo_occ, fock)

    mf.mo_coeff  = mo_can
    mf.mo_energy = mo_e
    mf.mo_occ    = mo_occ
    mf.e_tot     = e_tot
    mf.converged = converged
    if converged:
        log.note('TRAH converged in %d macro iterations: E= %.15g',
                 imacro, e_tot)
    else:
        log.note('TRAH not converged after %d macro iterations: '
                 'E= %.15g  |g|= %.3e', max_cycle, e_tot, gnorm)

    return e_tot


# PySCF fixed the ROHF Hessian upstream in pyscf/pyscf#3448 (commit
# 460a2a94, after 2.12.1). The old gen_g_hop_rohf is built on
# gen_g_hop_uhf; the fixed one is not, so only correct the old one
_ROHF_HESS_INCOMPLETE = 'gen_g_hop_uhf' in \
    newton_ah.gen_g_hop_rohf.__code__.co_names


def _rohf_hop(h_op, mo, mo_occ, fock_ao):
    """
    PySCF's ROHF/ROKS Hessian is the UHF one with a common rotation for
    both spins. That drops the second-order terms in which the rotations
    that are occupied-occupied (core-open, alpha) or virtual-virtual
    (open-virtual, beta) for one spin act together with its vir-occ
    rotations:

       sum_s 1/2 Tr F^s [B_s, [A_s, D_s]]

    (A_s: vir-occ part of kappa for spin s, B_s the remainder). These
    multiply the spin-resolved Fock matrices, which do not vanish at
    ROHF convergence, so the error persists and is largest along soft
    open-shell modes. This adds them to h_op (identical to the upstream
    fix to machine precision).
    """
    mask = scf.hf.uniq_var_indices(mo_occ)
    fock = [mo.T @ f @ mo for f in (fock_ao.focka, fock_ao.fockb)]
    proj = [(mo_occ > 0).astype(float), (mo_occ == 2).astype(float)]

    def h_op_full(x):
        k = scf.hf.unpack_uniq_var(x, mo_occ)
        grad = np.zeros_like(k)
        for f, d in zip(fock, proj):
            # vir-occ (a) and remaining (b) blocks of kappa for this spin
            vo  = np.outer(1. - d, d) + np.outer(d, 1. - d)
            a   = k * vo
            b   = k - a
            ad  = a * d[None, :] - d[:, None] * a      # [A,D]
            fb  = f @ b - b @ f                        # [F,B]
            ga  = (d[:, None] * fb - fb * d[None, :]).T
            gb  = (ad @ f - f @ ad).T
            grad += 0.5 * (ga * vo + gb * (1. - vo))
        # gradient w.r.t. the antisymmetric parameters; N.B. PySCF's h_op
        # is half the true Hessian
        return h_op(x) + 0.5 * (grad - grad.T)[mask]

    return h_op_full


def _ah_step(g, h_op, h_diag, trust, gnorm, max_micro, space=None,
             x_guess=None):
    """
    Davidson solution of the augmented Hessian eigenvalue problem

       | 0        alpha g^T | |1|       |1|
       | alpha g  H         | |x| = mu  |x|

    with the step x = v / (alpha v0) held within the trust radius by
    choice of alpha. The first subspace vector is g/|g|, so alpha only
    enters the subspace matrix and can be re-tuned without new
    sigma vectors.

    space = (bvec, svec) continues from a previous call at the same
    orbitals; x_guess (e.g. the last accepted step) seeds a new space.
    Returns the step, the predicted energy change, the subspace
    dimension and the subspace.
    """
    nvar = g.size
    tol  = max(min(0.1, gnorm) * gnorm, 1e-12)

    if space is not None:
        bvec, svec = space
    else:
        bvec = [g / gnorm]
        svec = [h_op(bvec[0])]
        if x_guess is not None and x_guess.size == g.size:
            t  = x_guess - np.dot(bvec[0], x_guess) * bvec[0]
            tn = np.linalg.norm(t)
            if tn > 1e-8 * np.linalg.norm(x_guess):
                bvec.append(t / tn)
                svec.append(h_op(bvec[-1]))

    alpha = 1.
    while True:
        B  = np.array(bvec)
        S  = np.array(svec)
        hb = B @ S.T
        hb = 0.5 * (hb + hb.T)
        gb = B @ g

        alpha, mu, c0, c = _solve_alpha(hb, gb, trust)

        # residual of the full augmented problem (the first component
        # vanishes since g lies in the subspace)
        bc  = c @ B
        res = alpha * c0 * g + c @ S - mu * bc
        rnorm = np.linalg.norm(res)
        if rnorm < tol * abs(alpha * c0) or len(bvec) >= nvar \
                or len(bvec) >= max_micro:
            break

        denom = h_diag - mu
        denom[np.abs(denom) < 1e-8] = 1e-8
        t = res / denom
        for b in bvec:
            t -= np.dot(b, t) * b
        for b in bvec:
            t -= np.dot(b, t) * b
        tn = np.linalg.norm(t)
        if tn < 1e-12:
            break
        t /= tn
        bvec.append(t)
        svec.append(h_op(t))

    # predicted energy change of the quadratic model. N.B. PySCF's
    # gen_g_hop returns half the true orbital gradient and Hessian
    x    = bc / (alpha * c0)
    xc   = c / (alpha * c0)
    pred = 2. * (np.dot(gb, xc) + 0.5 * np.dot(xc, hb @ xc))

    return x, pred, len(bvec), (bvec, svec)


def _solve_alpha(hb, gb, trust):
    """
    lowest eigenvector of the subspace augmented Hessian, with alpha
    increased by bisection until the step lies within the trust radius
    """
    def solve(alpha):
        n  = hb.shape[0]
        am = np.zeros((n+1, n+1))
        am[0, 1:] = alpha * gb
        am[1:, 0] = alpha * gb
        am[1:, 1:] = hb
        w, v = np.linalg.eigh(am)
        vec  = v[:, 0]
        if abs(vec[0]) < 1e-12:
            return w[0], vec[0], vec[1:], np.inf
        xn = np.linalg.norm(vec[1:]) / abs(alpha * vec[0])
        return w[0], vec[0], vec[1:], xn

    mu, c0, c, xn = solve(1.)
    if xn <= trust:
        return 1., mu, c0, c

    lo, hi = 1., 2.
    while True:
        mu, c0, c, xn = solve(hi)
        if xn <= trust or hi > 1e8:
            break
        lo, hi = hi, 2. * hi

    for _ in range(50):
        if 0.9 * trust <= xn <= trust:
            break
        mid = 0.5 * (lo + hi)
        mu_m, c0_m, c_m, xn_m = solve(mid)
        if xn_m > trust:
            lo = mid
        else:
            hi = mid
            mu, c0, c, xn = mu_m, c0_m, c_m, xn_m

    return hi, mu, c0, c


class _Stagnated(Exception):
    """raised from the DIIS callback to hand over to TRAH"""


def auto_kernel(mf, dm0=None, conv_tol=1e-8, conv_tol_grad=None, max_cycle=50,
                diis_min_iter=10, window=8, verbose=None):
    """
    AutoTRAH: run the default (DIIS) SCF of mf and switch to TRAH,
    starting from the lowest-energy density seen, if the orbital
    gradient stagnates or DIIS runs out of iterations. Stagnation:
    after diis_min_iter cycles, |g| over the last window cycles has
    not dropped below half of its best value before that
    """
    log  = logger.new_logger(mf, verbose)
    hist = {'g': [], 'e': np.inf, 'dm': None, 'cycle': 0}

    def callback(envs):
        hist['cycle'] = envs['cycle'] + 1
        hist['g'].append(envs['norm_gorb'])
        if envs['e_tot'] < hist['e']:
            hist['e']  = envs['e_tot']
            hist['dm'] = envs['dm']
        g = hist['g']
        if len(g) >= max(diis_min_iter, window + 1) and \
                min(g[-window:]) > 0.5 * min(g[:-window]):
            raise _Stagnated

    mf.callback = callback
    mf.conv_tol = conv_tol
    if conv_tol_grad is not None:
        mf.conv_tol_grad = conv_tol_grad
    mf.max_cycle = max_cycle
    try:
        e_tot = mf.kernel(dm0=dm0)
        if mf.converged:
            log.note('AutoTRAH: DIIS converged, TRAH not needed')
            return e_tot
        log.note('AutoTRAH: DIIS not converged in %d cycles, '
                 'switching to TRAH', hist['cycle'])
    except _Stagnated:
        log.note('AutoTRAH: DIIS stagnated at cycle %d (|g|= %.3e), '
                 'switching to TRAH', hist['cycle'], hist['g'][-1])
    finally:
        mf.callback = None

    return kernel(mf, dm0=hist['dm'], conv_tol=conv_tol, conv_tol_grad=conv_tol_grad,
                  max_cycle=max_cycle, verbose=verbose)
