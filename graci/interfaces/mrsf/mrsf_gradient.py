"""
Analytic nuclear gradients of RO-MRSF-TDDFT states: the driver that
combines the libmrsf gradient kernels (MO-space density-fitted
contractions), PySCF (XC kernel, AO derivative integrals, reference
gradient) and NumPy (Fock couplings, Z-vector solver, pairing terms).

Conventions (see the plan and ~/calculations/mrsf_dev/mrsf_grad_ref.py):
  generalised Fock  F_pq = dE/dT_pq  (C -> C T)
  rotation space    R = {C->O, C->V, O->V},  (p,q) = (higher, lower) class
  ROKS gradient     g_pq = F^E_pq - F^E_qp,  F^E_pq = 2 sum_s n^s_q f^s_pq
  Z-vector          H z = -R,  R_pq = F^omega_pq - F^omega_qp,  H = dg/dkappa
  relaxed density   P^s = T^s + p^{z,s},  p^{z,s}_pq = p^{z,s}_qp = (n^s_q - n^s_p) z_pq
  energy-weighted   W = 1/2 C (F^omega + F^z) C^T   (reference part inside PySCF's ROKS gradient)
  dE_I/dR = dE_ROKS/dR + sum h^R P + 2e(DF families) + XC probe - sum W S^R
"""
import sys
import time
import numpy as np
import scipy.linalg
from pyscf import df
from pyscf.grad import rhf as rhf_grad
from pyscf.grad import tduks as tduks_grad

import graci.io.output as output
import graci.interfaces.mrsf.mrsf_init as mrsf_init
import graci.interfaces.mrsf.mrsf_grad as mrsf_grad
import graci.interfaces.mrsf.mrsf_dfgrad as mrsf_dfgrad

SQ2 = np.sqrt(2.0)


def _sym(a):
    return 0.5*(a + a.T)


class _XCShim:
    """duck-typed object for pyscf.grad.tduks._contract_xc_kernel"""
    class _Base:
        pass

    def __init__(self, mf, C, occ_a, occ_b):
        self.mol = mf.mol
        self.verbose = 0
        self.stdout = mf.stdout
        self.base = _XCShim._Base()
        self.base._scf = _XCShim._Base()
        self.base._scf.grids = mf.grids
        self.base._scf._numint = mf._numint
        self.base._scf.xc = mf.xc
        self.base._scf.mo_coeff = (C, C)
        self.base._scf.mo_occ = (occ_a, occ_b)
        self.base._scf.do_nlc = lambda: False
        self.base.exclude_nlc = True


class GradientDriver:

    def __init__(self, ci, grad_obj):
        """ci: converged Mrsftddft object; grad_obj: the Mrsfgradient section"""
        self.ci, self.scf, self.opt = ci, ci.scf, grad_obj
        self.verbose = grad_obj.verbose
        scf = self.scf
        self.mol = scf.mol.pymol()
        self.mf = scf.pyscf_obj()
        mf = self.mf
        C = np.asarray(ci.mos)
        occ = np.asarray(ci.occ_ref)
        self.C, self.nmo, self.nao = C, C.shape[1], C.shape[0]
        if self.nmo != self.mol.nao:
            sys.exit('MRSF gradients require the full MO space (nmo = nao): '
                     'set mo_cutoff large enough in the $mrsftddft section')
        if ci.precision != 'double':
            sys.exit('MRSF gradients require precision = double in the $mrsftddft section')
        self.is_dft = (scf.xc.lower() != 'hf')
        omega, alpha, hyb = scf.hyb
        if abs(omega) > 1e-12:
            sys.exit('MRSF gradients: range-separated functionals are not yet supported')
        self.x_ref = hyb
        self.x_resp = ci.chf_eff
        self.kappa = np.array(ci.spc_eff, float)
        self.s = 1.0 if ci.mult == 1 else -1.0
        self.mult = ci.mult
        # orbital classes (same rule as the library: by occupation, index order)
        self.Cidx = np.where(np.abs(occ - 2.0) < 1e-6)[0]
        self.Oidx = np.where(np.abs(occ - 1.0) < 1e-6)[0]
        self.Vidx = np.where(np.abs(occ) < 1e-6)[0]
        self.H = np.concatenate([self.Cidx, self.Oidx])
        self.P = np.concatenate([self.Oidx, self.Vidx])
        self.nC, self.nV = len(self.Cidx), len(self.Vidx)
        self.nocca, self.nvirb = len(self.H), len(self.P)
        self.dims = (self.nocca, self.nvirb, None)
        self.occ = occ
        self.occ_a = (occ > 0.5).astype(float)
        self.occ_b = (occ > 1.5).astype(float)
        self.occH = occ[self.H]
        self.Hloc = {p: k for k, p in enumerate(self.H)}
        self.Ploc = {p: k for k, p in enumerate(self.P)}
        # MO spin Fock matrices and Roothaan orbital energies
        fock_ao = scf.fock_ao if scf.fock_ao is not None else scf.build_fock()
        self.fa = C.T @ fock_ao[0] @ C
        self.fb = C.T @ fock_ao[1] @ C
        # rotation space
        R = []
        for q in self.Cidx:
            for p in self.Oidx:
                R.append((p, q))
        for q in self.Cidx:
            for p in self.Vidx:
                R.append((p, q))
        for q in self.Oidx:
            for p in self.Vidx:
                R.append((p, q))
        self.R = R
        self.lzdim = len(R)
        self.Rp = np.array([p for p, q in R]); self.Rq = np.array([q for p, q in R])
        self.wa = self.occ_a[self.Rq] - self.occ_a[self.Rp]
        self.wb = self.occ_b[self.Rq] - self.occ_b[self.Rp]
        self.Rp_loc = np.array([self.Ploc[p] for p in self.Rp])
        self.Rq_loc = np.array([self.Hloc[q] for q in self.Rq])
        self.Rq_isC = np.array([q in set(self.Cidx) for q in self.Rq])
        # DF metric
        self.auxmol = mf.with_df.auxmol
        if self.auxmol is None:
            self.auxmol = df.addons.make_auxmol(self.mol, mf.with_df.auxbasis)
        self.naux = self.auxmol.nao_nr()
        if self.naux != ci.naux:
            sys.exit('MRSF gradients: auxiliary basis dimension mismatch (%d vs %d); '
                     'the Cholesky decomposition of the DF metric must not fall back '
                     'to an eigen-decomposition' % (self.naux, ci.naux))
        self.dims = (self.nocca, self.nvirb, self.naux)
        self.Lchol = scipy.linalg.cholesky(self.auxmol.intor('int2c2e', hermi=1), lower=True)
        # library: make sure the integrals of this CI object are loaded
        fock_mo = [self.fa, self.fb]
        mrsf_init.init(ci, fock_mo)
        mrsf_init.init_ints(ci, ci.eri_file())
        ierr = mrsf_grad.grad_init(ci.eri_file())
        if ierr != 0:
            sys.exit('mrsf_grad_init failed with error code %d' % ierr)
        self.Bso = mrsf_grad.bso(self.naux, self.nmo)
        # XC kernel cache
        if self.is_dft:
            ni = mf._numint
            self.ni = ni
            self._fxc = ni.cache_xc_kernel(self.mol, mf.grids, mf.xc, (C, C),
                                           (self.occ_a, self.occ_b), spin=1)[2]
        # reference gradient and one-electron derivative objects
        self.g0 = mf.nuc_grad_method()
        if self.is_dft:
            self.g0.grid_response = grad_obj.grid_response
        self.g0.verbose = 0
        t0 = time.time()
        self.grad_ref = self.g0.kernel()
        self.time_ref = time.time() - t0
        self.hcore_deriv = self.g0.hcore_generator(self.mol)
        self.s1 = self.g0.get_ovlp(self.mol)
        self.aoslices = self.mol.aoslice_by_atom()

    # ------------------------------------------------------------------
    # amplitudes and densities
    # ------------------------------------------------------------------
    def expanded_amplitudes(self, istate):
        """X~ (nvirb, nocca) of adiabatic state istate (0-based)"""
        irr, st = self.ci.state_sym(istate)
        x = self.ci.amps['adiabatic'][irr][:, st]
        X = np.array(x.reshape(self.nvirb, self.nocca, order='F'))
        pO1, pO2, hO1, hO2 = 0, 1, self.nC, self.nC + 1
        xoo = X[pO1, hO1]
        X[pO1, hO1] = xoo/SQ2
        if self.mult == 1:
            X[pO2, hO2] = -xoo/SQ2
        else:
            X[pO2, hO2] = xoo/SQ2
            X[pO1, hO2] = 0.0
            X[pO2, hO1] = 0.0
        return X

    def T_densities(self, Xt):
        nmo, H, P = self.nmo, self.H, self.P
        Ta = np.zeros((nmo, nmo)); Tb = np.zeros((nmo, nmo))
        Ta[np.ix_(H, H)] = -Xt.T @ Xt
        Tb[np.ix_(P, P)] = Xt @ Xt.T
        return Ta, Tb

    def z_to_local(self, z):
        """z (lzdim[, nvec]) -> Za, Zb (nvirb, nocca[, nvec]) with the
        occupation weights (n^s_q - n^s_p)"""
        z = np.asarray(z)
        single = (z.ndim == 1)
        if single:
            z = z[:, None]
        nvec = z.shape[1]
        Za = np.zeros((self.nvirb, self.nocca, nvec)); Zb = np.zeros((self.nvirb, self.nocca, nvec))
        Za[self.Rp_loc, self.Rq_loc, :] = self.wa[:, None]*z
        Zb[self.Rp_loc, self.Rq_loc, :] = self.wb[:, None]*z
        if single:
            return Za[:, :, 0], Zb[:, :, 0]
        return Za, Zb

    def pz_densities(self, z):
        pa = np.zeros((self.nmo, self.nmo)); pb = np.zeros((self.nmo, self.nmo))
        pa[self.Rp, self.Rq] = self.wa*z; pa[self.Rq, self.Rp] = self.wa*z
        pb[self.Rp, self.Rq] = self.wb*z; pb[self.Rq, self.Rp] = self.wb*z
        return pa, pb

    def vec_to_mat(self, z):
        Z = np.zeros((self.nmo, self.nmo))
        Z[self.Rp, self.Rq] = z
        return Z

    # ------------------------------------------------------------------
    # spin Fock response G^s[D] on the blocks needed
    # ------------------------------------------------------------------
    def fxc_blocks(self, dms):
        """f_xc[D] for AO densities dms (2, nvec, nao, nao) -> (2, nvec, nao, nao) AO"""
        v = self.ni.nr_uks_fxc(self.mol, self.mf.grids, self.mf.xc, None, dms, hermi=1, fxc=self._fxc)
        return np.asarray(v).reshape(2, -1, self.nao, self.nao)

    def G_full(self, Ahh, Za, Zb, jq_add=None, KbHP=None, KbHH=None, Pa_mo=None, Pb_mo=None):
        """G^a (all rows, columns H) and G^b (columns C) as nmo x nmo matrices"""
        Ghpa, Ghpb, Ghha, Ghhb = mrsf_grad.gfock(self.x_ref, self.dims, Ahh, Za, Zb, jq_add, True, True)
        Ghpa, Ghpb, Ghha, Ghhb = Ghpa[:, :, 0], Ghpb[:, :, 0], Ghha[:, :, 0], Ghhb[:, :, 0]
        if KbHP is not None:
            Ghpb = Ghpb - self.x_ref*KbHP
            Ghhb = Ghhb - self.x_ref*KbHH
        H, P = self.H, self.P
        Ga = np.zeros((self.nmo, self.nmo)); Gb = np.zeros((self.nmo, self.nmo))
        Ga[np.ix_(H, H)] = Ghha; Ga[np.ix_(P, H)] = Ghpa.T
        Gb[np.ix_(H, H)] = Ghhb; Gb[np.ix_(P, H)] = Ghpb.T
        if self.is_dft:
            C = self.C
            dm = np.array([[C @ Pa_mo @ C.T], [C @ Pb_mo @ C.T]])
            v = self.fxc_blocks(dm)
            Ga[:, H] += (C.T @ v[0, 0] @ C)[:, H]
            Gb[:, H] += (C.T @ v[1, 0] @ C)[:, H]
        return Ga, Gb

    # ------------------------------------------------------------------
    # pairing (spin-pair coupling) terms as rank-1 J-type pairs
    # ------------------------------------------------------------------
    def spc_terms(self, Xt):
        s, (kc, ko, kv) = self.s, self.kappa
        nC = self.nC
        O1, O2 = self.Oidx
        Cg, Vg = self.Cidx, self.Vidx
        xCO = [Xt[0, :nC], Xt[1, :nC]]
        xOV = [Xt[2:, nC], Xt[2:, nC + 1]]
        T = []
        T.append((s*kc,  (xCO[0], Cg, O2, 'col'), (xCO[0], Cg, O2, 'col')))
        T.append((s*kc,  (xCO[1], Cg, O1, 'col'), (xCO[1], Cg, O1, 'col')))
        T.append((-s*kc, (xCO[0], Cg, O1, 'col'), (xCO[1], Cg, O2, 'col')))
        T.append((-s*kc, (xCO[1], Cg, O2, 'col'), (xCO[0], Cg, O1, 'col')))
        T.append((s*ko,  (xOV[1], Vg, O1, 'col'), (xOV[1], Vg, O1, 'col')))
        T.append((s*ko,  (xOV[0], Vg, O2, 'col'), (xOV[0], Vg, O2, 'col')))
        T.append((-s*ko, (xOV[1], Vg, O2, 'col'), (xOV[0], Vg, O1, 'col')))
        T.append((-s*ko, (xOV[0], Vg, O1, 'col'), (xOV[1], Vg, O2, 'col')))
        T.append((2*s*kv,  (xCO[0], Cg, O1, 'col'), (xOV[1], Vg, O2, 'row')))
        T.append((-2*s*kv, (xCO[0], Cg, O2, 'col'), (xOV[1], Vg, O1, 'row')))
        T.append((-2*s*kv, (xCO[1], Cg, O1, 'col'), (xOV[0], Vg, O2, 'row')))
        T.append((2*s*kv,  (xCO[1], Cg, O2, 'col'), (xOV[0], Vg, O1, 'row')))
        return T

    def spec_matrix(self, spec):
        w, idx, o, orient = spec
        A = np.zeros((self.nmo, self.nmo))
        if orient == 'col':
            A[idx, o] = w
        else:
            A[o, idx] = w
        return A

    def spc_aux_vectors(self, Xt):
        out = []
        for t, specA, specB in self.spc_terms(Xt):
            vecs = []
            for (w, idx, o, orient) in (specA, specB):
                x = 0 if o == self.Oidx[0] else 1
                vecs.append(self.Bso[:, idx, x] @ w)
            out.append((t, specA, specB, vecs[0], vecs[1]))
        return out

    def full_from_blocks(self, MHH, MHP, MPP):
        M = np.zeros((self.nmo, self.nmo))
        H, P = self.H, self.P
        M[np.ix_(H, H)] = MHH
        M[np.ix_(H, P)] = MHP
        M[np.ix_(P, H)] = MHP.T
        M[np.ix_(P, P)] = MPP
        return M

    def gfock_spc(self, Xt, terms):
        """F^spc = sum_k t [J[B](A+A^T) + J[A](B+B^T)]"""
        jq = np.array([v for (_, _, _, a, b) in terms for v in (a, b)]).T
        JHH, JHP, JPP = mrsf_grad.jblocks(jq, self.dims)
        F = np.zeros((self.nmo, self.nmo))
        for k, (t, specA, specB, aQ, bQ) in enumerate(terms):
            JA = self.full_from_blocks(JHH[:, :, 2*k], JHP[:, :, 2*k], JPP[:, :, 2*k])
            JB = self.full_from_blocks(JHH[:, :, 2*k+1], JHP[:, :, 2*k+1], JPP[:, :, 2*k+1])
            A = self.spec_matrix(specA); B = self.spec_matrix(specB)
            F += t*(JB @ (A + A.T) + JA @ (B + B.T))
        return F

    def add_spc_to_families(self, Ghh, Fhp, terms):
        Hloc, Ploc, Cset = self.Hloc, self.Ploc, set(self.Cidx)
        for t, specA, specB, aQ, bQ in terms:
            for spec, cQ in ((specA, bQ), (specB, aQ)):
                w, idx, o, orient = spec
                for m, wm in zip(idx, w):
                    if m in Cset:
                        if orient == 'col':
                            Ghh[Hloc[m], Hloc[o], :] += 0.5*t*wm*cQ
                        else:
                            Ghh[Hloc[o], Hloc[m], :] += 0.5*t*wm*cQ
                    else:
                        Fhp[Hloc[o], Ploc[m], :] += 0.5*t*wm*cQ

    # ------------------------------------------------------------------
    # generalised Fock matrices
    # ------------------------------------------------------------------
    def cols_H(self, MH, MP):
        out = np.zeros((self.nmo, MH.shape[1]))
        out[self.H] = MH
        out[self.Vidx] = MP[2:]
        return out

    def gfock_omega(self, Xt, st, terms):
        fa, fb = self.fa, self.fb
        Ta, Tb = self.T_densities(Xt)
        Ta_loc = Ta[np.ix_(self.H, self.H)]
        F = 2.0*(fa @ Ta + fb @ Tb)
        Ga, Gb = self.G_full(Ta_loc, None, None, jq_add=st['jT'], KbHP=st['KbT_HP'], KbHH=st['KbT_HH'],
                             Pa_mo=Ta, Pb_mo=Tb)
        F[:, self.H] += 2.0*Ga[:, self.H]
        F[:, self.Cidx] += 2.0*Gb[:, self.Cidx]
        F[:, self.H] += self.cols_H(st['LaH'], st['LaP'])
        LbFull = np.zeros((self.nmo, self.nvirb))
        LbFull[self.H] = st['LbH']; LbFull[self.Vidx] = st['LbP'][2:]
        F[:, self.P] += LbFull
        F += self.gfock_spc(Xt, terms)
        return F

    def rhs(self, F):
        return F[self.Rp, self.Rq] - F[self.Rq, self.Rp]

    def gfock_z(self, z):
        fa, fb = self.fa, self.fb
        na, nb = self.occ_a, self.occ_b
        Z = self.vec_to_mat(z)
        Za, Zb = self.z_to_local(z)
        pa, pb = self.pz_densities(z)
        Ga, Gb = self.G_full(None, Za, Zb, Pa_mo=pa, Pb_mo=pb)
        F = np.zeros((self.nmo, self.nmo))
        F[:, self.H] += 2.0*Ga[:, self.H]
        F[:, self.Cidx] += 2.0*Gb[:, self.Cidx]
        for f, n in ((fa, na), (fb, nb)):
            F += 2.0*np.einsum('qs,s,ps->pq', Z, n, f)
            F += 2.0*np.einsum('rq,q,rp->pq', Z, n, f)
            F -= 2.0*np.einsum('rq,r,pr->pq', Z, n, f)
            F -= 2.0*np.einsum('qs,q,sp->pq', Z, n, f)
        return F

    # ------------------------------------------------------------------
    # Z-vector: operator and solvers
    # ------------------------------------------------------------------
    def hessian_apply(self, zmat):
        """H z for a batch of vectors zmat (lzdim, nvec)"""
        fa, fb, nmo = self.fa, self.fb, self.nmo
        nvec = zmat.shape[1]
        Za, Zb = self.z_to_local(zmat)
        Ghpa, Ghpb, _, _ = mrsf_grad.gfock(self.x_ref, self.dims, None, Za, Zb, None, True, False)
        if self.is_dft:
            C = self.C
            dms = np.empty((2, nvec, self.nao, self.nao))
            for v in range(nvec):
                pa, pb = self.pz_densities(zmat[:, v])
                dms[0, v] = C @ pa @ C.T
                dms[1, v] = C @ pb @ C.T
            vxc = self.fxc_blocks(dms)
            CH, CP = C[:, self.H], C[:, self.P]
            for v in range(nvec):
                Ghpa[:, :, v] += CH.T @ vxc[0, v] @ CP
                Ghpb[:, :, v] += CH.T @ vxc[1, v] @ CP
        out = np.zeros((self.lzdim, nvec))
        for v in range(nvec):
            K = np.zeros((nmo, nmo))
            K[self.Rp, self.Rq] += zmat[:, v]
            K[self.Rq, self.Rp] -= zmat[:, v]
            dfa = -K @ fa + fa @ K
            dfb = -K @ fb + fb @ K
            da = dfa[self.Rp, self.Rq] + Ghpa[self.Rq_loc, self.Rp_loc, v]
            db = dfb[self.Rp, self.Rq] + np.where(self.Rq_isC, Ghpb[self.Rq_loc, self.Rp_loc, v], 0.0)
            out[:, v] = 2.0*(self.wa*da + self.wb*db)
        return out

    def hessian_diag(self):
        fa, fb = self.fa, self.fb
        p, q = self.Rp, self.Rq
        return 2.0*(self.wa*(np.diag(fa)[p] - np.diag(fa)[q]) + self.wb*(np.diag(fb)[p] - np.diag(fb)[q]))

    def solve_pcg(self, Rmat, tol, maxiter):
        """block preconditioned conjugate gradients for H z = -R (lockstep,
        one operator call per iteration for the active vectors)"""
        n, nvec = Rmat.shape
        b = -Rmat
        d = self.hessian_diag()
        d = np.where(np.abs(d) < 1e-4, np.sign(d + 1e-300)*1e-4, d)
        x = np.zeros((n, nvec))
        r = b.copy()
        zp = r/d[:, None]
        p = zp.copy()
        rz = np.einsum('ik,ik->k', r, zp)
        active = np.arange(nvec)
        niter = np.zeros(nvec, int)
        resid = np.linalg.norm(r, axis=0)
        info = {'converged': resid < tol}
        for it in range(1, maxiter + 1):
            active = np.where(~info['converged'])[0]
            if len(active) == 0:
                break
            Ap = self.hessian_apply(p[:, active])
            for j, k in enumerate(active):
                pAp = p[:, k] @ Ap[:, j]
                alpha = rz[k]/pAp
                x[:, k] += alpha*p[:, k]
                r[:, k] -= alpha*Ap[:, j]
                resid[k] = np.linalg.norm(r[:, k])
                niter[k] = it
                if resid[k] < tol:
                    info['converged'][k] = True
                    continue
                znew = r[:, k]/d
                rznew = r[:, k] @ znew
                p[:, k] = znew + (rznew/rz[k])*p[:, k]
                rz[k] = rznew
            if self.verbose:
                with output.output_file(output.file_names['out_file'], 'a+') as f:
                    f.write('\n  Z-vector iteration %3d: max residual %.3e (%d active)' %
                            (it, resid[active].max(), len(active)))
        info['niter'] = niter
        info['resid'] = resid
        return x, info

    def solve_dense(self, Rmat):
        n = self.lzdim
        Hm = self.hessian_apply(np.eye(n))
        Hm = _sym(Hm)
        x = np.linalg.solve(Hm, -Rmat)
        resid = np.linalg.norm(Hm @ x + Rmat, axis=0)
        return x, {'converged': np.ones(Rmat.shape[1], bool), 'niter': np.zeros(Rmat.shape[1], int),
                   'resid': resid}

    # ------------------------------------------------------------------
    # the gradient of a list of states
    # ------------------------------------------------------------------
    def run(self, states):
        """states: list of adiabatic state indices (0-based); returns dict"""
        mol, C = self.mol, self.C
        nst = len(states)
        times = {}
        t0 = time.time()
        # pass 1: per-state RHS pieces
        Xts, terms_all, Rmat = [], [], np.zeros((self.lzdim, nst))
        for k, ist in enumerate(states):
            Xt = self.expanded_amplitudes(ist)
            st = mrsf_grad.state(self.x_resp, Xt, self.dims)
            terms = self.spc_aux_vectors(Xt)
            Fo = self.gfock_omega(Xt, st, terms)
            Rmat[:, k] = self.rhs(Fo)
            Xts.append(Xt); terms_all.append(terms)
            del st
        times['rhs'] = time.time() - t0
        # Z-vectors for all states at once
        t0 = time.time()
        if self.opt.zvec_solver == 'dense':
            Z, info = self.solve_dense(Rmat)
        else:
            Z, info = self.solve_pcg(Rmat, self.opt.zvec_tol, self.opt.zvec_iter)
        times['zvec'] = time.time() - t0
        for k, ist in enumerate(states):
            if not info['converged'][k]:
                output.print_message('WARNING: Z-vector of state %d not converged (residual %.3e)' %
                                     (ist + 1, info['resid'][k]))
        # pass 2: families, energy-weighted density, AO assembly
        t0 = time.time()
        fam_states, W_list, P_list = [], [], []
        for k, ist in enumerate(states):
            Xt = Xts[k]
            z = Z[:, k]
            st = mrsf_grad.state(self.x_resp, Xt, self.dims)
            terms = terms_all[k]
            Fo = self.gfock_omega(Xt, st, terms)
            Fz = self.gfock_z(z)
            W = 0.5*C @ _sym(Fo + Fz) @ C.T
            Ta, Tb = self.T_densities(Xt)
            pa, pb = self.pz_densities(z)
            Pa, Pb = Ta + pa, Tb + pb
            Za, Zb = self.z_to_local(z)
            self.add_spc_to_families(st['Ghh'], st['Fhp'], terms)
            g, pq = mrsf_grad.finish(self.x_ref, st, Ta[np.ix_(self.H, self.H)], Xt, Za, Zb, self.occH, self.dims)
            fam_states.append({'Xt': Xt, 'Ghh': st['Ghh'], 'Fhp': st['Fhp'], 'Yf': st['Yf'], 'g': g})
            W_list.append(W); P_list.append((C @ Pa @ C.T, C @ Pb @ C.T))
            del st
        times['families'] = time.time() - t0
        t0 = time.time()
        de2 = mrsf_dfgrad.df_grad_factorised(mol, self.auxmol, self.Lchol, C[:, self.H], C[:, self.P],
                                             fam_states, max_memory=max(500, 1000*self.ci.mem_budget))
        times['2e_ao'] = time.time() - t0
        t0 = time.time()
        grads = np.zeros((nst, mol.natm, 3))
        parts = []
        if self.is_dft:
            dexc_all = self.grad_xc_probe_multi(P_list)
        else:
            dexc_all = [np.zeros((mol.natm, 3)) for _ in range(nst)]
        for k in range(nst):
            Pa_ao, Pb_ao = P_list[k]
            Pt = Pa_ao + Pb_ao
            de1 = np.zeros((mol.natm, 3))
            deS = np.zeros((mol.natm, 3))
            for A, (s0, s1, p0, p1) in enumerate(self.aoslices):
                de1[A] = np.einsum('xpq,pq->x', self.hcore_deriv(A), Pt)
                deS[A] = -2.0*np.einsum('xpq,pq->x', self.s1[:, p0:p1], W_list[k][p0:p1])
            dexc = dexc_all[k]
            total = self.grad_ref + de1 + de2[k] + deS + dexc
            if mol.symmetry:
                total = rhf_grad.symmetrize(mol, total)
            grads[k] = total
            parts.append({'ref': self.grad_ref, '1e': de1, '2e': de2[k], 'S': deS, 'xc': dexc})
        times['1e_xc'] = time.time() - t0
        times['ref'] = self.time_ref
        return {'grad': grads, 'grad_ref': self.grad_ref, 'zvec_niter': info['niter'],
                'zvec_resid': info['resid'], 'zvec_converged': info['converged'],
                'parts': parts, 'times': times, 'P_ao': P_list}

    def grad_xc_probe(self, Pa_ao, Pb_ao):
        """XC probe term of one state through PySCF's TD-UKS kernel routine
        (reference implementation; grad_xc_probe_multi is used in production)"""
        mf, mol = self.mf, self.mol
        shim = _XCShim(mf, self.C, self.occ_a, self.occ_b)
        f1vo, f1oo, vxc1, k1ao = tduks_grad._contract_xc_kernel(shim, mf.xc, (Pa_ao, Pb_ao), None, True, False,
                                                                 max(2000, 1000*self.ci.mem_budget))
        return self._xc_probe_contract([f1vo], vxc1, [(Pa_ao, Pb_ao)])[0]

    def _xc_probe_contract(self, f1vo_list, vxc1, P_list):
        C = self.C
        Da = (C*self.occ_a) @ C.T
        Db = (C*self.occ_b) @ C.T
        aoslices = self.aoslices[:, 2]
        out = []
        for f1vo, (Pa_ao, Pb_ao) in zip(f1vo_list, P_list):
            t = np.zeros((3, self.nao))
            for s, (P, D) in enumerate(((Pa_ao, Da), (Pb_ao, Db))):
                t += 2.0*np.einsum('xpq,pq->xp', vxc1[s, 1:], P)
                t += 2.0*np.einsum('xpq,pq->xp', f1vo[s, 1:], D)
            out.append(np.add.reduceat(t, aoslices, axis=1).T)
        return out

    def grad_xc_probe_multi(self, P_list):
        """XC probe terms of several states in one pass over the grid:
        Tr(P^s dV_xc^s/dR) at fixed grid = contraction of P with the
        AO-derivative XC potential matrices plus the kernel potential of P
        contracted with the AO derivative of the reference density
        (same algebra as pyscf.grad.tduks._contract_xc_kernel)."""
        from pyscf.grad import tdrks as tdrks_grad
        mf, mol, ni = self.mf, self.mol, self.ni
        grids = mf.grids
        xctype = ni._xc_type(mf.xc)
        nao, C = self.nao, self.C
        nst = len(P_list)
        shls_slice = (0, mol.nbas)
        ao_loc = mol.ao_loc_nr()
        if xctype == 'LDA':
            fmat_, ao_deriv = tdrks_grad._lda_eval_mat_, 1
        elif xctype == 'GGA':
            fmat_, ao_deriv = tdrks_grad._gga_eval_mat_, 2
        else:
            raise NotImplementedError('MRSF gradients: XC type %s' % xctype)
        dms = [(0.5*(Pa + Pa.T), 0.5*(Pb + Pb.T)) for Pa, Pb in P_list]
        f1vo = [np.zeros((2, 4, nao, nao)) for _ in range(nst)]
        vxc1 = np.zeros((2, 4, nao, nao))
        max_memory = max(2000, 1000*self.ci.mem_budget)
        for ao, mask, weight, coords in ni.block_loop(mol, grids, nao, ao_deriv, max_memory):
            ao0 = ao[0] if xctype == 'LDA' else ao
            rho = (ni.eval_rho2(mol, ao0, C, self.occ_a, mask, xctype, with_lapl=False),
                   ni.eval_rho2(mol, ao0, C, self.occ_b, mask, xctype, with_lapl=False))
            vxc, fxc = ni.eval_xc_eff(mf.xc, rho, 2, xctype=xctype)[1:3]
            wv = vxc*weight
            fmat_(mol, vxc1[0], ao, wv[0], mask, shls_slice, ao_loc)
            fmat_(mol, vxc1[1], ao, wv[1], mask, shls_slice, ao_loc)
            for k in range(nst):
                rho1 = np.asarray((ni.eval_rho(mol, ao0, dms[k][0], mask, xctype, hermi=1, with_lapl=False),
                                   ni.eval_rho(mol, ao0, dms[k][1], mask, xctype, hermi=1, with_lapl=False)))
                if xctype == 'LDA':
                    rho1 = rho1[:, np.newaxis]
                wv = np.einsum('axg,axbyg,g->byg', rho1, fxc, weight)
                fmat_(mol, f1vo[k][0], ao, wv[0], mask, shls_slice, ao_loc)
                fmat_(mol, f1vo[k][1], ao, wv[1], mask, shls_slice, ao_loc)
        vxc1[:, 1:] *= -1
        for k in range(nst):
            f1vo[k][:, 1:] *= -1
        return self._xc_probe_contract(f1vo, vxc1, P_list)

    def finalise(self):
        mrsf_grad.grad_free()
