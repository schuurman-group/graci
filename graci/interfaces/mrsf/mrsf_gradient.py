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
Frozen core: omega is not invariant under rotations between the frozen
(k) and active (m) doubly occupied MOs; the canonical condition
c_km = (f^a + f^b)_km = 0 fixes them with multipliers
  zeta_km = -Q_km/h_km,  Q_km = F^omega_km - F^omega_mk,  h_km = (f^a + f^b)_kk - (f^a + f^b)_mm,
which add F^zeta (generalised Fock of sum zeta c) to the Z-vector right-hand
side and to W, and p^zeta = Zhat + Zhat^T (Zhat_km = zeta_km/2, both spins)
to the relaxed densities.
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
import graci.interfaces.mrsf.mrsf_xc as mrsf_xc
import graci.interfaces.mrsf.mrsf_hessian as mrsf_hessian

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
        # frozen core of the CI object (library rule: the nfc doubly occupied
        # MOs lowest in orbital energy) and the canonical-condition pairs
        self.nfc = int(getattr(ci, 'frozen_core', 0) or 0)
        emo = np.asarray(ci.emo, float)[:self.nmo]
        self.Cf = sorted(int(self.Cidx[h]) for h in np.argsort(emo[self.Cidx])[:self.nfc]) \
            if self.nfc > 0 else []
        self.Ca = [int(c) for c in self.Cidx if int(c) not in set(self.Cf)]
        self.N = [(k, m) for k in self.Cf for m in self.Ca]
        fsum = self.fa + self.fb
        self.hN = np.array([fsum[k, k] - fsum[m, m] for (k, m) in self.N])
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
        # the gradient session of the library is a standard (non-extended)
        # one with the full core: the amplitudes of the state are supplied
        # from Python and the gradient kernels work on the MRSF layout
        mrsf_init.set_extended(False)
        mrsf_init.set_frozen_core(0)
        mrsf_init.init(ci, fock_mo)
        mrsf_init.set_exchange(getattr(ci, 'exchange', 'auto'))
        mrsf_init.init_ints(ci, ci.eri_file())
        # Z-vector operator = orbital Hessian of the reference (hole-particle
        # DF block, Fock couplings and diagonal in the library; the XC
        # kernel on the grid, where the library grid module caches the AO
        # values, weights and libxc derivatives once)
        try:
            self.hess = mrsf_hessian.RoksHessian(mf, C, occ, self.fa, self.fb, self.x_ref,
                                                 ci.eri_file(), self.dims, grad_obj.mem_budget)
        except RuntimeError as err:
            sys.exit(str(err))
        if self.hess.lz != self.lzdim:
            sys.exit('MRSF gradients: inconsistent rotation space')
        self.Bso = mrsf_grad.bso(self.naux, self.nmo)
        self.dq_ref = mrsf_grad.dq(self.occ, self.naux)        # Coulomb vector of the reference density
        self.hdiag = self.hess.diag()
        if self.is_dft:
            self.ni = mf._numint
            self.xcgrid = self.hess.xcgrid
        # one-electron derivative integrals and the reference pieces that do
        # not go through the DF families: nuclear repulsion, the reference
        # densities, W_ref = sum_s D^s F^s D^s and the reference XC gradient
        t0 = time.time()
        g0 = mf.nuc_grad_method()          # only used for its derivative-integral helpers
        g0.verbose = 0
        self.hcore_deriv = g0.hcore_generator(self.mol)
        self.s1 = g0.get_ovlp(self.mol)
        self.aoslices = self.mol.aoslice_by_atom()
        self.grad_nuc = g0.grad_nuc(self.mol)
        self.Da = (C*self.occ_a) @ C.T
        self.Db = (C*self.occ_b) @ C.T
        self.W_ref = self.Da @ fock_ao[0] @ self.Da + self.Db @ fock_ao[1] @ self.Db
        # reference XC gradient: from the grid probe pass (in run) unless
        # the quadrature-weight response is requested (PySCF only)
        self.dexc_ref = np.zeros((self.mol.natm, 3))
        if self.is_dft and grad_obj.grid_response:
            from pyscf.grad import uks as uks_grad
            excsum, vmat = uks_grad.get_vxc_full_response(self.ni, self.mol, mf.grids, mf.xc,
                                                          np.array([self.Da, self.Db]))
            t = np.zeros((3, self.nao))
            for sdm, D in ((0, self.Da), (1, self.Db)):
                t += 2.0*np.einsum('xpq,pq->xp', vmat[sdm], D)
            self.dexc_ref = excsum + np.add.reduceat(t, self.aoslices[:, 2], axis=1).T
        elif self.is_dft:
            self.dexc_ref = None
        self.time_ref = time.time() - t0

    # ------------------------------------------------------------------
    # amplitudes and densities
    # ------------------------------------------------------------------
    def expanded_amplitudes(self, istate):
        """X~ (nvirb, nocca) of adiabatic state istate (0-based)"""
        irr, st = self.ci.state_sym(istate)
        x = self.ci.amps['adiabatic'][irr][:, st]
        X = np.array(x.reshape(self.nvirb, self.nocca, order='F'))
        if self.Cf:
            fcols = [self.Hloc[k] for k in self.Cf]
            if np.abs(X[:, fcols]).max() > 1e-12:
                sys.exit('MRSF gradients: amplitudes on frozen holes (inconsistent frozen core)')
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
    def G_full(self, Ahh, Za, Zb, jq_add=None, KbHP=None, KbHH=None, xcfac=None):
        """G^a (all rows, columns H) and G^b (columns C) as nmo x nmo
        matrices; xcfac = (Lfac, Rfac, rch) are the factor pairs of the
        AO densities for the XC kernel (one vector, see mrsf_xc.XCGrid.kernel)"""
        Ghpa, Ghpb, Ghha, Ghhb = mrsf_grad.gfock(self.x_ref, self.dims, Ahh, Za, Zb, jq_add, True, True)
        Ghpa, Ghpb, Ghha, Ghhb = Ghpa[:, :, 0], Ghpb[:, :, 0], Ghha[:, :, 0], Ghhb[:, :, 0]
        if KbHP is not None:
            Ghpb = Ghpb - self.x_ref*KbHP
            Ghhb = Ghhb - self.x_ref*KbHH
        H, P, V = self.H, self.P, self.Vidx
        Ga = np.zeros((self.nmo, self.nmo)); Gb = np.zeros((self.nmo, self.nmo))
        Ga[np.ix_(H, H)] = Ghha; Ga[np.ix_(P, H)] = Ghpa.T
        Gb[np.ix_(H, H)] = Ghhb; Gb[np.ix_(P, H)] = Ghpb.T
        if self.is_dft and xcfac is not None:
            Lfac, Rfac, rch = xcfac
            VHP, VHH = self.xcgrid.kernel(Lfac, Rfac, rch, True, True)
            Ga[np.ix_(H, H)] += VHH[:, :, 0, 0]; Ga[np.ix_(V, H)] += VHP[:, 2:, 0, 0].T
            Gb[np.ix_(H, H)] += VHH[:, :, 1, 0]; Gb[np.ix_(V, H)] += VHP[:, 2:, 1, 0].T
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
        xcfac = None
        if self.is_dft:
            # T^a = C_H Ta C_H^T, T^b = (C_P X~)(C_P X~)^T as factor pairs
            CH, CP = self.C[:, self.H], self.C[:, self.P]
            Lf = np.zeros((self.nao, self.nocca, 2, 1)); Rf = np.zeros((self.nao, self.nocca, 2, 1))
            Lf[:, :, 0, 0] = 0.5*(CH @ Ta_loc); Rf[:, :, 0, 0] = CH
            Lf[:, :, 1, 0] = 0.5*(CP @ Xt);     Rf[:, :, 1, 0] = CP @ Xt
            xcfac = (Lf, Rf, [[1], [0]])
        Ga, Gb = self.G_full(Ta_loc, None, None, jq_add=st['jT'], KbHP=st['KbT_HP'], KbHH=st['KbT_HH'],
                             xcfac=xcfac)
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
        xcfac = None
        if self.is_dft:
            # p^z,s = C_P Z_s C_H^T + h.c. as factor pairs
            CH, CP = self.C[:, self.H], self.C[:, self.P]
            Lf = np.zeros((self.nao, self.nocca, 2, 1)); Rf = np.zeros((self.nao, self.nocca, 2, 1))
            Lf[:, :, 0, 0] = CP @ Za; Rf[:, :, 0, 0] = CH
            Lf[:, :, 1, 0] = CP @ Zb; Rf[:, :, 1, 0] = CH
            xcfac = (Lf, Rf, [[1], [1]])
        Ga, Gb = self.G_full(None, Za, Zb, xcfac=xcfac)
        F = np.zeros((self.nmo, self.nmo))
        F[:, self.H] += 2.0*Ga[:, self.H]
        F[:, self.Cidx] += 2.0*Gb[:, self.Cidx]
        for f, n in ((fa, na), (fb, nb)):
            # 2 [ sum_s f_ps n_s Z_qs + n_q sum_r f_rp Z_rq - sum_r f_pr n_r Z_rq - n_q sum_s f_sp Z_qs ]
            fn = f*n
            F += 2.0*(fn @ Z.T + (f.T @ Z)*n - fn @ Z - (Z @ f).T*n)
        return F

    # ------------------------------------------------------------------
    # canonical multipliers of the frozen core
    # ------------------------------------------------------------------
    def zeta(self, Fo):
        """zeta_km = -Q_km/h_km on the frozen-active core pairs"""
        if not self.N:
            return np.zeros(0)
        Q = np.array([Fo[k, m] - Fo[m, k] for (k, m) in self.N])
        return -Q/self.hN

    def zeta_density(self, zeta):
        """p^zeta = Zhat + Zhat^T, Zhat_km = zeta_km/2 (MO, both spins)"""
        p = np.zeros((self.nmo, self.nmo))
        for (k, m), zt in zip(self.N, zeta):
            p[k, m] += 0.5*zt
            p[m, k] += 0.5*zt
        return p

    def gfock_zeta(self, zeta):
        """generalised Fock of sum_km zeta_km (f^a + f^b)_km:
        F = 2 sum_s [G^s[p] n^s + f^s p],  p = p^zeta (core-core, both spins),
        G^a = G^b = J[2p] - x K[p] + f_xc^s[p, p]"""
        p = self.zeta_density(zeta)
        F = 2.0*(self.fa + self.fb) @ p
        H, P, V = self.H, self.P, self.Vidx
        p_loc = p[np.ix_(H, H)]
        # gfock with D^a = p, D^b = 0: G^a = J[p] - x K[p], G^b = J[p]; their sum
        # is the spin Fock response of D^a = D^b = p
        Ghpa, Ghpb, Ghha, Ghhb = mrsf_grad.gfock(self.x_ref, self.dims, p_loc, None, None, None, True, True)
        G = np.zeros((self.nmo, self.nmo))
        G[np.ix_(H, H)] = Ghha[:, :, 0] + Ghhb[:, :, 0]
        G[np.ix_(P, H)] = (Ghpa[:, :, 0] + Ghpb[:, :, 0]).T
        Ga, Gb = G, G.copy()
        if self.is_dft:
            CH = self.C[:, H]
            Lf = np.zeros((self.nao, self.nocca, 2, 1)); Rf = np.zeros((self.nao, self.nocca, 2, 1))
            for s in (0, 1):
                Lf[:, :, s, 0] = 0.5*(CH @ p_loc); Rf[:, :, s, 0] = CH
            VHP, VHH = self.xcgrid.kernel(Lf, Rf, [[1], [1]], True, True)
            Ga[np.ix_(H, H)] += VHH[:, :, 0, 0]; Ga[np.ix_(V, H)] += VHP[:, 2:, 0, 0].T
            Gb[np.ix_(H, H)] += VHH[:, :, 1, 0]; Gb[np.ix_(V, H)] += VHP[:, 2:, 1, 0].T
        F[:, H] += 2.0*Ga[:, H]
        F[:, self.Cidx] += 2.0*Gb[:, self.Cidx]
        return F

    def zeta_exchange_fix(self, Ghh, p):
        """the families are finished with T^a + 2 p (Coulomb of both spins and
        alpha exchange of 2p, i.e. with D^a = 1_H); the beta exchange pair
        (p, D^b = 1_C) only has core columns: add back x/2 (p B_{.,O} + h.c.)
        (p is a hole-hole matrix: frozen-active core and O1-O2 pairs)"""
        H = self.H
        pHH = p[np.ix_(H, H)]
        for x, o in enumerate(self.Oidx):
            v = 0.5*self.x_ref*(pHH @ self.Bso[:, H, x].T)              # (nocca, naux)
            ho = self.Hloc[o]
            Ghh[:, ho, :] += v
            Ghh[ho, :, :] += v

    # ------------------------------------------------------------------
    # Z-vector: operator and solvers
    # ------------------------------------------------------------------
    def hessian_apply(self, zmat):
        """H z for a batch of vectors zmat (lzdim, nvec): the orbital Hessian
        of the reference (mrsf_hessian.RoksHessian)"""
        return self.hess.apply(zmat)

    def hessian_diag(self):
        return self.hdiag

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
    # per-state stages (the extended-method driver overrides these)
    # ------------------------------------------------------------------
    def exchange_pass(self, sd):
        """exchange-channel pass of the state (the large families)"""
        sd['st'] = mrsf_grad.state(self.x_resp, sd['Xt'], self.dims, self.dq_ref)

    def prepare_state(self, ist):
        """pass 1 of a state: amplitudes, exchange pass, generalised Fock,
        canonical multipliers of the frozen core; returns the state dict"""
        sd = {'ist': ist, 'Xt': self.expanded_amplitudes(ist)}
        self.exchange_pass(sd)
        sd['terms'] = self.spc_aux_vectors(sd['Xt'])
        sd['Fo'] = self.gfock_omega(sd['Xt'], sd['st'], sd['terms'])
        sd['zeta'] = self.zeta(sd['Fo'])
        sd['Fzeta'] = self.gfock_zeta(sd['zeta']) if self.N else np.zeros((self.nmo, self.nmo))
        return sd

    def extra_families(self, sd):
        """hook: additional family pieces (in place in sd['st']); returns
        the list of extra particle-particle family entries of the state"""
        return []

    def finish_state(self, sd, z):
        """pass 2 of a state: energy-weighted density, relaxed densities,
        families and probe factors; returns (entries, W, (Pa_ao, Pb_ao), fac)"""
        C = self.C
        Xt, st, terms, Fo = sd['Xt'], sd['st'], sd['terms'], sd['Fo']
        Fz = self.gfock_z(z)
        W = 0.5*C @ _sym(Fo + Fz + sd['Fzeta']) @ C.T
        Ta, Tb = self.T_densities(Xt)
        pa, pb = self.pz_densities(z)
        pzeta = self.zeta_density(sd['zeta'])
        Pa, Pb = Ta + pa + pzeta, Tb + pb + pzeta
        Za, Zb = self.z_to_local(z)
        self.add_spc_to_families(st['Ghh'], st['Fhp'], terms)
        Ta_fam = Ta[np.ix_(self.H, self.H)]
        if self.N:
            self.zeta_exchange_fix(st['Ghh'], pzeta)
            Ta_fam = Ta_fam + 2.0*pzeta[np.ix_(self.H, self.H)]
        extra = self.extra_families(sd)
        g, pq = mrsf_grad.finish(self.x_ref, st, Ta_fam, Xt, Za, Zb, self.occH, self.dims)
        entries = [{'Xt': Xt, 'Ghh': st['Ghh'], 'Fhp': st['Fhp'], 'Yf': st['Yf'], 'g': g}] + extra
        fac = self.P_factors(Xt, Za, Zb, pzeta)
        return entries, W, (C @ Pa @ C.T, C @ Pb @ C.T), fac

    def probe_g_factors(self, sd):
        """hook: factor pairs (LM, RM, LS) of the G-density probe terms of
        the extended method; None for standard MRSF"""
        return None

    def extra_density(self, sd):
        """hook: additional MO density contracted with the core-Hamiltonian
        derivative (the extended method's F^G density M); None otherwise"""
        return None

    # ------------------------------------------------------------------
    # the gradient of a list of states
    # ------------------------------------------------------------------
    def run(self, states):
        """states: list of adiabatic state indices (0-based); returns dict"""
        mol, C = self.mol, self.C
        nst = len(states)
        times = {}
        t0 = time.time()
        # pass 1: per-state RHS pieces; the generalised Fock matrices are
        # kept, the (large) families are kept only when they fit the budget
        fam_bytes = 8.0*self.naux*(self.nocca*self.nocca + 2*self.nocca*self.nvirb + self.nocca*self.nocca)
        keep_fam = nst*fam_bytes <= 0.5*self.ci.mem_budget*1.0e9
        sds = []
        Rmat = np.zeros((self.lzdim, nst))
        for k, ist in enumerate(states):
            sd = self.prepare_state(ist)
            Rmat[:, k] = self.rhs(sd['Fo'] + sd['Fzeta'])
            st_dq = np.array(sd['st']['dq'])
            if not keep_fam:
                sd['st'] = None
            sds.append(sd)
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
        fam_states, entry_state, W_list, P_list, fac_list, gfac_list, Mx_list = [], [], [], [], [], [], []
        for k, ist in enumerate(states):
            sd = sds[k]
            if sd['st'] is None:
                self.exchange_pass(sd)
            entries, W, P, fac = self.finish_state(sd, Z[:, k])
            fam_states += entries
            entry_state += [k]*len(entries)
            W_list.append(W); P_list.append(P); fac_list.append(fac)
            gfac_list.append(self.probe_g_factors(sd))
            Mx_list.append(self.extra_density(sd))
            sd['st'] = None
        times['families'] = time.time() - t0
        # reference two-electron families (first entry of the AO pass)
        Ghh_ref, g_ref = mrsf_grad.reffam(self.x_ref, st_dq, self.occH, self.dims)
        fam_states.insert(0, {'Xt': np.zeros((self.nvirb, self.nocca)), 'Ghh': Ghh_ref,
                              'Fhp': np.zeros((self.nocca, self.nvirb, self.naux), order='F'),
                              'Yf': np.zeros((self.nvirb, self.nocca, self.naux), order='F'), 'g': g_ref})
        t0 = time.time()
        de2_all = mrsf_dfgrad.df_grad_factorised(mol, self.auxmol, self.Lchol, C[:, self.H], C[:, self.P],
                                                 fam_states, max_memory=max(500, 1000*self.opt.mem_budget))
        de2_ref = de2_all[0]
        de2 = [np.zeros((mol.natm, 3)) for _ in range(nst)]
        for e, k in enumerate(entry_state):
            de2[k] += de2_all[e + 1]
        times['2e_ao'] = time.time() - t0
        # XC probe terms of the states and the fixed-grid reference XC
        # gradient in one pass over the grid
        t0 = time.time()
        if self.is_dft:
            dexc_ref_probe, dexc_all = self.grad_xc_probe_multi(fac_list, gfac_list)
            if self.dexc_ref is None:
                self.dexc_ref = dexc_ref_probe
        else:
            dexc_all = [np.zeros((mol.natm, 3)) for _ in range(nst)]
        # reference gradient
        de1_ref = np.zeros((mol.natm, 3)); deS_ref = np.zeros((mol.natm, 3))
        Dt = self.Da + self.Db
        for A, (s0, s1, p0, p1) in enumerate(self.aoslices):
            de1_ref[A] = np.einsum('xpq,pq->x', self.hcore_deriv(A), Dt)
            deS_ref[A] = -2.0*np.einsum('xpq,pq->x', self.s1[:, p0:p1], self.W_ref[p0:p1])
        self.grad_ref = self.grad_nuc + de1_ref + de2_ref + deS_ref + self.dexc_ref
        grads = np.zeros((nst, mol.natm, 3))
        parts = []
        for k in range(nst):
            Pa_ao, Pb_ao = P_list[k]
            Pt = Pa_ao + Pb_ao
            if Mx_list[k] is not None:
                Pt = Pt + C @ Mx_list[k] @ C.T
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
        if self.is_dft:
            times['xc_kernel'] = self.xcgrid.time_kernel
            times['xc_calls'] = self.xcgrid.ncalls_kernel
            times['xc_probe'] = self.xcgrid.time_probe
            times['xc_cached'] = self.xcgrid.cache
        # canonical multipliers: max |zeta| over the frozen-active core pairs per state and, for the
        # extended method (last pair of N), the multiplier of the O1-O2 rotation
        somo_gap = getattr(self, 'somo_gap', None)
        nsp = 1 if somo_gap is not None else 0
        zeta_max = np.array([np.abs(sd['zeta'][:len(sd['zeta'])-nsp]).max() if len(sd['zeta']) > nsp else 0.0
                             for sd in sds])
        zeta_12 = np.array([sd['zeta'][-1] for sd in sds]) if nsp else np.zeros(len(sds))
        self.debug = {'sds': sds, 'Z': Z, 'W': W_list, 'P': P_list}
        return {'grad': grads, 'grad_ref': self.grad_ref, 'zvec_niter': info['niter'],
                'zvec_resid': info['resid'], 'zvec_converged': info['converged'],
                'parts': parts, 'times': times, 'P_ao': P_list, 'zeta_max': zeta_max,
                'zeta_12': zeta_12, 'somo_gap': somo_gap}

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

    def P_factors(self, Xt, Za, Zb, pzeta=None):
        """relaxed densities P^s = T^s + p^{z,s} (+ p^zeta) as AO factor pairs
        (L, R), P = L R^T + R L^T with k = 2 nocca columns: alpha = [C_H Ta/2 | C_P Za]
        against [C_H | C_H], beta = [C_P X~/2 | C_P Zb] against [C_P X~ | C_H];
        the core-core p^zeta adds C_H p/2 to the columns paired with C_H"""
        CH, CP = self.C[:, self.H], self.C[:, self.P]
        nao, nocca = self.nao, self.nocca
        L = np.zeros((nao, 2*nocca, 2)); R = np.zeros((nao, 2*nocca, 2))
        L[:, :nocca, 0] = -0.5*(CH @ (Xt.T @ Xt)); R[:, :nocca, 0] = CH
        L[:, nocca:, 0] = CP @ Za;                 R[:, nocca:, 0] = CH
        CPX = CP @ Xt
        L[:, :nocca, 1] = 0.5*CPX;                 R[:, :nocca, 1] = CPX
        L[:, nocca:, 1] = CP @ Zb;                 R[:, nocca:, 1] = CH
        if pzeta is not None and self.N:
            add = 0.5*(CH @ pzeta[np.ix_(self.H, self.H)])
            L[:, :nocca, 0] += add
            L[:, nocca:, 1] += add
        return L, R

    def grad_xc_probe_multi(self, fac_list, gfac_list=None):
        """XC probe terms of all states and the reference XC gradient on the
        fixed grid through the library grid module; fac_list holds per state
        the factor pairs (L, R) of P^s (P_factors), gfac_list the G-density
        factor pairs (LM, RM, LS) of the extended method (or None entries).
        Returns (dexc_ref, [dexc_k])."""
        nst = len(fac_list)
        k = fac_list[0][0].shape[1]
        Lf = np.empty((self.nao, k, 2, nst), order='F'); Rf = np.empty((self.nao, k, 2, nst), order='F')
        for v, (L, R) in enumerate(fac_list):
            Lf[:, :, :, v] = L; Rf[:, :, :, v] = R
        gfac = None
        if gfac_list is not None and any(gf is not None for gf in gfac_list):
            kM = gfac_list[0][0].shape[1]; kS = gfac_list[0][2].shape[1]
            LM = np.empty((self.nao, kM, nst), order='F'); RM = np.empty((self.nao, kM, nst), order='F')
            LS = np.empty((self.nao, kS, nst), order='F')
            for v, gf in enumerate(gfac_list):
                LM[:, :, v], RM[:, :, v], LS[:, :, v] = gf
            gfac = (LM, RM, LS)
        if gfac is None:
            t = self.xcgrid.probe(Lf, Rf, np.zeros((2, nst), dtype=np.int32))
            tG = None
        else:
            t, tG = self.xcgrid.probe(Lf, Rf, np.zeros((2, nst), dtype=np.int32), gfac=gfac)
        de = -2.0*np.add.reduceat(t, self.aoslices[:, 2], axis=0)      # (natm, 3, nst+1)
        out = [de[:, :, v + 1] for v in range(nst)]
        if tG is not None:
            deG = -2.0*np.add.reduceat(tG, self.aoslices[:, 2], axis=0)
            out = [out[v] + deG[:, :, v] for v in range(nst)]
        return de[:, :, 0], out

    def grad_xc_probe_multi_pyscf(self, P_list):
        """reference implementation with PySCF's grid routines (unused):
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
        mrsf_grad.aograd_free()
        self.hess.free()
