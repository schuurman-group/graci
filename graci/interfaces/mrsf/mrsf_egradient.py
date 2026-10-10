"""
Analytic nuclear gradients of extended MRSF-TDDFT (EMRSF) states: the
driver of the standard method (mrsf_gradient.GradientDriver) with the
additional terms of the extended energy functional (plan Part XII; the
NumPy reference is ~/calculations/mrsf_dev/emrsf_grad_ref.py).

omega_E(v; C) = v^T A_E v, v = (x, y), A_E = [[A_M, c_cp Cpl], [c_cp Cpl^T, A_CV + A_G]]
  = omega_M(x; C)                                 (standard machinery, expanded x)
  + Tr(f^a T^a_G) + Tr(f^b T^b_G)                 (A_G Fock part: T^b_{O1O1} += |y|^2, T^a_{O2O2} -= |y|^2)
  + Tr(F^G M),  M = T^cv + sym(N_F)               (closed-shell KS matrix of G; T^cv_VV = Y Y^T,
                                                   T^cv_CC = -Y^T Y; N_F the Fock-type couplings)
  + two-electron pairs                            (CV Coulomb and exchange, covariant F' correction,
                                                   -c_H (O2O2|O1O1)|y|^2, two-electron couplings)
  + E_K = 1/2 int rho(D_s) f^(0)[rho_G] rho(D_s)  (singlet DFT; D_s = C (Yh + Yh^T) C^T)
Canonical multipliers (frozen-active core, O1-O2): mrsf_gradient.GradientDriver.
The generalised Fock, the DF families, the probe terms and the one-electron
density follow emrsf_grad_ref.EMRSFGrad term by term (its gates are the
finite differences of the energies); the closed-shell kernel set of the
grid module supplies v[f rho(M)], v[k rho(D_s)^2] and v[f rho(D_s)].
"""
import sys
import numpy as np

import graci.io.output as output
import graci.interfaces.mrsf.mrsf_grad as mrsf_grad
from graci.interfaces.mrsf.mrsf_gradient import GradientDriver, _sym, SQ2

# coupling coefficients of mrsf_extended.f90 (set_coefficients), per multiplicity
COEF = {1: dict(cJ1=-2.0, cK1=1.0, cF1=-1.0, cW2=-1.0, cW3=1.0, cJ4=2.0, cK4=-1.0, cF4=-1.0,
                cA5=-1.0, cB5=1.0, cK6=1.0, cJ6=-2.0, cG=SQ2),
        3: dict(cJ1=0.0, cK1=1.0, cF1=-1.0, cW2=-1.0, cW3=-1.0, cJ4=0.0, cK4=1.0, cF4=1.0,
                cA5=-1.0, cB5=-1.0, cK6=1.0, cJ6=0.0, cG=0.0)}


class EMRSFGradientDriver(GradientDriver):

    def __init__(self, ci, grad_obj):
        from graci.methods.mrsftddft import closed_shell_fock
        if str(getattr(ci, 'fprime', 'covariant')).lower() != 'covariant':
            sys.exit('$mrsfgradient: gradients of the extended method need the covariant '
                     "F' correction (fprime = covariant in the $mrsftddft section)")
        super().__init__(ci, grad_obj)
        C, nmo, nao = self.C, self.nmo, self.nao
        self.ncol = int(ci.ncol)
        self.singlet = (self.mult == 1)
        self.use_kernel = bool(self.is_dft and self.singlet)
        self.ccp = float(ci.ccp_eff)
        self.cH = float(ci.chf_eff)
        self.coef = COEF[self.mult]
        self.O1, self.O2 = int(self.Oidx[0]), int(self.Oidx[1])
        self.hO1, self.hO2 = self.Hloc[self.O1], self.Hloc[self.O2]
        self.Gidx = [int(c) for c in self.Cidx] + [self.O1]
        self.occG = np.zeros(nmo); self.occG[self.Gidx] = 2.0
        self.dq_G = mrsf_grad.dq(self.occG, self.naux)
        self.F_G = closed_shell_fock(self.scf, self.mf, C, self.Gidx)[1]
        self.BooQ = mrsf_grad.booq(self.dims)
        # canonical pair (O2, O1): the extended energy depends on the SOMO rotation
        self.N = list(self.N) + [(self.O2, self.O1)]
        fsum = self.fa + self.fb
        self.hN = np.array([fsum[k, k] - fsum[m, m] for (k, m) in self.N])
        self.somo_gap = 0.5*(fsum[self.O2, self.O2] - fsum[self.O1, self.O1])
        if self.somo_gap < 1.0e-2:
            output.print_message('WARNING: SOMO gap %.4f Ha: the extended-method gradient depends on '
                                 'the O1-O2 multiplier zeta_12 ~ 1/gap' % self.somo_gap)
        if self.is_dft:
            self.xcgrid.set_g_kernel(C, self.occG, self.use_kernel)
        self._cur = None

    # ------------------------------------------------------------------
    # amplitudes
    # ------------------------------------------------------------------
    def expanded_amplitudes(self, istate):
        """X~ (nvirb, nocca) of the MRSF part; the CV amplitudes Y (nV, nC)
        are kept in self._Y_last"""
        irr, st = self.ci.state_sym(istate)
        x = self.ci.amps['adiabatic'][irr][:, st]
        X = np.array(x.reshape(self.nvirb, self.ncol, order='F'))
        Xm = np.array(X[:, :self.nocca])
        self._Y_last = np.array(X[2:, self.nocca:])
        if self.Cf:
            fcols = [self.Hloc[k] for k in self.Cf]
            if max(np.abs(Xm[:, fcols]).max(), np.abs(self._Y_last[:, fcols]).max()) > 1e-12:
                sys.exit('MRSF gradients: amplitudes on frozen holes (inconsistent frozen core)')
        pO1, pO2, hO1, hO2 = 0, 1, self.nC, self.nC + 1
        xoo = Xm[pO1, hO1]
        Xm[pO1, hO1] = xoo/SQ2
        if self.mult == 1:
            Xm[pO2, hO2] = -xoo/SQ2
        else:
            Xm[pO2, hO2] = xoo/SQ2
            Xm[pO1, hO2] = 0.0
            Xm[pO2, hO1] = 0.0
        return Xm

    def compressed_slices(self, sd):
        """compressed MRSF amplitudes of the coupling classes from the
        expanded X~ (the OO slot folded back)"""
        Xt = sd['Xt']
        nC, hO1, hO2 = self.nC, self.nC, self.nC + 1
        return {'CO1': Xt[0, :nC].copy(), 'CO2': Xt[1, :nC].copy(),
                'O1V': Xt[2:, hO1].copy(), 'O2V': Xt[2:, hO2].copy(),
                'CV': Xt[2:, :nC].copy(), 'OO': Xt[0, hO1]*SQ2, 'G': Xt[0, hO2]}

    # ------------------------------------------------------------------
    # per-state stages
    # ------------------------------------------------------------------
    def T_densities(self, Xt):
        Ta, Tb = super().T_densities(Xt)
        y2 = self._cur['y2']
        Ta[self.O2, self.O2] -= y2
        Tb[self.O1, self.O1] += y2
        return Ta, Tb

    def exchange_pass(self, sd):
        """MRSF exchange pass with the A_G patches of the T densities, the
        CV exchange pass (mean-field partner D_G; the F' correction pairs are
        separate) and the DF products of the amplitudes used by the coupling terms"""
        super().exchange_pass(sd)
        st = sd['st']
        Y = sd['Y']
        nC, O1, O2, H, P = self.nC, self.O1, self.O2, self.H, self.P
        Yemb = np.zeros((self.nvirb, self.nocca)); Yemb[2:, :nC] = Y
        sd['Yemb'] = Yemb
        y2 = float(np.sum(Y*Y)); sd['y2'] = y2
        Bso = self.Bso
        st['jT'] = st['jT'] + y2*Bso[:, O1, 0]        # Coulomb vector of the T^b patch (the T^a patch is in Ahh)
        K11 = Bso[:, :, 0].T @ Bso[:, :, 0]                       # K[e_O1 e_O1^T]
        st['KbT_HP'] = st['KbT_HP'] + y2*K11[np.ix_(H, P)]
        st['KbT_HH'] = st['KbT_HH'] + y2*K11[np.ix_(H, H)]
        sd['stcv'] = mrsf_grad.state(self.cH, Yemb, self.dims, self.dq_G)
        xs = self.compressed_slices(sd)
        sd['xs'] = xs
        # DF products: sum_q B^Q_tq U_q for the amplitude vectors (naux, nmo, .)
        sd['bvY'] = mrsf_grad.bvec(Y, 'virt', self.nmo, self.naux)             # Y columns on the virtuals
        sd['bvXcv'] = mrsf_grad.bvec(xs['CV'], 'virt', self.nmo, self.naux)    # x_CV columns
        sd['bvxO2V'] = mrsf_grad.bvec(xs['O2V'], 'virt', self.nmo, self.naux)[:, :, 0]
        sd['bvxCO1'] = mrsf_grad.bvec(xs['CO1'], 'core', self.nmo, self.naux)[:, :, 0]
        t = 2.0*self.ccp; c = self.coef
        sd['w1'] = Y @ xs['CO1']                      # (nV) Y x_CO1
        sd['w2'] = Y @ xs['CO2']
        sd['u3'] = Y.T @ xs['O1V']                    # (nC)
        sd['u4'] = Y.T @ xs['O2V']
        sd['a1'] = t*c['cF1']; sd['a2'] = t*c['cF4']; sd['a3'] = t*c['cG']*xs['G']
        sd['bvw1'] = mrsf_grad.bvec(sd['w1'], 'virt', self.nmo, self.naux)[:, :, 0]
        sd['M'] = self.M_matrix(sd)
        sd['jy'] = mrsf_grad.jq_of(None, Yemb, self.dims)[:, 0]                 # Coulomb vector of Yh^T
        sd['Tcc'] = -Y.T @ Y
        sd['jTcc'] = np.einsum('ijQ,ij->Q', self.BooQ[:nC, :nC, :], sd['Tcc'])
        sd['jTpp'] = np.array(sd['stcv']['jT'])

    def M_matrix(self, sd):
        """M = T^cv + sym(N_F) (MO): the density contracted with F^G"""
        nmo, Y, Cidx, Vidx, O1, O2 = self.nmo, sd['Y'], self.Cidx, self.Vidx, self.O1, self.O2
        M = np.zeros((nmo, nmo))
        M[np.ix_(Vidx, Vidx)] = Y @ Y.T
        M[np.ix_(Cidx, Cidx)] = -Y.T @ Y
        N = np.zeros((nmo, nmo))
        N[O2, Vidx] += sd['a1']*sd['w1']
        N[Cidx, O1] += sd['a2']*sd['u4']
        N[np.ix_(Cidx, Vidx)] += sd['a3']*Y.T
        return M + _sym(N)

    def M_factors(self, sd):
        """factor pairs (LM, RM) of M = LM RM^T + RM LM^T (nao, 3 nC + 2)"""
        C, Y, nC = self.C, sd['Y'], self.nC
        CV, CC = C[:, self.Vidx], C[:, self.Cidx]
        LM = np.zeros((self.nao, 3*nC + 2)); RM = np.zeros((self.nao, 3*nC + 2))
        CVY = CV @ Y
        LM[:, :nC] = CVY;                      RM[:, :nC] = 0.5*CVY
        LM[:, nC:2*nC] = -0.5*(CC @ (Y.T @ Y)); RM[:, nC:2*nC] = CC
        LM[:, 2*nC] = 0.5*sd['a1']*C[:, self.O2]; RM[:, 2*nC] = CV @ sd['w1']
        LM[:, 2*nC + 1] = 0.5*sd['a2']*(CC @ sd['u4']); RM[:, 2*nC + 1] = C[:, self.O1]
        LM[:, 2*nC + 2:] = 0.5*sd['a3']*CC;    RM[:, 2*nC + 2:] = CVY
        return LM, RM

    def prepare_state(self, ist):
        sd = {'ist': ist, 'Xt': self.expanded_amplitudes(ist)}
        sd['Y'] = self._Y_last
        self._cur = sd
        self.exchange_pass(sd)
        sd['terms'] = self.spc_aux_vectors(sd['Xt'])
        sd['Fo'] = self.gfock_omega(sd['Xt'], sd['st'], sd['terms'])
        sd['zeta'] = self.zeta(sd['Fo'])
        sd['Fzeta'] = self.gfock_zeta(sd['zeta'])
        return sd

    def finish_state(self, sd, z):
        self._cur = sd
        return super().finish_state(sd, z)

    # ------------------------------------------------------------------
    # generalised Fock of omega_E
    # ------------------------------------------------------------------
    def full_matrix(self, MHH, MHP):
        """nmo x nmo symmetric matrix from its (H, H) and (H, P) blocks"""
        M = np.zeros((self.nmo, self.nmo))
        H, P = self.H, self.P
        M[np.ix_(H, H)] = MHH
        M[np.ix_(H, P)] = MHP
        M[np.ix_(P, H)] = MHP.T
        return M

    def gfock_omega(self, Xt, st, terms):
        sd = self._cur
        fa, fb, C, H, P, V = self.fa, self.fb, self.C, self.H, self.P, self.Vidx
        nocca, nao = self.nocca, self.nao
        Ta, Tb = self.T_densities(Xt)
        Ta_loc = Ta[np.ix_(H, H)]
        F = 2.0*(fa @ Ta + fb @ Tb)
        xcfac = None
        if self.is_dft:
            # T^a = C_H Ta C_H^T; T^b = (C_P X~)(C_P X~)^T + |y|^2 C_O1 C_O1^T
            CH, CP = C[:, H], C[:, P]
            k = nocca + 1
            Lf = np.zeros((nao, k, 2, 1)); Rf = np.zeros((nao, k, 2, 1))
            Lf[:, :nocca, 0, 0] = 0.5*(CH @ Ta_loc); Rf[:, :nocca, 0, 0] = CH
            Lf[:, :nocca, 1, 0] = 0.5*(CP @ Xt);     Rf[:, :nocca, 1, 0] = CP @ Xt
            Lf[:, nocca, 1, 0] = 0.5*sd['y2']*C[:, self.O1]; Rf[:, nocca, 1, 0] = C[:, self.O1]
            xcfac = (Lf, Rf, [[0], [0]])
        Ga, Gb = self.G_full(Ta_loc, None, None, jq_add=st['jT'], KbHP=st['KbT_HP'], KbHH=st['KbT_HH'],
                             xcfac=xcfac)
        F[:, H] += 2.0*Ga[:, H]
        F[:, self.Cidx] += 2.0*Gb[:, self.Cidx]
        F[:, H] += self.cols_H(st['LaH'], st['LaP'])
        LbFull = np.zeros((self.nmo, self.nvirb))
        LbFull[H] = st['LbH']; LbFull[V] = st['LbP'][2:]
        F[:, P] += LbFull
        F += self.gfock_spc(Xt, terms)
        # --- extended-method pieces ---
        stcv = sd['stcv']
        F[:, H] += self.cols_H(stcv['LaH'], stcv['LaP'])             # CV exchange channel
        LbFull = np.zeros((self.nmo, self.nvirb))
        LbFull[H] = stcv['LbH']; LbFull[V] = stcv['LbP'][2:]
        F[:, P] += LbFull
        M = sd['M']
        F += 2.0*self.F_G @ M                                         # index term of Tr(F^G M)
        G1 = self.G_closed(sd)                                        # D_G response of F^G
        F += 2.0*G1*self.occG[None, :]
        if self.use_kernel:
            V1m, K1 = sd['V1m'], sd['K1']
            S = np.zeros((self.nmo, self.nmo))
            S[np.ix_(V, self.Cidx)] = sd['Y']
            S = S + S.T
            F += 2.0*V1m @ S + K1*self.occG[None, :]
        F += self.pairs_gfock(sd)
        return F

    def G_closed(self, sd):
        """G_G[M] = J[M] - x/2 K[M] + f_xc^G M on the (all, H) columns, as an
        nmo x nmo matrix with the (H, all) rows filled symmetrically"""
        H, P, Cidx, Vidx = self.H, self.P, self.Cidx, self.Vidx
        M, stcv = sd['M'], sd['stcv']
        # partition of the symmetric M into the (H, H) block, the (V, H) block
        # (symmetrised by gfock) and the (V, V) block (Coulomb vector + K[Y Y^T])
        Mhh = np.asfortranarray(M[np.ix_(H, H)])
        Mph = np.array(M[np.ix_(P, H)]); Mph[:2, :] = 0.0
        Mph = np.asfortranarray(Mph)
        jq_pp = sd['jTpp']
        x = 0.5*self.x_ref
        Ghpa, Ghpb, Ghha, Ghhb = mrsf_grad.gfock(x, self.dims, Mhh, Mph, None, jq_pp, True, True)
        GHP = Ghpa[:, :, 0] - x*stcv['KbT_HP']
        GHH = Ghha[:, :, 0] - x*stcv['KbT_HH']
        if self.is_dft:
            LM, RM = self.M_factors(sd)
            LS = self.C[:, Vidx] @ sd['Y']
            MM, MK, V1 = self.xcgrid.potential_g(LM, RM, LS)
            VM = MM @ self.C                                          # v[f rho(M)] (H, all)
            GHH = GHH + VM[:, H]; GHP = GHP + VM[:, P]
            if self.use_kernel:
                sd['V1m'] = self.C.T @ V1 @ self.C
                K1H = MK @ self.C                                     # v[k rho(D_s)^2] (H, all)
                K1 = np.zeros((self.nmo, self.nmo))
                K1[H, :] = K1H; K1[:, H] = K1H.T
                sd['K1'] = K1
        return self.full_matrix(GHH, GHP)

    # ------------------------------------------------------------------
    # two-electron pairs: generalised Fock
    # ------------------------------------------------------------------
    def Jmat(self, jq):
        JHH, JHP, JPP = mrsf_grad.jblocks(jq, self.dims)
        return self.full_from_blocks(JHH[:, :, 0], JHP[:, :, 0], JPP[:, :, 0])

    def unit_pair(self, p, q, val=1.0):
        A = np.zeros((self.nmo, self.nmo)); A[p, q] = val
        return A

    def pairs_gfock(self, sd):
        """generalised Fock of the two-electron pairs (CV Coulomb, A_G, the
        covariant F' correction and the two-electron couplings):
        J-type (A, B, t): t [J[B](A + A^T) + J[A](B + B^T)]
        K-type (A, B, t): t [K[B] A^T + K[B]^T A + K[A] B^T + K[A]^T B]"""
        nmo, Cidx, Vidx, O1, O2 = self.nmo, self.Cidx, self.Vidx, self.O1, self.O2
        Bso, Y, xs, c = self.Bso, sd['Y'], sd['xs'], self.coef
        t = 2.0*self.ccp
        F = np.zeros((nmo, nmo))
        YhT = np.zeros((nmo, nmo)); YhT[np.ix_(Cidx, Vidx)] = Y.T
        Yh = YhT.T
        J_y = self.Jmat(sd['jy'])
        e11 = self.unit_pair(O1, O1); e22 = self.unit_pair(O2, O2); e21 = self.unit_pair(O2, O1)
        J_e11 = self.Jmat(Bso[:, O1, 0]); J_e22 = self.Jmat(Bso[:, O2, 1]); J_e21 = self.Jmat(Bso[:, O1, 1])
        def jpair(A, B, tt, JA, JB):
            return tt*(JB @ (A + A.T) + JA @ (B + B.T))
        # CO1-J: (Yh^T, row spec [O2, p] = x_CO1)
        if c['cJ1']:
            B = np.zeros((nmo, nmo)); B[O2, Cidx] = xs['CO1']
            F += jpair(YhT, B, t*c['cJ1'], J_y, self.Jmat(Bso[:, Cidx, 1] @ xs['CO1']))
        # CO2-W: (e_O2 e_O1^T, row spec [O2, b] = w2)
        if c['cW2']:
            B = np.zeros((nmo, nmo)); B[O2, Vidx] = sd['w2']
            F += jpair(e21, B, t*c['cW2'], J_e21, self.Jmat(Bso[:, Vidx, 1] @ sd['w2']))
        # O1V-W: (col spec [j, O1] = u3, e_O2 e_O1^T)
        if c['cW3']:
            A = np.zeros((nmo, nmo)); A[Cidx, O1] = sd['u3']
            F += jpair(A, e21, t*c['cW3'], self.Jmat(Bso[:, Cidx, 0] @ sd['u3']), J_e21)
        # O2V-J: (Yh^T, col spec [q, O1] = x_O2V)
        if c['cJ4']:
            B = np.zeros((nmo, nmo)); B[Vidx, O1] = xs['O2V']
            F += jpair(YhT, B, t*c['cJ4'], J_y, self.Jmat(Bso[:, Vidx, 0] @ xs['O2V']))
        # OO-J: (Yh^T, x_OO e_O2 e_O1^T)
        if c['cJ6'] and xs['OO']:
            F += jpair(YhT, xs['OO']*e21, t*c['cJ6'], J_y, xs['OO']*J_e21)
        # A_G two-electron part: (e_O2 e_O2^T, e_O1 e_O1^T) x (-c_H |y|^2)
        F += jpair(e22, e11, -self.cH*sd['y2'], J_e22, J_e11)
        # covariant F' correction: (T^cv_VV, e_x e_x^T) and (T^cv_CC, e_x e_x^T), +/-(1 - c_H)
        Tvv = np.zeros((nmo, nmo)); Tvv[np.ix_(Vidx, Vidx)] = Y @ Y.T
        Tcc = np.zeros((nmo, nmo)); Tcc[np.ix_(Cidx, Cidx)] = sd['Tcc']
        J_Tvv = self.Jmat(sd['jTpp']); J_Tcc = self.Jmat(sd['jTcc'])
        tc = 1.0 - self.cH
        F += jpair(Tvv, e22, tc, J_Tvv, J_e22) - jpair(Tvv, e11, tc, J_Tvv, J_e11)
        F += jpair(Tcc, e22, tc, J_Tcc, J_e22) - jpair(Tcc, e11, tc, J_Tcc, J_e11)
        # CV Coulomb (singlet): (Yh^T, Yh^T) x 2
        if self.singlet:
            F += jpair(YhT, YhT, 2.0, J_y, J_y)
        # ---- K-type pairs ----
        def kpair(A, B, tt, KA, KB):
            return tt*(KB @ A.T + KB.T @ A + KA @ B.T + KA.T @ B)
        K_Yh = None
        def K_Yh_get():
            # K[Yh]_{pr} = sum_Q sum_{q in V} B^Q_pq (sum_{s in C} y_qs B^Q_sr):
            # bdot over the core class of W(s, Q, w) = bvY(Q, w, s) gives K[Yh]^T
            W = np.ascontiguousarray(sd['bvY'].transpose(2, 0, 1))
            return mrsf_grad.bdot(W, 'core', nmo, self.naux).T
        # CO1-K: (Yh^T, col spec [p, O2] = x_CO1)
        if c['cK1']:
            K_Yh = K_Yh_get() if K_Yh is None else K_Yh
            B = np.zeros((nmo, nmo)); B[Cidx, O2] = xs['CO1']
            KB = sd['bvxCO1'].T @ Bso[:, :, 1]                      # K[x e_O2^T]
            F += kpair(YhT, B, t*c['cK1'], K_Yh.T, KB)
        # O2V-K: (col spec [q, O1] = x_O2V, Yh)
        if c['cK4']:
            K_Yh = K_Yh_get() if K_Yh is None else K_Yh
            A = np.zeros((nmo, nmo)); A[Vidx, O1] = xs['O2V']
            KA = sd['bvxO2V'].T @ Bso[:, :, 0]                      # K[x e_O1^T]
            F += kpair(A, Yh, t*c['cK4'], KA, K_Yh)
        # CV-A: (M12 = x_CV Y^T on V x V, e_O1 e_O2^T)
        e12 = self.unit_pair(O1, O2)
        K_e12 = Bso[:, :, 0].T @ Bso[:, :, 1]
        if c['cA5']:
            M12 = np.zeros((nmo, nmo)); M12[np.ix_(Vidx, Vidx)] = xs['CV'] @ Y.T
            KA = np.einsum('Qpj,Qrj->pr', sd['bvXcv'], sd['bvY'], optimize=True)
            F += kpair(M12, e12, t*c['cA5'], KA, K_e12)
        # CV-B: (M12c = Y^T x_CV on C x C, e_O1 e_O2^T)
        if c['cB5']:
            M12c_loc = Y.T @ xs['CV']
            M12c = np.zeros((nmo, nmo)); M12c[np.ix_(Cidx, Cidx)] = M12c_loc
            bvM = mrsf_grad.bvec(M12c_loc, 'core', nmo, self.naux)
            KA = mrsf_grad.bdot(np.ascontiguousarray(bvM.transpose(2, 0, 1)), 'core', nmo, self.naux).T
            F += kpair(M12c, e12, t*c['cB5'], KA, K_e12)
        # OO-K: (Yh^T, x_OO e_O1 e_O2^T)
        if c['cK6'] and xs['OO']:
            K_Yh = K_Yh_get() if K_Yh is None else K_Yh
            F += kpair(YhT, xs['OO']*e12, t*c['cK6'], K_Yh.T, xs['OO']*K_e12)
        return F

    # ------------------------------------------------------------------
    # families
    # ------------------------------------------------------------------
    def extra_families(self, sd):
        """CV exchange families, closed-shell mean-field pairs of M, the
        two-electron pair pieces and the A_G patches, added to the state's
        Ghh/Fhp/gpp in place; returns the CV particle-particle family entry"""
        st, stcv = sd['st'], sd['stcv']
        Y, xs, c = sd['Y'], sd['xs'], self.coef
        nC, nV, nocca = self.nC, self.nV, self.nocca
        H, P, Cidx, Vidx, O1, O2 = self.H, self.P, self.Cidx, self.Vidx, self.O1, self.O2
        hO1, hO2 = self.hO1, self.hO2
        Hloc, Ploc, Bso = self.Hloc, self.Ploc, self.Bso
        Ghh, Fhp = st['Ghh'], st['Fhp']
        hC = [Hloc[p] for p in Cidx]; pV = [Ploc[p] for p in Vidx]
        t = 2.0*self.ccp
        x = self.x_ref
        # (1) CV exchange channel: hole-hole family and metric part of the CV pass
        Ghh += stcv['Ghh']
        st['gpp'] += stcv['gpp']
        Yf_cv = np.array(stcv['Yf'])                                  # -c_H/4 Y^Q + d_G/4 Y, partner Yemb
        Yf_cv[:, hO1, :] = 0.0; Yf_cv[:, hO2, :] = 0.0                # the partner's SOMO-hole columns are zero:
                                                                      # free for the rank-1 pieces below
        # (2) A_G patches of T^b (|y|^2 at (O1, O1)): J pair with the reference
        # density and the beta exchange pair with the core
        y2 = sd['y2']
        Ghh[hO1, hO1, :] += 0.5*y2*self.dq_ref
        v = 0.5*x*y2*Bso[:, Cidx, 0].T                                # (nC, naux)
        Ghh[hC, hO1, :] -= v
        Ghh[hO1, hC, :] -= v
        # (3) closed-shell mean-field pairs of M: J pair (M, D_G), K pair (M, D_G) with x/2
        M = sd['M']
        dqG = self.dq_G
        Mhh = M[np.ix_(H, H)]
        Ghh += 0.5*np.einsum('Q,ij->ijQ', dqG, Mhh)
        MHV = M[np.ix_(H, Vidx)]
        Fhp[:, pV, :] += np.einsum('Q,ia->iaQ', dqG, MHV)
        pM = sd['jTpp'] + sd['jTcc'] + sd['a1']*(Bso[:, Vidx, 1] @ sd['w1']) \
            + sd['a2']*(Bso[:, Cidx, 0] @ sd['u4']) + sd['a3']*sd['jy']
        occG_H = self.occG[H]
        for i in range(nocca):
            Ghh[i, i, :] += 0.5*pM*occG_H[i]
        # K pair (M, D_G): Gamma -= x/2 (M B^Q 1_G + 1_G B^Q M)
        G_loc = [Hloc[p] for p in self.Gidx]
        isG = np.zeros(nocca); isG[G_loc] = 1.0
        MB = np.einsum('ij,jkQ->ikQ', Mhh, self.BooQ)                 # (M_HH B^Q_HH)
        bvY = sd['bvY']                                               # (naux, nmo, nC): sum_c B^Q_qc Y_cj
        # cross terms from the H x V blocks of M: sum_b M[p,b] B^Q_{bq}
        MBx = np.zeros_like(MB)
        MBx[hC, :, :] = 0.5*sd['a3']*np.transpose(bvY[:, H, :], (2, 1, 0))
        MBx[hO2, :, :] = 0.5*sd['a1']*sd['bvw1'][:, H].T
        MB += MBx
        Ghh -= 0.5*x*(MB*isG[None, :, None] + np.transpose(MB, (1, 0, 2))*isG[:, None, None])
        # particle rows: (M B^Q 1_G)_{aq}, a in V, q in G
        MaH = np.einsum('aj,jqQ->aqQ', M[np.ix_(Vidx, H)], self.BooQ)
        MaV = np.einsum('aj,Qqj->aqQ', Y, bvY[:, H, :])
        Fq = MaH + MaV                                                # (nV, nocca, naux)
        for q in G_loc:
            Fhp[q, pV, :] -= x*Fq[:, q, :]
        # (4) two-electron pairs (1/2 t A b^Q + 1/2 t a^Q B for J-type;
        #     1/2 t (A B^Q B^T + A^T B^Q B) for K-type)
        jy = sd['jy']
        def add_hp_YT(vecQ, fac):                                     # Fhp[j, b, Q] += fac y_bj vecQ
            Fhp[np.ix_(hC, pV)] += fac*np.einsum('bj,Q->jbQ', Y, vecQ)
        # CO1-J
        if c['cJ1']:
            tt = t*c['cJ1']
            add_hp_YT(Bso[:, Cidx, 1] @ xs['CO1'], 0.5*tt)
            Ghh[hO2, hC, :] += 0.5*tt*np.outer(xs['CO1'], jy)
        # CO2-W
        if c['cW2']:
            tt = t*c['cW2']
            Ghh[hO2, hO1, :] += 0.5*tt*(Bso[:, Vidx, 1] @ sd['w2'])
            Fhp[hO2, pV, :] += 0.5*tt*np.outer(sd['w2'], Bso[:, O1, 1])
        # O1V-W
        if c['cW3']:
            tt = t*c['cW3']
            Ghh[hC, hO1, :] += 0.5*tt*np.outer(sd['u3'], Bso[:, O1, 1])
            Ghh[hO2, hO1, :] += 0.5*tt*(Bso[:, Cidx, 0] @ sd['u3'])
        # O2V-J
        if c['cJ4']:
            tt = t*c['cJ4']
            add_hp_YT(Bso[:, Vidx, 0] @ xs['O2V'], 0.5*tt)
            Fhp[hO1, pV, :] += 0.5*tt*np.outer(xs['O2V'], jy)
        # OO-J
        if c['cJ6'] and xs['OO']:
            tt = t*c['cJ6']*xs['OO']
            add_hp_YT(Bso[:, O1, 1], 0.5*tt)
            Ghh[hO2, hO1, :] += 0.5*tt*jy
        # A_G two-electron: (e22, e11) x (-c_H y2)
        tt = -self.cH*y2
        Ghh[hO2, hO2, :] += 0.5*tt*Bso[:, O1, 0]
        Ghh[hO1, hO1, :] += 0.5*tt*Bso[:, O2, 1]
        # covariant F' correction, +/-(1 - c_H): PP via Yf_cv (partner Y), CC via Ghh
        tc = 1.0 - self.cH
        bdiff = Bso[:, O2, 1] - Bso[:, O1, 0]                         # b_2^Q - b_1^Q
        Yf_cv += 0.25*tc*np.einsum('Q,ai->aiQ', bdiff, sd['Yemb'])
        st['gpp'] += 0.5*tc*np.outer(sd['jTpp'], bdiff)
        Ghh[np.ix_(hC, hC)] += 0.5*tc*np.einsum('Q,ij->ijQ', bdiff, sd['Tcc'])
        Ghh[hO2, hO2, :] += 0.5*tc*(sd['jTpp'] + sd['jTcc'])
        Ghh[hO1, hO1, :] -= 0.5*tc*(sd['jTpp'] + sd['jTcc'])
        # CV Coulomb (singlet): (Yh^T, Yh^T) x 2
        if self.singlet:
            add_hp_YT(jy, 2.0)
        # K-type: CO1-K
        if c['cK1']:
            tt = t*c['cK1']
            vq = np.einsum('bp,Qb->pQ', Y, Bso[:, Vidx, 1])            # sum_b y_bp B^Q_{b O2}
            Ghh[np.ix_(hC, hC)] += 0.5*tt*np.einsum('pQ,q->pqQ', vq, xs['CO1'])
            Fhp[hO2, pV, :] += 0.5*tt*np.einsum('pr,Qr->pQ', Y, sd['bvxCO1'][:, Cidx])
        # O2V-K: (col spec [q, O1] = x_O2V, Yh)
        if c['cK4']:
            tt = t*c['cK4']
            # PP piece 1/2 t x (Y B^Q_{C,O1})^T: spare column hO1 of the CV partner
            Yf_cv[2:, hO1, :] += 0.25*tt*np.einsum('aj,Qj->aQ', Y, Bso[:, Cidx, 0])
            sd['Xp_hO1'] = xs['O2V']
            Ghh[hO1, hC, :] += 0.5*tt*np.einsum('Qr,rq->qQ', sd['bvxO2V'][:, Vidx], Y)
        # CV-A: (M12, e_O1 e_O2^T)
        if c['cA5']:
            tt = t*c['cA5']
            M12 = xs['CV'] @ Y.T
            Fhp[hO1, pV, :] += 0.5*tt*np.einsum('pr,Qr->pQ', M12, Bso[:, Vidx, 1])
            Fhp[hO2, pV, :] += 0.5*tt*np.einsum('rp,Qr->pQ', M12, Bso[:, Vidx, 0])
        # CV-B: (M12c, e_O1 e_O2^T)
        if c['cB5']:
            tt = t*c['cB5']
            M12c = Y.T @ xs['CV']
            Ghh[hC, hO1, :] += 0.5*tt*np.einsum('pr,Qr->pQ', M12c, Bso[:, Cidx, 1])
            Ghh[hC, hO2, :] += 0.5*tt*np.einsum('rp,Qr->pQ', M12c, Bso[:, Cidx, 0])
        # OO-K: (Yh^T, x_OO e_O1 e_O2^T)
        if c['cK6'] and xs['OO']:
            tt = t*c['cK6']*xs['OO']
            Ghh[hC, hO1, :] += 0.5*tt*np.einsum('bp,Qb->pQ', Y, Bso[:, Vidx, 1])
            Fhp[hO2, pV, :] += 0.5*tt*np.einsum('pr,Qr->pQ', Y, Bso[:, Cidx, 0])
        # metric contribution of the spare column: 2 sum_a (B^Q_PP x)_a Yf^Q'_a
        Xp = np.array(sd['Yemb'])
        if 'Xp_hO1' in sd:
            Xp[2:, hO1] = sd['Xp_hO1']
            st['gpp'] += 2.0*np.einsum('Qa,aR->QR', sd['bvxO2V'][:, P], Yf_cv[:, hO1, :])
        sd['Yf_cv'] = Yf_cv
        return [{'Xt': Xp, 'Ghh': np.zeros((nocca, nocca, self.naux), order='F'),
                 'Fhp': np.zeros((nocca, self.nvirb, self.naux), order='F'),
                 'Yf': np.asfortranarray(Yf_cv), 'g': np.zeros((self.naux, self.naux))}]

    # ------------------------------------------------------------------
    # probe factors and the one-electron density
    # ------------------------------------------------------------------
    def P_factors(self, Xt, Za, Zb, pzeta=None):
        """relaxed densities with the A_G patch of T^b (|y|^2 C_O1 C_O1^T) as an
        extra factor column (k = 2 nocca + 1)"""
        sd = self._cur
        CH, CP, C = self.C[:, self.H], self.C[:, self.P], self.C
        nao, nocca = self.nao, self.nocca
        Ta, Tb = self.T_densities(Xt)
        Ta_loc = Ta[np.ix_(self.H, self.H)]
        k = 2*nocca + 1
        L = np.zeros((nao, k, 2)); R = np.zeros((nao, k, 2))
        L[:, :nocca, 0] = 0.5*(CH @ Ta_loc);        R[:, :nocca, 0] = CH
        L[:, nocca:2*nocca, 0] = CP @ Za;           R[:, nocca:2*nocca, 0] = CH
        CPX = CP @ Xt
        L[:, :nocca, 1] = 0.5*CPX;                  R[:, :nocca, 1] = CPX
        L[:, nocca:2*nocca, 1] = CP @ Zb;           R[:, nocca:2*nocca, 1] = CH
        L[:, 2*nocca, 1] = 0.5*sd['y2']*C[:, self.O1]; R[:, 2*nocca, 1] = C[:, self.O1]
        if pzeta is not None and self.N:
            add = 0.5*(CH @ pzeta[np.ix_(self.H, self.H)])
            L[:, :nocca, 0] += add
            L[:, nocca:2*nocca, 1] += add
        return L, R

    def probe_g_factors(self, sd):
        if not self.is_dft:
            return None
        LM, RM = self.M_factors(sd)
        LS = self.C[:, self.Vidx] @ sd['Y']
        return LM, RM, LS

    def extra_density(self, sd):
        return sd['M']
