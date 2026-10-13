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
        with self.tm('setup: F_G'):
            self.F_G = self.session.fock_g(ci)                        # once per SCF (shared with the energy sections)
        with self.tm('setup: BooQ'):
            self.BooQ = mrsf_grad.booq(self.dims)
        with self.tm('setup: pair constants'):
            Bso = self.Bso
            self.BsoC = np.ascontiguousarray(Bso[:, self.Cidx, :])       # B^Q_{c O_x} (naux, nC, 2)
            self.BsoV = np.ascontiguousarray(Bso[:, self.Vidx, :])       # B^Q_{v O_x} (naux, nV, 2)
            self.K11 = Bso[:, :, 0].T @ Bso[:, :, 0]                      # K[e_O1 e_O1^T]
            self.K_e12 = Bso[:, :, 0].T @ Bso[:, :, 1]                    # K[e_O1 e_O2^T]
            JHH, JHP, JPP = mrsf_grad.jblocks(np.column_stack([Bso[:, self.O1, 0], Bso[:, self.O2, 1],
                                                               Bso[:, self.O1, 1]]), self.dims)
            self.J_e11, self.J_e22, self.J_e21 = [self.full_from_blocks(JHH[:, :, k], JHP[:, :, k], JPP[:, :, k])
                                                  for k in range(3)]
        # canonical pair (O2, O1): the extended energy depends on the SOMO rotation
        self.N = list(self.N) + [(self.O2, self.O1)]
        fsum = self.fa + self.fb
        self.hN = np.array([fsum[k, k] - fsum[m, m] for (k, m) in self.N])
        self.somo_gap = 0.5*(fsum[self.O2, self.O2] - fsum[self.O1, self.O1])
        if self.somo_gap < 1.0e-2:
            output.print_message('WARNING: SOMO gap %.4f Ha: the extended-method gradient depends on '
                                 'the O1-O2 multiplier zeta_12 ~ 1/gap' % self.somo_gap)
        if self.is_dft:
            with self.tm('setup: G kernel set'):
                self.xcgrid.ensure_g_kernel(self.Gidx, self.use_kernel, g_default=self.session.gidx)
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
        with self.tm('rhs: state (MRSF)'):
            super().exchange_pass(sd)
        st = sd['st']
        Y = sd['Y']
        nC, O1, O2, H, P = self.nC, self.O1, self.O2, self.H, self.P
        Yemb = np.zeros((self.nvirb, self.nocca)); Yemb[2:, :nC] = Y
        sd['Yemb'] = Yemb
        y2 = float(np.sum(Y*Y)); sd['y2'] = y2
        Bso = self.Bso
        st['jT'] = st['jT'] + y2*Bso[:, O1, 0]        # Coulomb vector of the T^b patch (the T^a patch is in Ahh)
        st['KbT_HP'] = st['KbT_HP'] + y2*self.K11[np.ix_(H, P)]
        st['KbT_HH'] = st['KbT_HH'] + y2*self.K11[np.ix_(H, H)]
        with self.tm('rhs: state (CV)'):
            sd['stcv'] = mrsf_grad.state(self.cH, Yemb, self.dims, self.dq_G)
        xs = self.compressed_slices(sd)
        sd['xs'] = xs
        t = 2.0*self.ccp; c = self.coef
        sd['w1'] = Y @ xs['CO1']                      # (nV) Y x_CO1
        sd['w2'] = Y @ xs['CO2']
        sd['u3'] = Y.T @ xs['O1V']                    # (nC)
        sd['u4'] = Y.T @ xs['O2V']
        sd['a1'] = t*c['cF1']; sd['a2'] = t*c['cF4']; sd['a3'] = t*c['cG']*xs['G']
        # DF products sum_q B^Q_tq U_q (naux, nmo, .) of the amplitude vectors on
        # the virtuals: one plane sweep for all of them
        with self.tm('rhs: bvec'):
            U = np.column_stack([Y, xs['CV'], xs['O2V'], sd['w1']])
            bv = mrsf_grad.bvec(U, 'virt', self.nmo, self.naux)
            sd['bvY'] = bv[:, :, :nC]                                # Y columns on the virtuals
            sd['bvXcv'] = bv[:, :, nC:2*nC]                          # x_CV columns
            sd['bvxO2V'] = bv[:, :, 2*nC]
            sd['bvw1'] = bv[:, :, 2*nC + 1]
            sd['bvxCO1'] = mrsf_grad.bvec(xs['CO1'], 'core', self.nmo, self.naux)[:, :, 0]
        with self.tm('rhs: M/jy/Tcc'):
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
        with self.tm('rhs: amplitudes'):
            sd = {'ist': ist, 'Xt': self.expanded_amplitudes(ist)}
        sd['Y'] = self._Y_last
        self._cur = sd
        with self.tm('rhs: exchange pass (all)'):
            self.exchange_pass(sd)
        with self.tm('rhs: spc vectors'):
            sd['terms'] = self.spc_aux_vectors(sd['Xt'])
        with self.tm('rhs: gfock_omega (all)'):
            sd['Fo'] = self.gfock_omega(sd['Xt'], sd['st'], sd['terms'])
        with self.tm('rhs: zeta'):
            sd['zeta'] = self.zeta(sd['Fo'])
        with self.tm('rhs: gfock_zeta'):
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
        with self.tm('rhs: G_full (gfock+kernel)'):
            Ga, Gb = self.G_full(Ta_loc, None, None, jq_add=st['jT'], KbHP=st['KbT_HP'], KbHH=st['KbT_HH'],
                                 xcfac=xcfac)
        F[:, H] += 2.0*Ga[:, H]
        F[:, self.Cidx] += 2.0*Gb[:, self.Cidx]
        F[:, H] += self.cols_H(st['LaH'], st['LaP'])
        LbFull = np.zeros((self.nmo, self.nvirb))
        LbFull[H] = st['LbH']; LbFull[V] = st['LbP'][2:]
        F[:, P] += LbFull
        # Coulomb matrices of the pairing terms and of the pair classes: one plane sweep
        with self.tm('rhs: jblocks'):
            xs, Y = sd['xs'], sd['Y']
            names = ['y', 'CO1', 'W2', 'W3', 'O2V', 'Tvv', 'Tcc']
            jq_pairs = [sd['jy'], self.BsoC[:, :, 1] @ xs['CO1'], self.BsoV[:, :, 1] @ sd['w2'],
                        self.BsoC[:, :, 0] @ sd['u3'], self.BsoV[:, :, 0] @ xs['O2V'], sd['jTpp'], sd['jTcc']]
            jq_spc = [v for (_, _, _, a, b) in terms for v in (a, b)]
            JHH, JHP, JPP = mrsf_grad.jblocks(np.array(jq_spc + jq_pairs).T, self.dims)
            n0 = len(jq_spc)
            sd['Jm'] = {nm: self.full_from_blocks(JHH[:, :, n0 + k], JHP[:, :, n0 + k], JPP[:, :, n0 + k])
                        for k, nm in enumerate(names)}
            sd['jv'] = dict(zip(names, jq_pairs))
        with self.tm('rhs: gfock_spc'):
            F += self.gfock_spc(Xt, terms, jbl=(JHH[:, :, :n0], JHP[:, :, :n0], JPP[:, :, :n0]))
        # --- extended-method pieces ---
        stcv = sd['stcv']
        F[:, H] += self.cols_H(stcv['LaH'], stcv['LaP'])             # CV exchange channel
        LbFull = np.zeros((self.nmo, self.nvirb))
        LbFull[H] = stcv['LbH']; LbFull[V] = stcv['LbP'][2:]
        F[:, P] += LbFull
        M = sd['M']
        F += 2.0*self.F_G @ M                                         # index term of Tr(F^G M)
        with self.tm('rhs: G_closed (all)'):
            G1 = self.G_closed(sd)                                        # D_G response of F^G
        F += 2.0*G1*self.occG[None, :]
        if self.use_kernel:
            # 2 V1m S with S = Y on (V, C) + h.c.: columns C <- 2 C^T V1 C_V Y, columns V <- 2 (C^T V1 C_C) Y^T
            F[:, self.Cidx] += 2.0*sd['VY']
            F[:, V] += 2.0*(sd['VC'] @ sd['Y'].T)
            F += sd['K1']*self.occG[None, :]
        with self.tm('rhs: pairs_gfock'):
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
        with self.tm('rhs: G_closed gfock'):
            Ghpa, Ghpb, Ghha, Ghhb = mrsf_grad.gfock(x, self.dims, Mhh, Mph, None, jq_pp, True, True)
        GHP = Ghpa[:, :, 0] - x*stcv['KbT_HP']
        GHH = Ghha[:, :, 0] - x*stcv['KbT_HH']
        if self.is_dft:
            # M and D_s as quadratic forms in [C_H | C_V Y | C_V w1] (shared with the probe)
            Bexp_g, gterms = self.g_quadratic(sd)
            with self.tm('rhs: potential_g'):
                MM, MK, VE, VH = self.xcgrid.potential_g_q(Bexp_g, gterms['M'], gterms['S'], self.use_kernel)
            VM = MM @ self.C                                          # v[f rho(M)] (H, all)
            GHH = GHH + VM[:, H]; GHP = GHP + VM[:, P]
            if self.use_kernel:
                nC = self.nC
                hC = [self.Hloc[p] for p in Cidx]
                sd['VY'] = self.C.T @ VE[:, :nC]                      # (v[f rho(D_s)] C_V Y) in the MO basis
                sd['VC'] = self.C.T @ VH[:, hC]                       # (v[f rho(D_s)] C_C) in the MO basis
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

    def pairs_gfock(self, sd):
        """generalised Fock of the two-electron pairs (CV Coulomb, A_G, the
        covariant F' correction and the two-electron couplings):
        J-type (A, B, t): t [J[B](A + A^T) + J[A](B + B^T)]
        K-type (A, B, t): t [K[B] A^T + K[B]^T A + K[A] B^T + K[A]^T B]
        with the sparse/low-rank operands handled as such: the J matrices come
        from the batched plane sweep of gfock_omega, the exchange matrices are
        never formed (only the rows/columns the operands select, each an
        O(naux nmo nC) contraction of the DF products of the amplitudes)"""
        nmo, naux, Cidx, Vidx, O1, O2 = self.nmo, self.naux, self.Cidx, self.Vidx, self.O1, self.O2
        Bso, BsoC, BsoV = self.Bso, self.BsoC, self.BsoV
        Y, xs, c, Jm = sd['Y'], sd['xs'], self.coef, sd['Jm']
        J_e11, J_e22, J_e21, K_e12 = self.J_e11, self.J_e22, self.J_e21, self.K_e12
        bvY, bvXcv = sd['bvY'], sd['bvXcv']
        t = 2.0*self.ccp
        js = self.jsparse
        F = np.zeros((nmo, nmo))

        def jpair(specA, specB, tt, JA, JB):
            return tt*(js(JB, *specA) + js(JA, *specB))
        # ---- J-type pairs ----
        if c['cJ1']:                                                        # CO1-J: (Yh^T, row [O2, C] = x_CO1)
            F += jpair(('Yh', Y), ('rc', O2, Cidx, xs['CO1']), t*c['cJ1'], Jm['y'], Jm['CO1'])
        if c['cW2']:                                                        # CO2-W: (e_O2 e_O1^T, row [O2, V] = w2)
            F += jpair(('unit', O2, O1, 1.0), ('rc', O2, Vidx, sd['w2']), t*c['cW2'], J_e21, Jm['W2'])
        if c['cW3']:                                                        # O1V-W: (col [C, O1] = u3, e_O2 e_O1^T)
            F += jpair(('rc', O1, Cidx, sd['u3']), ('unit', O2, O1, 1.0), t*c['cW3'], Jm['W3'], J_e21)
        if c['cJ4']:                                                        # O2V-J: (Yh^T, col [V, O1] = x_O2V)
            F += jpair(('Yh', Y), ('rc', O1, Vidx, xs['O2V']), t*c['cJ4'], Jm['y'], Jm['O2V'])
        if c['cJ6'] and xs['OO']:                                           # OO-J: (Yh^T, x_OO e_O2 e_O1^T)
            F += jpair(('Yh', Y), ('unit', O2, O1, xs['OO']), t*c['cJ6'], Jm['y'], xs['OO']*J_e21)
        # A_G two-electron part: (e_O2 e_O2^T, e_O1 e_O1^T) x (-c_H |y|^2)
        F += jpair(('unit', O2, O2, 1.0), ('unit', O1, O1, 1.0), -self.cH*sd['y2'], J_e22, J_e11)
        # covariant F' correction: (T^cv_VV, e_x e_x^T) and (T^cv_CC, e_x e_x^T), +/-(1 - c_H)
        tc = 1.0 - self.cH
        F += jpair(('Tvv', Y), ('unit', O2, O2, 1.0), tc, Jm['Tvv'], J_e22) \
            - jpair(('Tvv', Y), ('unit', O1, O1, 1.0), tc, Jm['Tvv'], J_e11)
        F += jpair(('Tcc', sd['Tcc']), ('unit', O2, O2, 1.0), tc, Jm['Tcc'], J_e22) \
            - jpair(('Tcc', sd['Tcc']), ('unit', O1, O1, 1.0), tc, Jm['Tcc'], J_e11)
        # CV Coulomb (singlet): (Yh^T, Yh^T) x 2
        if self.singlet:
            F += jpair(('Yh', Y), ('Yh', Y), 2.0, Jm['y'], Jm['y'])

        # ---- K-type pairs: K_Yh[p,r] = sum_Q sum_{s in C} bvY[Q,p,s] B^Q_{s r} ----
        cache = {}
        def K_Yh_rowO2():                                                   # K_Yh[O2, :]
            if 'rowO2' not in cache:
                cache['rowO2'] = mrsf_grad.bdot(np.ascontiguousarray(bvY[:, O2, :].T), 'core', nmo, naux)[:, 0]
            return cache['rowO2']
        def K_Yh_colO1():                                                   # K_Yh[:, O1]
            if 'colO1' not in cache:
                cache['colO1'] = np.einsum('Qps,Qs->p', bvY, BsoC[:, :, 0])
            return cache['colO1']
        # CO1-K: (Yh^T, col [C, O2] = x_CO1), K[x e_O2^T]_pq = sum_Q bvxCO1[Q,p] B^Q_{q O2}
        if c['cK1']:
            tt = t*c['cK1']; x = xs['CO1']; bv = sd['bvxCO1']
            F[:, Cidx] += tt*(bv.T @ (BsoV[:, :, 1] @ Y))                  # K[B] A^T
            F[:, Vidx] += tt*((bv[:, Cidx].T @ Bso[:, :, 1]).T @ Y.T)      # K[B]^T A
            F[:, Cidx] += tt*np.outer(K_Yh_rowO2(), x)                      # K[A] B^T
            u = np.einsum('scQ,c->Qs', self.BooQ[:self.nC, :self.nC, :], x)
            F[:, O2] += tt*np.einsum('Qps,Qs->p', bvY, u)                   # K[A]^T B
        # O2V-K: (col [V, O1] = x_O2V, Yh), K[x e_O1^T]_pq = sum_Q bvxO2V[Q,p] B^Q_{q O1}
        if c['cK4']:
            tt = t*c['cK4']; x = xs['O2V']; bv = sd['bvxO2V']
            F[:, Vidx] += tt*np.outer(K_Yh_colO1(), x)                      # K[B] A^T
            w = np.einsum('b,Qbs->Qs', x, bvY[:, Vidx, :])
            F[:, O1] += tt*mrsf_grad.bdot(np.ascontiguousarray(w.T), 'core', nmo, naux)[:, 0]   # K[B]^T A
            F[:, Vidx] += tt*((bv.T @ BsoC[:, :, 0]) @ Y.T)                 # K[A] B^T
            F[:, Cidx] += tt*(Bso[:, :, 0].T @ (bv[:, Vidx] @ Y))           # K[A]^T B
        # CV-A: (M12 = x_CV Y^T on V x V, e_O1 e_O2^T), K[M12]_pr = sum_Q sum_j bvXcv[Q,p,j] bvY[Q,r,j]
        if c['cA5']:
            tt = t*c['cA5']; xcv = xs['CV']
            F[:, Vidx] += tt*((K_e12[:, Vidx] @ Y) @ xcv.T)                 # K[B] A^T
            F[:, Vidx] += tt*((K_e12[Vidx, :].T @ xcv) @ Y.T)               # K[B]^T A
            F[:, O1] += tt*np.einsum('Qpj,Qj->p', bvXcv, bvY[:, O2, :])    # K[A] B^T: K[A][:, O2]
            F[:, O2] += tt*np.einsum('Qj,Qrj->r', bvXcv[:, O1, :], bvY)    # K[A]^T B: K[A][O1, :]
        # CV-B: (M12c = Y^T x_CV on C x C, e_O1 e_O2^T)
        if c['cB5']:
            tt = t*c['cB5']; M12c = Y.T @ xs['CV']
            F[:, Cidx] += tt*(K_e12[:, Cidx] @ M12c.T)                      # K[B] A^T
            F[:, Cidx] += tt*(K_e12[Cidx, :].T @ M12c)                      # K[B]^T A
            up = M12c @ BsoC[:, :, 1].T                                      # (nC, naux): sum_s' M12c[s,s'] B^Q_{s' O2}
            F[:, O1] += tt*mrsf_grad.bdot(np.ascontiguousarray(up), 'core', nmo, naux)[:, 0]        # K[A][:, O2]
            w = (BsoC[:, :, 0] @ M12c).T                                     # (nC, naux): sum_s B^Q_{O1 s} M12c[s,s']
            F[:, O2] += tt*mrsf_grad.bdot(np.ascontiguousarray(w), 'core', nmo, naux)[:, 0]         # K[A][O1, :]
        # OO-K: (Yh^T, x_OO e_O1 e_O2^T)
        if c['cK6'] and xs['OO']:
            tt = t*c['cK6']*xs['OO']
            F[:, Cidx] += tt*(K_e12[:, Vidx] @ Y)                           # K[B] A^T
            F[:, Vidx] += tt*(K_e12[Cidx, :].T @ Y.T)                       # K[B]^T A
            F[:, O1] += tt*K_Yh_rowO2()                                     # K[A] B^T
            F[:, O2] += tt*K_Yh_colO1()                                     # K[A]^T B
        return F

    # ------------------------------------------------------------------
    # families
    # ------------------------------------------------------------------
    def fam_bytes_state(self):
        base = super().fam_bytes_state()
        nC = self.nC
        return 2*base + 8.0*self.naux*self.nmo*(2*nC + 3) + 8.0*self.nvirb*self.nocca*self.naux

    def extra_families(self, sd):
        """CV exchange families, closed-shell mean-field pairs of M, the
        two-electron pair pieces and the A_G patches, added to the state's
        Ghh/Fhp/gpp in place (rank-1-per-Q pieces through grad_fam_outer, no
        family-size temporaries); returns the CV particle-particle family entry"""
        st, stcv = sd['st'], sd['stcv']
        Y, xs, c = sd['Y'], sd['xs'], self.coef
        nC, nV, nocca, nvirb, naux = self.nC, self.nV, self.nocca, self.nvirb, self.naux
        H, P, Cidx, Vidx, O1, O2 = self.H, self.P, self.Cidx, self.Vidx, self.O1, self.O2
        hO1, hO2 = self.hO1, self.hO2
        Hloc, Ploc, Bso, BsoC, BsoV = self.Hloc, self.Ploc, self.Bso, self.BsoC, self.BsoV
        Ghh, Fhp = st['Ghh'], st['Fhp']
        hC = [Hloc[p] for p in Cidx]
        t = 2.0*self.ccp
        x = self.x_ref
        jv, jy = sd['jv'], sd['jy']
        outer = mrsf_grad.fam_outer
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
        v = 0.5*x*y2*BsoC[:, :, 0].T                                  # (nC, naux)
        Ghh[hC, hO1, :] -= v
        Ghh[hO1, hC, :] -= v
        # (3) closed-shell mean-field pairs of M: J pair (M, D_G), K pair (M, D_G) with x/2
        M = sd['M']
        dqG = self.dq_G
        Mhh = M[np.ix_(H, H)]
        MHV = M[np.ix_(H, Vidx)]
        pM = sd['jTpp'] + sd['jTcc'] + sd['a1']*(BsoV[:, :, 1] @ sd['w1']) \
            + sd['a2']*(BsoC[:, :, 0] @ sd['u4']) + sd['a3']*jy
        occG_H = self.occG[H]
        outer(Ghh, 0, nocca, 0, nocca, [0.5*Mhh, 0.5*np.diag(occG_H)], [dqG, pM])
        outer(Fhp, 0, nocca, 2, nvirb, [MHV], [dqG])
        # K pair (M, D_G): Gamma -= x/2 (M B^Q 1_G + 1_G B^Q M)
        G_loc = [Hloc[p] for p in self.Gidx]
        isG = np.zeros(nocca); isG[G_loc] = 1.0
        BooQ2 = self.BooQ.reshape(nocca, nocca*naux, order='F')
        MB = (Mhh @ BooQ2).reshape(nocca, nocca, naux, order='F')     # (M_HH B^Q_HH)
        bvY = sd['bvY']                                               # (naux, nmo, nC): sum_c B^Q_qc Y_cj
        # cross terms from the H x V blocks of M: sum_b M[p,b] B^Q_{bq}
        MB[hC, :, :] += 0.5*sd['a3']*np.transpose(bvY[:, H, :], (2, 1, 0))
        MB[hO2, :, :] += 0.5*sd['a1']*sd['bvw1'][:, H].T
        Ghh -= 0.5*x*(MB*isG[None, :, None] + np.transpose(MB, (1, 0, 2))*isG[:, None, None])
        # particle rows: (M B^Q 1_G)_{aq}, a in V, q in G
        MaH = (M[np.ix_(Vidx, H)] @ BooQ2).reshape(nV, nocca, naux, order='F')
        Bt = np.ascontiguousarray(np.transpose(bvY[:, H, :], (2, 1, 0))).reshape(nC, nocca*naux)
        MaV = (Y @ Bt).reshape(nV, nocca, naux)
        Fq = MaH + MaV                                                # (nV, nocca, naux)
        Fhp[G_loc, 2:, :] -= x*np.transpose(Fq[:, G_loc, :], (1, 0, 2))
        # (4) two-electron pairs (1/2 t A b^Q + 1/2 t a^Q B for J-type;
        #     1/2 t (A B^Q B^T + A^T B^Q B) for K-type)
        # the (C, V) outer products Fhp[j, b, Q] += fac y_bj vec(Q) of the J-type pairs: one call
        mats, vecs = [], []
        if c['cJ1']:
            mats.append(0.5*t*c['cJ1']*Y.T); vecs.append(jv['CO1'])
            Ghh[hO2, hC, :] += 0.5*t*c['cJ1']*np.outer(xs['CO1'], jy)
        if c['cW2']:
            tt = t*c['cW2']
            Ghh[hO2, hO1, :] += 0.5*tt*jv['W2']
            Fhp[hO2, 2:, :] += 0.5*tt*np.outer(sd['w2'], Bso[:, O1, 1])
        if c['cW3']:
            tt = t*c['cW3']
            Ghh[hC, hO1, :] += 0.5*tt*np.outer(sd['u3'], Bso[:, O1, 1])
            Ghh[hO2, hO1, :] += 0.5*tt*jv['W3']
        if c['cJ4']:
            mats.append(0.5*t*c['cJ4']*Y.T); vecs.append(jv['O2V'])
            Fhp[hO1, 2:, :] += 0.5*t*c['cJ4']*np.outer(xs['O2V'], jy)
        if c['cJ6'] and xs['OO']:
            tt = t*c['cJ6']*xs['OO']
            mats.append(0.5*tt*Y.T); vecs.append(Bso[:, O1, 1])
            Ghh[hO2, hO1, :] += 0.5*tt*jy
        if self.singlet:                                              # CV Coulomb: (Yh^T, Yh^T) x 2
            mats.append(2.0*Y.T); vecs.append(jy)
        outer(Fhp, 0, nC, 2, nvirb, mats, vecs)
        # A_G two-electron: (e22, e11) x (-c_H y2)
        tt = -self.cH*y2
        Ghh[hO2, hO2, :] += 0.5*tt*Bso[:, O1, 0]
        Ghh[hO1, hO1, :] += 0.5*tt*Bso[:, O2, 1]
        # covariant F' correction, +/-(1 - c_H): PP via Yf_cv (partner Y), CC via Ghh
        tc = 1.0 - self.cH
        bdiff = Bso[:, O2, 1] - Bso[:, O1, 0]                         # b_2^Q - b_1^Q
        outer(Yf_cv, 0, nvirb, 0, nocca, [0.25*tc*sd['Yemb']], [bdiff])
        st['gpp'] += 0.5*tc*np.outer(sd['jTpp'], bdiff)
        outer(Ghh, 0, nC, 0, nC, [0.5*tc*sd['Tcc']], [bdiff])
        Ghh[hO2, hO2, :] += 0.5*tc*(sd['jTpp'] + sd['jTcc'])
        Ghh[hO1, hO1, :] -= 0.5*tc*(sd['jTpp'] + sd['jTcc'])
        # K-type: CO1-K
        if c['cK1']:
            tt = t*c['cK1']
            vq = Y.T @ BsoV[:, :, 1].T                                 # (nC, naux): sum_b y_bp B^Q_{b O2}
            Ghh[np.ix_(hC, hC)] += 0.5*tt*np.einsum('pQ,q->pqQ', vq, xs['CO1'])
            Fhp[hO2, 2:, :] += 0.5*tt*(Y @ sd['bvxCO1'][:, Cidx].T)
        # O2V-K: (col spec [q, O1] = x_O2V, Yh)
        if c['cK4']:
            tt = t*c['cK4']
            # PP piece 1/2 t x (Y B^Q_{C,O1})^T: spare column hO1 of the CV partner
            Yf_cv[2:, hO1, :] += 0.25*tt*(Y @ BsoC[:, :, 0].T)
            sd['Xp_hO1'] = xs['O2V']
            Ghh[hO1, hC, :] += 0.5*tt*(sd['bvxO2V'][:, Vidx] @ Y).T
        # CV-A: (M12, e_O1 e_O2^T)
        if c['cA5']:
            tt = t*c['cA5']
            M12 = xs['CV'] @ Y.T
            Fhp[hO1, 2:, :] += 0.5*tt*(M12 @ BsoV[:, :, 1].T)
            Fhp[hO2, 2:, :] += 0.5*tt*(M12.T @ BsoV[:, :, 0].T)
        # CV-B: (M12c, e_O1 e_O2^T)
        if c['cB5']:
            tt = t*c['cB5']
            M12c = Y.T @ xs['CV']
            Ghh[hC, hO1, :] += 0.5*tt*(M12c @ BsoC[:, :, 1].T)
            Ghh[hC, hO2, :] += 0.5*tt*(M12c.T @ BsoC[:, :, 0].T)
        # OO-K: (Yh^T, x_OO e_O1 e_O2^T)
        if c['cK6'] and xs['OO']:
            tt = t*c['cK6']*xs['OO']
            Ghh[hC, hO1, :] += 0.5*tt*(Y.T @ BsoV[:, :, 1].T)
            Fhp[hO2, 2:, :] += 0.5*tt*(Y @ BsoC[:, :, 0].T)
        # metric contribution of the spare column: 2 sum_a (B^Q_PP x)_a Yf^Q'_a
        Xp = np.array(sd['Yemb'])
        if 'Xp_hO1' in sd:
            Xp[2:, hO1] = sd['Xp_hO1']
            st['gpp'] += 2.0*(sd['bvxO2V'][:, P] @ Yf_cv[:, hO1, :])
        sd['Yf_cv'] = Yf_cv
        return [{'Xt': Xp, 'Yf': Yf_cv, 'g': np.zeros((naux, naux)), 'pp_only': True}]

    # ------------------------------------------------------------------
    # probe factors and the one-electron density
    # ------------------------------------------------------------------
    def P_quadratic(self, Xt, Za, Zb, pzeta=None):
        """relaxed densities with the A_G patch of T^b (|y|^2 C_O1 C_O1^T) on
        the hole-MO block of the beta form"""
        Bexp, terms = super().P_quadratic(Xt, Za, Zb, pzeta)
        nbb, idxb, Sb = terms['b']
        Sb[self.hO1, self.hO1] += self._cur['y2']
        return Bexp, terms

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
        return self.g_quadratic(sd)

    def g_quadratic(self, sd):
        """G-density terms of the extended method as quadratic forms in the
        basis [C_H | C_V Y | C_V w1]: M = T^cv + sym(N_F) and the CV trial
        density D_s = C_V Y C_C^T + h.c. (see M_matrix / M_factors); computed
        once per state"""
        if 'gq' in sd:
            return sd['gq']
        Y, nC, h = sd['Y'], self.nC, self.nocca
        CV = self.C[:, self.Vidx]
        Bexp = np.zeros((self.nao, nC + 1))
        Bexp[:, :nC] = CV @ Y
        Bexp[:, nC] = CV @ sd['w1']
        hC = np.array([self.Hloc[p] for p in self.Cidx])
        M = sd['M']
        nbM = h + nC + 1
        SM = np.zeros((nbM, nbM))
        SM[:h, :h] = M[np.ix_(self.H, self.H)]
        SM[h:h + nC, h:h + nC] = np.eye(nC)                               # T^cv_VV = Y Y^T
        for i in range(nC):                                               # M[C_i, V] = a3/2 Y^T_i
            SM[hC[i], h + i] = 0.5*sd['a3']; SM[h + i, hC[i]] = 0.5*sd['a3']
        SM[self.hO2, h + nC] = 0.5*sd['a1']; SM[h + nC, self.hO2] = 0.5*sd['a1']   # M[O2, V] = a1/2 w1
        idxM = np.concatenate([np.arange(1, h + 1), h + 1 + np.arange(nC + 1)])
        nbS = 2*nC
        SS = np.zeros((nbS, nbS))
        SS[:nC, nC:] = np.eye(nC); SS[nC:, :nC] = np.eye(nC)
        idxS = np.concatenate([hC + 1, h + 1 + np.arange(nC)])
        sd['gq'] = (Bexp, {'M': (nbM, idxM, SM), 'S': (nbS, idxS, SS)})
        return sd['gq']

    def extra_density(self, sd):
        return sd['M']
