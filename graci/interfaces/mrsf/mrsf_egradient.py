"""
Analytic gradients of the extended MRSF-TDDFT states (EMRSF): a
table-driven engine shared by the published parent G (ext_space =
published) and the O1-O2 rotation-invariant parent Pi = w_G G + w_D D
+ w_S S (ext_space = pi) of the extended method.

The energy of a state (x, y) above omega_M(x) is
  Tr(F_base M) + CV correction (T^cv, Dgamma) + E_W + 2 kappa_s (y|y)
  - c (jl|bd) y y + E_K + Tr(B H_pair) + sum_k tau_k E_k(x, y) + T4
with M = T^cv + sym(N), N the Fock-type coupling densities of the
parent tables (library getter mrsf_init.parent_tables, weighted by
w), the two-electron couplings E_k as J/K-type pairs of sparse/low-rank
operands read off the same tables, B the pair-block density (|y|^2
w w^T plus the w-response z_w of the Pi parent; |y|^2 e_G e_G^T for
the published parent) and T4 the HF-type SOMO tensor of the per-parent
coupling operators (Pi only). The published model is the special case
w = (1, 0, 0) with D_base = D_G, the covariant F' correction and the
canonical multiplier zeta_12 of the O1-O2 rotation; every term of the
published driver of Part XII is reproduced by the generic code.
"""
import sys
import numpy as np
import graci.io.output as output
import graci.interfaces.mrsf.mrsf_grad as mrsf_grad
import graci.interfaces.mrsf.mrsf_init as mrsf_init
import graci.interfaces.mrsf.mrsf_pi as mrsf_pi
from graci.interfaces.mrsf.mrsf_gradient import GradientDriver, _sym, SQ2

PAIRS_XY = ((0, 0), (1, 1), (0, 1))      # xy index of the W tables: (O1,O1), (O2,O2), (O1,O2)
KOO_HOLE = (0, 1, 0, 1)                  # (hole, particle) SOMO of the K_oo classes k = 0..3
KOO_PART = (1, 0, 0, 1)

# the G-parent table in the names of the Part XII driver (documentation and
# the table gate; the engine reads the library tables)
COEF = {1: dict(cJ1=-2.0, cK1=1.0, cF1=-1.0, cW2=-1.0, cW3=1.0, cJ4=2.0, cK4=-1.0, cF4=-1.0,
                cA5=-1.0, cB5=1.0, cK6=1.0, cJ6=-2.0, cG=SQ2),
        3: dict(cJ1=0.0, cK1=1.0, cF1=-1.0, cW2=-1.0, cW3=-1.0, cJ4=0.0, cK4=1.0, cF4=1.0,
                cA5=-1.0, cB5=-1.0, cK6=1.0, cJ6=0.0, cG=0.0)}


def _nz(v):
    return v is not None and np.abs(v).max() > 0.0


# ----------------------------------------------------------------------
# sparse / low-rank MO operands ("pieces") of the two-electron pairs
#   ('row', o, cls, v)   A[O_o, idx(cls)] = v
#   ('col', cls, o, v)   A[idx(cls), O_o] = v
#   ('cv', a, Y)         A[c, b] = a y_bc           (C x V)
#   ('vv', L, R)         A[a, b] = sum_j L_aj R_bj  (V x V)
#   ('cc', M)            A[c, d] = M_cd             (C x C)
#   ('unit', x, y, val)  A = val e_{O_x} e_{O_y}^T
# J-type pair (A, B, t): energy t sum_pqrs A_pq (pq|rs) B_rs
# K-type pair (A, B, t): energy t sum_pqrs A_pq B_rs (pr|qs)
# ----------------------------------------------------------------------
class EMRSFGradientDriver(GradientDriver):

    def __init__(self, ci):
        self.space = str(getattr(ci, 'ext_space', 'published')).lower()
        if self.space == 'published' and str(getattr(ci, 'fprime', 'covariant')).lower() != 'covariant':
            sys.exit('$mrsftddft: gradients of the extended method need the covariant '
                     "F' correction (fprime = covariant in the $mrsftddft section)")
        super().__init__(ci)
        # explicit amplitude vectors (derivative couplings): stacked (nvirb, nocca + nC)
        # arrays [X~ | Y with zero SOMO rows], see expanded_amplitudes
        self.supports_xts = True
        C, nmo, nao = self.C, self.nmo, self.nao
        self.ncol = int(ci.ncol)
        self.singlet = (self.mult == 1)
        self.ext_model = str(getattr(ci, 'ext_model', 'kernel')).lower()
        self.use_kernel = bool(self.is_dft and self.singlet and self.ext_model == 'kernel')
        self.kappa_s = 1.0 if self.ext_model == 'kernel' else float(self.kappa[2])
        self.ccp = float(ci.ccp_eff)
        self.cH = float(ci.chf_eff)
        self.t2 = 2.0*self.ccp
        self.O1, self.O2 = int(self.Oidx[0]), int(self.Oidx[1])
        self.Ol = [self.O1, self.O2]
        self.hO1, self.hO2 = self.Hloc[self.O1], self.Hloc[self.O2]
        self.hO = [self.hO1, self.hO2]
        self.hC = [self.Hloc[p] for p in self.Cidx]
        self.Cl = [int(c) for c in self.Cidx]; self.Vl = [int(v) for v in self.Vidx]
        self.tabs = [mrsf_init.parent_tables(self.mult, X) for X in (1, 2, 3)]
        with self.tm('setup: BooQ'):
            self.BooQ = mrsf_grad.booq(self.dims)
        with self.tm('setup: pair constants'):
            Bso = self.Bso
            self.BsoC = np.ascontiguousarray(Bso[:, self.Cidx, :])       # B^Q_{c O_x} (naux, nC, 2)
            self.BsoV = np.ascontiguousarray(Bso[:, self.Vidx, :])       # B^Q_{v O_x} (naux, nV, 2)
            self.K_e = {}
            for x in (0, 1):
                for y in (0, 1):
                    self.K_e[(x, y)] = Bso[:, :, x].T @ Bso[:, :, y]      # (p O_x|q O_y)
            JHH, JHP, JPP = mrsf_grad.jblocks(np.column_stack([Bso[:, self.O1, 0], Bso[:, self.O2, 1],
                                                               Bso[:, self.O2, 0]]), self.dims)
            Je = [self.full_from_blocks(JHH[:, :, k], JHP[:, :, k], JPP[:, :, k]) for k in range(3)]
            self.J_e = {(0, 0): Je[0], (1, 1): Je[1], (0, 1): Je[2], (1, 0): Je[2]}
        # ---- model pieces ----
        self.n_base = np.zeros(nmo)
        self.n_base[self.Cidx] = 2.0
        self.pi = None
        if self.space == 'published':
            self.w = np.array([1.0, 0.0, 0.0])
            self.n_base[self.O1] = 2.0
            self.Gidx = self.Cl + [self.O1]
            self.cv_corr = (np.array([[-1.0, 0.0], [0.0, 1.0]]), 1.0 - self.cH, 0.0, 0.0)
            self.has_T4 = False
            extraN = [(self.O2, self.O1)]
            with self.tm('setup: F_G'):
                self.F_base = self.session.fock_g(ci)                  # once per SCF (shared with the energies)
            self.with_kxc = self.use_kernel
            if self.is_dft:
                with self.tm('setup: G kernel set'):
                    self.xcgrid.ensure_g_kernel(self.Gidx, self.use_kernel, g_default=self.session.gidx)
        else:
            self.n_base[[self.O1, self.O2]] = 1.0
            with self.tm('setup: pair block'):
                ints = mrsf_pi.somo_integrals(self.Bso, self.O1, self.O2)
                self.Hp = mrsf_pi.pair_block(self.fa, self.fb, self.O1, self.O2, ints, self.cH)
                w, lam = mrsf_pi.pi_weights(self.Hp)
                self.lam_pair = np.asarray(lam, float)
                self.U = np.linalg.eigh(0.5*(self.Hp + self.Hp.T))[1]
                if self.U[0, 0] < 0.0:
                    self.U[:, 0] = -self.U[:, 0]
                pinfo = getattr(ci, 'pi_info', None)
                if pinfo is not None and np.abs(np.asarray(pinfo['w'], float) - w).max() > 1e-7:
                    sys.exit('$mrsftddft: the parent weights of the energy section differ from the '
                             'recomputed pair block (max difference %.2e)' %
                             np.abs(np.asarray(pinfo['w'], float) - w).max())
                self.w = w
            self.pair_gap = float(self.lam_pair[1] - self.lam_pair[0])
            if self.pair_gap < 1.0e-2:
                output.print_message('WARNING: pair-block gap %.4f Ha: the EMRSF-Pi gradient carries the '
                                     'w-response z_w ~ 1/gap' % self.pair_gap)
            self.dgam = mrsf_pi.delta_gamma_pi(self.w)
            self.cv_corr = (self.dgam, self.cH, -0.5*self.cH, 1.0 if self.is_dft else 0.0)
            self.has_T4 = True
            extraN = []
            with self.tm('setup: F_ens'):
                self.F_base = self.session.fock_ens()
            self.with_kxc = bool(self.is_dft)
            if self.is_dft:
                with self.tm('setup: ens kernel set'):
                    self.session.ensure_kernel_ens(True)
                    vdg = getattr(self.session, 'pi_vdg', None)
                    if vdg is not None and np.abs(vdg[0] - self.dgam).max() < 1e-12:
                        self.V_dg = vdg[1]                                 # from the energy section's pi_setup
                    else:
                        dm = np.zeros((nmo, nmo))
                        dm[np.ix_(self.Ol, self.Ol)] = self.dgam
                        self.V_dg = C.T @ self.session.fxc_potential_ao(C @ dm @ C.T) @ C
                        self.session.pi_vdg = (np.array(self.dgam), self.V_dg)
                    self.xcgrid.set_dgamma(self.dgam)
            self.pi = {'w': self.w, 'pair': self.lam_pair, 'gap': self.pair_gap}
        self.lam = self.n_base[self.H]/2.0
        self.dq_base = mrsf_grad.dq(self.n_base, self.naux)
        self.tab = {k: sum(self.w[X]*self.tabs[X][k] for X in range(3)) for k in self.tabs[0]}
        if self.is_dft:
            self.xcgrid.set_hole_weights(self.lam)
        # probe flags: the kernel-model terms (f[D_s], D_s), (k[D_s, D_s], D_base) only with the CV
        # kernel; the E_W terms of the Pi parent (which need the third derivatives installed) for
        # every DFT Pi case
        self.probe_with_kxc = self.use_kernel
        self.probe_with_w = bool(self.space == 'pi' and self.is_dft)
        # canonical pairs: frozen-active core (base) + the O1-O2 rotation of the published parent
        self.N = list(self.N) + extraN
        fsum = self.fa + self.fb
        self.hN = np.array([fsum[k, k] - fsum[m, m] for (k, m) in self.N])
        self.somo_gap = None
        if self.space == 'published':
            self.somo_gap = 0.5*(fsum[self.O2, self.O2] - fsum[self.O1, self.O1])
            if self.somo_gap < 1.0e-2:
                output.print_message('WARNING: SOMO gap %.4f Ha: the extended-method gradient depends on '
                                     'the O1-O2 multiplier zeta_12 ~ 1/gap' % self.somo_gap)
        self._cur = None
        self.pi_diag = []

    # ------------------------------------------------------------------
    # amplitudes
    # ------------------------------------------------------------------
    def amplitudes(self, istate):
        """expanded MRSF amplitudes X~ (nvirb, nocca) and the CV amplitudes Y (nV, nC)"""
        irr, st = self.ci.state_sym(istate)
        x = self.ci.amps['adiabatic'][irr][:, st]
        X = np.array(x.reshape(self.nvirb, self.ncol, order='F'))
        Xm = np.array(X[:, :self.nocca])
        Y = np.array(X[2:, self.nocca:])
        if self.Cf:
            fcols = [self.Hloc[k] for k in self.Cf]
            if max(np.abs(Xm[:, fcols]).max(), np.abs(Y[:, fcols]).max()) > 1e-12:
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
        return Xm, Y

    def expanded_amplitudes(self, istate):
        """the stacked expanded amplitudes (nvirb, nocca + nC) of a state:
        [X~ | Y] with the CV amplitudes Y in the virtual rows of the last nC
        columns (SOMO rows zero); every term of the extended response
        functional is a quadratic form in this array (no excitation-energy
        term), so sums and differences of states define the couplings"""
        Xm, Y = self.amplitudes(istate)
        V = np.zeros((self.nvirb, self.nocca + self.nC))
        V[:, :self.nocca] = Xm
        V[2:, self.nocca:] = Y
        return V

    def split_amplitudes(self, V):
        """stacked array -> (X~, Y)"""
        V = np.asarray(V)
        return np.array(V[:, :self.nocca]), np.array(V[2:, self.nocca:])

    def slices(self, Xt):
        """compressed MRSF amplitudes of the coupling classes from X~"""
        nC, hO1, hO2 = self.nC, self.nC, self.nC + 1
        return {'CO': [Xt[0, :nC].copy(), Xt[1, :nC].copy()],
                'OV': [Xt[2:, hO1].copy(), Xt[2:, hO2].copy()],
                'CV': Xt[2:, :nC].copy(),
                'slot': [float(Xt[0, hO1])*SQ2, float(Xt[0, hO2]), float(Xt[1, hO1])]}

    # ------------------------------------------------------------------
    # Fock-type coupling densities and the pair lists of a table
    # ------------------------------------------------------------------
    def N_pieces(self, tb, xs, Y):
        """N (scaled by 2 c_cp) as rows [O_x, V], columns [C, O_x] and the
        (C, V) coefficient a of a Y^T"""
        t = self.t2
        rows, cols = [None, None], [None, None]
        for x in (0, 1):
            v = np.zeros(self.nV); u = np.zeros(self.nC)
            for s in (0, 1):
                if tb['Fco'][s, x] != 0.0:
                    v += t*tb['Fco'][s, x]*(Y @ xs['CO'][s])
                if tb['Fov'][s, x] != 0.0:
                    u += t*tb['Fov'][s, x]*(Y.T @ xs['OV'][s])
            rows[x] = v if _nz(v) else None
            cols[x] = u if _nz(u) else None
        a = t*sum(tb['Fblk'][k]*xs['slot'][k] for k in range(3))
        return {'rows': rows, 'cols': cols, 'a': float(a)}

    def N_matrix(self, N):
        M = np.zeros((self.nmo, self.nmo))
        for x in (0, 1):
            if N['rows'][x] is not None:
                M[self.Ol[x], self.Vl] += N['rows'][x]
            if N['cols'][x] is not None:
                M[self.Cl, self.Ol[x]] += N['cols'][x]
        if N['a'] != 0.0:
            M[np.ix_(self.Cl, self.Vl)] += N['a']*self._cur['Y'].T
        return M

    def N_as_pieces(self, N, Y):
        out = []
        for x in (0, 1):
            if N['rows'][x] is not None:
                out.append(('row', x, 'V', N['rows'][x]))
            if N['cols'][x] is not None:
                out.append(('col', 'C', x, N['cols'][x]))
        if N['a'] != 0.0:
            out.append(('cv', N['a'], Y))
        return out

    def table_pairs(self, tb, xs, Y, scale=1.0):
        """J/K-type pairs of the two-electron coupling classes of a table"""
        t = self.t2*scale
        cv = ('cv', 1.0, Y)
        xcv = xs['CV']
        pairs = []
        for s in (0, 1):
            for x in (0, 1):
                c = tb['Jco'][s, x]
                if c:
                    pairs.append(('J', ('row', x, 'C', xs['CO'][s]), cv, t*c))
                c = tb['Kco'][s, x]
                if c:
                    pairs.append(('K', ('col', 'C', x, xs['CO'][s]), cv, t*c))
                c = tb['Jov'][s, x]
                if c:
                    pairs.append(('J', ('row', x, 'V', xs['OV'][s]), cv, t*c))
                c = tb['Kov'][s, x]
                if c:
                    pairs.append(('K', ('row', x, 'V', xs['OV'][s]), cv, t*c))
            for xy in range(3):
                a, b = PAIRS_XY[xy]
                for z in (0, 1):
                    c = tb['Wco'][s, xy, z]
                    if c:
                        pairs.append(('J', ('unit', a, b, 1.0), ('row', z, 'V', Y @ xs['CO'][s]), t*c))
                    c = tb['Wov'][s, xy, z]
                    if c:
                        pairs.append(('J', ('unit', a, b, 1.0), ('col', 'C', z, Y.T @ xs['OV'][s]), t*c))
        for x in (0, 1):
            for y in (0, 1):
                c = tb['A'][x, y]
                if c:
                    pairs.append(('K', ('vv', xcv, Y), ('unit', x, y, 1.0), t*c))
                c = tb['B'][x, y]
                if c:
                    pairs.append(('K', ('cc', Y.T @ xcv), ('unit', x, y, 1.0), t*c))
        for slot in range(3):
            xv = xs['slot'][slot]
            if xv == 0.0:
                continue
            c = tb['Joo'][slot]
            if c:
                pairs.append(('J', ('unit', 0, 1, xv), cv, t*c))
            for k in range(4):
                c = tb['Koo'][slot, k]
                if c:                                      # (A, B) symmetric: the SOMO unit as the B operand
                    pairs.append(('K', cv, ('unit', KOO_HOLE[k], KOO_PART[k], xv), t*c))
        return pairs

    def model_pairs(self, sd):
        """the fixed two-electron pairs of the extended energy: CV Coulomb,
        the pair-block two-electron part E2(B), the CV correction on T^cv
        and (Pi) the SOMO tensor T4 of the coupling operators"""
        Y, c = sd['Y'], self.cH
        B = sd['B']
        cv = ('cv', 1.0, Y)
        Tvv = ('vv', Y, Y); Tcc = ('cc', sd['Tcc'])
        pairs = []
        if self.singlet:
            pairs.append(('J', cv, cv, 2.0*self.kappa_s))
        s_ = 0.5*B[2, 2]
        t_ = (B[0, 2] + B[1, 2])/SQ2
        e11, e22, e12 = ('unit', 0, 0, 1.0), ('unit', 1, 1, 1.0), ('unit', 0, 1, 1.0)
        for A, Bp, val in ((e11, e22, -c*(B[0, 0] + B[1, 1])), (e11, e11, -c*s_), (e22, e22, -c*s_),
                           (e12, e12, -c*(2.0*B[0, 1] - B[2, 2])), (e12, e11, -2.0*c*t_), (e12, e22, 2.0*c*t_)):
            if val != 0.0:
                pairs.append(('J', A, Bp, val))
        dg, cJ, cK, _ = self.cv_corr
        for x in (0, 1):
            for y in (0, 1):
                if dg[x, y] == 0.0:
                    continue
                if cJ:
                    pairs.append(('J', Tvv, ('unit', x, y, 1.0), cJ*dg[x, y]))
                    pairs.append(('J', Tcc, ('unit', x, y, 1.0), cJ*dg[x, y]))
                if cK:
                    pairs.append(('K', Tvv, ('unit', x, y, 1.0), cK*dg[x, y]))
                    pairs.append(('K', Tcc, ('unit', x, y, 1.0), cK*dg[x, y]))
        if self.has_T4:
            for P in self.N_as_pieces(sd['Nt'], Y):
                pairs.append(('J', P, e11, 1.0)); pairs.append(('J', P, e22, -1.0))
                pairs.append(('K', P, e11, -0.5)); pairs.append(('K', P, e22, 0.5))
        return pairs

    # ------------------------------------------------------------------
    # primitives on the pieces
    # ------------------------------------------------------------------
    def idx(self, cls):
        return self.Cl if cls == 'C' else self.Vl

    def Bs(self, cls):
        return self.BsoC if cls == 'C' else self.BsoV

    def bv(self, P):
        """DF products sum_q B^Q_tq v_q of a row/col piece (from the batch)"""
        sd = self._cur
        key = id(P[3]) if P[0] in ('row', 'col') else None
        if key not in sd['bvmap']:
            cls = P[2] if P[0] == 'row' else P[1]
            sd['bvmap'][key] = mrsf_grad.bvec(P[3], 'core' if cls == 'C' else 'virt', self.nmo, self.naux)[:, :, 0]
        return sd['bvmap'][key]

    def bv_of(self, L):
        """DF products of a (nV, k) matrix on the virtuals (Y and x_CV cached)"""
        sd = self._cur
        key = id(L)
        if key not in sd['bvmap']:
            sd['bvmap'][key] = mrsf_grad.bvec(L, 'virt', self.nmo, self.naux)
        return sd['bvmap'][key]

    def jvec(self, P):
        """Coulomb vector a^Q = sum_pq B^Q_pq A_pq"""
        sd = self._cur
        kind = P[0]
        if kind == 'row':
            return self.Bs(P[2])[:, :, P[1]] @ P[3]
        if kind == 'col':
            return self.Bs(P[1])[:, :, P[2]] @ P[3]
        if kind == 'cv':
            return P[1]*sd['jy']
        if kind == 'vv':
            if P[1] is sd['Y'] and P[2] is sd['Y']:
                return sd['jTpp']
            return np.einsum('aj,Qaj->Q', P[1], self.bv_of(P[2])[:, self.Vidx, :])
        if kind == 'cc':
            if P[1] is sd['Tcc']:
                return sd['jTcc']
            return np.einsum('ijQ,ij->Q', self.BooQ[:self.nC, :self.nC, :], P[1])
        if kind == 'unit':
            return P[3]*self.Bso[:, self.Ol[P[2]], P[1]]
        raise ValueError(kind)

    def Japply(self, J, P):
        """J (A + A^T)"""
        out = np.zeros((self.nmo, self.nmo))
        kind = P[0]
        Cl, Vl, Ol = self.Cl, self.Vl, self.Ol
        if kind in ('row', 'col'):
            if kind == 'row':
                o, idx, v = Ol[P[1]], self.idx(P[2]), P[3]
            else:
                o, idx, v = Ol[P[2]], self.idx(P[1]), P[3]
            out[:, idx] += np.outer(J[:, o], v)
            out[:, o] += J[:, idx] @ v
        elif kind == 'cv':
            a, Y = P[1], P[2]
            out[:, Vl] += a*(J[:, Cl] @ Y.T)
            out[:, Cl] += a*(J[:, Vl] @ Y)
        elif kind == 'vv':
            L, R = P[1], P[2]
            out[:, Vl] += (J[:, Vl] @ R) @ L.T + (J[:, Vl] @ L) @ R.T
        elif kind == 'cc':
            M = P[1]
            out[:, Cl] += J[:, Cl] @ (M + M.T)
        elif kind == 'unit':
            x, y, val = P[1], P[2], P[3]
            out[:, Ol[y]] += val*J[:, Ol[x]]
            out[:, Ol[x]] += val*J[:, Ol[y]]
        else:
            raise ValueError(kind)
        return out

    def Kapply(self, K, P):
        """K A^T + K^T A for a full matrix K"""
        out = np.zeros((self.nmo, self.nmo))
        kind = P[0]
        Cl, Vl, Ol = self.Cl, self.Vl, self.Ol
        if kind == 'cv':
            a, Y = P[1], P[2]
            out[:, Cl] += a*(K[:, Vl] @ Y)
            out[:, Vl] += a*(K[Cl, :].T @ Y.T)
        elif kind == 'vv':
            L, R = P[1], P[2]
            out[:, Vl] += (K[:, Vl] @ R) @ L.T + (K[Vl, :].T @ L) @ R.T
        elif kind == 'cc':
            M = P[1]
            out[:, Cl] += K[:, Cl] @ M.T + K[Cl, :].T @ M
        elif kind == 'row':
            o, idx, v = Ol[P[1]], self.idx(P[2]), P[3]
            out[:, o] += K[:, idx] @ v
            out[:, idx] += np.outer(K[o, :], v)
        elif kind == 'col':
            o, idx, v = Ol[P[2]], self.idx(P[1]), P[3]
            out[:, idx] += np.outer(K[:, o], v)
            out[:, o] += K[idx, :].T @ v
        else:
            raise ValueError(kind)
        return out

    def pkey(self, P, x):
        """memoisation key of a piece and a SOMO index (per state)"""
        kind = P[0]
        if kind == 'cv':
            return ('cv', id(P[2]), float(P[1]), x)
        if kind == 'vv':
            return ('vv', id(P[1]), id(P[2]), x)
        if kind == 'cc':
            return ('cc', id(P[1]), x)
        if kind == 'unit':
            return ('unit', P[1], P[2], float(P[3]), x)
        return (kind, id(P[3]), P[1], P[2], x)

    def K_cols(self, P, x):
        """(K[A][:, O_x], K[A^T][:, O_x]) with K[A]_pq = sum_rs A_rs (pr|qs) (memoised per state)"""
        sd = self._cur
        if P[0] == 'cv' and P[1] != 1.0:
            k1, k2 = self.K_cols(('cv', 1.0, P[2]), x)
            return P[1]*k1, P[1]*k2
        key = ('K',) + self.pkey(P, x)
        if key not in sd['kcache']:
            sd['kcache'][key] = self._K_cols(P, x)
        return sd['kcache'][key]

    def _K_cols(self, P, x):
        sd = self._cur
        kind = P[0]
        Bso, BsoC, BsoV = self.Bso, self.BsoC, self.BsoV
        nmo, naux = self.nmo, self.naux
        bdot = mrsf_grad.bdot
        if kind == 'cv':
            a, Y = P[1], P[2]
            bvY = self.bv_of(Y)
            k1 = a*bdot(np.ascontiguousarray(bvY[:, self.Ol[x], :].T), 'core', nmo, naux)[:, 0]
            k2 = a*np.einsum('Qpc,Qc->p', bvY, BsoC[:, :, x])
            return k1, k2
        if kind == 'vv':
            L, R = P[1], P[2]
            bvL, bvR = self.bv_of(L), self.bv_of(R)
            return (np.einsum('Qpj,Qj->p', bvL, bvR[:, self.Ol[x], :]),
                    np.einsum('Qpj,Qj->p', bvR, bvL[:, self.Ol[x], :]))
        if kind == 'cc':
            M = P[1]
            return (bdot(np.ascontiguousarray(M @ BsoC[:, :, x].T), 'core', nmo, naux)[:, 0],
                    bdot(np.ascontiguousarray(M.T @ BsoC[:, :, x].T), 'core', nmo, naux)[:, 0])
        if kind in ('row', 'col'):
            if kind == 'row':
                o, cls, v = P[1], P[2], P[3]
            else:
                o, cls, v = P[2], P[1], P[3]
            bvv = self.bv(P)
            direct = Bso[:, :, o].T @ (self.Bs(cls)[:, :, x] @ v)     # K[row] column, K[col^T] column
            trans = bvv.T @ Bso[:, self.Ol[o], x]
            return (direct, trans) if kind == 'row' else (trans, direct)
        raise ValueError(kind)

    def AB_cols(self, P, x):
        """((A B^Q)[:, O_x], (A^T B^Q)[:, O_x]) as placement lists
        [(rowclass, array)], rowclass 'C' (nC, naux), 'V' (nV, naux) or
        ('O', o) (naux) (memoised per state)"""
        sd = self._cur
        if P[0] == 'cv' and P[1] != 1.0:
            AB, AtB = self.AB_cols(('cv', 1.0, P[2]), x)
            return [(c, P[1]*a) for c, a in AB], [(c, P[1]*a) for c, a in AtB]
        key = ('AB',) + self.pkey(P, x)
        if key not in sd['kcache']:
            sd['kcache'][key] = self._AB_cols(P, x)
        return sd['kcache'][key]

    def _AB_cols(self, P, x):
        sd = self._cur
        kind = P[0]
        Bso, BsoC, BsoV = self.Bso, self.BsoC, self.BsoV
        if kind == 'cv':
            a, Y = P[1], P[2]
            bvY = self.bv_of(Y)
            return [('C', a*bvY[:, self.Ol[x], :].T)], [('V', a*(Y @ BsoC[:, :, x].T))]
        if kind == 'vv':
            L, R = P[1], P[2]
            bvL, bvR = self.bv_of(L), self.bv_of(R)
            return [('V', L @ bvR[:, self.Ol[x], :].T)], [('V', R @ bvL[:, self.Ol[x], :].T)]
        if kind == 'cc':
            M = P[1]
            return [('C', M @ BsoC[:, :, x].T)], [('C', M.T @ BsoC[:, :, x].T)]
        if kind in ('row', 'col'):
            if kind == 'row':
                o, cls, v = P[1], P[2], P[3]
            else:
                o, cls, v = P[2], P[1], P[3]
            vecQ = self.Bs(cls)[:, :, x] @ v                           # sum_q B^Q_{q O_x} v_q
            outer = np.outer(v, Bso[:, self.Ol[x], o])                 # v_q B^Q_{O_o O_x}
            if kind == 'row':
                return [(('O', o), vecQ)], [(cls, outer)]
            return [(cls, outer)], [(('O', o), vecQ)]
        raise ValueError(kind)

    def place_col(self, Ghh, Fhp, target_x, entries, scale):
        """family placement of column O_target pieces"""
        ht = self.hO[target_x]
        for cls, arr in entries:
            if cls == 'C':
                Ghh[self.hC, ht, :] += scale*arr
            elif cls == 'V':
                Fhp[ht, 2:, :] += scale*arr
            else:
                Ghh[self.hO[cls[1]], ht, :] += scale*arr

    def famJ(self, Ghh, Fhp, P, vec, t):
        """1/2 t A (x) vec into the families (PP pieces of T^cv_VV into the
        CV partner family of the state)"""
        sd = self._cur
        kind = P[0]
        s = 0.5*t
        if kind == 'row':
            o, cls, v = self.hO[P[1]], P[2], P[3]
            if cls == 'V':
                Fhp[o, 2:, :] += s*np.outer(v, vec)
            else:
                Ghh[o, self.hC, :] += s*np.outer(v, vec)
        elif kind == 'col':
            cls, o, v = P[1], self.hO[P[2]], P[3]
            if cls == 'V':
                Fhp[o, 2:, :] += s*np.outer(v, vec)
            else:
                Ghh[self.hC, o, :] += s*np.outer(v, vec)
        elif kind == 'cv':
            mrsf_grad.fam_outer(Fhp, 0, self.nC, 2, self.nvirb, [s*P[1]*P[2].T], [vec])
        elif kind == 'vv':
            if not (P[1] is sd['Y'] and P[2] is sd['Y']):
                raise NotImplementedError('J-type pair with a general V x V operand')
            sd['pp_mats'].append(0.5*s*sd['Yemb']); sd['pp_vecs'].append(vec)
            sd['st']['gpp'] += s*np.outer(sd['jTpp'], vec)
        elif kind == 'cc':
            mrsf_grad.fam_outer(Ghh, 0, self.nC, 0, self.nC, [s*P[1]], [vec])
        elif kind == 'unit':
            Ghh[self.hO[P[1]], self.hO[P[2]], :] += s*P[3]*vec
        else:
            raise ValueError(kind)

    # ------------------------------------------------------------------
    # evaluation of the pair lists: energies, generalised Fock, families
    # ------------------------------------------------------------------
    def pair_energy(self, pair):
        kind, A, B, t = pair
        if kind == 'J':
            return t*float(np.dot(self.jvec(A), self.jvec(B)))
        if B[0] == 'unit':
            x, y, val = B[1], B[2], B[3]
            _, AtB = self.AB_cols(A, x)
            e = 0.0
            for cls, arr in AtB:
                if cls == 'C':
                    e += float(np.einsum('cQ,Qc->', arr, self.BsoC[:, :, y]))
                elif cls == 'V':
                    e += float(np.einsum('bQ,Qb->', arr, self.BsoV[:, :, y]))
                else:
                    e += float(np.dot(arr, self.Bso[:, self.Ol[y], cls[1]]))
            return t*val*e
        # B = cv: A = col(C, O_x, v) or row(O_x, V, v)
        k1, k2 = self.K_cols(B, A[2] if A[0] == 'col' else A[1])
        if A[0] == 'col':
            return t*float(np.dot(A[3], k1[self.Cl]))
        return t*float(np.dot(A[3], k2[self.Vl]))

    def pairs_gfock(self, pairs, Jm):
        nmo = self.nmo
        F = np.zeros((nmo, nmo))
        Ol = self.Ol
        for pair in pairs:
            kind, A, B, t = pair
            if kind == 'J':
                F += t*(self.Japply(Jm[self.jkey(B)], A) + self.Japply(Jm[self.jkey(A)], B))
                continue
            if B[0] == 'unit':
                x, y, val = B[1], B[2], B[3]
                F += t*val*self.Kapply(self.K_e[(x, y)], A)
                kA, kAt = self.K_cols(A, y)
                F[:, Ol[x]] += t*val*kA
                kA, kAt = self.K_cols(A, x)
                F[:, Ol[y]] += t*val*kAt
                continue
            # B = cv(1, Y)
            Y = B[2]
            bvY = self.bv_of(Y)
            Bso, BsoC = self.Bso, self.BsoC
            if A[0] == 'col':                                  # K_co: A[c, O_x] = v_c
                x, v = A[2], A[3]
                bvx = self.bv(A)
                kcv, kcvt = self.K_cols(B, x)
                F[:, self.Cl] += t*np.outer(kcv, v)                                   # K[B] A^T
                u = np.einsum('scQ,c->Qs', self.BooQ[:self.nC, :self.nC, :], v)
                F[:, Ol[x]] += t*np.einsum('Qpc,Qc->p', bvY, u)                       # K[B]^T A
                F[:, self.Cl] += t*(bvx.T @ bvY[:, Ol[x], :])                         # K[A] B^T
                F[:, self.Vl] += t*(Bso[:, :, x].T @ (bvx[:, self.Cl] @ Y.T))        # K[A]^T B
            else:                                              # K_ov: A[O_x, b] = v_b
                x, v = A[1], A[3]
                bvx = self.bv(A)
                kcv, kcvt = self.K_cols(B, x)
                w = np.einsum('b,Qbs->sQ', v, bvY[:, self.Vidx, :])
                F[:, Ol[x]] += t*mrsf_grad.bdot(np.ascontiguousarray(w), 'core', nmo, self.naux)[:, 0]  # K[B] A^T
                F[:, self.Vl] += t*np.outer(kcvt, v)                                  # K[B]^T A
                F[:, self.Cl] += t*(Bso[:, :, x].T @ (bvx[:, self.Vl] @ Y))          # K[A] B^T
                F[:, self.Vl] += t*(bvx.T @ (BsoC[:, :, x] @ Y.T))                    # K[A]^T B
        return F

    def pairs_families(self, pairs, Ghh, Fhp):
        sd = self._cur
        for pair in pairs:
            kind, A, B, t = pair
            if kind == 'J':
                self.famJ(Ghh, Fhp, A, self.jvec(B), t)
                self.famJ(Ghh, Fhp, B, self.jvec(A), t)
                continue
            if B[0] == 'unit':
                x, y, val = B[1], B[2], B[3]
                AB, _ = self.AB_cols(A, y)
                self.place_col(Ghh, Fhp, x, AB, 0.5*t*val)
                _, AtB = self.AB_cols(A, x)
                self.place_col(Ghh, Fhp, y, AtB, 0.5*t*val)
                continue
            Y = B[2]
            bvY = self.bv_of(Y)
            if A[0] == 'col':                                  # K_co
                x, v = A[2], A[3]
                bvx = self.bv(A)
                Ghh[np.ix_(self.hC, self.hC)] += 0.5*t*np.einsum('c,Qd->cdQ', v, bvY[:, self.Ol[x], :])
                Fhp[self.hO[x], 2:, :] += 0.5*t*(Y @ bvx[:, self.Cl].T)
            else:                                              # K_ov: PP piece through the spare column hO_x
                x, v = A[1], A[3]
                bvx = self.bv(A)
                Ghh[self.hO[x], self.hC, :] += 0.5*t*(bvx[:, self.Vl] @ Y).T
                sd['Xp'][2:, self.hO[x]] += t*v
                sd['spare'][x] = True

    def jkey(self, P):
        kind = P[0]
        if kind in ('row', 'col'):
            return (kind, id(P[3]), P[1] if kind == 'row' else P[2], P[2] if kind == 'row' else P[1])
        if kind == 'cv':
            return ('cv', float(P[1]))
        if kind == 'vv':
            return ('vv', id(P[1]), id(P[2]))
        if kind == 'cc':
            return ('cc', id(P[1]))
        return ('unit', P[1], P[2], float(P[3]))

    def Jm_batch(self, pairs, extra_vecs):
        """Coulomb matrices of every J-pair operand (one plane sweep) and of
        the extra (name, vector) entries; unit operands use the stored J_e"""
        keys, vecs = [], []
        for pair in pairs:
            if pair[0] != 'J':
                continue
            for P in (pair[1], pair[2]):
                k = self.jkey(P)
                if k in keys or P[0] == 'unit':
                    continue
                keys.append(k); vecs.append(self.jvec(P))
        names = [nm for nm, _ in extra_vecs]
        allv = vecs + [v for _, v in extra_vecs]
        Jm = {}
        if allv:
            JHH, JHP, JPP = mrsf_grad.jblocks(np.array(allv).T, self.dims)
            for k, key in enumerate(keys + names):
                Jm[key] = self.full_from_blocks(JHH[:, :, k], JHP[:, :, k], JPP[:, :, k])
        for pair in pairs:
            if pair[0] != 'J':
                continue
            for P in (pair[1], pair[2]):
                if P[0] == 'unit':
                    Jm[self.jkey(P)] = P[3]*self.J_e[(P[1], P[2])]
        return Jm

    # ------------------------------------------------------------------
    # per-state stages
    # ------------------------------------------------------------------
    def T_densities(self, Xt):
        Ta, Tb = super().T_densities(Xt)
        sd = self._cur
        Ol = self.Ol
        Ta[np.ix_(Ol, Ol)] += sd['Ta_p']
        Tb[np.ix_(Ol, Ol)] += sd['Tb_p']
        return Ta, Tb

    def pair_patches(self, B):
        """T^a/T^b SOMO patches of Tr(B H_pair)"""
        s_ = 0.5*B[2, 2]
        d = (B[1, 2] - B[0, 2])/SQ2
        Tb_p = np.array([[B[0, 0] + s_, d], [d, B[1, 1] + s_]])
        Ta_p = np.array([[-(B[1, 1] + s_), d], [d, -(B[0, 0] + s_)]])
        return Ta_p, Tb_p

    def exchange_pass(self, sd):
        """MRSF exchange pass (unpatched), CV exchange pass (mean-field
        partner D_base), the amplitude pieces and DF products"""
        with self.tm('rhs: state (MRSF)'):
            super().exchange_pass(sd)
        Y, Xt = sd['Y'], sd['Xt']
        nC, O1, O2, H, P = self.nC, self.O1, self.O2, self.H, self.P
        Yemb = np.zeros((self.nvirb, self.nocca)); Yemb[2:, :nC] = Y
        sd['Yemb'] = Yemb
        sd['y2'] = float(np.sum(Y*Y))
        with self.tm('rhs: state (CV)'):
            sd['stcv'] = mrsf_grad.state(self.cH, Yemb, self.dims, self.dq_base)
        self._cur = sd
        if 'xs' not in sd:
            # amplitude pieces (kept across the two passes: the pair lists hold these objects)
            xs = self.slices(Xt)
            sd['xs'] = xs
            sd['Tcc'] = -Y.T @ Y
            # Fock-type densities of the parents (unit weight) and their weighted sums
            sd['N_X'] = [self.N_pieces(tb, xs, Y) for tb in self.tabs]
            sd['Nw'] = self.N_pieces(self.tab, xs, Y)
            if self.has_T4:
                tbt = {k: self.w[0]*self.tabs[0][k] - self.w[1]*self.tabs[1][k] for k in self.tab}
                sd['Nt'] = self.N_pieces(tbt, xs, Y)
        xs = sd['xs']
        sd['bvmap'] = {}
        sd['kcache'] = {}
        # DF products on the virtuals / core of the vectors used by the pairs: one sweep each
        with self.tm('rhs: bvec'):
            vcols = [Y, xs['CV'], xs['OV'][0][:, None], xs['OV'][1][:, None]]
            vkeys = [id(Y), id(xs['CV']), id(xs['OV'][0]), id(xs['OV'][1])]
            # partner columns of the K_ov classes (spare SOMO-hole columns of the CV family)
            sd['Xp_cols'] = [None, None]
            for x in (0, 1):
                v = sum(self.t2*self.tab['Kov'][y, x]*xs['OV'][y] for y in (0, 1) if self.tab['Kov'][y, x] != 0.0)
                if np.ndim(v) and _nz(v):
                    sd['Xp_cols'][x] = v
                    vcols.append(v[:, None]); vkeys.append(('Xp', x))
            ccols = [xs['CO'][0][:, None], xs['CO'][1][:, None]]
            ckeys = [id(xs['CO'][0]), id(xs['CO'][1])]
            for N in ([sd['Nw']] + ([sd['Nt']] if self.has_T4 else [])):
                for x in (0, 1):
                    if N['rows'][x] is not None:
                        vcols.append(N['rows'][x][:, None]); vkeys.append(id(N['rows'][x]))
                    if N['cols'][x] is not None:
                        ccols.append(N['cols'][x][:, None]); ckeys.append(id(N['cols'][x]))
            bv = mrsf_grad.bvec(np.column_stack(vcols), 'virt', self.nmo, self.naux)
            k0 = 0
            for key, col in zip(vkeys, vcols):
                n = col.shape[1]
                sd['bvmap'][key] = np.ascontiguousarray(bv[:, :, k0:k0 + n]) if n > 1 else np.ascontiguousarray(bv[:, :, k0])
                k0 += n
            bc = mrsf_grad.bvec(np.column_stack(ccols), 'core', self.nmo, self.naux)
            for k, key in enumerate(ckeys):
                sd['bvmap'][key] = np.ascontiguousarray(bc[:, :, k])
        with self.tm('rhs: M/jy/Tcc'):
            sd['jy'] = mrsf_grad.jq_of(None, Yemb, self.dims)[:, 0]                 # Coulomb vector of Yh^T
            sd['jTcc'] = np.einsum('ijQ,ij->Q', self.BooQ[:nC, :nC, :], sd['Tcc'])
            sd['jTpp'] = np.array(sd['stcv']['jT'])
            sd['M'] = self.M_matrix(sd)
        if 'Tb_p' in sd:
            # second pass (families recomputed): the pair-block patches of the state data
            self.apply_pair_patches(sd)

    def M_matrix(self, sd):
        """M = T^cv + sym(N_w) (MO): the density contracted with F_base"""
        nmo, Y, Cl, Vl = self.nmo, sd['Y'], self.Cl, self.Vl
        M = np.zeros((nmo, nmo))
        M[np.ix_(Vl, Vl)] = Y @ Y.T
        M[np.ix_(Cl, Cl)] = sd['Tcc']
        return M + _sym(self.N_matrix(sd['Nw']))

    def set_pair_density(self, sd, zw):
        """pair-block density B = |y|^2 w w^T + 1/2 (z_w w^T + w z_w^T), the
        T patches and the Coulomb/exchange patches of the MRSF state data"""
        w = self.w
        B = sd['y2']*np.outer(w, w)
        if zw is not None:
            B = B + 0.5*(np.outer(zw, w) + np.outer(w, zw))
        sd['B'] = B
        sd['Ta_p'], sd['Tb_p'] = self.pair_patches(B)
        self.apply_pair_patches(sd)

    def apply_pair_patches(self, sd):
        """Coulomb vector and exchange matrices of the beta density of the
        MRSF state data with the pair-block patch T^b_p added (from the
        unpatched values saved in pass 1)"""
        st, Bso, H, P = sd['st'], self.Bso, self.H, self.P
        Tb = sd['Tb_p']
        jT = np.array(sd['jT0']); KHP = np.array(sd['KbT_HP0']); KHH = np.array(sd['KbT_HH0'])
        for x in (0, 1):
            for y in (0, 1):
                if Tb[x, y] == 0.0:
                    continue
                jT += Tb[x, y]*Bso[:, self.Ol[y], x]
                KHP += Tb[x, y]*self.K_e[(x, y)][np.ix_(H, P)]
                KHH += Tb[x, y]*self.K_e[(x, y)][np.ix_(H, H)]
        st['jT'], st['KbT_HP'], st['KbT_HH'] = jT, KHP, KHH

    def V_somo(self, sd):
        """SOMO block of c (J - 1/2 K)[T^cv] + V_T"""
        Ol, Bso, BsoC = self.Ol, self.Bso, self.BsoC
        jT = sd['jTpp'] + sd['jTcc']
        V = np.zeros((2, 2))
        KHH = sd['stcv']['KbT_HH']
        for x in (0, 1):
            for y in (0, 1):
                J = float(np.dot(jT, Bso[:, Ol[y], x]))
                K = KHH[self.hO[x], self.hO[y]] + float(np.einsum('Qc,cd,Qd->', BsoC[:, :, x], sd['Tcc'], BsoC[:, :, y]))
                V[x, y] = self.cH*(J - 0.5*K)
        if 'VT' in sd:
            V += sd['VT'][np.ix_(self.hO, Ol)]
        return V

    def g_w(self, sd):
        """d omega / d w at fixed orbitals and amplitudes (unit energies of the
        parents and the Dgamma_Pi dependence of the CV correction)"""
        e = np.zeros(3)
        Fb = self.F_base
        for X in range(3):
            pairs = self.table_pairs(self.tabs[X], sd['xs'], sd['Y'])
            e[X] = sum(self.pair_energy(p) for p in pairs)
            N = sd['N_X'][X]
            Fcpl = Fb
            if X < 2:
                sign = 1.0 if X == 0 else -1.0
                Fcpl = Fb + sign*(self.J_e[(0, 0)] - self.J_e[(1, 1)] - 0.5*(self.K_e[(0, 0)] - self.K_e[(1, 1)]))
            e[X] += float(np.sum(Fcpl*self.N_matrix(N)))
        V = self.V_somo(sd)
        wG, wD, wS = self.w
        dV = V[0, 0] - V[1, 1]
        g = e + np.array([2*wG*dV - 2*SQ2*wS*V[0, 1], -2*wD*dV + 2*SQ2*wS*V[0, 1], 2*SQ2*(wD - wG)*V[0, 1]])
        return g, e

    def solve_zw(self, g):
        w, lam, U = self.w, self.lam_pair, self.U
        gp = g - np.dot(w, g)*w
        zw = np.zeros(3)
        for k in (1, 2):
            uk = U[:, k]
            zw -= uk*np.dot(uk, gp)/(lam[k] - lam[0])
        return zw

    def prepare_state(self, ist, Xt=None, extra=None, Y=None):
        with self.tm('rhs: amplitudes'):
            if Xt is None:
                Xt, Y = self.amplitudes(ist)
            elif np.shape(Xt)[1] == self.nocca + self.nC and Y is None:
                Xt, Y = self.split_amplitudes(Xt)          # stacked explicit vector
            sd = {'ist': ist, 'Xt': np.array(Xt), 'Y': np.array(Y)}
        self._cur = sd
        with self.tm('rhs: exchange pass (all)'):
            self.exchange_pass(sd)
        st = sd['st']
        sd['jT0'] = np.array(st['jT']); sd['KbT_HP0'] = np.array(st['KbT_HP']); sd['KbT_HH0'] = np.array(st['KbT_HH'])
        # kernel potentials of M, D_s, T (the G-type set at the base density)
        if self.is_dft:
            with self.tm('rhs: potential_g'):
                Bexp_g, gterms = self.g_quadratic(sd)
                want_t = self.probe_with_w
                MM, MK, MT, MK2, VE, VH = self.xcgrid.potential_g_q2(Bexp_g, gterms['M'], gterms['S'],
                                                                    gterms['T'], self.use_kernel, want_t)
                sd['VM'] = MM @ self.C                                        # v[f rho(M)] (H, all)
                if self.use_kernel:
                    nC = self.nC
                    sd['VY'] = self.C.T @ VE[:, :nC]                          # v[f rho(D_s)] C_V Y (MO)
                    sd['VC'] = self.C.T @ VH[:, self.hC]                      # v[f rho(D_s)] C_C (MO)
                    sd['K1'] = MK @ self.C                                    # v[k rho(D_s)^2] (H, all)
                if want_t:
                    sd['VT'] = MT @ self.C                                    # v[f rho(T)] (H, all)
                    sd['Kmix'] = MK2 @ self.C                                 # v[k rho(T) rho(Dg)] (H, all)
        # parent-weight response (Pi) and the pair-block density
        zw = None
        if self.space == 'pi':
            with self.tm('rhs: w-response'):
                g, e = self.g_w(sd)
                zw = self.solve_zw(g)
                sd['gw'], sd['eX'], sd['zw'] = g, e, zw
        self.set_pair_density(sd, zw)
        with self.tm('rhs: spc vectors'):
            sd['terms'] = self.spc_aux_vectors(sd['Xt'])
        with self.tm('rhs: gfock_omega (all)'):
            sd['Fo'] = self.gfock_omega(sd['Xt'], sd['st'], sd['terms'])
        if self.space == 'pi':
            Fo = sd['Fo']
            sd['q12'] = float(Fo[self.O1, self.O2] - Fo[self.O2, self.O1])
            self.pi_diag.append({'ist': ist, 'zw': np.linalg.norm(zw), 'q12': sd['q12']})
        with self.tm('rhs: zeta'):
            sd['zeta'] = self.zeta(sd['Fo'] if extra is None else sd['Fo'] + extra)
        with self.tm('rhs: gfock_zeta'):
            sd['Fzeta'] = self.gfock_zeta(sd['zeta']) if self.N else np.zeros((self.nmo, self.nmo))
        return sd

    def finish_state(self, sd, z):
        self._cur = sd
        return super().finish_state(sd, z)

    # ------------------------------------------------------------------
    # generalised Fock of omega_E
    # ------------------------------------------------------------------
    def full_matrix(self, MHH, MHP):
        M = np.zeros((self.nmo, self.nmo))
        H, P = self.H, self.P
        M[np.ix_(H, H)] = MHH
        M[np.ix_(H, P)] = MHP
        M[np.ix_(P, H)] = MHP.T
        return M

    def gfock_omega(self, Xt, st, terms):
        sd = self._cur
        fa, fb, C, H, P, V = self.fa, self.fb, self.C, self.H, self.P, self.Vidx
        nocca, nao, Ol = self.nocca, self.nao, self.Ol
        Ta, Tb = self.T_densities(Xt)
        Ta_loc = Ta[np.ix_(H, H)]
        F = 2.0*(fa @ Ta + fb @ Tb)
        xcfac = None
        if self.is_dft:
            # T^a = C_H Ta C_H^T; T^b = (C_P X~)(C_P X~)^T + C_O Tb_p C_O^T
            CH, CP = C[:, H], C[:, P]
            k = nocca + 2
            Lf = np.zeros((nao, k, 2, 1)); Rf = np.zeros((nao, k, 2, 1))
            Lf[:, :nocca, 0, 0] = 0.5*(CH @ Ta_loc); Rf[:, :nocca, 0, 0] = CH
            Lf[:, :nocca, 1, 0] = 0.5*(CP @ Xt);     Rf[:, :nocca, 1, 0] = CP @ Xt
            Lf[:, nocca:, 1, 0] = 0.5*(C[:, Ol] @ sd['Tb_p']); Rf[:, nocca:, 1, 0] = C[:, Ol]
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
        # pair lists of the state and the Coulomb matrices of every J operand (one plane sweep
        # with the pairing-term vectors)
        with self.tm('rhs: pair lists'):
            sd['pairs'] = self.table_pairs(self.tab, sd['xs'], sd['Y']) + self.model_pairs(sd)
        with self.tm('rhs: jblocks'):
            jq_spc = [v for (_, _, _, a, b) in terms for v in (a, b)]
            Jm = self.Jm_batch(sd['pairs'], [('spc%d' % k, v) for k, v in enumerate(jq_spc)])
            n0 = len(jq_spc)
            sd['Jm'] = Jm
        with self.tm('rhs: gfock_spc'):
            JHH = np.zeros((self.nocca, self.nocca, n0)); JHP = np.zeros((self.nocca, self.nvirb, n0))
            JPP = np.zeros((self.nvirb, self.nvirb, n0))
            for k in range(n0):
                Jk = Jm['spc%d' % k]
                JHH[:, :, k] = Jk[np.ix_(H, H)]; JHP[:, :, k] = Jk[np.ix_(H, P)]; JPP[:, :, k] = Jk[np.ix_(P, P)]
            F += self.gfock_spc(Xt, terms, jbl=(JHH, JHP, JPP))
        # --- extended-method pieces ---
        stcv = sd['stcv']
        F[:, H] += self.cols_H(stcv['LaH'], stcv['LaP'])             # CV exchange channel
        LbFull = np.zeros((self.nmo, self.nvirb))
        LbFull[H] = stcv['LbH']; LbFull[V] = stcv['LbP'][2:]
        F[:, P] += LbFull
        M = sd['M']
        F += 2.0*self.F_base @ M                                      # index term of Tr(F_base M)
        with self.tm('rhs: G_closed (all)'):
            G1 = self.G_closed(sd)                                    # D_base response of F_base
        F += 2.0*G1*self.n_base[None, :]
        if self.use_kernel:
            # 2 V1 S with S = Y on (V, C) + h.c.: columns C <- 2 C^T V1 C_V Y, columns V <- 2 (C^T V1 C_C) Y^T
            F[:, self.Cidx] += 2.0*sd['VY']
            F[:, V] += 2.0*(sd['VC'] @ sd['Y'].T)
            F[:, H] += sd['K1'].T*self.n_base[H][None, :]              # K1 symmetric: columns H
        if self.probe_with_w:
            # E_W = Tr(V_dg T^cv): 2 V_dg T + 2 V_T Dgamma + 2 K_mix n_base
            Tcv = np.zeros((self.nmo, self.nmo))
            Tcv[np.ix_(self.Vl, self.Vl)] = sd['Y'] @ sd['Y'].T
            Tcv[np.ix_(self.Cl, self.Cl)] = sd['Tcc']
            F += 2.0*self.V_dg @ Tcv
            VTO = sd['VT'][self.hO, :].T                               # V_T[:, O] (nmo, 2)
            F[:, Ol] += 2.0*(VTO @ self.dgam)
            F[:, H] += 2.0*sd['Kmix'].T*self.n_base[H][None, :]
        with self.tm('rhs: pairs_gfock'):
            F += self.pairs_gfock(sd['pairs'], Jm)
        return F

    def G_closed(self, sd):
        """G_base[M] = J[M] - x/2 K[M] + f_xc M on the (all, H) columns, as an
        nmo x nmo matrix with the (H, all) rows filled symmetrically"""
        H, P = self.H, self.P
        M, stcv = sd['M'], sd['stcv']
        Mhh = np.asfortranarray(M[np.ix_(H, H)])
        Mph = np.array(M[np.ix_(P, H)]); Mph[:2, :] = 0.0
        Mph = np.asfortranarray(Mph)
        x = 0.5*self.x_ref
        with self.tm('rhs: G_closed gfock'):
            Ghpa, Ghpb, Ghha, Ghhb = mrsf_grad.gfock(x, self.dims, Mhh, Mph, None, sd['jTpp'], True, True)
        GHP = Ghpa[:, :, 0] - x*stcv['KbT_HP']
        GHH = Ghha[:, :, 0] - x*stcv['KbT_HH']
        if self.is_dft:
            VM = sd['VM']
            GHH = GHH + VM[:, H]; GHP = GHP + VM[:, P]
        return self.full_matrix(GHH, GHP)

    # ------------------------------------------------------------------
    # families
    # ------------------------------------------------------------------
    def fam_bytes_state(self):
        base = super().fam_bytes_state()
        nC = self.nC
        return 2*base + 8.0*self.naux*self.nmo*(2*nC + 8) + 8.0*self.nvirb*self.nocca*self.naux

    def extra_families(self, sd):
        """CV exchange families, closed-shell mean-field pairs of M, the
        pair-block patches and every two-electron pair of the pair lists,
        added to the state's Ghh/Fhp/gpp in place; returns the CV
        particle-particle family entry (partner Y with the spare SOMO-hole
        columns of the K_ov classes)"""
        self._cur = sd
        st, stcv = sd['st'], sd['stcv']
        Y, Tb = sd['Y'], sd['Tb_p']
        nC, nV, nocca, nvirb, naux = self.nC, self.nV, self.nocca, self.nvirb, self.naux
        H, P, Ol, hO, hC = self.H, self.P, self.Ol, self.hO, self.hC
        Hloc, Bso, BsoC, BsoV = self.Hloc, self.Bso, self.BsoC, self.BsoV
        Ghh, Fhp = st['Ghh'], st['Fhp']
        x = self.x_ref
        outer = mrsf_grad.fam_outer
        # (1) CV exchange channel: hole-hole family and metric part of the CV pass
        Ghh += stcv['Ghh']
        st['gpp'] += stcv['gpp']
        Yf_cv = np.array(stcv['Yf'])                                  # -c_H/4 Y^Q + d/4 Y, partner Yemb
        Yf_cv[:, hO, :] = 0.0                                         # the partner's SOMO-hole columns are free
        sd['Xp'] = np.array(sd['Yemb'])
        sd['spare'] = [False, False]
        sd['pp_mats'], sd['pp_vecs'] = [], []
        # (2) pair-block patch of T^b: J pair with the reference density and the beta
        # exchange pair with the core
        for a in (0, 1):
            for b in (0, 1):
                if Tb[a, b] == 0.0:
                    continue
                Ghh[hO[a], hO[b], :] += 0.5*Tb[a, b]*self.dq_ref
                v = 0.5*x*Tb[a, b]*BsoC[:, :, b].T                      # (nC, naux)
                Ghh[hC, hO[a], :] -= v
                Ghh[hO[a], hC, :] -= v
        # (3) closed-shell mean-field pairs of M with D_base: J pair and K pair (x/2, weights n/2)
        M = sd['M']
        Nw = sd['Nw']
        Mhh = M[np.ix_(H, H)]
        MHV = M[np.ix_(H, self.Vidx)]
        pM = sd['jTpp'] + sd['jTcc'] + Nw['a']*sd['jy']
        for xx in (0, 1):
            if Nw['rows'][xx] is not None:
                pM = pM + BsoV[:, :, xx] @ Nw['rows'][xx]
            if Nw['cols'][xx] is not None:
                pM = pM + BsoC[:, :, xx] @ Nw['cols'][xx]
        occ_H = self.n_base[H]
        outer(Ghh, 0, nocca, 0, nocca, [0.5*Mhh, 0.5*np.diag(occ_H)], [self.dq_base, pM])
        outer(Fhp, 0, nocca, 2, nvirb, [MHV], [self.dq_base])
        lam = self.lam
        BooQ2 = self.BooQ.reshape(nocca, nocca*naux, order='F')
        MB = (Mhh @ BooQ2).reshape(nocca, nocca, naux, order='F')     # (M_HH B^Q_HH)
        bvY = self.bv_of(Y)                                           # (naux, nmo, nC): sum_b B^Q_qb Y_bj
        # cross terms from the H x V blocks of M: sum_b M[p,b] B^Q_{bq}
        if Nw['a'] != 0.0:
            MB[hC, :, :] += 0.5*Nw['a']*np.transpose(bvY[:, H, :], (2, 1, 0))
        for xx in (0, 1):
            if Nw['rows'][xx] is not None:
                MB[hO[xx], :, :] += 0.5*self.bv(('row', xx, 'V', Nw['rows'][xx]))[:, H].T
        Ghh -= 0.5*x*(MB*lam[None, :, None] + np.transpose(MB, (1, 0, 2))*lam[:, None, None])
        # particle rows: (M B^Q Lambda)_{aq}, a in V
        MaH = (M[np.ix_(self.Vidx, H)] @ BooQ2).reshape(nV, nocca, naux, order='F')
        Bt = np.ascontiguousarray(np.transpose(bvY[:, H, :], (2, 1, 0))).reshape(nC, nocca*naux)
        MaV = (Y @ Bt).reshape(nV, nocca, naux)
        Fq = MaH + MaV                                                # (nV, nocca, naux)
        Fhp[:, 2:, :] -= x*lam[:, None, None]*np.transpose(Fq, (1, 0, 2))
        # (4) the two-electron pairs
        self.pairs_families(sd['pairs'], Ghh, Fhp)
        if sd['pp_mats']:
            outer(Yf_cv, 0, nvirb, 0, nocca, sd['pp_mats'], sd['pp_vecs'])
        # K_ov pieces: the spare SOMO-hole columns of the CV partner and their metric part
        Xp = sd['Xp']
        for xx in (0, 1):
            if not sd['spare'][xx]:
                continue
            Yf_cv[2:, hO[xx], :] += 0.25*(Y @ BsoC[:, :, xx].T)
            if ('Xp', xx) in sd['bvmap'] and np.abs(sd['Xp_cols'][xx] - Xp[2:, hO[xx]]).max() < 1e-12:
                bvp = sd['bvmap'][('Xp', xx)]
            else:
                bvp = mrsf_grad.bvec(Xp[:, hO[xx]][2:], 'virt', self.nmo, naux)[:, :, 0]
            st['gpp'] += 2.0*(bvp[:, P] @ Yf_cv[:, hO[xx], :])
        sd['Yf_cv'] = Yf_cv
        return [{'Xt': Xp, 'Yf': Yf_cv, 'g': np.zeros((naux, naux)), 'pp_only': True}]

    # ------------------------------------------------------------------
    # probe factors and the one-electron density
    # ------------------------------------------------------------------
    def P_quadratic(self, Xt, Za, Zb, pzeta=None):
        """relaxed densities with the pair-block patch of T^b on the hole-MO
        block of the beta form"""
        Bexp, terms = super().P_quadratic(Xt, Za, Zb, pzeta)
        nbb, idxb, Sb = terms['b']
        Tb = self._cur['Tb_p']
        for a in (0, 1):
            for b in (0, 1):
                Sb[self.hO[a], self.hO[b]] += Tb[a, b]
        return Bexp, terms

    def probe_g_factors(self, sd):
        if not self.is_dft:
            return None
        return self.g_quadratic(sd)

    def g_quadratic(self, sd):
        """quadratic forms (basis [C_H | C_V Y | C_V r_1 | C_V r_2]) of the
        G-type terms: M = T^cv + sym(N_w), the CV trial density D_s = C_V Y
        C_C^T + h.c. and T^cv (Pi); computed once per state"""
        if 'gq' in sd:
            return sd['gq']
        Y, nC, h = sd['Y'], self.nC, self.nocca
        CV = self.C[:, self.Vidx]
        Nw = sd['Nw']
        rows = [xx for xx in (0, 1) if Nw['rows'][xx] is not None]
        Bexp = np.zeros((self.nao, nC + len(rows)))
        Bexp[:, :nC] = CV @ Y
        for k, xx in enumerate(rows):
            Bexp[:, nC + k] = CV @ Nw['rows'][xx]
        hC = np.array(self.hC)
        M = sd['M']
        nbM = h + nC + len(rows)
        SM = np.zeros((nbM, nbM))
        SM[:h, :h] = M[np.ix_(self.H, self.H)]
        SM[h:h + nC, h:h + nC] = np.eye(nC)                               # T^cv_VV = Y Y^T
        for i in range(nC):                                               # M[C_i, V] = a/2 Y^T_i
            SM[hC[i], h + i] = 0.5*Nw['a']; SM[h + i, hC[i]] = 0.5*Nw['a']
        for k, xx in enumerate(rows):                                     # M[O_x, V] = r_x/2
            SM[self.hO[xx], h + nC + k] = 0.5; SM[h + nC + k, self.hO[xx]] = 0.5
        idxM = np.concatenate([np.arange(1, h + 1), h + 1 + np.arange(nC + len(rows))])
        nbS = 2*nC
        SS = np.zeros((nbS, nbS))
        SS[:nC, nC:] = np.eye(nC); SS[nC:, :nC] = np.eye(nC)
        idxS = np.concatenate([hC + 1, h + 1 + np.arange(nC)])
        ST = np.zeros((nbS, nbS))
        ST[:nC, :nC] = sd['Tcc']; ST[nC:, nC:] = np.eye(nC)
        gq = {'M': (nbM, idxM, SM), 'S': (nbS, idxS, SS)}
        if self.probe_with_w:
            gq['T'] = (nbS, idxS, ST)
        else:
            gq['T'] = (0, np.zeros(0, int), np.zeros((0, 0)))
        sd['gq'] = (Bexp, gq)
        return sd['gq']

    def extra_density(self, sd):
        return sd['M']

    def run(self, states, **kw):
        self.pi_diag = []
        res = super().run(states, **kw)
        if self.pi is not None:
            res['pi'] = {'w': np.array(self.w), 'pair': np.array(self.lam_pair), 'gap': self.pair_gap,
                         'zw_norm': np.array([d['zw'] for d in self.pi_diag]),
                         'q12': np.array([d['q12'] for d in self.pi_diag])}
        return res

    def gradients_and_nacs(self, states, pairs, **kw):
        """gradients and derivative couplings (see GradientDriver): the
        per-vector diagnostics of the parent Pi (|z_w|, Q_12) reported for
        the state vectors only (Q_12 of a pair vector is not zero: its
        extra right-hand side is not invariant)"""
        from graci.interfaces.mrsf.mrsf_gradient import nac_vectors
        gres, out = super().gradients_and_nacs(states, pairs, **kw)
        if gres.get('pi', None) is not None:
            need, pos = nac_vectors(states, pairs)
            sel = [pos[s] for s in states]
            pi = dict(gres['pi'])
            pi['zw_norm'] = np.asarray(pi['zw_norm'])[sel]
            pi['q12'] = np.asarray(pi['q12'])[sel]
            gres['pi'] = pi
        return gres, out

    def finalise(self):
        if self.is_dft:
            self.xcgrid.clear_hole_weights()
            self.xcgrid.clear_dgamma()
        super().finalise()
