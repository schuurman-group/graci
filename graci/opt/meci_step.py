"""
Constrained step of the minimum-energy conical intersection (MECI)
optimiser ($meciopt). Pure numpy: no GRaCI objects, so the step can be
tested on model surfaces.

Notation (Cartesian, Bohr and Hartree, vectors of length 3N): states
I < J in energy order, Ebar = (E_I + E_J)/2 with gradient gbar,
dE = E_J - E_I with gradient g, and h the interstate coupling vector
<I|dH/dR|J> (finite at the intersection, sign arbitrary). Translations
and rotations are projected out of every vector.

The Lagrangian L = Ebar - lam^T r is minimised subject to two
constraints: r1 = dE (target 0) and r2, whose value is 0 at every
geometry and whose gradient is h. The scheme is the reduced-restricted
quasi-Newton method of Anglada and Bofill (J. Comput. Chem. 18, 992
(1997)) as used in OpenMolcas (De Vico, Olivucci and Lindh, JCTC 1,
1029 (2005)):

  - lam: least-squares multipliers, lam = (R^T R)^-1 R^T gbar, R = [g h]
  - T = [Tb, Tti]: Tb spans the branching plane {g, h}, Tti its
    complement in the internal (translation/rotation-free) space
  - gap step: dxb = Tb dy, dy = -(R^T Tb)^-1 r, i.e. the Newton step of
    the linearised constraints (g.dxb = -dE, h.dxb = 0)
  - seam step: restricted-step RFO on the reduced gradient
    Tti^T (gbar + W dxb) and the reduced Hessian Tti^T W Tti
  - W: model Hessian at the current geometry plus BFGS updates with the
    Lagrangian gradients gbar_k - R_k lam_k of the last few steps, each
    with its own multipliers (gbar projected off its own branching
    plane). Close to the intersection the adiabatic g and h rotate
    within the plane from step to step and lam with them, so products
    R_k lam with the current lam would be inconsistent; the projected
    gradients depend only on the plane.

Steps are accepted and the trust radius is adapted with a merit
function chosen by the kind of step. When the gap step dominates, the
exact penalty (L1 merit) function phi = Ebar + mu dE, mu > |lam_1|: far
from the seam the gap step raises Ebar while it lowers dE, which a test
on the Lagrangian alone would reject. When the seam step dominates, the
Lagrangian L = Ebar - lam_1 dE of the current multipliers: along a
curved seam the gap re-opens to second order, which the linear gap term
of the L1 merit does not predict (the Maratos effect: good steps get
ratios well below 1, the trust radius stalls and longer steps are
rejected), while the quasi-Newton Lagrangian Hessian contains that
curvature. The gap enters the convergence test explicitly.
"""

import time
import numpy as np

# relative rank tolerance of the constraint gradients
RANK_TOL = 1.e-6
# norm below which a constraint gradient is treated as null (Eh/Bohr)
NULL_GRAD = 1.e-10
# |predicted change| below which the trust radius is not updated (Eh)
PRED_FLOOR = 1.e-6
# a step is rejected if the merit function rises by more than this (Eh)
# against the prediction
REJECT_RISE = 1.e-4
# penalty weight of the merit function: mu >= MU_SCALE |lam_1| + MU_MIN,
# never decreased
MU_SCALE = 1.25
MU_MIN = 0.05
# maximum number of consecutive rejected steps / failed evaluations
MAX_REJECT = 3
MAX_FAIL = 2


def tr_basis(x, tol=1.e-6):
    """
    orthonormal basis (3N, ntr) of the infinitesimal translations and
    rotations (about the centroid) of the geometry x (3N,); ntr = 6, 5
    for a linear geometry
    """
    X = np.asarray(x, dtype=float).reshape(-1, 3)
    n = X.shape[0]
    c = X - X.mean(axis=0)
    vecs = []
    for a in range(3):
        t = np.zeros((n, 3))
        t[:, a] = 1.0
        vecs.append(t.ravel())
    for a in range(3):
        e = np.zeros(3)
        e[a] = 1.0
        vecs.append(np.cross(e, c).ravel())
    U, s, _ = np.linalg.svd(np.array(vecs).T, full_matrices=False)
    return U[:, s > tol*s[0]]


def complement(B, n):
    """orthonormal basis of the orthogonal complement of the
    orthonormal columns of B (n, m) in R^n"""
    if B.shape[1] == 0:
        return np.eye(n)
    w, V = np.linalg.eigh(np.eye(n) - B @ B.T)
    return V[:, w > 0.5]


def rsrfo(H, g, trust, tol=1.e-12, maxit=200):
    """
    restricted-step rational-function step for the quadratic model
    g.s + 1/2 s.H.s (Besalu and Bofill, TCA 100, 265 (1998)): the RFO
    step s = -(H - mu0)^-1 g, mu0 the lowest eigenvalue of the augmented
    Hessian [[H, g], [g^T, 0]], if |s| <= trust; otherwise the shifted
    step with |s| = trust and mu < mu0, which is the exact trust-region
    solution. Includes the hard case (no gradient component along the
    lowest, negative-curvature mode). Returns (s, predicted change,
    on_boundary).
    """
    n = g.size
    if n == 0:
        return np.zeros(0), 0.0, False
    b, U = np.linalg.eigh(0.5*(H + H.T))
    gh = U.T @ g
    gnorm = np.linalg.norm(gh)
    nz = np.abs(gh) > 1.e-12*max(gnorm, 1.e-300)

    def shifted(mu):
        st = np.zeros(n)
        st[nz] = -gh[nz]/(b[nz] - mu)
        return st

    def result(st, boundary):
        pred = float(gh @ st + 0.5*np.sum(b*st**2))
        return U @ st, pred, boundary

    if not nz.any():
        # no gradient: a step only along negative curvature
        if b[0] >= 0.:
            return np.zeros(n), 0.0, False
        st = np.zeros(n)
        st[0] = trust
        return result(st, True)

    # lowest eigenvalue of the augmented Hessian of the coupled modes
    bn, gn = b[nz], gh[nz]
    m = bn.size
    A = np.zeros((m + 1, m + 1))
    A[:m, :m] = np.diag(bn)
    A[:m, m] = gn
    A[m, :m] = gn
    mu0 = np.linalg.eigvalsh(A)[0]

    # admissible shifts lie below every eigenvalue of H - mu, including
    # those of modes without a gradient component
    hard = (not nz[0]) and b[0] < mu0
    mu_hi = b[0] if hard else mu0
    s_hi = shifted(mu_hi)
    if np.linalg.norm(s_hi) <= trust:
        if not hard:
            return result(s_hi, False)
        # hard case: fill the trust radius along the lowest mode
        s_hi[0] = np.sqrt(max(trust**2 - np.sum(s_hi**2), 0.))
        return result(s_hi, True)

    # |s(mu)| = trust for mu in (mu_lo, mu_hi): Newton on 1/|s| - 1/trust,
    # safeguarded by bisection (|s| increases monotonically with mu)
    mu_lo = min(bn[0], b[0]) - gnorm/trust
    mu = mu_lo
    for it in range(maxit):
        st = shifted(mu)
        sn = np.linalg.norm(st)
        phi = 1.0/sn - 1.0/trust
        if abs(sn - trust) <= tol*trust:
            break
        if phi > 0.:
            mu_lo = mu
        else:
            mu_hi = mu
        dsn = np.sum(gn**2/(bn - mu)**3)/sn
        mu_new = mu + phi*sn**2/dsn
        if not (mu_lo < mu_new < mu_hi):
            mu_new = 0.5*(mu_lo + mu_hi)
        mu = mu_new
    return result(shifted(mu), True)


def bfgs(W, s, y, damp=0.2):
    """BFGS update of W with the step s and gradient change y, with
    Powell damping (s.y >= damp s.W.s) so that W stays positive
    definite on the span of the steps"""
    if s @ s < 1.e-12:
        return W
    Ws = W @ s
    sWs = s @ Ws
    if sWs <= 1.e-12*(s @ s):
        return W
    sy = s @ y
    if sy < damp*sWs:
        theta = (1.0 - damp)*sWs/(sWs - sy)
        y = theta*y + (1.0 - theta)*Ws
        sy = s @ y
    return W - np.outer(Ws, Ws)/sWs + np.outer(y, y)/sy


class SeamPoint:
    """
    one geometry of the MECI search: energies, gradients and coupling of
    the states I < J and the derived branching-plane quantities. h_ref:
    the coupling of the previous point, to which the sign of h is
    aligned; tr = False for model spaces without rigid-body motion
    """
    def __init__(self, x, e, grads, h, h_ref=None, tr=True, info=None):
        self.x = np.asarray(x, dtype=float).ravel()
        n = self.x.size
        self.info = {} if info is None else info
        self.P = tr_basis(self.x) if tr else np.zeros((n, 0))

        gI = self.proj(np.asarray(grads[0], dtype=float).ravel())
        gJ = self.proj(np.asarray(grads[1], dtype=float).ravel())
        hv = self.proj(np.asarray(h, dtype=float).ravel())
        if h_ref is not None and hv @ h_ref < 0.:
            hv = -hv

        self.e = (float(e[0]), float(e[1]))
        self.ebar = 0.5*(self.e[0] + self.e[1])
        self.de = self.e[1] - self.e[0]
        self.gbar = 0.5*(gI + gJ)
        self.g = gJ - gI
        self.h = hv

        # active constraints: the gap if its gradient is not null, the
        # coupling if it is not null and not parallel to g
        cols, vals, act = [], [], []
        gn = np.linalg.norm(self.g)
        if gn > NULL_GRAD:
            cols.append(self.g)
            vals.append(self.de)
            act.append(0)
        hn = np.linalg.norm(hv)
        if hn > NULL_GRAD:
            hperp = hv if gn <= NULL_GRAD else hv - (hv @ self.g)/gn**2*self.g
            if np.linalg.norm(hperp) > RANK_TOL*hn:
                cols.append(hv)
                vals.append(0.)
                act.append(1)
        self.active = act
        self.R = np.array(cols).T if cols else np.zeros((n, 0))
        self.r = np.array(vals)

        # least-squares multipliers (0 for inactive constraints) and the
        # gradient of the Lagrangian
        self.lam = np.zeros(2)
        if act:
            self.lam[act] = np.linalg.lstsq(self.R, self.gbar, rcond=None)[0]
        self.gl = self.gbar - self.R2 @ self.lam

        # branching plane and seam space
        if act:
            U, _, _ = np.linalg.svd(self.R, full_matrices=False)
            self.Tb = U
        else:
            self.Tb = np.zeros((n, 0))
        self.Tti = complement(np.hstack([self.P, self.Tb]), n)
        self.gseam = self.Tti @ (self.Tti.T @ self.gbar)

    @property
    def R2(self):
        """[g h], including inactive columns"""
        return np.column_stack([self.g, self.h])

    def proj(self, v):
        """v without its translation/rotation components"""
        return v - self.P @ (self.P.T @ v)


class Step:
    """a proposed step: dx, its gap part dxb, the predicted change of
    the merit function, whether it was truncated, the trust radius, the
    penalty weight mu and the merit function ('lagrangian' when the seam
    part of the step dominates, else 'l1')"""
    def __init__(self, dx, dxb, pred, trunc, trust, mu, merit):
        self.dx = dx
        self.dxb = dxb
        self.pred = pred
        self.trunc = trunc
        self.trust = trust
        self.mu = mu
        self.merit = merit


def characterise(pt):
    """
    first-order characterisation of the intersection at pt (Fdez. Galvan
    et al., JCTC 12, 3636 (2016)): the half gap gradient g/2 and h are
    rotated into an orthogonal pair x, y (|x| >= |y|); pitch
    d = sqrt((x^2 + y^2)/2), asymmetry D = (x^2 - y^2)/(x^2 + y^2), tilt
    s = gbar in the plane (sx, sy >= 0 fix the signs of x and y),
    relative tilt sigma = |s|/d, heading theta = atan2(sy, sx);
    P < 1 peaked, > 1 sloped; B < 1 bifurcating, > 1 single path
    """
    gt = 0.5*pt.g
    h = pt.h
    gg, hh, gh = gt @ gt, h @ h, gt @ h
    beta = 0.5*np.arctan2(2.0*gh, gg - hh)
    x = np.cos(beta)*gt + np.sin(beta)*h
    y = np.cos(beta)*h - np.sin(beta)*gt
    xx, yy = x @ x, y @ y
    if yy > xx:
        x, y, xx, yy = y, x, yy, xx
    xn, yn = np.sqrt(xx), np.sqrt(yy)
    ux = x/xn if xn > 0. else x
    uy = y/yn if yn > 0. else y
    sx, sy = pt.gbar @ ux, pt.gbar @ uy
    if sx < 0.:
        sx, ux = -sx, -ux
    if sy < 0.:
        sy, uy = -sy, -uy
    pitch = np.sqrt(0.5*(xx + yy))
    asym = (xx - yy)/(xx + yy) if xx + yy > 0. else 0.
    tilt = np.sqrt(sx**2 + sy**2)
    sigma = tilt/pitch if pitch > 0. else 0.
    theta = np.arctan2(sy, sx)
    peaked = sigma**2/(1.0 - asym**2)*(1.0 - asym*np.cos(2.0*theta)) \
        if abs(asym) < 1. else np.inf
    if asym > 0.:
        bif = (sigma/(2.0*asym))**(2./3.) * \
            (((1.0 + asym)*np.cos(theta)**2)**(1./3.) +
             ((1.0 - asym)*np.sin(theta)**2)**(1./3.))
    else:
        bif = np.inf
    return {'pitch': pitch, 'asymmetry': asym, 'tilt': tilt,
            'rel_tilt': sigma, 'heading': theta, 'P': peaked, 'B': bif,
            'x': ux, 'y': uy}


class MeciSolver:
    """
    the MECI iterations on top of an evaluator x -> (x, (E_I, E_J),
    (g_I, g_J), h, info) or None on failure. hess_guess: x -> model
    Hessian (3N, 3N); conv: thresholds gap, grad_max, grad_rms,
    step_max, step_rms, energy
    """
    def __init__(self, hess_guess, conv, trust=0.3, trust_min=0.01,
                 trust_max=0.5, hess_update='bfgs', hess_window=5,
                 max_iter=50, eig_floor=5.e-3, gap_cap=0.866, tr=True):
        self.hess_guess = hess_guess
        self.conv = conv
        self.trust = trust
        self.trust_min = trust_min
        self.trust_max = trust_max
        self.hess_update = hess_update
        self.hess_window = hess_window
        self.max_iter = max_iter
        self.eig_floor = eig_floor
        self.gap_cap = gap_cap
        self.tr = tr
        # secant pairs (older point, newer point) of the BFGS updates
        self.pairs = []
        # penalty weight of the merit function
        self.mu = 0.
        # wall time spent in the optimiser itself (s)
        self.wall = 0.

    def point(self, raw, ref=None):
        """SeamPoint from the evaluator output, h aligned to ref"""
        x, e, grads, h, info = raw
        return SeamPoint(x, e, grads, h,
                         h_ref=None if ref is None else ref.h,
                         tr=self.tr, info=info)

    def hessian(self, pt):
        """model Hessian at pt (eigenvalues in the internal space floored
        at eig_floor) plus the BFGS updates of the stored pairs with the
        Lagrangian gradients of the points' own multipliers"""
        n = pt.x.size
        H = np.asarray(self.hess_guess(pt.x), dtype=float)
        C = complement(pt.P, n)
        w, V = np.linalg.eigh(C.T @ (0.5*(H + H.T)) @ C)
        V = C @ V
        W = (V*np.maximum(w, self.eig_floor)) @ V.T
        if self.hess_update == 'bfgs':
            for a, b in self.pairs[-self.hess_window:]:
                s = pt.proj(b.x - a.x)
                y = pt.proj(b.gl - a.gl)
                W = bfgs(W, s, y)
        return W

    def step(self, pt):
        """the step from pt within the current trust radius"""
        trust = self.trust
        W = self.hessian(pt)

        # Newton step of the linearised constraints in the branching
        # plane; it has priority and is capped at gap_cap*trust
        trunc = False
        dxb = np.zeros(pt.x.size)
        if pt.active:
            dy = -np.linalg.solve(pt.R.T @ pt.Tb, pt.r)
            dxb = pt.Tb @ dy
            nb = np.linalg.norm(dxb)
            if nb > self.gap_cap*trust:
                dxb *= self.gap_cap*trust/nb
                trunc = True

        # RS-RFO step in the seam space, total step within the trust radius
        tx = np.sqrt(max(trust**2 - dxb @ dxb, 0.))
        gx = pt.Tti.T @ (pt.gbar + W @ dxb)
        Wx = pt.Tti.T @ W @ pt.Tti
        dz, _, boundary = rsrfo(Wx, gx, tx)
        dx = dxb + pt.Tti @ dz

        # predicted change of the merit function: the Lagrangian when the
        # seam step dominates, else Ebar + mu dE (the coupling constraint
        # is satisfied to first order by construction)
        self.mu = max(self.mu, MU_SCALE*abs(pt.lam[0]) + MU_MIN)
        if np.linalg.norm(dxb) < np.linalg.norm(dz):
            merit = 'lagrangian'
            pred = float(pt.gl @ dx + 0.5*dx @ W @ dx)
        else:
            merit = 'l1'
            pred = float(pt.gbar @ dx + 0.5*dx @ W @ dx +
                         self.mu*(abs(pt.de + pt.g @ dx) - abs(pt.de)))
        return Step(dx, dxb, pred, trunc or boundary, trust, self.mu, merit)

    def converged(self, pt, stp, ebar_prev):
        """convergence flags at pt with the proposed step stp"""
        c = self.conv
        n = pt.x.size
        gs = pt.gseam
        dx = stp.dx
        flags = {
            'gap': abs(pt.de) <= c['gap'],
            'grad_max': np.max(np.abs(gs)) <= c['grad_max'],
            'grad_rms': np.linalg.norm(gs)/np.sqrt(n) <= c['grad_rms'],
            'step_max': np.max(np.abs(dx)) <= c['step_max'] and not stp.trunc,
            'step_rms': np.linalg.norm(dx)/np.sqrt(n) <= c['step_rms'] and not stp.trunc,
            'energy': ebar_prev is not None and abs(pt.ebar - ebar_prev) <= c['energy']}
        conv = flags['gap'] and flags['grad_max'] and flags['grad_rms'] and \
            ((flags['step_max'] and flags['step_rms']) or flags['energy'])
        return conv, flags

    def add_pair(self, a, b):
        """store the secant pair (a, b); the window is cleared when the
        branching plane itself has turned (smallest principal cosine
        between the planes below 0.5). Rotations of g and h within the
        plane, which are large close to the intersection, do not matter
        because the secant pairs use the Lagrangian gradients of each
        point's own multipliers (gbar projected off its own plane)"""
        if a.Tb.shape[1] == b.Tb.shape[1] and a.Tb.shape[1] > 0:
            cos = np.linalg.svd(a.Tb.T @ b.Tb, compute_uv=False)
            if cos[-1] < 0.5:
                self.pairs = []
                return
        self.pairs.append((a, b))
        self.pairs = self.pairs[-self.hess_window:]

    def assess(self, cur, new, stp):
        """accept or reject the trial point new reached from cur with stp
        (merit function of the step), and update the trust radius;
        returns (accepted, rho)"""
        if stp.merit == 'lagrangian':
            dphi = (new.ebar - cur.ebar) - cur.lam[0]*(new.de - cur.de)
        else:
            dphi = (new.ebar - cur.ebar) + stp.mu*(abs(new.de) - abs(cur.de))
        rho = dphi/stp.pred if abs(stp.pred) > PRED_FLOOR else None
        self.add_pair(cur, new)
        dxn = np.linalg.norm(stp.dx)
        if (rho is None or rho < 0.) and dphi > REJECT_RISE:
            self.trust = max(self.trust_min, 0.25*dxn)
            return False, rho
        if rho is not None:
            if rho < 0.25:
                self.trust = max(self.trust_min, 0.5*min(self.trust, dxn))
            elif rho > 0.75 and dxn >= 0.8*self.trust:
                self.trust = min(self.trust_max, 1.5*self.trust)
        return True, rho

    def run(self, evaluate, raw0, callback=None, on_accept=None,
            on_reject=None):
        """
        iterate from the evaluator output raw0 at the start geometry.
        callback(it, pt, stp, flags, rho, nrej) is called once per
        iteration (before the next evaluation); on_accept/on_reject are
        called for every trial point. Returns a dict with the status
        ('converged', 'max_iter', 'rejected', 'failed'), the last point
        and the history
        """
        t0 = time.perf_counter()
        cur = self.point(raw0)
        hist = []
        ebar_prev = None
        rho = None
        nrej = 0
        status = 'max_iter'
        it = 0
        while True:
            stp = self.step(cur)
            conv, flags = self.converged(cur, stp, ebar_prev)
            self.wall += time.perf_counter() - t0
            hist.append({'x': cur.x.copy(), 'e': cur.e, 'lam': cur.lam.copy(),
                         'trust': stp.trust, 'rho': rho,
                         'grad_max': float(np.max(np.abs(cur.gseam))),
                         'grad_rms': float(np.linalg.norm(cur.gseam)/np.sqrt(cur.x.size)),
                         'step_max': float(np.max(np.abs(stp.dx))),
                         'step_rms': float(np.linalg.norm(stp.dx)/np.sqrt(cur.x.size)),
                         'info': cur.info})
            if callback is not None:
                callback(it, cur, stp, flags, rho, nrej)
            t0 = time.perf_counter()
            if conv:
                status = 'converged'
                break
            if it >= self.max_iter:
                break

            # trial points until one is accepted
            nrej = 0
            nfail = 0
            new = None
            while True:
                self.wall += time.perf_counter() - t0
                raw = evaluate(cur.x + stp.dx)
                t0 = time.perf_counter()
                if raw is None:
                    nfail += 1
                    if nfail > MAX_FAIL:
                        status = 'failed'
                        break
                    self.trust = max(self.trust_min, 0.5*np.linalg.norm(stp.dx))
                    stp = self.step(cur)
                    continue
                new = self.point(raw, ref=cur)
                ok, rho = self.assess(cur, new, stp)
                if ok:
                    if on_accept is not None:
                        on_accept(new)
                    break
                if on_reject is not None:
                    on_reject(new)
                nrej += 1
                if nrej > MAX_REJECT:
                    status = 'rejected'
                    break
                stp = self.step(cur)
            if status in ('failed', 'rejected'):
                break
            ebar_prev = cur.ebar
            cur = new
            it += 1

        return {'status': status, 'point': cur, 'niter': it, 'history': hist}
