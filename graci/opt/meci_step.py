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

Far from the seam this composite step is a poor guide: the gap step
moves along g whatever that costs. When g is a stiff coordinate along
which both states rise (lam_1 < -1/2) the seam is reached far above the
MECI, although a soft coordinate may close the gap at second order. With
approach = 'penalty' the step is there one restricted-step RFO step on
the exact penalty function phi = Ebar + mu dE in the space without
rigid-body motion and without h (the Sl1QP step of Fletcher: the
minimiser of the model of phi within the trust radius):

  - Hessian of phi: the model Hessian plus Bofill updates with the
    gradients gbar + mu g of the last few steps. They keep the negative
    curvature of the gap, which the damped BFGS update of W discards
  - mu starts at MU_START and is doubled until the step closes at least
    the fraction FEAS_FRAC of what a pure gap step would close in the
    linear model, so the seam comes closer at every step and strongly
    sloped intersections (|lam_1| > mu) are reached
  - as soon as the step would cross the linearised seam, the composite
    step takes over

The gap step closes the gap to first order only. On a curved seam the
seam step dz leaves the seam by a drift of about 1/2 |dz|^2 kappa in the
branching plane, so every long seam step re-opens the gap and the next
step has to close it again. With gap_correction the drift of each
composite step is measured at the next point (the Newton gap step there,
plus the correction that was applied) and stored as a curvature kappa
along the direction of its seam step; the next composite step adds the
predicted drift of its own seam step, -1/2 (dz.z_prev)^2 kappa, to the
gap step (a second-order correction, no extra evaluations).

algorithm = 'projected' is the two-phase projected-gradient algorithm
of Sicilia, Blancafort, Bearpark and Robb (J. Chem. Theory Comput. 4,
257 (2008)), which builds on Bearpark, Robb and Schlegel (Chem. Phys.
Lett. 223, 269 (1994)):

  - far phase, dE above gap_switch: a quasi-Newton step on the composite
    gradient 2 dE g/|g| + P gbar (P projects on the seam space), with a
    BFGS Hessian updated with that gradient. It is not the gradient of a
    function, so there is no merit function: every step is accepted and
    the trust radius follows the norm of the composite gradient
  - near phase: the Newton step of the linearised constraints (the
    published step along g alone misses the tip of a cone whose g and h
    are not orthogonal) plus the RS-RFO step in the seam space on the
    projected gradient, without the coupling W dxb to the gap step;
    accepted and limited like the composite step

The steps are built in a working space given by a coordinate object:
Cartesian coordinates without rigid-body motion, or redundant internal
coordinates (internal_coords.py). Every point carries a twin with its
vectors in that space (the point itself in Cartesian coordinates);
energies, the seam gradient of the convergence test, the step
statistics and the characterisation of the intersection are Cartesian
in both cases.

Steps are accepted and the trust radius is adapted with a merit
function chosen by the kind of step. For approach steps and when the gap
step dominates a composite step, the exact penalty (L1 merit) function
phi = Ebar + mu dE, mu > |lam_1|: far from the seam the gap step raises
Ebar while it lowers dE, which a test on the Lagrangian alone would
reject. When the seam step dominates, the Lagrangian L = Ebar - lam_1 dE
of the current multipliers: along a curved seam the gap re-opens to
second order, which the linear gap term of the L1 merit does not predict
(the Maratos effect: good steps get ratios well below 1, the trust
radius stalls and longer steps are rejected), while the quasi-Newton
Lagrangian Hessian contains that curvature. The gap enters the
convergence test explicitly.
"""

import time
import numpy as np

# relative rank tolerance of the constraint gradients
RANK_TOL = 1.e-6
# norm below which a constraint gradient is treated as null (Eh/Bohr)
NULL_GRAD = 1.e-10
# |predicted change| below which the trust radius is not updated (Eh):
# the larger of PRED_FLOOR and NOISE_FAC times the energy threshold
# (below it the ratio of the merit function is noise)
PRED_FLOOR = 1.e-6
NOISE_FAC = 10.0
# a step is rejected if the merit function rises by more than this (Eh)
# against the prediction
REJECT_RISE = 1.e-4
# penalty weight of the merit function: mu >= MU_SCALE |lam_1| + MU_MIN
# in composite steps, never decreased
MU_SCALE = 1.25
MU_MIN = 0.05
# approach steps: start value of mu, its growth factor and the maximum
# number of passes of the steering loop, and the fraction of the gap
# closure of a pure gap step that the step must reach in the linear model
MU_START = 1.0
MU_STEER = 2.0
MAX_STEER = 12
FEAS_FRAC = 0.1
# a secant pair enters the Hessian of the penalty function only if the
# gap gradients of its two points agree in direction (cosine above this)
PAIR_COS = 0.5
# trust radius: growth factor, and the fraction of the step length it is
# set to after the first and after further consecutive rejected steps
TRUST_GROW = 2.0
REJECT_SHRINK = (0.5, 0.25)
# far phase of the projected algorithm: the radius is set to FAR_SHRINK
# times the step length when the composite gradient grows, and
# multiplied by FAR_GROW when it falls below FAR_DROP times its old norm
# on a step at the boundary (FAR_DROP = 0.75 held the radius down for
# eleven steps of an ethylene run in which the gradient fell by 10 % a
# step)
FAR_SHRINK = 0.5
FAR_GROW = 1.5
FAR_DROP = 1.0
# second-order gap correction: the curvature is measured only on seam
# steps longer than SOC_MIN_STEP (Bohr), and the correction is at most
# SOC_MAX_FRAC times the length of the seam step
SOC_MIN_STEP = 0.02
SOC_MAX_FRAC = 0.25
# maximum number of consecutive rejected steps / failed evaluations /
# shortened steps after a failed back-transformation
MAX_REJECT = 3
MAX_FAIL = 2
MAX_BACK = 4


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


def bofill(W, s, y):
    """Bofill update of W with the step s and gradient change y (J.
    Comput. Chem. 15, 1 (1994)): phi SR1 + (1 - phi) PSB with
    phi = (xi.s)^2/(xi.xi s.s), xi = y - W s. It meets the secant
    condition and, unlike the damped BFGS update, keeps negative
    curvature"""
    ss = s @ s
    if ss < 1.e-12:
        return W
    xi = y - W @ s
    xx = xi @ xi
    if xx < 1.e-24:
        return W
    sx = s @ xi
    phi = sx*sx/(xx*ss)
    W = W + (1.0 - phi)*((np.outer(xi, s) + np.outer(s, xi))/ss
                         - sx*np.outer(s, s)/ss**2)
    if abs(sx) > 1.e-10*np.sqrt(xx*ss):
        W = W + phi*np.outer(xi, xi)/sx
    return W


class SeamPoint:
    """
    one geometry of the MECI search: energies, gradients and coupling of
    the states I < J and the derived branching-plane quantities. h_ref:
    the coupling of the previous point, to which the sign of h is
    aligned; tr = False for model spaces without rigid-body motion. P:
    orthonormal basis of the directions excluded from the optimisation,
    given for the working-space twin of a point in internal coordinates
    (the redundancies) in place of the translations and rotations.
    w is the twin whose vectors the steps are built from: the point
    itself in Cartesian coordinates
    """
    def __init__(self, x, e, grads, h, h_ref=None, tr=True, info=None, P=None):
        self.x = np.asarray(x, dtype=float).ravel()
        n = self.x.size
        self.info = {} if info is None else info
        if P is not None:
            self.P = P
        else:
            self.P = tr_basis(self.x) if tr else np.zeros((n, 0))
        self.w = self

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

        # kept by the solver: the floored model Hessian and the basis of
        # the approach steps
        self.w0 = None
        self.th = None

    @property
    def R2(self):
        """[g h], including inactive columns"""
        return np.column_stack([self.g, self.h])

    def approach_basis(self):
        """orthonormal basis of the space of the approach steps: no
        rigid-body motion and, if the coupling constraint is active, no
        component along h"""
        if self.th is None:
            cols = [self.P]
            if 1 in self.active:
                cols.append((self.h/np.linalg.norm(self.h))[:, None])
            self.th = complement(np.hstack(cols), self.x.size)
        return self.th

    def proj(self, v):
        """v without its translation/rotation components"""
        return v - self.P @ (self.P.T @ v)


class Step:
    """a proposed step: dx, its gap part dxb, the predicted change of
    the merit function, whether it was truncated, the trust radius, the
    penalty weight mu, the merit function ('lagrangian' when the seam
    part of a composite step dominates, else 'l1') and the kind of step
    ('penalty': approach step, or 'composite'). Diagnostics: the
    Hessian W of the step gives the model curvature along dx; hmin is
    the lowest eigenvalue of the reduced Hessian and npairs the number
    of secant pairs in W. Composite steps also carry the second-order
    gap correction c (included in dx and dxb), the seam part zs of the
    step and whether the gap step was truncated"""
    def __init__(self, dx, dxb, pred, trunc, trust, mu, merit,
                 kind='composite', W=None, hmin=np.nan, npairs=0,
                 c=None, zs=None, gap_trunc=False):
        self.dx = dx
        self.dxb = dxb
        # set by MeciSolver.step: the geometry the step leads to, its
        # Cartesian displacement and the iterations of the
        # back-transformation (dx itself is in the working space)
        self.x_new = None
        self.dx_cart = dx
        self.bt_iter = 0
        self.c = np.zeros(dx.size) if c is None else c
        self.zs = np.zeros(dx.size) if zs is None else zs
        self.gap_trunc = gap_trunc
        self.pred = pred
        self.trunc = trunc
        self.trust = trust
        self.mu = mu
        self.merit = merit
        self.kind = kind
        dd = float(dx @ dx)
        self.curv = float(dx @ W @ dx)/dd if W is not None and dd > 0. else np.nan
        self.hmin = float(hmin)
        self.npairs = int(npairs)


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


class CoordinateFailure(Exception):
    """the coordinate system cannot describe the geometry"""


class CartesianCoords:
    """
    Cartesian working space of the optimiser: the twin of a point is the
    point itself, and differences and displacements are plain vector
    operations. Other coordinate systems (internal_coords.py) provide
    the same four methods; the steps are built from the twins (pt.w),
    while energies, the seam gradient of the convergence test and the
    characterisation stay with the Cartesian points
    """
    name = 'cartesian'
    nprim = 0
    nrebuild = 0

    def attach(self, pt, rebuild=True):
        """set the working-space twin pt.w; returns whether the
        coordinates were rebuilt for the geometry"""
        pt.w = pt
        return False

    def diff(self, a, b):
        """working-space displacement from the point a to the point b"""
        return b.x - a.x

    def displace(self, pt, dx):
        """the geometry for the working-space step dx from pt: (x,
        Cartesian displacement, converged, iterations)"""
        return pt.x + dx, dx, True, 0

    def hess_guess(self, solver, pt):
        """guess Hessian in the working space at pt"""
        return solver.hess_guess(pt.x)


class MeciSolver:
    """
    the MECI iterations on top of an evaluator x -> (x, (E_I, E_J),
    (g_I, g_J), h, info) or None on failure. hess_guess: x -> model
    Hessian (3N, 3N); conv: thresholds gap, grad_max, grad_rms,
    step_max, step_rms, energy; approach: 'penalty' (approach steps on
    the penalty function far from the seam) or 'newton' (composite
    steps only); gap_correction: second-order correction of the gap
    step from the measured curvature of the seam; coords: the
    coordinate system of the steps (Cartesian if None); algorithm:
    'lagrangian' or 'projected' (two-phase projected-gradient algorithm
    with its far phase above the gap gap_switch)
    """
    def __init__(self, hess_guess, conv, trust=0.3, trust_min=0.01,
                 trust_max=0.5, hess_update='bfgs', hess_window=5,
                 max_iter=50, eig_floor=5.e-3, gap_cap=0.866, tr=True,
                 approach='penalty', gap_correction=False, coords=None,
                 algorithm='lagrangian', gap_switch=5.e-3):
        self.hess_guess = hess_guess
        self.algorithm = algorithm
        self.gap_switch = gap_switch
        self.coords = CartesianCoords() if coords is None else coords
        # messages for the caller (rebuilt or abandoned coordinates) and
        # the number of failed back-transformations
        self.events = []
        self.nbt_fail = 0
        self.conv = conv
        self.trust = trust
        self.trust_min = trust_min
        self.trust_max = trust_max
        self.hess_update = hess_update
        # number of secant pairs kept for the Hessian updates; 0 = all
        # (since the last reset of the window)
        self.hess_window = hess_window
        self.max_iter = max_iter
        self.eig_floor = eig_floor
        self.gap_cap = gap_cap
        self.tr = tr
        self.approach = approach
        self.gap_correction = gap_correction
        # curvature of the seam measured on the last composite step:
        # (kappa (n,), unit direction of its seam step (n,)) or None
        self.curv = None
        # secant pairs (older point, newer point) of the BFGS updates of
        # the seam Hessian (cleared when the branching plane turns) and
        # of the Hessian of the penalty function (all pairs)
        self.pairs = []
        self.pairs_all = []
        # penalty weight of the merit function
        self.mu = 0.
        # number of consecutive rejected steps
        self.nrej_run = 0
        # secant curvature s.y/s.s of the last trial step (diagnostic)
        self.sec_curv = np.nan
        # wall time spent in the optimiser itself (s)
        self.wall = 0.

    def window(self, pairs):
        """the last hess_window secant pairs of the list (all for 0)"""
        return pairs[-self.hess_window:] if self.hess_window > 0 else pairs

    def point(self, raw, ref=None):
        """SeamPoint from the evaluator output, h aligned to ref, with
        its working-space twin. If the coordinates are rebuilt for the
        new geometry, the twin of ref is recomputed in the new set and
        the secant pairs and the seam curvature, which belong to the old
        set, are dropped. If the coordinates fail, the optimisation
        continues in Cartesian coordinates"""
        x, e, grads, h, info = raw
        pt = SeamPoint(x, e, grads, h,
                       h_ref=None if ref is None else ref.h,
                       tr=self.tr, info=info)
        try:
            rebuilt = self.coords.attach(pt)
            if rebuilt:
                self.events.append(self.coords.name + ' coordinates rebuilt')
                if ref is not None:
                    self.coords.attach(ref, rebuild=False)
        except CoordinateFailure as err:
            self.events.append(str(err) + ': continuing in Cartesian coordinates')
            self.coords = CartesianCoords()
            self.coords.attach(pt)
            if ref is not None:
                self.coords.attach(ref)
            rebuilt = True
        if rebuilt:
            self.pairs = []
            self.pairs_all = []
            self.curv = None
        return pt

    def model_hessian(self, pt):
        """model Hessian at pt in the working space, with the eigenvalues
        in the space of the allowed motions floored at eig_floor (kept
        with the twin of the point)"""
        w = pt.w
        if w.w0 is None:
            n = w.x.size
            H = np.asarray(self.coords.hess_guess(self, pt), dtype=float)
            C = getattr(w, 'V', None)
            if C is None:
                C = complement(w.P, n)
            ev, V = np.linalg.eigh(C.T @ (0.5*(H + H.T)) @ C)
            V = C @ V
            w.w0 = (V*np.maximum(ev, self.eig_floor)) @ V.T
        return w.w0

    def hessian(self, pt):
        """seam (Lagrangian) Hessian at pt: the model Hessian plus the
        BFGS updates of the stored pairs with the Lagrangian gradients of
        the points' own multipliers"""
        w = pt.w
        W = self.model_hessian(pt)
        if self.hess_update == 'bfgs':
            for a, b in self.window(self.pairs):
                s = w.proj(self.coords.diff(a, b))
                y = w.proj(b.w.gl - a.w.gl)
                W = bfgs(W, s, y)
        return W

    def penalty_hessian(self, pt, mu):
        """Hessian of the penalty function Ebar + mu dE at pt: the model
        Hessian plus Bofill updates with the gradients gbar + mu g of the
        stored pairs (none for hess_update = 'none'). A pair whose gap
        gradients differ in direction lies on two sides of the cone,
        where the gradient of dE jumps, and is left out. Returns the
        Hessian and the number of pairs used"""
        w = pt.w
        W = self.model_hessian(pt)
        npairs = 0
        if self.hess_update == 'none':
            return W, npairs
        for a, b in self.window(self.pairs_all):
            ag, bg = a.w.g, b.w.g
            if ag @ bg < PAIR_COS*np.linalg.norm(ag)*np.linalg.norm(bg):
                continue
            s = w.proj(self.coords.diff(a, b))
            y = w.proj((b.w.gbar + mu*bg) - (a.w.gbar + mu*ag))
            W = bofill(W, s, y)
            npairs += 1
        return W, npairs

    def penalty_step(self, pt):
        """
        approach step from pt: the restricted-step RFO step on the
        penalty function Ebar + mu dE in the space without rigid-body
        motion and h. mu is doubled until the step closes the fraction
        FEAS_FRAC of the gap closure of a pure gap step in the linear
        model; it is kept only if the step is taken. Returns None, i.e.
        the composite step is to be taken, if the gap constraint is
        inactive, the gap is within its threshold, the step would cross
        the linearised seam or the steering loop is exhausted
        """
        w = pt.w
        if 0 not in w.active or pt.de <= self.conv['gap']:
            return None
        trust = self.trust
        Th = w.approach_basis()
        gb = Th.T @ w.gbar
        gh = Th.T @ w.g
        best = min(pt.de, self.gap_cap*trust*np.linalg.norm(gh))
        mu = max(self.mu, MU_START)
        for ipass in range(MAX_STEER):
            W, npairs = self.penalty_hessian(pt, mu)
            Wh = Th.T @ W @ Th
            dz, _, _ = rsrfo(Wh, gb + mu*gh, trust)
            dx = Th @ dz
            lin = pt.de + w.g @ dx
            if lin <= 0.:
                return None
            if pt.de - lin >= FEAS_FRAC*best:
                self.mu = mu
                # predicted change of Ebar + mu dE (dE stays positive in
                # the linear model): the value of the RS-RFO model
                pred = float(w.gbar @ dx + 0.5*dx @ W @ dx + mu*(lin - pt.de))
                return Step(dx, np.zeros(w.x.size), pred, True, trust, mu,
                            'l1', kind='penalty', W=W,
                            hmin=np.linalg.eigvalsh(0.5*(Wh + Wh.T))[0],
                            npairs=npairs)
            mu *= MU_STEER
        return None

    # ------------------------------------------------------------------
    # projected-gradient algorithm
    # ------------------------------------------------------------------
    def composite_gradient(self, w):
        """gradient of the far phase at the twin w: 2 dE g/|g| plus the
        seam gradient (Bearpark, Robb and Schlegel)"""
        if 0 not in w.active:
            return w.gseam.copy()
        return 2.0*w.de*w.g/np.linalg.norm(w.g) + w.gseam

    def far_hessian(self, pt):
        """Hessian of the far phase at pt: the model Hessian plus BFGS
        updates with the composite gradients of the stored pairs; and the
        number of pairs"""
        w = pt.w
        W = self.model_hessian(pt)
        npairs = 0
        if self.hess_update == 'bfgs':
            for a, b in self.window(self.pairs_all):
                s = w.proj(self.coords.diff(a, b))
                y = w.proj(self.composite_gradient(b.w) - self.composite_gradient(a.w))
                W = bfgs(W, s, y)
                npairs += 1
        return W, npairs

    def far_step(self, pt):
        """far-phase step of the projected algorithm: the RS-RFO step on
        the composite gradient within the trust radius"""
        w = pt.w
        n = w.x.size
        W, npairs = self.far_hessian(pt)
        C = getattr(w, 'V', None)
        if C is None:
            C = complement(w.P, n)
        Wc = C.T @ W @ C
        dz, pred, _ = rsrfo(Wc, C.T @ self.composite_gradient(w), self.trust)
        return Step(C @ dz, np.zeros(n), pred, True, self.trust, self.mu, 'none',
                    kind='far', W=W, hmin=np.linalg.eigvalsh(0.5*(Wc + Wc.T))[0],
                    npairs=npairs)

    def projected_step(self, pt):
        """near-phase step of the projected algorithm: the Newton step of
        the linearised constraints in the branching plane and the RS-RFO
        step in the seam space on the projected gradient (without the
        coupling W dxb of the reduced gradient of the composite step)"""
        W = self.hessian(pt)
        w = pt.w
        trust = self.trust
        dxb = self.newton_gap_step(w)
        nb = np.linalg.norm(dxb)
        trunc = nb > self.gap_cap*trust
        if trunc:
            dxb = dxb*self.gap_cap*trust/nb
        tx = np.sqrt(max(trust**2 - dxb @ dxb, 0.))
        Wx = w.Tti.T @ W @ w.Tti
        dz, _, boundary = rsrfo(Wx, w.Tti.T @ w.gbar, tx)
        zs = w.Tti @ dz
        dx = dxb + zs
        self.mu = max(self.mu, MU_SCALE*abs(w.lam[0]) + MU_MIN)
        if np.linalg.norm(dxb) < np.linalg.norm(dz):
            merit = 'lagrangian'
            pred = float(w.gl @ dx + 0.5*dx @ W @ dx)
        else:
            merit = 'l1'
            pred = float(w.gbar @ dx + 0.5*dx @ W @ dx +
                         self.mu*(abs(w.de + w.g @ dx) - abs(w.de)))
        hmin = np.linalg.eigvalsh(0.5*(Wx + Wx.T))[0] if Wx.size else np.nan
        npairs = len(self.window(self.pairs)) if self.hess_update == 'bfgs' else 0
        return Step(dx, dxb, pred, trunc or boundary, trust, self.mu, merit,
                    kind='composite', W=W, hmin=hmin, npairs=npairs,
                    zs=zs, gap_trunc=trunc)

    def build_step(self, pt):
        """the working-space step from pt within the current trust
        radius. algorithm = 'projected': the far-phase step while the gap
        is above gap_switch, else the near-phase step. algorithm =
        'lagrangian': an approach step on the penalty function while the
        linearised seam is out of reach (approach = 'penalty'), else the
        composite step"""
        if self.algorithm == 'projected':
            if 0 in pt.w.active and pt.de > self.gap_switch:
                return self.far_step(pt)
            return self.projected_step(pt)
        if self.approach == 'penalty':
            stp = self.penalty_step(pt)
            if stp is not None:
                return stp
        return self.composite_step(pt)

    def step(self, pt):
        """the step from pt with the geometry it leads to (x_new) and its
        Cartesian displacement (dx_cart), which the convergence test and
        the step statistics use. If the back-transformation from internal
        coordinates fails, the trust radius is halved and the step built
        again, at most MAX_BACK times; after that its first, linear
        iterate is taken"""
        for itry in range(MAX_BACK + 1):
            stp = self.build_step(pt)
            stp.x_new, stp.dx_cart, ok, stp.bt_iter = self.coords.displace(pt, stp.dx)
            if ok:
                break
            self.nbt_fail += 1
            if itry < MAX_BACK:
                self.trust = max(self.trust_min, 0.5*np.linalg.norm(stp.dx))
        return stp

    @staticmethod
    def newton_gap_step(pt):
        """Newton step of the linearised constraints in the branching
        plane (g.dx = -dE, h.dx = 0), untruncated"""
        if not pt.active:
            return np.zeros(pt.x.size)
        dy = -np.linalg.solve(pt.R.T @ pt.Tb, pt.r)
        return pt.Tb @ dy

    def gap_correction_vector(self, pt, zs):
        """second-order gap correction of a composite step from pt with
        the seam part zs: minus the drift 1/2 (zs.z_prev)^2 kappa predicted
        from the curvature measured on the last composite step, in the
        branching plane of pt, at most SOC_MAX_FRAC |zs| long"""
        n = pt.x.size
        s = np.linalg.norm(zs)
        if self.curv is None or s == 0.:
            return np.zeros(n)
        kappa, zprev = self.curv
        drift = 0.5*(zs @ zprev)**2*kappa
        c = -pt.Tb @ (pt.Tb.T @ drift)
        cn = np.linalg.norm(c)
        if cn > SOC_MAX_FRAC*s:
            c *= SOC_MAX_FRAC*s/cn
        return c

    def update_curvature(self, cur, new, stp):
        """measure the curvature of the seam on the accepted step stp from
        cur to new: the drift of its seam step out of the seam is minus
        the Newton gap step at new minus the correction applied with the
        step (projected on the branching plane of new). Only composite
        steps with an untruncated gap step, no turn of the branching
        plane and a seam step longer than SOC_MIN_STEP are used; the
        estimate is discarded after approach steps, truncated gap steps
        and turns of the plane, and kept otherwise"""
        if not self.gap_correction:
            return
        cur, new = cur.w, new.w
        # a step built in coordinates that were rebuilt since is of no use
        ok = stp.kind == 'composite' and not stp.gap_trunc and \
            stp.zs.size == new.x.size == cur.x.size and \
            0 in cur.active and 0 in new.active and \
            cur.Tb.shape[1] == new.Tb.shape[1]
        if ok:
            cos = np.linalg.svd(cur.Tb.T @ new.Tb, compute_uv=False)
            ok = cos[-1] >= 0.5
        if not ok:
            self.curv = None
            return
        s = np.linalg.norm(stp.zs)
        if s < SOC_MIN_STEP:
            return
        drift = -new.Tb @ (new.Tb.T @ (self.newton_gap_step(new) + stp.c))
        self.curv = (2.0*drift/s**2, stp.zs/s)

    def composite_step(self, pt):
        """the composite step from pt: Newton step of the linearised
        constraints in the branching plane (with the second-order gap
        correction) and RS-RFO step in the seam space"""
        W = self.hessian(pt)
        pt = pt.w
        trust = self.trust
        Wx = pt.Tti.T @ W @ pt.Tti
        dxb0 = self.newton_gap_step(pt)

        def build(c):
            # gap step (priority, capped at gap_cap*trust) and RS-RFO step
            # in the seam space, total step within the trust radius
            dxb = dxb0 + c
            nb = np.linalg.norm(dxb)
            trunc = nb > self.gap_cap*trust
            if trunc:
                dxb = dxb*self.gap_cap*trust/nb
            tx = np.sqrt(max(trust**2 - dxb @ dxb, 0.))
            gx = pt.Tti.T @ (pt.gbar + W @ dxb)
            dz, _, boundary = rsrfo(Wx, gx, tx)
            return dxb, dz, trunc, boundary

        c = np.zeros(pt.x.size)
        dxb, dz, trunc, boundary = build(c)
        if self.gap_correction and pt.active and not trunc:
            c = self.gap_correction_vector(pt, pt.Tti @ dz)
            if c @ c > 0.:
                step2 = build(c)
                if step2[2]:
                    c = np.zeros(pt.x.size)
                else:
                    dxb, dz, trunc, boundary = step2
        zs = pt.Tti @ dz
        dx = dxb + zs

        # predicted change of the merit function: the Lagrangian when the
        # seam step dominates, else Ebar + mu dE (the coupling constraint
        # is satisfied to first order by construction; the correction c
        # only cancels the predicted drift, so it is left out of the
        # linear gap)
        self.mu = max(self.mu, MU_SCALE*abs(pt.lam[0]) + MU_MIN)
        if np.linalg.norm(dxb) < np.linalg.norm(dz):
            merit = 'lagrangian'
            pred = float(pt.gl @ dx + 0.5*dx @ W @ dx)
        else:
            merit = 'l1'
            pred = float(pt.gbar @ dx + 0.5*dx @ W @ dx +
                         self.mu*(abs(pt.de + pt.g @ (dx - c)) - abs(pt.de)))
        hmin = np.linalg.eigvalsh(0.5*(Wx + Wx.T))[0] if Wx.size else np.nan
        npairs = len(self.window(self.pairs)) if self.hess_update == 'bfgs' else 0
        return Step(dx, dxb, pred, trunc or boundary, trust, self.mu, merit,
                    kind='composite', W=W, hmin=hmin, npairs=npairs,
                    c=c, zs=zs, gap_trunc=trunc)

    def converged(self, pt, stp, ebar_prev):
        """convergence flags at pt with the proposed step stp"""
        c = self.conv
        n = pt.x.size
        gs = pt.gseam
        dx = stp.dx_cart
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
        point's own multipliers (gbar projected off its own plane). The
        pairs of the Hessian of the penalty function are kept through
        such turns"""
        self.pairs_all.append((a, b))
        self.pairs_all = self.window(self.pairs_all)
        if a.w.Tb.shape[1] == b.w.Tb.shape[1] and a.w.Tb.shape[1] > 0:
            cos = np.linalg.svd(a.w.Tb.T @ b.w.Tb, compute_uv=False)
            if cos[-1] < 0.5:
                self.pairs = []
                return
        self.pairs.append((a, b))
        self.pairs = self.window(self.pairs)

    def assess(self, cur, new, stp):
        """accept or reject the trial point new reached from cur with stp
        (merit function of the step), and update the trust radius;
        returns (accepted, rho). No ratio, and so no update of the trust
        radius, when the predicted change is below the noise floor"""
        cw, nw = cur.w, new.w
        if stp.kind == 'far':
            # far phase of the projected algorithm: no merit function,
            # every step is accepted; the radius follows the norm of the
            # composite gradient
            self.add_pair(cur, new)
            gc, gn = self.composite_gradient(cw), self.composite_gradient(nw)
            s = cw.proj(self.coords.diff(cur, new))
            ss = s @ s
            self.sec_curv = float(s @ (gn - gc))/ss if ss > 0. else np.nan
            q = np.linalg.norm(gn)/max(np.linalg.norm(gc), 1.e-30)
            dxn = np.linalg.norm(stp.dx)
            if q > 1.0:
                self.trust = max(self.trust_min, FAR_SHRINK*dxn)
            elif q < FAR_DROP and dxn >= 0.8*self.trust:
                self.trust = min(self.trust_max, FAR_GROW*self.trust)
            self.nrej_run = 0
            return True, None
        if stp.merit == 'lagrangian':
            dphi = (new.ebar - cur.ebar) - cw.lam[0]*(new.de - cur.de)
        else:
            dphi = (new.ebar - cur.ebar) + stp.mu*(abs(new.de) - abs(cur.de))
        floor = max(PRED_FLOOR, NOISE_FAC*self.conv['energy'])
        rho = dphi/stp.pred if abs(stp.pred) > floor else None
        self.add_pair(cur, new)

        # diagnostic: secant curvature along the step, of the function
        # whose Hessian the step used
        s = cw.proj(self.coords.diff(cur, new))
        if stp.kind == 'penalty':
            y = (nw.gbar + stp.mu*nw.g) - (cw.gbar + stp.mu*cw.g)
        else:
            y = nw.gl - cw.gl
        ss = s @ s
        self.sec_curv = float(s @ y)/ss if ss > 0. else np.nan

        dxn = np.linalg.norm(stp.dx)
        if (rho is None or rho < 0.) and dphi > REJECT_RISE:
            shrink = REJECT_SHRINK[min(self.nrej_run, len(REJECT_SHRINK) - 1)]
            self.trust = max(self.trust_min, shrink*dxn)
            self.nrej_run += 1
            return False, rho
        self.nrej_run = 0
        if rho is not None:
            if rho < 0.25:
                self.trust = max(self.trust_min, 0.5*min(self.trust, dxn))
            elif rho > 0.75 and dxn >= 0.8*self.trust:
                self.trust = min(self.trust_max, TRUST_GROW*self.trust)
        return True, rho

    def run(self, evaluate, raw0, callback=None, on_accept=None,
            on_reject=None):
        """
        iterate from the evaluator output raw0 at the start geometry.
        callback(it, pt, stp, flags, rho, nrej, acc) is called once per
        iteration (before the next evaluation): stp is the step proposed
        at pt; rho, nrej and acc belong to the accepted step that led to
        pt (acc: its kind, merit function, mu, model and secant
        curvature along the step, lowest eigenvalue of the reduced
        Hessian, number of secant pairs and length of the second-order
        gap correction; None at the start).
        on_accept/on_reject are called for every trial point. Returns a
        dict with the status ('converged', 'max_iter', 'rejected',
        'failed'), the last point and the history
        """
        t0 = time.perf_counter()
        cur = self.point(raw0)
        hist = []
        ebar_prev = None
        rho = None
        acc = None
        nrej = 0
        status = 'max_iter'
        it = 0
        while True:
            stp = self.step(cur)
            conv, flags = self.converged(cur, stp, ebar_prev)
            self.wall += time.perf_counter() - t0
            hist.append({'x': cur.x.copy(), 'e': cur.e, 'lam': cur.lam.copy(),
                         'trust': stp.trust, 'rho': rho, 'acc': acc,
                         'grad_max': float(np.max(np.abs(cur.gseam))),
                         'grad_rms': float(np.linalg.norm(cur.gseam)/np.sqrt(cur.x.size)),
                         'step_max': float(np.max(np.abs(stp.dx_cart))),
                         'step_rms': float(np.linalg.norm(stp.dx_cart)/np.sqrt(cur.x.size)),
                         'info': cur.info})
            if callback is not None:
                callback(it, cur, stp, flags, rho, nrej, acc)
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
                raw = evaluate(stp.x_new)
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
                    acc = {'kind': stp.kind, 'merit': stp.merit, 'mu': stp.mu,
                           'curv_model': stp.curv, 'curv_secant': self.sec_curv,
                           'hess_min': stp.hmin, 'npairs': stp.npairs,
                           'gap_corr': float(np.linalg.norm(stp.c)),
                           'bt_iter': stp.bt_iter}
                    self.update_curvature(cur, new, stp)
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
