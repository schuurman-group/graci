"""
Model Hessian of Lindh, Bernhardsson, Karlstrom and Malmqvist (Chem.
Phys. Lett. 241, 423 (1995)) in Cartesian coordinates, the starting
Hessian of the MECI optimiser. Pure numpy/scipy.

H = sum_t k_t b_t b_t^T over all stretches, bends and torsions, b_t the
Wilson B-vector of the primitive t, with force constants

  stretch i-j      k = k_r   rho_ij
  bend    i-j-k    k = k_phi rho_ij rho_jk
  torsion i-j-k-l  k = k_tau rho_ij rho_jk rho_kl

rho_ij = exp(alpha_ij (r_ref,ij^2 - r_ij^2)), alpha and r_ref depending
on the periodic-table rows of i and j (H-He, Li-Ne, Na and beyond).
Terms whose rho products fall below the screening thresholds are left
out. Near-linear bends are replaced by the two linear bends
perpendicular to the bond axis; torsions over a near-linear angle are
left out. The result is translation/rotation free and positive
semi-definite (Bohr, Hartree).
"""

import numpy as np
import scipy.sparse as sparse

from graci.opt.meci_step import tr_basis

# parameters of the 1995 paper (rows H-He, Li-Ne, Na-)
ALPHA = np.array([[1.0000, 0.3949, 0.3949],
                  [0.3949, 0.2800, 0.2800],
                  [0.3949, 0.2800, 0.2800]])
R_REF = np.array([[1.35, 2.10, 2.53],
                  [2.10, 2.87, 3.40],
                  [2.53, 3.40, 3.40]])
K_STR = 0.45
K_BEND = 0.15
K_TORS = 0.005

# screening thresholds on the rho products of each term type
RHO_STR = 1.e-8
RHO_BEND = 1.e-6
RHO_TORS = 1.e-5
# bends with sin(angle) below this are treated as linear; torsions over
# such an angle are skipped
SIN_LINEAR = 0.0872    # sin(175 deg)
SIN_TORS = 0.1736      # sin(170 deg)


def period_row(Z):
    """0, 1, 2 for H-He, Li-Ne and Na onwards"""
    Z = np.asarray(Z)
    return np.where(Z <= 2, 0, np.where(Z <= 10, 1, 2))


def rho_matrix(X, Z):
    """rho_ij for all atom pairs (zero diagonal)"""
    row = period_row(Z)
    a = ALPHA[row[:, None], row[None, :]]
    r0 = R_REF[row[:, None], row[None, :]]
    d2 = np.sum((X[:, None, :] - X[None, :, :])**2, axis=2)
    rho = np.exp(a*(r0**2 - d2))
    np.fill_diagonal(rho, 0.)
    return rho


def b_stretch(X, i, j):
    """B-vectors (nt, 2, 3) of the distances i-j"""
    u = X[i] - X[j]
    u /= np.linalg.norm(u, axis=1)[:, None]
    return np.stack([u, -u], axis=1)


def b_bend(X, i, j, k):
    """B-vectors (nt, 3, 3) of the angles i-j-k (j central) and sin"""
    u = X[i] - X[j]
    v = X[k] - X[j]
    lu = np.linalg.norm(u, axis=1)
    lv = np.linalg.norm(v, axis=1)
    eu = u/lu[:, None]
    ev = v/lv[:, None]
    cos = np.clip(np.sum(eu*ev, axis=1), -1., 1.)
    sin = np.sqrt(1. - cos**2)
    sin_safe = np.where(sin > 0., sin, 1.)
    bi = (cos[:, None]*eu - ev)/(lu*sin_safe)[:, None]
    bk = (cos[:, None]*ev - eu)/(lv*sin_safe)[:, None]
    return np.stack([bi, -bi - bk, bk], axis=1), sin


def b_linear_bend(X, i, j, k):
    """B-vectors (2 nt, 3, 3) of the two linear bends perpendicular to
    the axis of the near-collinear triples i-j-k (angle near 180 deg, or
    near 0 deg with i and k on the same side of j, where the sign of the
    k term keeps the coordinate invariant under rotations)"""
    u = X[i] - X[j]
    v = X[k] - X[j]
    lu = np.linalg.norm(u, axis=1)
    lv = np.linalg.norm(v, axis=1)
    eu = u/lu[:, None]
    side = -np.sign(np.sum(u*v, axis=1))
    lv = lv*np.where(side == 0., 1., side)
    # first perpendicular direction: cross product with the Cartesian
    # axis least aligned with the bond axis
    ax = np.zeros_like(eu)
    ax[np.arange(len(eu)), np.argmin(np.abs(eu), axis=1)] = 1.
    w1 = np.cross(eu, ax)
    w1 /= np.linalg.norm(w1, axis=1)[:, None]
    w2 = np.cross(eu, w1)
    out = []
    for w in (w1, w2):
        bi = w/lu[:, None]
        bk = w/lv[:, None]
        out.append(np.stack([bi, -bi - bk, bk], axis=1))
    return np.concatenate(out, axis=0)


def b_torsion(X, i, j, k, l):
    """B-vectors (nt, 4, 3) of the dihedral angles i-j-k-l"""
    b1 = X[j] - X[i]
    b2 = X[k] - X[j]
    b3 = X[l] - X[k]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    l2 = np.linalg.norm(b2, axis=1)
    bi = -(l2/np.sum(n1**2, axis=1))[:, None]*n1
    bl = (l2/np.sum(n2**2, axis=1))[:, None]*n2
    f1 = (np.sum(b1*b2, axis=1)/l2**2)[:, None]
    f3 = (np.sum(b3*b2, axis=1)/l2**2)[:, None]
    bj = f3*bl - (1. + f1)*bi
    bk = f1*bi - (1. + f3)*bl
    return np.stack([bi, bj, bk, bl], axis=1)


def primitives(X, Z):
    """
    the screened primitives: list of (atoms (nt, k), B-vectors
    (nt, k, 3), force constants (nt,)) for stretches, bends (incl.
    linear bends) and torsions
    """
    X = np.asarray(X, dtype=float).reshape(-1, 3)
    rho = rho_matrix(X, Z)
    terms = []

    # stretches
    i, j = np.nonzero(np.triu(rho > RHO_STR, 1))
    if i.size:
        terms.append((np.stack([i, j], axis=1), b_stretch(X, i, j),
                      K_STR*rho[i, j]))

    # bends i-j-k, i < k, j central
    P3 = rho[:, :, None]*rho[:, None, :]          # (j, i, k)
    j, i, k = np.nonzero(P3 > RHO_BEND)
    keep = i < k
    i, j, k = i[keep], j[keep], k[keep]
    if i.size:
        f = K_BEND*rho[i, j]*rho[j, k]
        B, sin = b_bend(X, i, j, k)
        lin = sin < SIN_LINEAR
        if np.any(~lin):
            terms.append((np.stack([i, j, k], axis=1)[~lin], B[~lin], f[~lin]))
        if np.any(lin):
            Bl = b_linear_bend(X, i[lin], j[lin], k[lin])
            at = np.stack([i[lin], j[lin], k[lin]], axis=1)
            terms.append((np.concatenate([at, at]), Bl,
                          np.concatenate([f[lin], f[lin]])))

    # torsions i-j-k-l over the central pairs j < k: triples i-j-k first,
    # then the atoms l bonded to k
    tors = []
    A = rho[:, :, None]*rho[None, :, :]           # (i, j, k)
    ti, tj, tk = np.nonzero(A > RHO_TORS)
    keep = (tj < tk) & (ti != tk)
    ti, tj, tk = ti[keep], tj[keep], tk[keep]
    if ti.size:
        prod = A[ti, tj, tk][:, None]*rho[tk, :]
        m, l = np.nonzero(prod > RHO_TORS)
        keep = (l != tj[m]) & (l != ti[m])
        m, l = m[keep], l[keep]
        if m.size:
            tors.append(np.stack([ti[m], tj[m], tk[m], l], axis=1))
    if tors:
        at = np.concatenate(tors)
        i, j, k, l = at.T
        # skip torsions over near-linear angles
        _, s1 = b_bend(X, i, j, k)
        _, s2 = b_bend(X, j, k, l)
        ok = (s1 > SIN_TORS) & (s2 > SIN_TORS)
        if np.any(ok):
            at = at[ok]
            i, j, k, l = at.T
            terms.append((at, b_torsion(X, i, j, k, l),
                          K_TORS*rho[i, j]*rho[j, k]*rho[k, l]))
    return terms


def lindh(x, Z):
    """Lindh model Hessian (3N, 3N) at the geometry x (3N,) in Bohr for
    the nuclear charges Z"""
    X = np.asarray(x, dtype=float).reshape(-1, 3)
    n3 = X.size
    rows, cols, vals, fc = [], [], [], []
    nt = 0
    for at, B, f in primitives(X, Z):
        m, k = at.shape
        r = np.repeat(np.arange(nt, nt + m), 3*k)
        c = (3*at[:, :, None] + np.arange(3)[None, None, :]).ravel()
        rows.append(r)
        cols.append(c)
        vals.append(B.ravel())
        fc.append(f)
        nt += m
    if nt == 0:
        return np.zeros((n3, n3))
    Bm = sparse.csr_matrix((np.concatenate(vals),
                            (np.concatenate(rows), np.concatenate(cols))),
                           shape=(nt, n3))
    F = sparse.diags(np.concatenate(fc))
    H = (Bm.T @ (F @ Bm)).toarray()
    # linear bends of near-collinear (not exactly collinear) triples are
    # only approximately rotation invariant: project the rigid-body
    # motions out
    P = tr_basis(X.ravel())
    HP = H @ P
    H = H - HP @ P.T - P @ HP.T + P @ (P.T @ HP) @ P.T
    return 0.5*(H + H.T)
