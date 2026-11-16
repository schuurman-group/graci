"""
Redundant internal coordinates for the MECI optimiser ($meciopt,
coordinates = internal), after Peng, Ayala, Schlegel and Frisch (J.
Comput. Chem. 17, 49 (1996)). numpy, and the covalent radii of PySCF.

Primitives, from the geometry at build time:
  bonds       pairs closer than BOND_FAC (r_cov,i + r_cov,j); fragments
              are joined through their shortest interatomic distance
              (joining bonds)
  bends       i-j-k over two bonds at j; above 165 degrees two linear
              bends replace the angle: the sums of the projections of
              the two bond unit vectors on two lab-fixed axes
              perpendicular to the chain
  torsions    i-j-k-l over three bonds; over a collinear run j...k the
              torsions i-j...k-l of the atoms bonded to its two ends
              (a cumulene has no twist coordinate otherwise)
  impropers   one dihedral i-j-k-l (j the centre) at every atom with
              exactly three bonds (a planar centre has no out-of-plane
              coordinate otherwise)

At a geometry x the values q(x) and the Wilson matrix B = dq/dx define
the working space of the optimiser. Rigid-body motion is projected out
of B: the lab-fixed linear bends are not rotation invariant once the
chain bends, and the projection keeps the Cartesian steps free of
rotations.
  - the eigenvectors W of B^T B with non-zero eigenvalues L give the
    orthonormal basis V = B W L^-1/2 of the range of B in the space of
    the primitives; the rank must be 3N - 6 (3N - 5 for a linear
    molecule)
  - gradients transform as g_q = B (B^T B)^+ g_x, so the scalar product
    of two transformed vectors is the Cartesian one in the metric
    (B^T B)^+
  - the complement of V, the redundancies, is excluded from the
    optimisation like the translations and rotations in Cartesian
    coordinates
A displacement dq is transformed back iteratively,
x <- x + (B^T B)^+ B^T (q0 + dq - q(x)), until x no longer changes. The
target q0 + dq lies on the tangent plane and is in general not reached
exactly, so the residual in q is not a criterion.
"""

import numpy as np
from pyscf.data import radii

import graci.opt.model_hessian as model_hessian
from graci.opt.meci_step import SeamPoint, CoordinateFailure, tr_basis, complement

# bonds: closer than BOND_FAC times the sum of the covalent radii at
# build time; the primitives are rebuilt when a bond grows beyond
# BREAK_FAC times that sum (joining bonds excepted) or a non-bonded pair
# comes closer than CONTACT_FAC times it
BOND_FAC = 1.3
BREAK_FAC = 2.0
CONTACT_FAC = 1.1
# an angle becomes two linear bends above 165 degrees and an ordinary
# bend again below 160 degrees (cosines); dihedrals need both inner
# angles away from 0 and 180 degrees by more than 15 degrees. (With a
# limit of 175 degrees a cumulene bent to 174 degrees has ill-conditioned
# torsions and, with one angle on each side of the limit, no twist
# coordinate.)
COS_LIN_ON = np.cos(np.radians(165.))
COS_LIN_OFF = np.cos(np.radians(160.))
COS_DEGEN = np.cos(np.radians(15.))
# guess force constant of the improper dihedrals (Eh/rad^2; the
# out-of-plane value of Schlegel, Theor. Chim. Acta 66, 333 (1984)),
# times the rho factors of the three bonds
K_OOP = 0.045
# back-transformation: maximum number of iterations, rms change of x
# (Bohr) at convergence, growth of the change that counts as divergence
BT_MAXIT = 25
BT_TOL = 1.e-10
BT_DIVERGE = 1.e2
# relative threshold of the eigenvalues of B^T B
RANK_TOL = 1.e-8


class RankError(CoordinateFailure):
    """the primitives do not span the internal motions of the molecule"""


def cos_angle(X, i, j, k):
    """cosines of the angles i-j-k (j central)"""
    u = X[i] - X[j]
    v = X[k] - X[j]
    c = np.sum(u*v, axis=-1)/(np.linalg.norm(u, axis=-1)*np.linalg.norm(v, axis=-1))
    return np.clip(c, -1., 1.)


def dihedral(X, i, j, k, l):
    """dihedral angles i-j-k-l in (-pi, pi], in the sign convention of
    model_hessian.b_torsion"""
    b1 = X[j] - X[i]
    b2 = X[k] - X[j]
    b3 = X[l] - X[k]
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    e2 = b2/np.linalg.norm(b2, axis=-1)[..., None]
    return np.arctan2(np.sum(np.cross(n1, n2)*e2, axis=-1), np.sum(n1*n2, axis=-1))


def fragments(conn):
    """fragment label of every atom for the connectivity matrix conn"""
    n = conn.shape[0]
    label = -np.ones(n, dtype=int)
    nfrag = 0
    for a in range(n):
        if label[a] >= 0:
            continue
        stack = [a]
        label[a] = nfrag
        while stack:
            b = stack.pop()
            for c in np.nonzero(conn[b])[0]:
                if label[c] < 0:
                    label[c] = nfrag
                    stack.append(c)
        nfrag += 1
    return label


class RedundantInternals:
    """
    redundant internal coordinates of a molecule with the nuclear
    charges Z; guess: 'lindh' (force constants of the Lindh model for
    each primitive) or 'unit'. build(x) sets the primitives up at the
    geometry x (Bohr); the optimiser uses attach, diff, displace and
    hess_guess
    """
    name = 'internal'

    def __init__(self, Z, guess='lindh'):
        self.Z = np.asarray(Z, dtype=int)
        self.rcov = np.asarray(radii.COVALENT)[self.Z]
        self.rsum = self.rcov[:, None] + self.rcov[None, :]
        self.guess = guess
        # number of rebuilds of the primitives after the first build
        self.nrebuild = 0
        self.built = False

    # ------------------------------------------------------------------
    # primitives
    # ------------------------------------------------------------------
    def build(self, x):
        """the primitives at the geometry x; RankError if they do not
        span the 3N - 6 (3N - 5) internal motions"""
        X = np.asarray(x, dtype=float).reshape(-1, 3)
        N = X.shape[0]
        d = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        eye = np.eye(N, dtype=bool)
        conn = (d < BOND_FAC*self.rsum) & ~eye
        join = np.zeros((N, N), dtype=bool)
        while True:
            label = fragments(conn | join)
            if label.max() == 0:
                break
            dd = np.where(label[:, None] != label[None, :], d, np.inf)
            a, b = np.unravel_index(np.argmin(dd), dd.shape)
            join[a, b] = join[b, a] = True
        conn = conn | join
        bi, bj = np.nonzero(np.triu(conn, 1))
        self.conn = conn
        self.bonds = np.stack([bi, bj], axis=1)
        self.joining = join[bi, bj]
        nbrs = [np.nonzero(conn[a])[0] for a in range(N)]

        # bends and linear bends
        bends, lin = [], []
        for j in range(N):
            nb = nbrs[j]
            for a in range(len(nb)):
                for b in range(a + 1, len(nb)):
                    i, k = int(nb[a]), int(nb[b])
                    if cos_angle(X, i, j, k) < COS_LIN_ON:
                        lin.append((i, j, k))
                    else:
                        bends.append((i, j, k))
        self.bends = np.array(bends, dtype=int).reshape(-1, 3)
        self.lin3 = np.array(lin, dtype=int).reshape(-1, 3)
        # two lab-fixed axes perpendicular to each linear chain
        axes = []
        for (i, j, k) in lin:
            eu = (X[i] - X[j])/np.linalg.norm(X[i] - X[j])
            ax = np.zeros(3)
            ax[np.argmin(np.abs(eu))] = 1.
            w1 = np.cross(eu, ax)
            w1 /= np.linalg.norm(w1)
            axes += [w1, np.cross(eu, w1)]
        self.lin = np.repeat(self.lin3, 2, axis=0)
        self.lin_axis = np.array(axes, dtype=float).reshape(-1, 3)

        # torsions over three bonds, both inner angles ordinary bends
        linset = set(lin) | set((k, j, i) for (i, j, k) in lin)
        tors = []
        for (j, k) in self.bonds:
            j, k = int(j), int(k)
            for i in nbrs[j]:
                i = int(i)
                if i == k or (i, j, k) in linset:
                    continue
                for l in nbrs[k]:
                    l = int(l)
                    if l == j or l == i or (j, k, l) in linset:
                        continue
                    tors.append((i, j, k, l))
        self.ntors = len(tors)

        # torsions across the collinear runs
        runs = set()
        for (i, j, k) in lin:
            chain = [i, j, k]
            for end in (-1, 0):
                grown = True
                while grown:
                    grown = False
                    a, b = (chain[-2], chain[-1]) if end == -1 else (chain[1], chain[0])
                    for m in nbrs[b]:
                        m = int(m)
                        if m not in chain and (a, b, m) in linset:
                            chain = chain + [m] if end == -1 else [m] + chain
                            grown = True
                            break
            runs.add(tuple(chain) if chain[0] < chain[-1] else tuple(reversed(chain)))
        # a run is also taken with an end atom added whose angle to the far
        # end is still nearly linear (a bent chain): the torsion over the
        # shorter run would be degenerate there
        work = sorted(runs)
        seen = set(work)
        bridge = set()
        while work:
            chain = work.pop()
            j, k = chain[0], chain[-1]
            for (far, end, front) in ((k, j, True), (j, k, False)):
                for m in nbrs[end]:
                    m = int(m)
                    if m in chain or abs(cos_angle(X, far, end, m)) <= COS_DEGEN:
                        continue
                    new = (m,) + chain if front else chain + (m,)
                    new = new if new[0] < new[-1] else tuple(reversed(new))
                    if new not in seen:
                        seen.add(new)
                        work.append(new)
            for i in nbrs[j]:
                i = int(i)
                if i in chain or abs(cos_angle(X, i, j, k)) > COS_DEGEN:
                    continue
                for l in nbrs[k]:
                    l = int(l)
                    if l in chain or l == i or abs(cos_angle(X, j, k, l)) > COS_DEGEN:
                        continue
                    bridge.add((i, j, k, l))
        bridge = sorted(bridge)
        self.nbridge = len(bridge)

        # impropers at the three-coordinate centres (centre second)
        improper = []
        for a in range(N):
            if len(nbrs[a]) != 3:
                continue
            b, c, e = (int(v) for v in nbrs[a])
            for (p, q, r) in ((b, c, e), (c, e, b), (e, b, c)):
                if abs(cos_angle(X, p, a, q)) < COS_DEGEN and \
                        abs(cos_angle(X, a, q, r)) < COS_DEGEN:
                    improper.append((p, a, q, r))
                    break
        self.nimproper = len(improper)
        self.dih = np.array(tors + bridge + improper, dtype=int).reshape(-1, 4)

        nb, na, nl, nd = len(self.bonds), len(self.bends), len(self.lin), len(self.dih)
        self.nprim = nb + na + nl + nd
        self.is_dih = np.zeros(self.nprim, dtype=bool)
        self.is_dih[nb + na + nl:] = True
        self.built = True

        ntr = tr_basis(X.ravel()).shape[1]
        rank = self.frame(X.ravel())['rank']
        if rank != X.size - ntr:
            raise RankError('redundant internal coordinates: rank %d of %d '
                            '(%d primitives)' % (rank, X.size - ntr, self.nprim))
        return

    def summary(self):
        """the numbers of primitives of each type"""
        return ('{:d} bonds, {:d} bends, {:d} linear bends, {:d} torsions, '
                '{:d} impropers'.format(len(self.bonds), len(self.bends), len(self.lin),
                                        self.ntors + self.nbridge, self.nimproper))

    def needs_rebuild(self, x):
        """whether the primitives no longer suit the geometry x: a broken
        bond (joining bonds excepted), a new contact, an angle that has
        become linear or a linear one that has bent, a dihedral with a
        nearly collinear inner angle"""
        X = np.asarray(x, dtype=float).reshape(-1, 3)
        d = np.linalg.norm(X[:, None, :] - X[None, :, :], axis=2)
        bi, bj = self.bonds.T
        if np.any((d[bi, bj] > BREAK_FAC*self.rsum[bi, bj]) & ~self.joining):
            return True
        if np.any((d < CONTACT_FAC*self.rsum) & ~self.conn &
                  ~np.eye(len(X), dtype=bool)):
            return True
        if len(self.bends) and np.any(cos_angle(X, *self.bends.T) < COS_LIN_ON):
            return True
        if len(self.lin3) and np.any(cos_angle(X, *self.lin3.T) > COS_LIN_OFF):
            return True
        if len(self.dih):
            i, j, k, l = self.dih.T
            if np.any(np.abs(cos_angle(X, i, j, k)) > COS_DEGEN) or \
                    np.any(np.abs(cos_angle(X, j, k, l)) > COS_DEGEN):
                return True
        return False

    # ------------------------------------------------------------------
    # values, B matrix and frame
    # ------------------------------------------------------------------
    def values(self, x):
        """q(x): distances (Bohr), angles, linear bends and dihedrals
        (rad)"""
        X = np.asarray(x, dtype=float).reshape(-1, 3)
        out = [np.linalg.norm(X[self.bonds[:, 0]] - X[self.bonds[:, 1]], axis=1)]
        if len(self.bends):
            out.append(np.arccos(cos_angle(X, *self.bends.T)))
        if len(self.lin):
            i, j, k = self.lin.T
            eu = X[i] - X[j]
            ev = X[k] - X[j]
            eu = eu/np.linalg.norm(eu, axis=1)[:, None]
            ev = ev/np.linalg.norm(ev, axis=1)[:, None]
            out.append(np.sum(self.lin_axis*(eu + ev), axis=1))
        if len(self.dih):
            out.append(dihedral(X, *self.dih.T))
        return np.concatenate(out)

    def bmatrix(self, x):
        """Wilson matrix dq/dx (nprim, 3N), rigid-body motion projected
        out"""
        X = np.asarray(x, dtype=float).reshape(-1, 3)
        B = np.zeros((self.nprim, X.size))
        row = 0

        def scatter(atoms, vecs):
            nonlocal row
            nt, na = atoms.shape
            r = np.repeat(np.arange(row, row + nt), 3*na)
            c = (3*atoms[:, :, None] + np.arange(3)[None, None, :]).ravel()
            np.add.at(B, (r, c), vecs.ravel())
            row += nt

        scatter(self.bonds, model_hessian.b_stretch(X, self.bonds[:, 0], self.bonds[:, 1]))
        if len(self.bends):
            scatter(self.bends, model_hessian.b_bend(X, *self.bends.T)[0])
        if len(self.lin):
            i, j, k = self.lin.T
            u = X[i] - X[j]
            v = X[k] - X[j]
            lu = np.linalg.norm(u, axis=1)[:, None]
            lv = np.linalg.norm(v, axis=1)[:, None]
            eu, ev = u/lu, v/lv
            w = self.lin_axis
            bi = (w - np.sum(w*eu, axis=1)[:, None]*eu)/lu
            bk = (w - np.sum(w*ev, axis=1)[:, None]*ev)/lv
            scatter(self.lin, np.stack([bi, -bi - bk, bk], axis=1))
        if len(self.dih):
            scatter(self.dih, model_hessian.b_torsion(X, *self.dih.T))
        P = tr_basis(X.ravel())
        return B - (B @ P) @ P.T

    def frame(self, x):
        """q, B and the working-space frame at x: V (orthonormal basis of
        the range of B), agrad = B (B^T B)^+ (gradients to the working
        space; its transpose takes dq to dx) and the rank"""
        q = self.values(x)
        B = self.bmatrix(x)
        lam, W = np.linalg.eigh(B.T @ B)
        keep = lam > RANK_TOL*lam[-1]
        Wk = W[:, keep]/np.sqrt(lam[keep])
        V = B @ Wk
        return {'q': q, 'B': B, 'V': V, 'agrad': V @ Wk.T, 'rank': int(np.sum(keep))}

    # ------------------------------------------------------------------
    # interface of the optimiser
    # ------------------------------------------------------------------
    def attach(self, pt, rebuild=True):
        """give the point pt its working-space twin pt.w; the primitives
        are rebuilt first if they no longer suit the geometry (and
        rebuild is set). Returns whether they were rebuilt"""
        rebuilt = False
        if not self.built:
            self.build(pt.x)
        elif rebuild and self.needs_rebuild(pt.x):
            self.build(pt.x)
            self.nrebuild += 1
            rebuilt = True
        f = self.frame(pt.x)
        ntr = pt.P.shape[1]
        if f['rank'] != pt.x.size - ntr:
            if not rebuild or rebuilt:
                raise RankError('redundant internal coordinates: rank %d of %d'
                                % (f['rank'], pt.x.size - ntr))
            self.build(pt.x)
            self.nrebuild += 1
            rebuilt = True
            f = self.frame(pt.x)
        A = f['agrad']
        w = SeamPoint(f['q'], pt.e, (A @ (pt.gbar - 0.5*pt.g), A @ (pt.gbar + 0.5*pt.g)),
                      A @ pt.h, tr=False, info=pt.info,
                      P=complement(f['V'], self.nprim))
        w.V = f['V']
        pt.w = w
        return rebuilt

    def diff(self, a, b):
        """q(b) - q(a), dihedral differences wrapped to (-pi, pi]"""
        return self.wrap(b.w.x - a.w.x)

    def wrap(self, dq):
        dq = np.array(dq, dtype=float)
        m = self.is_dih
        dq[m] = np.pi - np.mod(np.pi - dq[m], 2.0*np.pi)
        return dq

    def displace(self, pt, dq):
        """Cartesian geometry for the working-space step dq from pt:
        returns (x, x - pt.x, converged, iterations). If the iteration
        fails, x is its first, linear iterate"""
        target = pt.w.x + dq
        x = pt.x.copy()
        first, rms0 = None, None
        for it in range(BT_MAXIT):
            f = self.frame(x)
            dx = f['agrad'].T @ self.wrap(target - f['q'])
            x = x + dx
            rms = np.sqrt(np.mean(dx*dx))
            if it == 0:
                first, rms0 = x.copy(), rms
            if rms < BT_TOL:
                return x, x - pt.x, True, it + 1
            if rms > BT_DIVERGE*max(rms0, BT_TOL):
                break
        return first, first - pt.x, False, BT_MAXIT

    def hess_guess(self, solver, pt):
        """diagonal guess Hessian in the primitives at pt: the force
        constants of the Lindh model, or 0.5 for guess = 'unit'"""
        if self.guess != 'lindh':
            return 0.5*np.eye(self.nprim)
        X = pt.x.reshape(-1, 3)
        rho = model_hessian.rho_matrix(X, self.Z)
        k = [model_hessian.K_STR*rho[self.bonds[:, 0], self.bonds[:, 1]]]
        for tri in (self.bends, self.lin):
            if len(tri):
                k.append(model_hessian.K_BEND*rho[tri[:, 0], tri[:, 1]]*rho[tri[:, 1], tri[:, 2]])
        if len(self.dih):
            i, j, m, l = self.dih.T
            kd = np.empty(len(self.dih))
            t, b = self.ntors, self.ntors + self.nbridge
            kd[:t] = model_hessian.K_TORS*rho[i[:t], j[:t]]*rho[j[:t], m[:t]]*rho[m[:t], l[:t]]
            kd[t:b] = model_hessian.K_TORS*rho[i[t:b], j[t:b]]*rho[m[t:b], l[t:b]]
            kd[b:] = K_OOP*rho[j[b:], i[b:]]*rho[j[b:], m[b:]]*rho[j[b:], l[b:]]
            k.append(kd)
        return np.diag(np.concatenate(k))
