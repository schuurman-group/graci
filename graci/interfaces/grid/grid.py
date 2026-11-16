"""
Molecular quadrature grids from libgrid.

GraciGrids is a pyscf.dft.gen_grid.Grids whose build() uses libgrid for
the atomic grids (Treutler-Ahlrichs radial grids, Lebedev angular grids,
NWChem pruning, PySCF's level tables and point order), the partition
weights (Becke, or the local Delaunay cells and weights of Gill &
Tsoukalas, JCP 165, 174104 (2026)) and PySCF's spatial box ordering.
The remaining attributes (atm_idx, quadrature_weights, padding, non0tab)
are set exactly as Grids.build sets them, so every consumer of mf.grids
(libscf, libmrsf, PySCF numint) is unchanged.

weight_response(grids, q) contracts the nuclear derivatives of the
weights with per-point integrands (grid response of gradients).
"""
import sys
import numpy as np
from pyscf import gto
from pyscf.dft import gen_grid, radi
from pyscf.lib import logger
import graci.core.libs as libs

RADIAL = {'treutler': 1, 'gauss_chebyshev': 2, 'mura_knowles': 3, 'becke': 4, 'delley': 5}
PRUNE  = {'none': 0, 'nwchem': 1, 'sg1': 2, 'treutler': 3}
ADJUST = {'none': 0, 'treutler': 1, 'becke': 2}
STEP   = {'b3': 1, 'g8': 2, 'gmu': 2}
MODE   = {'becke': 1, 'delaunay': 2}

_RADI_FUNCS = {radi.treutler_ahlrichs: 'treutler', radi.gauss_chebyshev: 'gauss_chebyshev',
               radi.mura_knowles: 'mura_knowles', radi.becke: 'becke', radi.delley: 'delley'}
_PRUNE_FUNCS = {gen_grid.nwchem_prune: 'nwchem', gen_grid.sg1_prune: 'sg1',
                gen_grid.treutler_prune: 'treutler', None: 'none'}
_ADJUST_FUNCS = {radi.treutler_atomic_radii_adjust: 'treutler',
                 radi.becke_atomic_radii_adjust: 'becke', None: 'none'}


def _lib(name, args):
    """libgrid call with error check"""
    libs.lib_load('grid')
    out = libs.lib_func(name, args)
    ierr = out[-1] if isinstance(out, tuple) else out
    if int(ierr) != 0:
        sys.exit('\n ' + name + ': libgrid error ' + str(int(ierr)) + '\n')
    return out


def nuclear_charges(mol):
    """nuclear charges from the atom symbols (as PySCF's grid builder)"""
    return np.array([gto.charge(mol.atom_symbol(ia)) for ia in range(mol.natm)], dtype=np.int32)


class GraciGrids(gen_grid.Grids):
    """
    Grids built by libgrid. PySCF attributes (level, prune, radi_method,
    radii_adjust, atom_grid) keep their meaning; new attributes:
      partition     'becke' | 'delaunay': cells and weights (Gill & Tsoukalas)
      cells         override of the cell definition ('becke' | 'delaunay')
      weights_mode  override of the weight definition
      step          'b3' (Becke) | 'g8' | 'gmu' (with step_mu)
      step_mu       exponent of the G_mu step
      neighbour_tol tolerance (bohr) of the Delaunay contact test
      neighbour_shells  1: Delaunay neighbours (the paper's protocol); 2
                    (default): neighbours of neighbours as well. One shell
                    gives mEh quadrature errors for planar molecules on the
                    PySCF grids (the excluded cells are poorly suppressed with
                    only 3-6 step factors); two shells reproduce Becke's
                    accuracy at a cost that is still O(1) per point
      builder       'graci' | 'pyscf' (PySCF's own build, for comparisons)
    After build(): nnb, nb (neighbour lists, 1-based), knear (nearest
    nucleus of every point, 1-based, 0 = padding) and the arrays used by
    weight_response.
    """
    partition = 'becke'
    cells = None
    weights_mode = None
    step = 'b3'
    step_mu = 8.0
    neighbour_tol = 1.0e-6
    neighbour_shells = 2
    builder = 'graci'

    _keys = gen_grid.Grids._keys.union({'partition', 'cells', 'weights_mode', 'step', 'step_mu',
                                        'neighbour_tol', 'neighbour_shells', 'builder', 'nnb', 'nb', 'knear'})

    def __init__(self, mol, partition='becke', step='b3', adjust='treutler', step_mu=8.0, shells=2):
        gen_grid.Grids.__init__(self, mol)
        self.partition = partition
        self.step = step
        self.step_mu = step_mu
        self.neighbour_shells = int(shells)
        self.radii_adjust = {'treutler': radi.treutler_atomic_radii_adjust,
                             'becke': radi.becke_atomic_radii_adjust, 'none': None}[adjust]
        self.nnb = None
        self.nb = None
        self.knear = None
        self._resp = None

    def __setattr__(self, key, val):
        if key in ('partition', 'cells', 'weights_mode', 'step', 'step_mu', 'neighbour_tol', 'neighbour_shells',
                   'builder'):
            self.reset()
        gen_grid.Grids.__setattr__(self, key, val)

    def dump_flags(self, verbose=None):
        gen_grid.Grids.dump_flags(self, verbose)
        logger.info(self, 'libgrid partition: %s cells, %s weights, %d neighbour shell(s), %s step, builder %s',
                    self.cells or self.partition, self.weights_mode or self.partition, self.neighbour_shells,
                    self.step, self.builder)
        return self

    # ------------------------------------------------------------------
    def _codes(self):
        """libgrid scheme codes of the current settings"""
        if self.radi_method not in _RADI_FUNCS:
            raise NotImplementedError('GraciGrids: radial scheme ' + str(self.radi_method))
        if self.prune not in _PRUNE_FUNCS:
            raise NotImplementedError('GraciGrids: pruning scheme ' + str(self.prune))
        if self.radii_adjust not in _ADJUST_FUNCS:
            raise NotImplementedError('GraciGrids: radii adjustment ' + str(self.radii_adjust))
        if self.radii_adjust is not None and self.atomic_radii is not radi.BRAGG_RADII:
            raise NotImplementedError('GraciGrids: only Bragg radii are tabulated')
        if self.becke_scheme is not gen_grid.original_becke:
            raise NotImplementedError('GraciGrids: becke_scheme ' + str(self.becke_scheme))
        cells = MODE[(self.cells or self.partition).lower()]
        weights = MODE[(self.weights_mode or self.partition).lower()]
        step = STEP[self.step.lower()]
        mu = 8.0 if self.step.lower() == 'g8' else float(self.step_mu)
        return dict(radial=RADIAL[_RADI_FUNCS[self.radi_method]], prune=PRUNE[_PRUNE_FUNCS[self.prune]],
                    adjust=ADJUST[_ADJUST_FUNCS[self.radii_adjust]], cells=cells, weights=weights,
                    step=step, mu=mu)

    def _user_sizes(self, mol):
        """per-atom (nrad, nang) overrides from atom_grid (0 = default)"""
        natm = mol.natm
        nrad = np.zeros(natm, dtype=np.int32)
        nang = np.zeros(natm, dtype=np.int32)
        atom_grid = self.atom_grid
        if isinstance(atom_grid, (list, tuple)):
            nrad[:] = atom_grid[0]
            nang[:] = atom_grid[1]
        elif atom_grid:
            default = atom_grid.get('default', None)
            for ia in range(natm):
                cfg = atom_grid.get(mol.atom_symbol(ia), default)
                if cfg is not None:
                    nrad[ia], nang[ia] = cfg
        return nrad, nang

    def build(self, mol=None, with_non0tab=False, sort_grids=True, **kwargs):
        if self.builder == 'pyscf':
            return gen_grid.Grids.build(self, mol, with_non0tab=with_non0tab, sort_grids=sort_grids, **kwargs)
        if mol is None:
            mol = self.mol
        if self.verbose >= logger.WARN:
            self.check_sanity()
        codes = self._codes()
        natm = mol.natm
        R = np.asfortranarray(mol.atom_coords().T, dtype=float)
        Z = nuclear_charges(mol)
        nrad_user, nang_user = self._user_sizes(mol)

        # atomic grids
        npts, _ = _lib('grid_size', [natm, Z, int(self.level), codes['radial'], codes['prune'],
                                     nrad_user, nang_user, 0, 0])
        npts = int(npts)
        xyz = np.zeros((3, npts), order='F')
        iatom = np.zeros(npts, dtype=np.int32)
        vol = np.zeros(npts)
        _lib('grid_atomic', [natm, R, Z, int(self.level), codes['radial'], codes['prune'], nrad_user,
                             nang_user, npts, xyz, iatom, vol, 0])

        # neighbour lists (Delaunay cells or weights)
        ld = max(natm - 1, 1)
        nnb = np.zeros(natm, dtype=np.int32)
        nb = np.zeros((ld, natm), dtype=np.int32, order='F')
        if codes['cells'] == MODE['delaunay'] or codes['weights'] == MODE['delaunay']:
            _lib('grid_neighbours', [natm, R, float(self.neighbour_tol), nnb, nb, 0])
            if int(self.neighbour_shells) > 1:
                extend_shells(nnb, nb, int(self.neighbour_shells))

        # partition weights
        w = np.zeros(npts)
        knear = np.zeros(npts, dtype=np.int32)
        _lib('grid_partition', [natm, R, Z, codes['adjust'], codes['step'], codes['mu'], codes['cells'],
                                codes['weights'], npts, xyz, iatom, vol, nnb, nb, w, knear, 0])

        coords = np.ascontiguousarray(xyz.T)
        if sort_grids:
            perm = np.zeros(npts, dtype=np.int32)
            _lib('grid_sort', [natm, R, npts, xyz, perm, 0])
            idx = perm - 1
            coords = coords[idx]
            w = w[idx]
            iatom = iatom[idx]
            vol = vol[idx]
            knear = knear[idx]

        if self.alignment > 1:
            padding = gen_grid._padding_size(npts, self.alignment)
            logger.debug(self, 'Padding %d grids', padding)
            if padding > 0:
                coords = np.vstack([coords, np.repeat([[1e-4]*3], padding, axis=0)])
                w = np.hstack([w, np.zeros(padding)])
                iatom = np.hstack([iatom, np.zeros(padding, dtype=np.int32)])
                vol = np.hstack([vol, np.zeros(padding)])
                knear = np.hstack([knear, np.zeros(padding, dtype=np.int32)])

        self.coords = coords
        self.weights = w
        self.atm_idx = (iatom - 1).astype(np.int32)
        self.quadrature_weights = vol
        self.knear = knear
        self.nnb = nnb
        self.nb = nb
        # arrays for the weight response (grid order, padding = atom 0)
        self._resp = dict(R=R, Z=Z, codes=codes, xyz=np.asfortranarray(coords.T), iatom=iatom, vol=vol,
                          nnb=nnb, nb=nb)

        if with_non0tab:
            self.non0tab = self.make_mask(mol, self.coords)
            self.screen_index = self.non0tab
        else:
            self.screen_index = self.non0tab = None
        logger.info(self, 'tot grids = %d', len(self.weights))
        return self

    def reset(self, mol=None):
        gen_grid.Grids.reset(self, mol)
        self._resp = None
        return self

    def neighbour_lists(self):
        """list of 0-based neighbour index arrays of every atom"""
        if self.nnb is None:
            self.build()
        return [np.sort(self.nb[:self.nnb[i], i] - 1) for i in range(self.nnb.size)]

    def sigma(self):
        """sigma(r) = sum over all cells - sum over the local cells (eq. 22)
        at every grid point (grid order; diagnostics, O(natm N))"""
        if self._resp is None:
            self.build()
        d = self._resp
        natm = d['R'].shape[1]
        npts = d['iatom'].size
        if d['nnb'].sum() == 0 and natm > 1:
            nnb = np.zeros(natm, dtype=np.int32)
            nb = np.zeros((max(natm - 1, 1), natm), dtype=np.int32, order='F')
            _lib('grid_neighbours', [natm, d['R'], float(self.neighbour_tol), nnb, nb, 0])
        else:
            nnb, nb = d['nnb'], d['nb']
        sig = np.zeros(npts)
        _lib('grid_sigma', [natm, d['R'], d['Z'], d['codes']['adjust'], d['codes']['step'], d['codes']['mu'],
                            d['codes']['cells'], npts, d['xyz'], d['iatom'], nnb, nb, sig, 0])
        return sig


def extend_shells(nnb, nb, shells):
    """replace the neighbour lists (1-based, in place) by their closure over
    `shells` hops (shells = 2: neighbours of neighbours)"""
    natm = nnb.size
    lists = [set(nb[:nnb[i], i] - 1) for i in range(natm)]
    for _ in range(shells - 1):
        new = []
        for i in range(natm):
            s = set(lists[i])
            for j in lists[i]:
                s |= lists[j]
            s.discard(i)
            new.append(s)
        lists = new
    for i in range(natm):
        l = np.array(sorted(lists[i]), dtype=np.int32)
        nnb[i] = l.size
        nb[:l.size, i] = l + 1


def xc_grid_response(ni, mol, grids, xc, dms, blksize=4000):
    """
    Grid response of the XC energy E = sum_p w_p F(rho(r_p)) for the spin
    densities dms = (D_a, D_b) (ROKS/UKS) or a single D (RKS): the
    derivative of the weights, sum_p F_p dw_p/dR_A (libgrid), plus the
    translation of the grid points of atom A with their atom,
    sum_{p in A} w_p grad_r F(r_p) = sum_{p in A} w_p sum_s [v_s grad rho_s
    + sum_k v_s^k grad(d_k rho_s)] (LDA and GGA). This is PySCF's `excsum`
    of get_vxc_full_response. Returns (natm, 3).
    """
    if not isinstance(grids, GraciGrids) or grids._resp is None:
        raise RuntimeError('xc_grid_response: the grid must be a built GraciGrids')
    xctype = ni._xc_type(xc)
    if xctype not in ('LDA', 'GGA'):
        raise NotImplementedError('xc_grid_response: ' + xctype + ' functionals')
    dms = np.asarray(dms)
    if dms.ndim == 2:
        dms = np.array([0.5*dms, 0.5*dms])
    natm = mol.natm
    ngrid = grids.coords.shape[0]
    atm_idx = grids.atm_idx
    weights = grids.weights
    q = np.zeros(ngrid)
    trans = np.zeros((natm, 3))
    deriv = 1 if xctype == 'LDA' else 2
    hidx = {(0, 0): 4, (0, 1): 5, (0, 2): 6, (1, 1): 7, (1, 2): 8, (2, 2): 9}
    for p0 in range(0, ngrid, blksize):
        p1 = min(p0 + blksize, ngrid)
        ao = ni.eval_ao(mol, grids.coords[p0:p1], deriv=deriv)      # (4|10, nb, nao)
        nb = p1 - p0
        # densities and their gradients: c0 = phi D, rho = c0.phi, d_j rho = 2 c0.d_j phi
        c0 = [ao[0] @ D for D in dms]
        r0 = [np.einsum('pi,pi->p', c, ao[0]) for c in c0]
        r1 = [2.0*np.einsum('pi,xpi->xp', c, ao[1:4]) for c in c0]
        if xctype == 'LDA':
            rho = np.array(r0)
        else:
            rho = np.array([np.vstack([r0[s][None, :], r1[s]]) for s in range(2)])
        exc, vxc = ni.eval_xc_eff(xc, rho, deriv=1, xctype=xctype)[:2]
        exc = np.asarray(exc).reshape(-1)
        vxc = np.asarray(vxc).reshape(2, -1, nb)
        q[p0:p1] = (r0[0] + r0[1])*exc
        gF = np.zeros((3, nb))
        for s in range(2):
            gF += vxc[s, 0]*r1[s]
            if xctype == 'GGA':
                # sum_k v_k d_j d_k rho = 2 sum_mu,nu D_mu,nu [d_j phi_mu psi_nu + phi_mu chi_j,nu],
                # psi = sum_k v_k d_k phi, chi_j = sum_k v_k d_j d_k phi; with eta = psi D
                v = vxc[s, 1:4]
                psi = np.einsum('kp,kpi->pi', v, ao[1:4])
                eta = psi @ dms[s]
                for j in range(3):
                    chi = sum(v[k][:, None]*ao[hidx[(min(j, k), max(j, k))]] for k in range(3))
                    gF[j] += 2.0*(np.einsum('pi,pi->p', eta, ao[j+1]) + np.einsum('pi,pi->p', c0[s], chi))
        ia = atm_idx[p0:p1]
        live = ia >= 0
        wf = weights[p0:p1][live, None]*gF[:, live].T                 # (nlive, 3)
        np.add.at(trans, ia[live], wf)
    return weight_response(grids, q) + trans


def weight_response(grids, q):
    """
    grad[G, :, c] = sum_p q[p, c] dw_p/dR_G for the integrands q (npts,)
    or (npts, nq) in grid order (padding rows ignored; grid points move
    with their parent atoms). Returns (natm, 3) or (natm, 3, nq).
    """
    if not isinstance(grids, GraciGrids) or grids._resp is None:
        raise RuntimeError('weight_response: the grid must be a built GraciGrids')
    d = grids._resp
    natm = d['R'].shape[1]
    npts = d['iatom'].size
    q = np.asarray(q, dtype=float)
    squeeze = (q.ndim == 1)
    q = np.asfortranarray(q.reshape(npts, -1))
    nq = q.shape[1]
    grad = np.zeros((3, natm, nq), order='F')
    _lib('grid_wderiv_contract', [natm, d['R'], d['Z'], d['codes']['adjust'], d['codes']['step'],
                                  d['codes']['mu'], d['codes']['cells'], d['codes']['weights'], npts,
                                  d['xyz'], d['iatom'], d['vol'], d['nnb'], d['nb'], nq, q, grad, 0])
    grad = np.ascontiguousarray(grad.transpose(1, 0, 2))
    return grad[:, :, 0] if squeeze else grad
