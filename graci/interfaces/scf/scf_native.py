"""
GRaCI's native SCF backend ($scf backend = native): the SCF iterations and
the Fock builds (DF-J, DF-K with occupied-orbital factors, XC quadrature,
DIIS, symmetry-blocked diagonalisation) run in the Fortran library libscf
(graci/dep/scf). PySCF supplies the basis, the DF tensor (built by the DF
object of Scf.pyscf_obj(), so that every post-SCF consumer sees the same
operator), the grid, the one-electron integrals, the parsed functional and
the symmetry-adapted orbitals.
"""

import os
import sys
import time
import numpy as np
from pyscf import gto
from pyscf import scf as pyscf_scf
from pyscf.dft import libxc
import graci.core.libs as libs

# cycles with exact spin Fock matrices before the open-row shortcut
EXACT_UNTIL = 3


def check_input(scf_obj, mf):
    """hard errors for the settings the native backend does not support"""
    unsupported = []
    if scf_obj.cosmo:
        unsupported.append('cosmo')
    if scf_obj.x2c:
        unsupported.append('x2c')
    if not scf_obj.mol.use_df:
        unsupported.append('conventional (non-density-fitted) integrals')
    if np.any(np.atleast_1d(scf_obj.chk_stable)):
        unsupported.append('chk_stable')
    if scf_obj.diag_method.lower() != 'cdiis':
        unsupported.append('diag_method = ' + str(scf_obj.diag_method))
    if scf_obj.lvl_shift is not None:
        unsupported.append('lvl_shift')
    if scf_obj.damp_fac is not None:
        unsupported.append('damp_fac')
    if scf_obj.xc != 'hf':
        ni = mf._numint
        xctype = ni._xc_type(scf_obj.xc)
        if xctype not in ('LDA', 'GGA'):
            unsupported.append('functional type ' + str(xctype))
        if ni.libxc.is_nlc(scf_obj.xc):
            unsupported.append('non-local correlation')
        omega = ni.rsh_and_hybrid_coeff(scf_obj.xc, spin=mf.mol.spin)[0]
        if abs(omega) > 1e-12:
            unsupported.append('range-separated functionals')
    if len(unsupported) > 0:
        sys.exit('\n $scf backend = native does not support: ' +
                 ', '.join(unsupported) + '\n')


def functional(scf_obj, mf):
    """libxc components (ids, families 1 LDA / 2 GGA, coefficients) and the
    exact-exchange fraction of the functional, as PySCF parses it"""
    if scf_obj.xc == 'hf':
        return [], [], [], 1.0
    _, fn = libxc.parse_xc(scf_obj.xc)
    ids = [int(f) for f, _ in fn]
    coefs = [float(c) for _, c in fn]
    fam = [1 if libxc.xc_type(f) == 'LDA' else 2 for f in ids]
    cx = float(mf._numint.rsh_and_hybrid_coeff(scf_obj.xc, spin=mf.mol.spin)[2])
    return ids, fam, coefs, cx


def df_tensor(mf):
    """the packed DF tensor (naux, nao(nao+1)/2) of the DF object of mf,
    built in memory (mandatory for the native backend)"""
    wdf = mf.with_df
    cderi = getattr(wdf, '_cderi', None)
    if isinstance(cderi, np.ndarray):
        return cderi
    total_mb = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / 1e6
    wdf.max_memory = max(wdf.max_memory, 0.9 * total_mb)
    wdf.build()
    cderi = getattr(wdf, '_cderi', None)
    if not isinstance(cderi, np.ndarray):
        sys.exit('\n $scf backend = native: the DF tensor does not fit in '
                 'memory\n')
    return cderi


def guess_factors(scf_obj, mf, pymol, guess, nalpha, nbeta):
    """start factor X0 (nao, k0a) and k0b: alpha density X0 X0^T, beta
    density from the first k0b columns (RKS: D = 2 X0 X0^T)"""
    S = pymol.intor_symmetric('int1e_ovlp')
    open_shell = nalpha > nbeta

    # a previous Scf object (another geometry, or a re-optimisation):
    # its MOs projected onto the present basis, as Scf.guess_dm
    if guess is not None:
        mol0 = guess.mol.mol_obj
        mo = pyscf_scf.addons.project_mo_nr2nr(mol0, guess.orbs, pymol)
        norm = np.einsum('pi,pi->i', mo, S @ mo)
        mo = mo / np.sqrt(norm)
        occ0 = np.asarray(guess.orb_occ)
        closed = mo[:, occ0 > 1.5]
        single = mo[:, (occ0 > 0.5) & (occ0 < 1.5)]
        if open_shell:
            return np.hstack([closed, single]), closed.shape[1]
        X = np.hstack([closed, single / np.sqrt(2.0)])
        return X, X.shape[1]

    # a restart file, else PySCF's guess density: alpha = beta = D/2,
    # natural orbitals weighted by sqrt(n/2)
    init_mos = scf_obj.make_chkfile_name(scf_obj.init_guess)
    if os.path.exists(init_mos):
        D = pyscf_scf.hf.init_guess_by_chkfile(pymol, init_mos)
    else:
        D = mf.get_init_guess(pymol, key=scf_obj.init_guess)
    D = np.asarray(D)
    if D.ndim == 3:
        D = D[0] + D[1]
    # natural orbitals: (S D S) c = n S c, i.e. S^1/2 D S^1/2 u = n u,
    # c = S^-1/2 u
    w, v = np.linalg.eigh(S)
    n, U = np.linalg.eigh((v * np.sqrt(w)) @ v.T @ D @ (v * np.sqrt(w)) @ v.T)
    C = (v / np.sqrt(w)) @ v.T @ U
    keep = n > 1e-10
    X = C[:, keep] * np.sqrt(0.5 * n[keep])
    return X, X.shape[1]


def run(scf_obj, pymol, guess):
    """
    one native SCF for the Scf object scf_obj (molecule pymol, optional
    guess Scf object). Returns a dict with 'converged' and 'info' (cycles,
    gradient norm, timings) and, if converged, mo_coeff, mo_occ,
    mo_energy, e_tot, fock_ao, hyb, naux and the DF tensor cderi
    """
    t_start = time.time()
    mf = scf_obj.pyscf_obj()
    check_input(scf_obj, mf)
    nao = pymol.nao_nr()
    nalpha, nbeta = pymol.nelec

    # one-electron operators and the DF tensor
    h = np.asfortranarray(mf.get_hcore(), dtype=float)
    S = np.asfortranarray(pymol.intor_symmetric('int1e_ovlp'), dtype=float)
    enuc = float(pymol.energy_nuc())
    cderi = df_tensor(mf)
    naux = int(mf.with_df.auxmol.nao_nr())
    t_data = time.time() - t_start

    # near-linear dependencies (all nao orbitals are kept downstream)
    smin = np.linalg.eigvalsh(S)[0]

    # functional and grid
    ids, fam, coefs, cx = functional(scf_obj, mf)
    nfunc = len(ids)
    nspin = 2 if nalpha > nbeta else 1
    if nfunc > 0:
        ncomp = 4 if 2 in fam else 1
        grids = mf.grids
        if grids.coords is None:
            grids.build()
        coords = np.ascontiguousarray(grids.coords, dtype=float)
        weights = np.ascontiguousarray(grids.weights, dtype=float)

    # symmetry-adapted orbitals
    if pymol.symmetry:
        so = np.asfortranarray(np.hstack(pymol.symm_orb), dtype=float)
        nso = np.array([b.shape[1] for b in pymol.symm_orb], dtype=np.int32)
        irrid = np.array(pymol.irrep_id, dtype=np.int32)
    else:
        so = np.asfortranarray(np.eye(nao))
        nso = np.array([nao], dtype=np.int32)
        irrid = np.array([0], dtype=np.int32)

    # basis in libcint's layout and the cart->sph tables
    atm = np.ascontiguousarray(pymol._atm, dtype=np.int32)
    bas = np.ascontiguousarray(pymol._bas, dtype=np.int32)
    env = np.ascontiguousarray(pymol._env, dtype=float)
    aoloc = np.ascontiguousarray(pymol.ao_loc_nr(), dtype=np.int32)
    lmax = int(bas[:, 1].max())
    tabs, off, pos = [], np.zeros(lmax + 1, dtype=np.int32), 0
    for l in range(lmax + 1):
        T = np.asarray(gto.cart2sph(l), dtype=float)
        off[l] = pos
        tabs.append(T.ravel(order='F'))
        pos += T.size
    c2s = np.concatenate(tabs)

    # start factor
    X0, k0b = guess_factors(scf_obj, mf, pymol, guess, nalpha, nbeta)
    X0 = np.asfortranarray(X0, dtype=float)

    # library session (freed at the end: the grid and the DF pointer are
    # specific to this geometry)
    ierr = libs.lib_func('scf_set_basis',
                         [pymol.natm, pymol.nbas, env.size, atm.ravel(), bas.ravel(), env, nao,
                          aoloc, bool(pymol.cart), lmax, c2s.size, c2s, off, '', 0])
    if ierr != 0:
        sys.exit('\n scf_set_basis: error ' + str(ierr) + '\n')
    if nfunc > 0:
        ierr = libs.lib_func('scf_set_grid',
                             [weights.size, coords.T, weights, ncomp, float(scf_obj.mem_budget), '', 0])
        if ierr != 0:
            sys.exit('\n scf_set_grid: error ' + str(ierr) + '\n')
    ids_a = np.array(ids if nfunc > 0 else [0], dtype=np.int32)
    fam_a = np.array(fam if nfunc > 0 else [0], dtype=np.int32)
    coef_a = np.array(coefs if nfunc > 0 else [0.0], dtype=float)
    ierr = libs.lib_func('scf_set_operator',
                         [nao, cderi.shape[0], h, S, enuc, cderi.T, cx, nfunc, ids_a, fam_a, coef_a,
                          nspin, nso.size, nso, irrid, so, 0])
    if ierr != 0:
        sys.exit('\n scf_set_operator: error ' + str(ierr) + '\n')

    conv_tol_grad = scf_obj.conv_tol_grad
    if conv_tol_grad is None:
        conv_tol_grad = np.sqrt(scf_obj.conv_tol)
    diis_space = scf_obj.diis_space if scf_obj.diis_space > 0 else 14
    diis_start = scf_obj.diis_start if scf_obj.diis_start is not None else 1
    opts = np.array([scf_obj.conv_tol, conv_tol_grad, scf_obj.max_iter, diis_space,
                     diis_start, 0.0, EXACT_UNTIL, 1.0 if scf_obj.verbose else 0.0])
    C = np.zeros((nao, nao), order='F')
    e = np.zeros(nao)
    occ = np.zeros(nao)
    irr = np.zeros(nao, dtype=np.int32)
    Fa = np.zeros((nao, nao), order='F')
    Fb = np.zeros((nao, nao), order='F')
    ener = np.zeros(8)
    info = np.zeros(8, dtype=np.int32)
    libs.lib_func('scf_reset_timings', [])
    t0 = time.time()
    (C, e, occ, irr, Fa, Fb, ener, info, ierr) = libs.lib_func(
        'scf_run', [nao, nalpha, nbeta, X0.shape[1], k0b, X0, opts, C, e, occ, irr, Fa, Fb,
                    ener, info, 0])
    t_run = time.time() - t0
    t = libs.lib_func('scf_timings', [np.zeros(16)])
    libs.lib_func('scf_free', [])
    if ierr != 0:
        sys.exit('\n scf_run: error ' + str(ierr) + '\n')

    info_d = {'converged'      : bool(info[0] == 1),
              'cycles'         : int(info[1]),
              'shortcut_cycles': int(info[2]),
              'occ_rule_ok'    : bool(info[3] == 1),
              'gradient_norm'  : float(ener[6]),
              'delta_e'        : float(ener[7]),
              'time_run'       : float(t_run),
              'time_data'      : float(t_data),
              'time_j'         : float(t[0]),
              'time_k'         : float(t[1]),
              'time_xc'        : float(t[3]),
              'ao_cache'       : bool(t[13] > 0.5),
              'min_overlap_eig': float(smin)}
    if not info_d['converged']:
        return {'converged': False, 'info': info_d}

    if scf_obj.xc == 'hf':
        hyb = (0., 0., 1.)
    else:
        hyb = tuple(float(x) for x in
                    mf._numint.rsh_and_hybrid_coeff(scf_obj.xc, spin=pymol.spin))

    return {'converged': True,
            'mo_coeff' : np.array(C),
            'mo_occ'   : np.array(occ),
            'mo_energy': np.array(e),
            'e_tot'    : float(ener[0]),
            'fock_ao'  : np.array([Fa, Fb]),
            'hyb'      : hyb,
            'naux'     : naux,
            'cderi'    : cderi,
            'info'     : info_d}
