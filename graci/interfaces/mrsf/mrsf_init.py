"""
Module for initialising/finalising the mrsf (MRSF-TDDFT) library
"""

import numpy as np
import graci.core.libs as libs

#
def init(ci_method, fock_mo):
    """
    Initialise the mrsf library for the Mrsftddft object ci_method.
    fock_mo: (2, nmo, nmo) array of the alpha/beta spin Fock matrices
    of the reference in the MO basis of ci_method.mos
    """

    nmo   = ci_method.nmo
    nel   = ci_method.nel
    mosym = np.array(ci_method.mosym, dtype=np.int64)
    occ   = np.array(ci_method.occ_ref, dtype=float)
    moen  = np.array(ci_method.emo, dtype=float)
    fa    = np.reshape(np.asarray(fock_mo[0], dtype=float), (nmo*nmo), order='F')
    fb    = np.reshape(np.asarray(fock_mo[1], dtype=float), (nmo*nmo), order='F')
    chf   = float(ci_method.chf_eff)
    spc   = np.array(ci_method.spc_eff, dtype=float)
    escf  = float(ci_method.scf.energy)

    if ci_method.scf.mol.sym_indx <= 0:
        pgrp = 1
    else:
        pgrp = ci_method.scf.mol.sym_indx + 1

    args = (nmo, nel, int(ci_method.scf.mult), mosym, occ, moen, fa, fb,
            chf, spc, escf, pgrp, ci_method.label, ci_method.verbose)
    libs.lib_func('mrsf_initialise', args)

    return

#
def init_ints(ci_method, eri_file):
    """load the density-fitted MO integrals from eri_file"""

    args = ('df', ci_method.precision, eri_file, ci_method.vv_storage,
            float(ci_method.mem_budget))
    libs.lib_func('mrsf_int_initialise', args)

    return

#
def set_extended(flag):
    """select the extended (EMRSF) response space; must precede init()"""

    libs.lib_func('mrsf_set_extended', (bool(flag),))

    return

#
def set_frozen_core(nfc):
    """number of doubly occupied MOs (the lowest by orbital energy)
    excluded from the response space; must precede init()"""

    libs.lib_func('mrsf_set_frozen_core', (int(nfc),))

    return

#
def set_exchange(mode):
    """exchange term of the sigma vectors: 'auto' (explicit (ij|ab) tensor
    when it fits mem_budget, else the DF plane sweep), 'tensor', 'df';
    must precede init_ints()"""

    m = {'auto': 0, 'tensor': 1, 'df': 2}[str(mode).lower()]
    libs.lib_func('mrsf_set_exchange', (m,))

    return

#
def init_ext(ci_method, fdft_mo, eri_file):
    """extended-method data: fdft_mo (nmo, nmo) closed-shell KS matrix of
    the configuration G in the MO basis, coupling scale ci_method.ccp_eff,
    kernel flag ci_method.use_kernel; after init_ints() and, when the
    kernel is used, after the grid initialisation (mrsf_xc.XCGrid)"""

    nmo = int(ci_method.nmo)
    fdft = np.asfortranarray(np.asarray(fdft_mo, dtype=float).reshape(nmo, nmo))
    args = (nmo, fdft, float(ci_method.ccp_eff), bool(ci_method.use_kernel), eri_file)
    libs.lib_func('mrsf_ext_initialise', args)

    return

#
def ext_dims():
    """returns (ncol, ncv): columns of the response vector and CV slots"""

    return libs.lib_func('mrsf_get_ext_dims', (0, 0))

#
def dims():
    """returns (nocca, nvirb, xdim, naux)"""

    return libs.lib_func('mrsf_get_dims', (0, 0, 0, 0))

#
def finalize():
    """finalise the mrsf library"""

    libs.lib_func('mrsf_finalise', ())

    return

#
def report_timings():
    """print the library timings of the current section and reset the
    counters (the loaded integrals are kept)"""

    libs.lib_func('mrsf_report_timings', ())

    return
