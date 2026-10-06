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
