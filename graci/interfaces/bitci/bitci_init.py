"""
Module for loading/finalizing bitci library
"""

import numpy as np
import graci.core.libs as libs
import sys as sys

#
def init(ci_method):
    """Initialize the bitci library"""

    # (note that the pgrp variable uses Fortran indexing)
    imult = ci_method.mult
    nel   = ci_method.nel 
    nmo   = ci_method.nmo
    mosym = np.array(ci_method.mosym)
    moen  = np.array(ci_method.emo)

    if ci_method.scf.mol.sym_indx <= 0:
        pgrp = 1
    else:
        pgrp = ci_method.scf.mol.sym_indx + 1

    escf  = ci_method.scf.energy

    # CVS core MO flags
    cvs_flag = np.zeros(nmo, dtype=int)
    for i in ci_method.icvs:
        cvs_flag[i-1] = 1

    # Degenerate partner index array (1-indexed, Fortran convention).
    # degen_group uses group IDs (both partners share the same positive int, 0=singleton).
    # degen_orbs[i] = 1-indexed canonical MO index of the partner of MO i, or 0.
    degen_grp = ci_method.scf.degen_group
    degen_orbs = np.zeros(nmo, dtype=np.int32)
    if degen_grp is not None:
        for gid in np.unique(degen_grp):
            if gid == 0:
                continue
            members = np.where(degen_grp == gid)[0]
            if len(members) == 2:
                i, j = members
                if i >= nmo or j >= nmo:
                    continue
                degen_orbs[i] = j + 1
                degen_orbs[j] = i + 1

    ham   = ci_method.hamiltonian

    label = ci_method.label

    verbose = ci_method.verbose

    # call to bitci_initialise
    args = (imult, nel, nmo, mosym, moen, pgrp, escf, cvs_flag,
            degen_orbs, ham, label, verbose)
    libs.lib_func('bitci_initialise', args)

    # optional overriding of the Hamiltonian parameters
    if ci_method.hparam is not None:
        nparam = ci_method.hparam.size
        args   = (nparam, ci_method.hparam)
        libs.lib_func('override_hparam', args)

    return

def finalize():
    """finalize the bitci library"""

    args = ()
    libs.lib_func('bitci_finalise', args)

    return
