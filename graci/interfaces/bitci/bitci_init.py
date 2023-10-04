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
        
    ham   = ci_method.hamiltonian
        
    label = ci_method.label

    verbose = ci_method.verbose
    
    # call to bitci_initialise
    args = (imult, nel, nmo, mosym, moen, pgrp, escf, cvs_flag,
            ham, label, verbose)
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
