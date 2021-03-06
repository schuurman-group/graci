"""
Module for loading/finalizing bitci library
"""

import numpy as np
import graci.core.libs as libs
import sys as sys

#
def init(ci_method):
    """Initialize the bitci library"""

    # (note that the pgrp and iham variables use Fortran indexing)
    imult = ci_method.mult
    nel   = ci_method.nel 
    nmo   = ci_method.scf.nmo
    mosym = np.array(ci_method.scf.orb_sym)
    moen  = np.array(ci_method.scf.orb_ener)

    if ci_method.scf.mol.sym_indx <= 0:
        pgrp = 1
    else:
        pgrp = ci_method.scf.mol.sym_indx + 1

    escf  = ci_method.scf.energy

    ham   = ci_method.hamiltonian
        
    label = ci_method.label

    # call to bitci_initialise
    args = (imult, nel, nmo, mosym, moen, pgrp, escf, ham, label)
    libs.lib_func('bitci_initialise', args)

    return

def finalize():
    """finalize the bitci library"""

    args = ()
    libs.lib_func('bitci_finalise', args)

    return
