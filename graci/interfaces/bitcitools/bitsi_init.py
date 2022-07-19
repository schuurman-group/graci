"""
Module for loading/finalizing bitci library
"""

import sys as sys
import numpy as np
import graci.core.libs as libs

#
def init(bra, ket, calctype, verbose):
    """Initialize the bitsi library"""

    # if the number of mos is different between bra and ket, end
    if bra.scf.nmo != ket.scf.nmo:
        sys.exit('bra.scf.nmo != ket.scf.nmo init_bitsi')

    # set all variables that have to be passed to bitsi_initialise
    # (note that the pgrp uses Fortran indexing)
    multBra = bra.mult
    multKet = ket.mult
    nelBra  = bra.nel
    nelKet  = ket.nel
    nmo     = bra.scf.nmo

    # catch instances of C1 symmetry
    if bra.scf.mol.sym_indx <= 0:
        pgrp = 1
    else:
        pgrp = bra.scf.mol.sym_indx + 1

    # call to bitsi_initialise
    args = (multBra, multKet, nelBra, nelKet, nmo, pgrp, calctype,
            verbose)
    libs.lib_func('bitsi_initialise', args)

    return

def finalize():
    """finalize the bitci library"""

    args = ()
    libs.lib_func('bitsi_finalise', args)

    return

