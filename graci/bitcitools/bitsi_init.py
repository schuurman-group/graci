"""
Module for loading/finalizing bitci library
"""

import numpy as np
import graci.core.libs as libs

#
def init(si_method):
    """Initialize the bitci library"""

    bra = si_method.bra_obj
    ket = si_method.ket_obj

    # if the number of mos is different between bra and ket, end
    if bra.scf.nmo != ket.scf.nmo:
        sys.exit('bra.scf.nmo != ket.scf.nmo init_bitsi')

    # set all variables that have to be passed to bitsi_initialise
    # (note that the pgrp uses Fortran indexing)
    multBra = bra.scf.mol.mult
    multKet = ket.scf.mol.mult
    nelBra   = bra.scf.mol.nel
    nelKet   = ket.scf.mol.nel
    nmo      = bra.scf.nmo

    # not sure what this is about...
    if bra.scf.mol.sym_indx <= 0:
        pgrp = 1
    else:
        pgrp = bra.scf.mol.sym_indx + 1

    # call to bitsi_initialise
    args = (multBra, multKet, nelBra, nelKet, nmo, pgrp)
    libs.lib_func('bitsi_initialise', args)

    return

def finalize():
    """finalize the bitci library"""

    args = ()
    libs.lib_func('bitsi_finalise', args)

    return

