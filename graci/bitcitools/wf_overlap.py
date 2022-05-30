"""
Module for the calculation of wave function overlaps
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing

@timing.timed
def extract(bra, ket):
    """
    Extraction of the determinant representation of the bra and
    ket MRCI wave functions
    """

    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()
    
    # bitci wave function objects
    wfn_bra = bra.bitci_mrci()
    wfn_ket = bra.bitci_mrci()

    # Extract the bra wave functions
    for irr in range(nirr_bra):

        # number of states for this irrep
        nstates = bra.n_states_sym(irr)
        if nstates == 0:
            continue

        # bitci configuration and eigenvectos scratch files
        conf_file = wfn_bra.conf_name[irr]
        vec_file  = wfn_bra.ci_name[irr]

        # bitwf determinant wave function file number
        wf_scr = 0
        
        # extract the determinant representation wave functions
        # for this irrep
        args   = (conf_file, vec_file, nstates, 'bra', wf_scr)
        wf_scr = libs.lib_func('detwf', args)
        
    return
    
@timing.timed
def overlap(bra, ket, overlap_list):
    """
    Calculation of the overlaps between all pairs of states
    in overlap_list using the determinant representation
    of the wave functions
    """

    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()

    # bitci bra mrci wfn object
    bra_wfn = bra.bitci_mrci()
   
    # bitci ket mrci wfn object
    ket_wfn = ket.bitci_mrci()
    
    # overlaps for all irreps
    overlap = [[[] for i in range(nirr_bra)] for j in range(nirr_ket)]
    
    return overlap
