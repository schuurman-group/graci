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
    wfn_ket = ket.bitci_mrci()

    # bitwf wave function scratch file numbers
    wfunits_bra = []
    wfunits_ket = []
    
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
        args   = (irr, conf_file, vec_file, nstates, 'bra', wf_scr)
        wf_scr = libs.lib_func('detwf', args)
        wfunits_bra.append(wf_scr)
        
    # Extract the ket wave functions
    for irr in range(nirr_ket):
    
        # number of states for this irrep
        nstates = ket.n_states_sym(irr)
        if nstates == 0:
            continue
    
        # bitci configuration and eigenvectos scratch files
        conf_file = wfn_ket.conf_name[irr]
        vec_file  = wfn_ket.ci_name[irr]
    
        # bitwf determinant wave function file number
        wf_scr = 0
        
        # extract the determinant representation wave functions
        # for this irrep
        args   = (irr, conf_file, vec_file, nstates, 'ket', wf_scr)
        wf_scr = libs.lib_func('detwf', args)
        wfunits_ket.append(wf_scr)

    return wfunits_bra, wfunits_ket
    
@timing.timed
def overlap(bra, ket, bra_wfunit, ket_wfunit, overlap_list):
    """
    Calculation of the overlaps between all pairs of states
    in overlap_list using the determinant representation
    of the wave functions
    """

    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()

    # bitci bra and ket mrci wfn object
    wfn_bra = bra.bitci_mrci()
   
    # bitci ket mrci wfn object
    wfn_ket = ket.bitci_mrci()    
    
    # overlaps for all irreps
    overlap = [[[] for i in range(nirr_bra)] for j in range(nirr_ket)]

    # loop over irreps: note that we are enforcing equal bra and ket
    # point groups and that the wave function overlaps will be zero
    # unless the bra and ket irreps are the same
    for irr in range(bra.n_irrep()):

        # pairs of states for this bra irrep and ket irrep
        # bitwf uses Fortran indexing for these, hence the +1
        npairs = len(overlap_list[irr][irr])

        if npairs == 0:
            continue

        overlap_pairs = 1 + np.reshape(
            np.array(overlap_list[irr][irr], dtype=int),
            (2*npairs), order='F')

        # total number of bra and ket roots for this irrep
        bra_tot = bra.n_states_sym(irr)
        ket_tot = ket.n_states_sym(irr)

        # wave function overlap array
        sij = np.zeros((2*npairs), dtype=np.float64)

        # bitwf wave function file numbers
        bra_unit = bra_wfunit[irr]
        ket_unit = ket_wfunit[irr]

        # compute the overlaps for all requested pairs of states
        args = (irr, bra_tot, ket_tot, npairs, overlap_pairs,
                bra_unit, ket_unit, sij)
        sij  = libs.lib_func('detoverlap', args)

        # Reshape sij to shape (2,npairs)?
        
    sys.exit()
        
    return overlap