"""
Module for the calculation of wave function overlaps
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.core.molecule as molecule

@timing.timed
def extract(bra, ket, rep='adiabatic'):
    """
    Extraction of the determinant representation of the bra and
    ket MRCI wave functions
    """

    # number of irreps
    nirr_bra = bra.n_irrep()
    nirr_ket = ket.n_irrep()
    
    # bitci wave function objects
    wfn_bra = bra.bitci_mrci()
    wfn_ket = ket.bitci_mrci()

    # bitwf wave function scratch file numbers
    wfunits_bra = [0 for i in range(nirr_bra)]
    wfunits_ket = [0 for i in range(nirr_ket)]
    
    # Extract the bra wave functions
    for irr in bra.irreps_nonzero():
    
        # number of states for this irrep
        nstates = bra.n_states_sym(irr)
        if nstates == 0:
            continue
    
        # bitci configuration and eigenvector scratch files
        conf_file = wfn_bra.conf_name[rep][irr]
        vec_file  = wfn_bra.ci_name[rep][irr]
    
        # bitwf determinant wave function file number
        wf_scr = 0
        
        # extract the determinant representation wave functions
        # for this irrep
        args   = (irr, conf_file, vec_file, nstates, 'bra', wf_scr)
        wf_scr = libs.lib_func('detwf', args)

        wfunits_bra[irr] = wf_scr
        
    # Extract the ket wave functions
    for irr in ket.irreps_nonzero():
    
        # number of states for this irrep
        nstates = ket.n_states_sym(irr)
        if nstates == 0:
            continue
    
        # bitci configuration and eigenvectos scratch files
        conf_file = wfn_ket.conf_name[rep][irr]
        vec_file  = wfn_ket.ci_name[rep][irr]
    
        # bitwf determinant wave function file number
        wf_scr = 0
        
        # extract the determinant representation wave functions
        # for this irrep
        args   = (irr, conf_file, vec_file, nstates, 'ket', wf_scr)
        wf_scr = libs.lib_func('detwf', args)
        wfunits_ket[irr] = wf_scr

    return wfunits_bra, wfunits_ket

@timing.timed
def overlap(bra, ket, bra_wfunit, ket_wfunit, overlap_list, norm_thresh,
            det_thresh):

    if bra.scf.mol.sym_indx == ket.scf.mol.sym_indx:
        overlap = overlap_same_sym(bra, ket, bra_wfunit, ket_wfunit,
                                   overlap_list, norm_thresh,
                                   det_thresh)
    else:
        overlap = overlap_diff_sym(bra, ket, bra_wfunit, ket_wfunit,
                                   overlap_list, norm_thresh,
                                   det_thresh)
        
    
    return overlap

def overlap_same_sym(bra, ket, bra_wfunit, ket_wfunit, overlap_list,
                     norm_thresh, det_thresh):
    """
    Calculation of the overlaps between all pairs of states
    in overlap_list using the determinant representation
    of the wave functions. Makes use of the bra and ket point
    groups being the same
    """

    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()

    # bitci bra and ket mrci wfn object
    wfn_bra = bra.bitci_mrci()
   
    # bitci ket mrci wfn object
    wfn_ket = ket.bitci_mrci()    

    # frozen core orbitals
    ncore_el = np.array([molecule.atom_ncore[n] for n in
                         [molecule.atom_name.index(m)
                          for m in bra.scf.mol.asym]])
    ncore_mo = int(np.sum(ncore_el/2))

    # fill in the core MO indices
    # (bitX libraries have the MOs in energy order and
    # use Fortran indexing)
    icore = np.arange(0, ncore_mo, 1) + 1
    ncore = icore.size

    # deleted core orbital flag
    # (hard-wired to True for now)
    delete_core = True
    
    # overlaps for all irreps
    overlap = [[[] for i in range(nirr_bra)] for j in range(nirr_ket)]

    # fill in the overlaps that are non-zero by symmetry
    for irr in bra.irreps_nonzero():

        # number of pairs of bra and ket states
        npairs = len(overlap_list[irr][irr])
        if npairs == 0:
            continue

        # pairs of states for this bra irrep and ket irrep
        # bitwf uses Fortran indexing for these, hence the +1
        overlap_pairs = 1 + np.reshape(
            np.array(overlap_list[irr][irr], dtype=int),
            (2*npairs), order='F')

        # total number of bra and ket roots for this irrep
        bra_tot = bra.n_states_sym(irr)
        ket_tot = ket.n_states_sym(irr)

        # wave function overlap array
        sij = np.zeros((npairs), dtype=np.float64)

        # bitwf wave function file numbers
        bra_unit = bra_wfunit[irr]
        ket_unit = ket_wfunit[irr]

        # compute the overlaps for all requested pairs of states
        args = (irr, irr, bra_tot, ket_tot, npairs, overlap_pairs,
                bra_unit, ket_unit, norm_thresh, det_thresh,
                ncore, icore, delete_core, sij)

        sij  = libs.lib_func('detoverlap', args)

        overlap[irr][irr] = sij
    
    # fill in the overlaps that are zero by symmetry
    for bra_irr in bra.irreps_nonzero():
        for ket_irr in ket.irreps_nonzero():
            if bra_irr != ket_irr:
                npairs = len(overlap_list[bra_irr][ket_irr])
                if npairs == 0:
                    continue
                overlap[bra_irr][ket_irr] = np.zeros((npairs),
                                                     dtype=np.float64)
                
    return overlap

def overlap_diff_sym(bra, ket, bra_wfunit, ket_wfunit, overlap_list,
                     norm_thresh, det_thresh):
    """
    Calculation of the overlaps between all pairs of states
    in overlap_list using the determinant representation
    of the wave functions.
    """

    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()

    # bitci bra and ket mrci wfn object
    wfn_bra = bra.bitci_mrci()
   
    # bitci ket mrci wfn object
    wfn_ket = ket.bitci_mrci()    

    # frozen core orbitals
    ncore_el = np.array([molecule.atom_ncore[n] for n in
                         [molecule.atom_name.index(m)
                          for m in bra.scf.mol.asym]])
    ncore_mo = int(np.sum(ncore_el/2))

    # fill in the core MO indices
    # (bitX libraries have the MOs in energy order and
    # use Fortran indexing)
    icore = np.arange(0, ncore_mo, 1) + 1
    ncore = icore.size

    # deleted core orbital flag
    # (hard-wired to True for now)
    delete_core = True
    
    # overlaps for all irreps
    overlap = [[[] for i in range(nirr_ket)] for j in range(nirr_bra)]

    # fill in the overlaps
    for bra_irr in bra.irreps_nonzero():

        for ket_irr in ket.irreps_nonzero():

            # number of pairs of bra and ket states
            npairs = len(overlap_list[bra_irr][ket_irr])
            if npairs == 0:
                continue

            # pairs of states for this bra irrep and ket irrep
            # bitwf uses Fortran indexing for these, hence the +1
            overlap_pairs = 1 + np.reshape(
                np.array(overlap_list[bra_irr][ket_irr], dtype=int),
                (2*npairs), order='F')

            # total number of bra and ket roots for this irrep
            bra_tot = bra.n_states_sym(bra_irr)
            ket_tot = ket.n_states_sym(ket_irr)

            # wave function overlap array
            sij = np.zeros((npairs), dtype=np.float64)

            # bitwf wave function file numbers
            bra_unit = bra_wfunit[bra_irr]
            ket_unit = ket_wfunit[ket_irr]

            # compute the overlaps for all requested pairs of states
            args = (bra_irr, ket_irr, bra_tot, ket_tot, npairs, overlap_pairs,
                    bra_unit, ket_unit, norm_thresh, det_thresh,
                    ncore, icore, delete_core, sij)

            sij  = libs.lib_func('detoverlap', args)

            overlap[bra_irr][ket_irr] = sij
            
    return overlap
