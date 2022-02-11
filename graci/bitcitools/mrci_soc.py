"""
Module for the calculation of MRCI triplet transition density
matrices <psi_m|T_ij^(1,k)|psi'_n>, k=0,+1
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing

@timing.timed
def redmat(bra, ket, trans_list):
    """Calculation of the MRCI reduced matrix elements 
       for all pairs of states in trans_list"""
    
    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()

    # number of molecular orbitals
    nmo = bra.scf.nmo

     # bitci bra mrci wfn object
    bra_wfn = bra.bitci_mrci()
   
    # bitci ket mrci wfn object
    ket_wfn = ket.bitci_mrci()

    # Triplet TDMs for all irreps
    umat = [[[] for i in range(nirr_bra)] for j in range(nirr_ket)]

    # Loop over pairs of irreps for the initial and final manifolds
    for ket_irr in range(nirr_ket):
        for bra_irr in range(nirr_bra):

            # pairs of states for this bra irrep and ket irrep
            # bitsi uses Fortran indexing for these, hence the +1
            npairs = len(trans_list[bra_irr][ket_irr])

            if npairs == 0:
                continue

            tdm_pairs = 1 + np.reshape(
                np.array(trans_list[bra_irr][ket_irr], dtype=int),
                (2*npairs), order='F')

            # total number of bra and ket roots for this irrep
            bra_tot = bra.n_state_sym(bra_irr)
            ket_tot = ket.n_state_sym(ket_irr)

            # Triplet TDM array
            uij = np.zeros((nmo*nmo*npairs), dtype=np.float64)

            # configuration and eigenvector files: bra wfn
            bra_conf = bra_wfn.conf_name[bra_irr]
            bra_vec  = bra_wfn.ci_name[bra_irr]

            # configuration and eigenvector files: ket wfn
            ket_conf = ket_wfn.conf_name[ket_irr]
            ket_vec  = ket_wfn.ci_name[ket_irr]

            # Compute the reduced matrix elements for all states in
            # this irrep
            args = (bra_irr, ket_irr, bra_tot, ket_tot, npairs, 
                    tdm_pairs, uij,  bra_conf, bra_vec, ket_conf, 
                    ket_vec)
            uij  = libs.lib_func('redmat_mrci', args)

            # Add the reduced matrix elements to the list
            umat[bra_irr][ket_irr] = \
                np.reshape(uij, (nmo, nmo, npairs), order='F')
            
    return umat

@timing.timed
def clebsch_gordan(bra, ket):
    """Retrieval of all the Clebsch-Gordan coefficients needed to compute
       the SOC matrix elements between a given pair of sets of bra and ket
       MRCI wave functions"""

    # bra and ket spin multiplicities
    mult_bra = bra.mult
    mult_ket = ket.mult
        
    # array of Clebsh-Gordan coefficients
    cgcoe = np.zeros((3*mult_ket, mult_bra), dtype=np.float64)
    
    # fetch the coefficients
    
    return cgcoe
