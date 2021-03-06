"""
Module for the calculation of MRCI one-electron reduced
transition density matrices
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing

@timing.timed
def tdm(bra, ket, trans_list):
    """Calculation of the MRCI 1-TDMs for all pairs of 
      states in trans_list"""

    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()

    # number of molecular orbitals
    nmo = bra.scf.nmo

    # bitci bra mrci wfn object
    bra_wfn = bra.bitci_mrci()
   
    # bitci ket mrci wfn object
    ket_wfn = ket.bitci_mrci()

    # 1-TDMs for all irreps
    rho = [[[] for i in range(nirr_bra)] for j in range(nirr_ket)]
    
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
            bra_tot = bra.n_states_sym(bra_irr)
            ket_tot = ket.n_states_sym(ket_irr)

            # 1-TDM array
            rhoij = np.zeros((nmo*nmo*npairs), dtype=np.float64)

            # configuration and eigenvector files: bra wfn
            bra_conf = bra_wfn.conf_name[bra_irr]
            bra_vec  = bra_wfn.ci_name[bra_irr]

            # configuration and eigenvector files: ket wfn
            ket_conf = ket_wfn.conf_name[ket_irr]
            ket_vec  = ket_wfn.ci_name[ket_irr]

            # Compute the 1-TDMs for all states in this irrep
            args = (bra_irr, ket_irr, bra_tot, ket_tot, npairs, 
                    tdm_pairs, rhoij,  bra_conf, bra_vec, ket_conf, 
                    ket_vec)
            rhoij = libs.lib_func('transition_density_mrci', args)

            # Add the 1-TDMs to the list
            rho[bra_irr][ket_irr] = np.reshape(rhoij, 
                                          (nmo, nmo, npairs),
                                           order='F')
            
            ## Temporary: save the 1-TDMs to disk
            #pairs = np.reshape(tdm_pairs, (npairs, 2), order='F')
            #for n in range(npairs):
            #    ibra = pairs[n][0]
            #    iket = pairs[n][1]
            #    if bra_irr == ket_irr and ibra == iket:
            #        continue
            #    f = 'xdm_'+str(bra_irr)+str(ibra)+'_' \
            #        +str(ket_irr)+str(iket)
            #    np.save(f, rho[bra_irr][ket_irr][:,:,n])

    return rho

