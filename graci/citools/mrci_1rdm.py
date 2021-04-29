"""
Module for the calculation of MRCI one-electron reduced
density matrices
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

def rdm(ci, scf):
    """Calculation of the MRCI 1-RDMs for all states"""

    # MRCI configuration file scratch file numbers
    ci_confunits = np.array(ci.mrci_wfn.conf_units, dtype=int)
    
    # MRCI eigenvector scratch file numbers
    ci_ciunits = np.array(ci.mrci_wfn.ci_units, dtype=int)
    
    # No. MOs
    nmo = scf.nmo

    # 1-RDMs for all irreps
    dmat_all = []
    
    # Loop over irreps
    for irr in range(len(ci.nstates)):
        states = [n for n in range(ci.nstates[irr])]
        
        # States for which the 1-RDMs are required (note that bitCI
        # uses Fortran indexing here)
        state_indx = np.array([n+1 for n in states], dtype=int)

        # No. states
        nstates = len(states)
        
        # 1-RDM array
        dmat = np.zeros((nmo*nmo*nstates), dtype=np.float64)
        
        # Compute the 1-RDM for all states in this irrep
        args = (irr, nstates, state_indx, dmat, ci_confunits,
                ci_ciunits)
        (dmat) = libs.lib_func('density_mrci', args)

        # Add the 1-RDMs to the list
        dmat_all.append(np.reshape(dmat, (nmo, nmo, nstates),
                                   order='F'))

    # Save the 1-RDMs
    ci.mrci_wfn.set_dmat(dmat_all)

    # Test
    #occ, vec = np.linalg.eigh(dmat_all[1][:,:,0])
    #indx = np.abs(occ).argsort()[::-1]
    #occ  = occ[indx]
    #vec  = vec[:, indx]
    #for i in range(nmo):
    #   print(i,occ[i])
    
    return

