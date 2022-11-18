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

@timing.timed
def rdm(ci_method):
    """Calculation of the MRCI 1-RDMs for all states"""

    # number of molecular orbitals
    nmo = ci_method.nmo

    # bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # MRCI configuration file scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)
    
    # MRCI eigenvector scratch file numbers
    ci_ciunits = np.array(mrci_wfn.ci_units['adiabatic'], dtype=int)
    
    # 1-RDMs for all irreps
    dmat_all = []
    
    # Loop over irreps
    for irr in range(ci_method.n_irrep()):
        states = [n for n in range(ci_method.n_states_sym(irr))]
        
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
        dmat = libs.lib_func('density_mrci', args)
 
        # Add the 1-RDMs to the list
        dmat_all.append(np.reshape(dmat, (nmo, nmo, nstates),
                                   order='F'))

        ## Temporary: save the 1-RDMs to disk
        #for i in range(nstates):
        #    f = 'xdm_'+str(irr)+str(i+1)+'_' \
        #        +str(irr)+str(i+1)
        #    np.save(f, dmat_all[irr][:,:,i])
        
    return dmat_all

