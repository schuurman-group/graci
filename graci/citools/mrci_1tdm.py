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

def tdm(si_method):
    """Calculation of the MRCI 1-RDMs for all states"""

    # number of molecular orbitals
    nmo = ci_method.scf.nmo

    # bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # MRCI configuration file scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)
    
    # MRCI eigenvector scratch file numbers
    ci_ciunits = np.array(mrci_wfn.ci_units, dtype=int)
    
    # 1-RDMs for all irreps
    rho_all = []
    
    # Loop over irreps
    for irr in range(ci_method.n_irrep()):
        states = [n for n in range(ci_method.n_states(irr))]
        
        # States for which the 1-RDMs are required (note that bitCI
        # uses Fortran indexing here)
        state_indx = np.array([n+1 for n in states], dtype=int)

        # No. states
        nstates = len(states)
        
        # 1-RDM array
        rho = np.zeros((nmo*nmo*n_bra*n_ket), dtype=np.float64)
        
        # Compute the 1-RDM for all states in this irrep
        args = (irr, nstates, state_indx, dmat, ci_confunits,
                ci_ciunits)
        rho = libs.lib_func('density_mrci', args)

        # Add the 1-RDMs to the list
        rho_all.append(np.reshape(rho, (nmo, nmo, n_bra, n_ket),
                                   order='F'))

    return rho_all

