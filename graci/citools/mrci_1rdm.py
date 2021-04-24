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

def rdm(ci, states, irr, scf):
    """Calculation of the MRCI 1-RDMs for a single irrep"""

    # States for which the 1-RDMs are required (note that bitCI
    # uses Fortran indexing here)
    state_indx = np.array([n+1 for n in states], dtype=int)
    
    # Name of the MRCI configuration file
    conf_name = ci.mrci_wfn.conf_name[irr]

    # Name of the MRCI eigenvector file
    ci_name = ci.mrci_wfn.ci_name[irr]

    # No. states
    nstates = len(states)

    # No. MOs
    nmo = scf.nmo
    
    # Density matrices for all states
    dmat = np.zeros((nmo*nmo*nstates), dtype=np.float64)
    
    # Compute the 1-RDM for all states in this irrep
    args = (irr, nstates, state_indx, dmat, conf_name, ci_name)
    (dmat) = libs.lib_func('density_mrci', args)
    
    return

