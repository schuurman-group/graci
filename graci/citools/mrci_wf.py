"""
Module for the writing of the determinant expansions of the
MRCI wave functions to disk
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

def extract_wf(ci_method):
    """Extraction of the determinant expansions of the 
    MRCI wave functions for all states"""

    # bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # MRCI configuration file scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)
    
    # MRCI eigenvector scratch file numbers
    ci_ciunits = np.array(mrci_wfn.ci_units, dtype=int)

    # Loop over irreps
    for irr in range(ci_method.n_irrep()):
        states = [n for n in range(ci_method.n_state_sym(irr))]
        
        # States for which the WFs are required (note that bitCI
        # uses Fortran indexing here)
        state_indx = np.array([n+1 for n in states], dtype=int)

        # No. states
        nstates = len(states)

        # Name of the HDF5 output file
        # Temporary hack, we will sort this out later
        outfile = 'wf.h5'
        
        # No. determinants
        ndet    = 0
        
        # Extract the WFs for all states in this irrep
        args = (irr, nstates, state_indx, ci_confunits,
                ci_ciunits, outfile, ndet)
        ndet = libs.lib_func('wf_mrci', args)
        
    return
