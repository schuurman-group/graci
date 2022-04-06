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

@timing.timed
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
        
        # States for which the WFs are required (note that bitCI
        # uses Fortran indexing here)
        states = np.array([n+1 for n in 
               range(ci_method.n_states_sym(irr))], dtype=int)

        # No. states
        nstates = states.size

        # state labeling of stored wfn (adiabatic label)
        state_lbls = np.array([ci_method.state_index(irr, states[n]-1) 
                                for n in range(nstates)], dtype=int)

        # Name of the HDF5 output file
        outfile = output.file_names['chkpt_file'] 
       
        # method label into which we store wfn datasets
        grp_name = type(ci_method).__name__ + '.' + str(ci_method.label)

        # No. determinants
        ndet    = 0
        
        # Extract the WFs for all states in this irrep
        args = (irr, nstates, states, ci_confunits,
                ci_ciunits, outfile, grp_name, state_lbls, ndet)
        
        ndet = libs.lib_func('wf_mrci', args)
        
    return
