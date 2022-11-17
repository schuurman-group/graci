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

    # Initialise the determinant bit string and eigenvector arrays
    ci_method.det_strings['adiabatic'] = []
    ci_method.vec_det['adiabatic']     = []
    
    # Loop over irreps
    for irr in range(ci_method.n_irrep()):
        
        # States for which the WFs are required (note that bitCI
        # uses Fortran indexing here)
        states = np.array([n+1 for n in 
               range(ci_method.n_states_sym(irr))], dtype=int)

        # No. states
        nstates = states.size

        # No. determinants
        ndet = 0

        # Determinant wave function scratch file numbers
        ci_wfunit = 0

        # Norm-based truncation threshold
        norm_thresh = ci_method.wf_thresh

        # Truncated no. determinants
        ndet_trunc = 0
        
        # Convert the WFs from the CSF to determinant basis for
        # all states in this irrep
        args              = (irr, nstates, states, ci_confunits,
                             ci_ciunits, ndet, ci_wfunit)
        (ndet, ci_wfunit) = libs.lib_func('wf_mrci', args)

        # Get the truncated no. determinants
        args       = (ci_wfunit, nstates, norm_thresh, ndet_trunc)
        ndet_trunc = libs.lib_func('ndet_truncated', args)

        # Retrieve the truncated list of determinants and
        # eigenvectors
        n_int      = 0
        args       = [n_int]
        n_int      = libs.lib_func('get_n_int', args)
        dets       = np.zeros((n_int*2*ndet_trunc), dtype=np.int64)
        vec        = np.zeros((ndet_trunc*nstates), dtype=float)
        args       = (ci_wfunit, ndet, ndet_trunc, nstates, dets,
                      vec, norm_thresh)
        (det, vec) = libs.lib_func('retrieve_det_truncated', args)

        # Reshaping
        det = np.reshape(det, (n_int,2,ndet_trunc), order='F')
        vec = np.reshape(vec, (ndet_trunc,nstates), order='F')
       
        # Save the determinant bit strings and eigenvectors
        # for this irrep
        ci_method.det_strings['adiabatic'].append(det)
        ci_method.vec_det['adiabatic'].append(vec)
        
    return
