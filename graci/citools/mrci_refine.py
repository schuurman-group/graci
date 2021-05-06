"""
Module for the refinement of the reference space in an MRCI calculation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert

def refine_ref_space(ci_method):
    """Refinement of the reference space"""

    # Start timing
    timing.start('mrci_refine')

    # number of irreps given by length of ci.nstates vector
    nirr = len(ci_method.nstates)

    # bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # bitci refernce space wfn object
    ref_wfn  = ci_method.bitci_ref()

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)

    # Bitci eigenvector scratch file numbers
    ci_ciunits = np.array(mrci_wfn.ci_units, dtype=int)

    # No. roots per irrep
    nstates = ci_method.nstates

    # Configuration selection threshold. We will just hardcode this
    # for now
    cthrsh = 0.055
    
    # Minimum reference space norm
    min_norm = 0.

    # Scratch file numbers for the updated reference configurations
    ref_confunits = np.array(ref_wfn.conf_units, dtype=int)

    # New number of reference space configurations per irrep
    ref_nconf = np.zeros(nirr, dtype=int)
    
    # Refine the reference space
    args = (ci_confunits, ref_confunits, ci_ciunits, nstates, cthrsh, 
            min_norm, ref_nconf)
    (confunits_ref, min_norm, ref_nconf) = \
            libs.lib_func('refine_ref_space', args)

    # Stop timing
    timing.stop('mrci_refine')
    
    #return min_norm.value
    return min_norm, ref_nconf, ref_confunits
