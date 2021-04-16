"""
Module for the refinement of the reference space in an MRCI calculation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert

def refine_ref_space(ci):
    """Refinement of the reference space"""

    # Start timing
    timing.start('mrci_refine')

    # number of irreps given by length of ci.nstates vector
    nirr = len(ci.nstates)

    # Bitci MRCI configuration scratch file numbers
    confscrM = np.array(ci.mrci_conf.confscr, dtype=int)

    # Bitci eigenvector scratch file numbers
    vecscr = np.array(ci.mrci_conf.vecscr, dtype=int)

    # No. roots per irrep
    nstates = ci.nstates

    # Configuration selection threshold. We will just hardcode this
    # for now
    cthrsh = 0.1

    # Minimum reference space norm
    min_norm = 0.

    # Scratch file numbers for the updated reference configurations
    confscrR = np.array(ci.ref_conf.confscr, dtype=int)

    # New number of reference space configurations per irrep
    nconf0 = np.zeros(nirr, dtype=int)
    
    # Refine the reference space
    args = (confscrM, confscrR, vecscr, nstates, cthrsh, 
            min_norm, nconf0)
    (confscrM, confscrR, vecscr, nstates, cthrsh, min_norm, nconf0) = \
            libs.lib_func('refine_ref_space', args)

    # Set the number of reference space configurations
    ci.ref_conf.set_nconf(nconf0)
    
    # Set the reference space configurations scratch file number
    ci.ref_conf.set_confscr(confscrR)
        
    # Stop timing
    timing.stop('mrci_refine')
    
    #return min_norm.value
    return min_norm
