"""
Module for the refinement of the reference space in an MRENPT2 calculation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert

@timing.timed
def refine_ref_space(ci_method):
    """Refinement of the reference space"""

    # number of irreps given by length of ci.nstates vector
    nirr = ci_method.n_irrep()

    # bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # bitci refernce space wfn object
    ref_wfn  = ci_method.bitci_ref()

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)

    # Bitci eigenvector scratch file numbers
    ci_ciunits = np.array(mrci_wfn.ci_units, dtype=int)

    # Bitci Q-space info scratch file numbers
    ci_qunits = np.array(ci_method.qunits, dtype=int)

    # No. roots per irrep
    nstates = ci_method.n_states_sym()

    # Configuration selection parameters
    alpha = 0.061
    beta  = 3.3
    
    # Minimum reference space norm
    min_norm = 0.

    # Scratch file numbers for the updated reference configurations
    ref_confunits = np.array(ref_wfn.conf_units, dtype=int)

    # New number of reference space configurations per irrep
    ref_nconf = np.zeros(nirr, dtype=int)
    
    # Refine the reference space
    args = (ci_confunits, ref_confunits, ci_ciunits, ci_qunits,
            nstates, alpha, beta, min_norm, ref_nconf)
    (confunits_ref, min_norm, ref_nconf) = \
            libs.lib_func('refine_ref_space_pt2', args)

    return min_norm, ref_nconf, ref_confunits
