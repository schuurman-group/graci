"""
Module for the truncation of DFT/MRCI(2) 1st-order corrected
wave functions
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing

@timing.timed
def truncate_wf(ci_method):
    """
    Removal of deadwood configurations from the DFT/MRCI(2)
    1st-order corrected wave functions
    """
    # the bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()
    
    # number of configurations per irrep
    nconf = np.array(mrci_wfn.nconf, dtype=int)

    # bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)

    # bitci MRCI eigenvector scratch file numbers
    ciunits = mrci_wfn.ci_units['adiabatic']
    
    for irrep in range(ci_method.n_irrep()):
        thresh       = ci_method.truncate_thresh
        nroots       = ci_method.n_states_sym(irrep)
        nconf_new    = 0
        args         = (irrep, nroots, ci_confunits[irrep],
                        ciunits[irrep], thresh, nconf_new)
        (nconf_new)  = libs.lib_func('truncate_mrci_wf', args)
        nconf[irrep] = nconf_new

    return nconf
