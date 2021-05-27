"""
Module for the pruning of the reference space
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.io.output as output
import graci.io.convert as convert

def prune(ci_method):
    """removes deadwood from the reference space"""

    # number of irreps
    nirr = ci_method.n_irrep()

    # number of ref space configurations after pruning
    nconf_new = []
    
    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci_method.n_states(irrep)

        # Bitci ref conf scratch file number
        confunit = ci_method.ref_wfn.conf_units[irrep]

        # Bitci ref eigenvector scratch file number
        ciunit   = ci_method.ref_wfn.ci_units[irrep]

        # Number of ref confs for the current irrep
        nconf    = ci_method.ref_wfn.nconf[irrep]

        # Call to the bitci ref space pruning routine
        args = (irrep, nroots, confunit, nconf, ciunit)
        (nconf) = libs.lib_func('prune_ref_space',args)
        
    
    sys.exit()

    # We need to return the new numbers of ref confs!
    # (to be saved to the bitciwfn class in the dftmrci module)

    # Should we be using nroots+nextra for *all* calculations if we
    # are going to be removing ref space confs?
    
    return
