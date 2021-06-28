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
    nconf_new = np.zeros(nirr, dtype=int)

    # Number of roots per irrep
    nroots = ci_method.nstates

    # Bitci ref conf scratch file numbers
    confunit = ci_method.ref_wfn.conf_units
    
    # Bitci ref eigenvector scratch file numbers
    ciunit   = np.array(ci_method.ref_wfn.ci_units, dtype=int)
    
    # Number of ref confs for each irrep
    nconf    = ci_method.ref_wfn.nconf
    
    # Call to the bitci ref space pruning routine
    args = (nroots, confunit, nconf, ciunit)
    (nconf_new) = libs.lib_func('prune_ref_space',args)

    return nconf_new
