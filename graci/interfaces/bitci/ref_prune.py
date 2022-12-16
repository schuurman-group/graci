"""
Module for the pruning of the reference space
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.core.libs as libs
import graci.io.convert as convert

@timing.timed
def prune(ci_method):
    """removes deadwood from the reference space"""

    # number of irreps
    nirr = ci_method.n_irrep()

    # number of ref space configurations after pruning
    nconf_new = np.zeros(nirr, dtype=int)

    # Number of roots per irrep
    nroots = ci_method.nstates

    # Number of extra roots per irrep
    nextra = np.array(ci_method.nextra['max'])

    # Bitci ref conf scratch file numbers
    confunit = ci_method.ref_wfn.conf_units['adiabatic']
    
    # Bitci ref eigenvector scratch file numbers
    ciunit   = np.array(ci_method.ref_wfn.ci_units['adiabatic'], dtype=int)
    
    # Number of ref confs for each irrep
    nconf    = ci_method.ref_wfn.nconf['adiabatic']

    # Call to the bitci ref space pruning routine
    args = (nroots, nextra, confunit, nconf, ciunit)
    (nconf_new) = libs.lib_func('prune_ref_space',args)

    return nconf_new
