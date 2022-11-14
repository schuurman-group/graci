"""
Module for the calculation of QDPT Type-I and Type-II diabatic states
within the DFT/MRCI(2) framework
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.output as output
import graci.io.convert as convert
import graci.core.molecule as molecule

@timing.timed
def type_ii(ci_method0, ci_method1):
    """
    Type-II QDPT diabatisation within the DFT/MRCI(2) framework
    Here, ci_method0 corresponds to the previous geometry, R_n-1,
    and ci_method1 to the current geometry R_n
    """

    # initialise the list of ADT matrices, one per irrep
    adt_matrices = []

    # GVVPT2 regularizer index
    ireg = ci_method1.allowed_regularizer.index(ci_method1.regularizer)+1

    # nirr is given by the length of the nstates vector in ci obj
    nirr = ci_method1.n_irrep()
    
    # the bitci mrci wfn object
    mrci_wfn = ci_method1.bitci_mrci()

    # bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)

    # bitci ref space eigenvector scratch file numbers
    ref_ciunits = np.array(ci_method1.ref_wfn.ci_units, dtype=int)    

    # A-vector scratch file number
    Aunits = ci_method1.Aunits

    print('\n\n Aunits:', Aunits)
    sys.exit()
    
    # loop over irreps
    #for irr in range(nirr):
    
        
    
    sys.exit()
    
    return
