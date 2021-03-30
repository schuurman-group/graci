"""
Module for the refinement of the reference space in an MRCI calculation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.methods.molecule as molecule
import graci.io.convert as convert

def refine_ref_space(mol, ci, conf0, confsd, lib_bitci):
    """Refinement of the reference space"""

    # Start timing
    timing.start('mrci_refine')

    # Number of irreps
    if mol.comp_sym != 'c1':
        nirrep = molecule.nirrep[mol.sym_indx]
    else: 
        nirrep = 1

    # Bitci MRCI configuration scratch file numbers
    confscrM = convert.convert_ctypes(np.array(confsd.confscr, 
                                      dtype=int),dtype='int32')
        
    # Bitci eigenvector scratch file numbers
    vecscr = convert.convert_ctypes(np.array(confsd.vecscr, 
                                      dtype=int), dtype='int32')

    # No. roots per irrep
    nstates = convert.convert_ctypes(ci.nstates, dtype='int32')

    # Configuration selection threshold. We will just hardcode this
    # for now
    cthrsh = convert.convert_ctypes(0.1, dtype='double')

    # Minimum reference space norm
    min_norm = convert.convert_ctypes(0., dtype='double')
    
    # Scratch file numbers for the updated reference configurations
    confscrR = convert.convert_ctypes(np.array(conf0.confscr, 
                                      dtype=int), dtype='int32')

    # New number of reference space configurations per irrep
    nconf0 = np.zeros(nirrep, dtype=int)
    nconf0 = convert.convert_ctypes(nconf0, dtype='int32')
    
    # Refine the reference space
    lib_bitci.refine_ref_space(confscrM, confscrR, vecscr, nstates,
                               ctypes.byref(cthrsh),
                               ctypes.byref(min_norm),
                               nconf0)

    # Convert the nconf0 and confscr ctypes array to lists
    # (note that slicing a ctypes array will automatically
    # produce a list)
    nconf0=nconf0[:]
    confscrR=confscrR[:]

    # Set the number of reference space configurations
    conf0.set_nconf(nconf0)
    
    # Set the reference space configurations scratch file number
    conf0.set_confscr(confscrR)
        
    # Stop timing
    timing.stop('mrci_refine')
    
    return min_norm.value
