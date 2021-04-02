"""
Module for generating the MRCI configurations from single
and double excitations out of an arbitrary reference space
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

def generate(ci, lib_bitci):
    """generate the MRCI configurations"""

    # Construct the molecule object
    timing.start('mrci_space')

    # number of irreps
    nirr = len(ci.nstates)

    # Print the section header
    output.print_mrcispace_header()
    
    # Bitci reference configuration scratch file numbers
    conf0scr = convert.convert_ctypes(np.array(ci.ref_conf.confscr, 
                               dtype=int), dtype='int32')

    # Bitci MRCI configuration scratch file numbers
    confMscr = convert.convert_ctypes(np.zeros(nirr, dtype=int),
                               dtype='int32')

    # Number of MRCI configurations per irrep
    nconf = convert.convert_ctypes(np.zeros(nirr, dtype=int),
                               dtype='int32')
    
    # Energy of the highest-lying reference space
    # state of interest relative to the ground state
    emax = convert.convert_ctypes(ci.ref_conf.ener.max() - 
                                  ci.ref_conf.ener.min(), 
                                  dtype='double')
        
    # Number of roots for each irrep
    nroots = convert.convert_ctypes(np.array(ci.nstates, 
                               dtype=int),dtype='int32')

    # Reference space eigenvector scratch file numbers
    vec0scr = convert.convert_ctypes(np.array(ci.ref_conf.vecscr, 
                               dtype=int), dtype='int32')
    
    # Loop over irreps
    for i in range(nirr):

        # Irrep number
        irrep = ctypes.c_int32(i)
        
        # Construct the MRCI configurations for the i'th irrep
        lib_bitci.generate_mrci_confs(ctypes.byref(irrep),
                                      nroots,
                                      conf0scr,
                                      confMscr,
                                      nconf,
                                      ctypes.byref(emax))

        # Optional filtering based on the ASCI selection
        # criterion
        if ci.asci != 'off':
            thrsh = ci.asci_thresh[ci.asci]
            thrsh = convert.convert_ctypes(thrsh, dtype='double')
            lib_bitci.filter_asci(ctypes.byref(thrsh),
                                  ctypes.byref(irrep),
                                  nroots,
                                  confMscr,
                                  vec0scr,
                                  nconf)
            
    # Convert the nconf and confMscr ctypes array to lists
    # (note that slicing a ctypes array will automatically
    # produce a list)
    nconf=nconf[:]
    confMscr=confMscr[:]

    # Set the number of MRCI configurations
    ci.mrci_conf.set_nconf(nconf)
    
    # Set the reference space determinants scratch file number
    ci.mrci_conf.set_confscr(confMscr)

    # Stop timing
    timing.stop('mrci_space')
    
    return


