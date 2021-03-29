"""
Module for generating the MRCI configurations from single
and double excitations out of an arbitrary reference space
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.methods.molecule as molecule
import graci.methods.wavefunction as wavefunction
import graci.io.convert as convert
import graci.io.output as output

def generate(mol, ci, conf0, lib_bitci):
    """generate the MRCI configurations"""

    # Construct the molecule object
    timing.start('mrci_space')

    # Print the section header
    output.print_mrcispace_header()
    
    # Initialise the reference space wavefunction object
    confsd = wavefunction.Wavefunction()

    # Number of irreps
    if mol.comp_sym != 'c1':
        nirrep = molecule.nirrep[mol.sym_indx]
    else:
        nirrep = 1
    
    # Bitci reference configuration scratch file numbers
    conf0scr = convert.convert_ctypes(np.array(conf0.confscr, dtype=int),
                               dtype='int32')

    # Bitci MRCI configuration scratch file numbers
    confMscr = convert.convert_ctypes(np.zeros(nirrep, dtype=int),
                               dtype='int32')

    # Number of MRCI configurations per irrep
    nconf = convert.convert_ctypes(np.zeros(nirrep, dtype=int),
                               dtype='int32')
    
    # Energy of the highest-lying reference space
    # state of interest relative to the ground state
    emax = convert.convert_ctypes(conf0.ener.max()-conf0.ener.min(),
                               dtype='double')
        
    # Number of roots for each irrep
    nroots = convert.convert_ctypes(np.array(ci.nstates, 
                               dtype=int),dtype='int32')

    # Reference space eigenvector scratch file numbers
    vec0scr = convert.convert_ctypes(np.array(conf0.vecscr, dtype=int),
                               dtype='int32')
    
    # Loop over irreps
    for i in range(nirrep):

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
    confsd.set_nconf(nconf)
    
    # Set the reference space determinants scratch file number
    confsd.set_confscr(confMscr)

    # Stop timing
    timing.stop('mrci_space')
    
    return confsd


