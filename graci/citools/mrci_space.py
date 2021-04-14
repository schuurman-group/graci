"""
Module for generating the MRCI configurations from single
and double excitations out of an arbitrary reference space
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

def generate(scf, ci, lib_bitci):
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
    
    # Energy of the highest-lying reference space state
    nonzero = ci.ref_conf.ener[ci.ref_conf.ener != 0]
    emax = convert.convert_ctypes(nonzero.max(), dtype='double')

    # CVS core MO flags
    cvsflag = np.zeros(scf.nmo, dtype=int)
    for i in ci.icvs:
        cvsflag[i-1] = 1
    cvsflag = convert.convert_ctypes(cvsflag, dtype='int32')
    
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
                                      ctypes.byref(emax),
                                      cvsflag)

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
    
    # Set the MRCI configuration scratch file number
    ci.mrci_conf.set_confscr(confMscr)

    # Retrieve the MRCI configuration scratch file names
    name     = convert.convert_ctypes(' '*255, dtype='string')
    confname = []
    for i in range(nirr):
        scrnum = convert.convert_ctypes(ci.mrci_conf.confscr[i],
                                        dtype='int32')
        lib_bitci.retrieve_filename(ctypes.byref(scrnum), name)
        confname.append(bytes.decode(name.value))
    ci.mrci_conf.set_confname(confname)
    
    # Stop timing
    timing.stop('mrci_space')
    
    return


