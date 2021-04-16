"""
Module for generating the MRCI configurations from single
and double excitations out of an arbitrary reference space
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

def generate(scf, ci):
    """generate the MRCI configurations"""

    # Construct the molecule object
    timing.start('mrci_space')

    # number of irreps
    nirr = len(ci.nstates)

    # Print the section header
    output.print_mrcispace_header()
    
    # Bitci reference configuration scratch file numbers
    conf0scr = np.array(ci.ref_conf.confscr, dtype=int)

    # Bitci MRCI configuration scratch file numbers
    confMscr = np.zeros(nirr, dtype=int)

    # Number of MRCI configurations per irrep
    nconf = np.zeros(nirr, dtype=int)

    # Energy of the highest-lying reference space state
    emax = ci.ref_conf.ener[ci.ref_conf.ener != 0].max()

    # CVS core MO flags
    cvsflag = np.zeros(scf.nmo, dtype=int)
    for i in ci.icvs:
        cvsflag[i-1] = 1
    
    vec0scr = np.array(ci.ref_conf.vecscr, dtype=int)
    nroots  = ci.nstates

    # Loop over irreps
    for irrep in range(nirr):

        args = (irrep, nroots, conf0scr, confMscr, nconf, emax, cvsflag)
        (irrep, nroots, conf0scr, confMscr, nconf, emax, cvsflag) = \
                libs.lib_func('generate_mrci_confs',args)
        
        # Optional filtering based on the ASCI selection
        # criterion
        if ci.asci != 'off':
            thrsh = ci.asci_thresh[ci.asci]
            args = (thrsh, irrep, nroots, confMscr, vec0scr, nconf)
            (thrsh, irrep, nroots, confMscr, vec0scr, nconf) = \
                libs.lib_func('filter_asci', args)

    # Set the number of MRCI configurations
    ci.mrci_conf.set_nconf(nconf)
    
    # Set the MRCI configuration scratch file number
    ci.mrci_conf.set_confscr(confMscr)

    # Retrieve the MRCI configuration scratch file names
    confname = []
    name     = ' '*255
    for i in range(nirr):
        args = (ci.mrci_conf.confscr[irrep], name)
        (ci.mrci_conf.confscr[irrep], name) = \
                libs.lib_func('retrieve_filename', args)
        confname.append(name)
    ci.mrci_conf.set_confname(confname)
    
    # Stop timing
    timing.stop('mrci_space')
    
    return


