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
    ref_confunits = np.array(ci.ref_wfn.conf_units, dtype=int)

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.zeros(nirr, dtype=int)

    # Number of MRCI configurations per irrep
    nconf = np.zeros(nirr, dtype=int)

    # Energy of the highest-lying reference space state
    emax = ci.ref_wfn.ener[ci.ref_wfn.ener != 0].max()

    # CVS core MO flags
    cvsflag = np.zeros(scf.nmo, dtype=int)
    for i in ci.icvs:
        cvsflag[i-1] = 1
    
    ref_ciunits = np.array(ci.ref_wfn.ci_units, dtype=int)
    nroots      = ci.nstates

    # Loop over irreps
    for irrep in range(nirr):

        args = (irrep, nroots, ref_confunits, ci_confunits, \
                nconf, emax, cvsflag)
        (ci_confunits, nconf) = libs.lib_func('generate_mrci_confs',args)
        
        # Optional pruning of the configuration space
        if ci.prune != 'off':
            thrsh = ci.prune_thresh[ci.prune]
            args = (thrsh, irrep, nroots, ci_confunits, ref_ciunits, nconf)
            nconf = libs.lib_func('mrci_prune', args)

    # Set the number of MRCI configurations
    ci.mrci_wfn.set_nconf(nconf)
    
    # Set the MRCI configuration scratch file number
    ci.mrci_wfn.set_confunits(ci_confunits)

    # Retrieve the MRCI configuration scratch file names
    confname = []
    name     = ' '
    for i in range(nirr):
        args = (ci.mrci_wfn.conf_units[i], name)
        name = libs.lib_func('retrieve_filename', args)
        confname.append(name)
    ci.mrci_wfn.set_confname(confname)

    # Stop timing
    timing.stop('mrci_space')
    
    return

 
