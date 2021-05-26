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

def generate(ci_method):
    """generate the MRCI configurations"""

    # Construct the molecule object
    timing.start('mrci_space')

    # number of irreps
    nirr    = ci_method.n_irrep()

    # number of mos
    nmo     = ci_method.scf.nmo

    # bitci reference space wfn
    ref_wfn = ci_method.bitci_ref()

    # Print the section header
    output.print_mrcispace_header()
    
    # Bitci reference configuration scratch file numbers
    ref_confunits = np.array(ref_wfn.conf_units, dtype=int)

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.zeros(nirr, dtype=int)

    # Number of MRCI configurations per irrep
    nconf = np.zeros(nirr, dtype=int)

    # Energy of the highest-lying reference space state
    emax = ci_method.ref_ener[ci_method.ref_ener != 0].max()

    # CVS core MO flags
    cvsflag = np.zeros(nmo, dtype=int)
    for i in ci_method.icvs:
        cvsflag[i-1] = 1
    
    ref_ciunits = np.array(ref_wfn.ci_units, dtype=int)
    nroots      = ci_method.n_states()

    # Pruning: no. extra roots to include in the ENPT2 calculation
    nextra = ci_method.prune_extra
    
    # Loop over irreps
    for irrep in range(nirr):

        args = (irrep, nroots, ref_confunits, ci_confunits, \
                nconf, emax, cvsflag)
        (ci_confunits, nconf) = libs.lib_func('generate_mrci_confs',args)
        
        # Optional pruning of the configuration space
        if ci_method.prune != 'off':
            thrsh = ci_method.prune_thresh[ci_method.prune]
            args = (thrsh, irrep, nroots, nextra, ci_confunits,
                    ref_ciunits, nconf)
            nconf = libs.lib_func('mrci_prune', args)

    # Retrieve the MRCI configuration scratch file names
    confname = []
    name     = ' '
    for i in range(nirr):
        args = (ci_confunits[i], name)
        name = libs.lib_func('retrieve_filename', args)
        confname.append(name)

    # Stop timing
    timing.stop('mrci_space')
    
    return nconf, ci_confunits, confname

 
