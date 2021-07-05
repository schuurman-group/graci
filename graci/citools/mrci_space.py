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
    nroots      = ci_method.n_state_sym()

    # Pruning: bitci Q-space energy correction scracth file numbers
    eq_units = np.zeros(nirr, dtype=int)
    
    # Loop over irreps
    for irrep in range(nirr):

        args = (irrep, nroots, ref_confunits, ci_confunits, \
                nconf, emax, cvsflag)
        (ci_confunits, nconf) = libs.lib_func('generate_mrci_confs',args)
        
        # Optional pruning of the configuration space
        if ci_method.pmrci:
            thrsh = ci_method.prune_thresh
            nextra = ci_method.nextra['prune'][irrep]            
            args = (thrsh, irrep, nroots, nextra, ci_confunits,
                    ref_ciunits, nconf, eq_units)
            (nconf, eq_units) = libs.lib_func('mrci_prune', args)
            
    # Retrieve the MRCI configuration scratch file names
    confname = []
    name     = ' '
    for i in range(nirr):
        args = (ci_confunits[i], name)
        name = libs.lib_func('retrieve_filename', args)
        confname.append(name)

    # Stop timing
    timing.stop('mrci_space')
    
    return nconf, ci_confunits, confname, eq_units

def set_prune_vars(ci_method):
    """Handles the setting of the MRCI pruning logical flag
    and threshold variables"""

    # Pruning not requested: return False and 1.0
    if ci_method.prune == 'off' and ci_method.prune_thresh == 1.:
        return False, 1.

    # Pruning reqested: return True and the requested threshold
    # Note that the 'prune_thresh' keyword takes precedent over
    # the 'prune' keyword
    if ci_method.prune_thresh != 1.:
        return True, ci_method.prune_thresh
    else:
        return True, ci_method.prune_dict[ci_method.prune]
