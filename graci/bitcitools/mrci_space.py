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

@timing.timed
def generate(ci_method):
    """generate the MRCI configurations"""

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
    cvs_flag = np.zeros(nmo, dtype=int)
    for i in ci_method.icvs:
        cvs_flag[i-1] = 1

    # DDCI flag
    ddci_flag = ci_method.ddci
    
    ref_ciunits = np.array(ref_wfn.ci_units, dtype=int)
    nroots      = ci_method.n_states_sym()

    # Pruning: bitci Q-space energy correction scracth file numbers
    eq_units = np.zeros(nirr, dtype=int)

    # Generate the MRCI configurations for all irreps
    args = (nroots, ref_confunits, ci_confunits, nconf, emax, cvs_flag,
            ddci_flag)
    (ci_confunits, nconf) = libs.lib_func('generate_mrci_confs',args)
    
    # Loop over irreps
    for irrep in range(nirr):
        
        # Optional pruning of the configuration space
        # (note that not all mrci-type methods will have pruning
        # attributes)
        try:
            if ci_method.prune:
                thrsh = ci_method.prune_thresh
                nextra = ci_method.nextra['prune'][irrep]            
                args = (thrsh, irrep, nroots, nextra, ci_confunits,
                        ref_ciunits, nconf, eq_units)
                (nconf, eq_units) = libs.lib_func('mrci_prune', args)
        except AttributeError:
            pass
    
    # Retrieve the MRCI configuration scratch file names
    confname = []
    name     = ' '
    for i in range(nirr):
        args = (ci_confunits[i], name)
        name = libs.lib_func('retrieve_filename', args)
        confname.append(name)

    return nconf, ci_confunits, confname, eq_units

