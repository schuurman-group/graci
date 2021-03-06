"""
Module for calculating the ENPT2 energy and wave function corrections
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.output as output
import graci.io.convert as convert

@timing.timed
def corrections(ci_method):
    """Calculation of the ENPT2 corrections"""

    # is this a multistate MR-ENPT2 calculation?
    multistate = ci_method.multistate

    # ISA shift
    shift = ci_method.shift
    
    # nirr is given by the length of the nstates vector in ci obj
    nirr = ci_method.n_irrep()

    # the bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)

    # Bitci ref space eigenvector scratch file numbers
    ref_ciunits = np.array(ci_method.ref_wfn.ci_units, dtype=int)

    # Initialise the list to hold the eigenvector and Q-space
    # scratch file numbers
    ciunit  = 0
    ciunits = []
    qunit   = 0
    qunits  = []
    dspunit  = 0
    dspunits = []

    # Number of configurations per irrep (this can change if wave
    # function truncation is being used)
    nconf = np.array(ci_method.mrci_wfn.nconf, dtype=int)
    
    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep)

        # Number of extra roots
        nextra = ci_method.nextra['enpt2'][irrep]
        
        args = (irrep, nroots, nextra, shift, multistate, ci_confunits,
                ciunit, ref_ciunits, qunit, dspunit)

        ciunit, qunit, dspunit = libs.lib_func('mrenpt2', args)

        # Bitci eigenvector, Q-space, and intruder state
        # scratch numbers
        ciunits.append(ciunit)
        qunits.append(qunit)
        dspunits.append(dspunit)

    # Retrieve the MR-ENPT2 energies
    maxroots = max(ci_method.n_states_sym())
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci_method.n_states_sym(irrep) > 0:

            # Number of roots for the current irrep
            nroots = ci_method.n_states_sym(irrep)

            args = (ciunits[irrep], nroots, ener[irrep,:nroots])
            (ener[irrep,:nroots]) = \
                    libs.lib_func('retrieve_energies', args)

    # Retrieve the MR-ENPT2 eigenvector scratch file names
    ciname = []
    name   = ''
    for irrep in range(nirr):
        args = (ciunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        ciname.append(name)

    # Optional truncation of the MR-ENPT2 1st-order corrected
    # wave functions
    if ci_method.truncate:
        for irrep in range(nirr):
            thresh       = ci_method.truncate_thresh
            nroots       = ci_method.n_states_sym(irrep)
            nconf_new    = 0
            args         = (irrep, nroots, ci_confunits[irrep],
                            ciunits[irrep], thresh, nconf_new)
            (nconf_new)  = libs.lib_func('truncate_mrci_wf', args)
            nconf[irrep] = nconf_new
    
    # Print the report of the MR-ENPT2 states 
    output.print_dftmrenpt2_states_header()
    ciunits = np.array(ciunits, dtype=int)
    nstates = ci_method.n_states_sym()
    args = (ci_confunits, ciunits, nstates)
    libs.lib_func('print_mrci_states', args)  
    
    return ciunits, ciname, ener, qunits, dspunits, nconf
