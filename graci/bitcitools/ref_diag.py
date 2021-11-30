"""
Module for performing the reference space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

@timing.timed
def diag(ci_method):
    """Diagonalisation of the reference space Hamiltonian"""

    # length of nstates vector is the number of irreps
    nirr    = ci_method.n_irrep()

    # bitci reference space wfn
    ref_wfn = ci_method.bitci_ref()

    # Print the section header
    output.print_refdiag_header()
    
    # Bitci reference configuration scratch file numbers
    confunits = np.array(ref_wfn.conf_units, dtype=int)

    # Numbers of configurations
    nconf = np.array(ref_wfn.nconf, dtype=int)

    # Bitci eigenvector scratch file numbers
    ciunit   = 0
    ciunits = []
    
    # Loop over irreps
    for irrep in range(nirr):
    
        # Number of roots for the current irrep
        nroots = ci_method.n_state_sym(irrep)

        # Number of extra roots
        nextra = ci_method.nextra['max'][irrep]

        # Call to the bitci reference space diagonalisation routine
        args = (irrep, nroots+nextra, confunits, nconf, ciunit)
        (nroots, ciunit) = libs.lib_func('ref_diag_mrci',args)

        # Bitci eigenvector scratch number
        ciunits.append(ciunit)
    
        # If the number of reference space configurations for the
        # current irrep is less than the requested number of roots
        # then reset nstates accordingly
        if nroots < ci_method.n_state_sym(irrep):
            ci_method.nstates[irrep] = nroots
    
    # Retrieve the reference space energies
    maxroots = max(ci_method.n_state_sym())
    ener = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci_method.n_state_sym(irrep) > 0:
    
            # Number of roots for the current irrep
            nroots = ci_method.n_state_sym(irrep)

            # Retrieve the energies
            args = (ciunits[irrep], nroots, ener[irrep, :nroots],
                    np.array([n+1 for n in range(nroots)], dtype=int))
            
            (ener[irrep, :nroots]) = \
                    libs.lib_func('retrieve_some_energies', args)
            
    return ciunits, ener 

@timing.timed
def n_extra(ci_method):
    """Determination of the number of extra reference space
    eigenvectors needed"""

    # Initialisation
    nextra = {'prune' : None,
              'guess' : None,
              'max'   : None}

    # length of nstates vector is the number of irreps
    nirr    = ci_method.n_irrep()
    
    # If pruning is off and ENPT2 guess vectors are not being used,
    # then nextra = 0
    if ci_method.prune == 'off' and ci_method.diag_guess != 'enpt2':
        for key in nextra:
            nextra[key] = [0 for n in range(nirr)]
        return nextra
            
    # Number of extra ref space roots needed for pruning
    if ci_method.prune != 'off':
        nextra['prune'] = [ci_method.prune_extra for n in range(nirr)]
    else:
        nextra['prune'] = [0 for n in range(nirr)]
    
    # Number of extra ref space roots needed for guess vector generation
    if ci_method.diag_guess == 'enpt2':
        if len(ci_method.diag_blocksize) == 0:
            nextra['guess'] = [int(round(max(n+5, n*1.3)) - n) 
                                        for n in ci_method.nstates]
        else:
            nextra['guess'] = ci_method.diag_blocksize - ci_method.nstates
    else:
        nextra['guess'] = [0 for n in range(nirr)]

    # Number of extra ref space roots to be calculated
    nextra['max'] = []
    for n in range(nirr):
        nextra['max'].append(max([nextra[key][n] for key in nextra
                                  if key != 'max']))
            
    return nextra
