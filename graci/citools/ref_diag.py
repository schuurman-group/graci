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

def diag(mol, ci):
    """Diagonalisation of the reference space Hamiltonian"""

    # length of nstates vector is the number of irreps
    nirr = len(ci.nstates)

    # Construct the molecule object
    timing.start('ref_diag')

    # Print the section header
    output.print_refdiag_header()
    
    # Bitci reference configuration scratch file numbers
    confscr = np.array(ci.ref_conf.confscr, dtype=int)

    # Numbers of configurations
    nconf = np.array(ci.ref_conf.nconf, dtype=int)

    # Bitci eigenvector scratch file numbers
    vecscr1 = 0
    vecscr  = []
    
    # Loop over irreps
    for irrep in range(nirr):
    
        # Number of roots for the current irrep
        nroots = ci.n_states(irrep)

        # Call to the bitci reference space diagonalisation routine
        args = (irrep, nroots, confscr, nconf, vecscr1)
        (irrep, nroots, confscr, nconf, vecscr1) = \
                libs.lib_func('ref_diag_mrci',args)

        # Bitci eigenvector scratch number
        vecscr.append(vecscr1)
    
        # If the number of reference space configurations for the
        # current irrep is less than the requested number of roots
        # then reset nstates accordingly

        if nroots < ci.nstates[irrep]:
            ci.nstates[irrep] = nroots
    
    # Retrieve the reference space energies
    maxroots = max(ci.nstates)
    ener = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci.n_states(irrep) > 0:
    
            # Number of roots for the current irrep
            nroots = ci.n_states(irrep)

            # Retrieve the energies
            ener1 = np.zeros(nroots, dtype=float)
            args = (vecscr[irrep], nroots, ener1)
            (vecscr[irrep], nroots, ener1) = \
                    libs.lib_func('retrieve_energies', args)
            ener[irrep, :nroots] = ener1
        
    # Save the list of bitci eigenvector scratch numbers
    ci.ref_conf.set_vecscr(vecscr)
    
    # Save the reference space state energies
    ci.ref_conf.set_ener(np.transpose(ener))
    
    # Print the section header
    output.print_refdiag_summary(mol, ci)
    
    # Stop timing
    timing.stop('ref_diag')
    
    return
