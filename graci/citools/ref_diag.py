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
    confunits = np.array(ci.ref_wfn.conf_units, dtype=int)

    # Numbers of configurations
    nconf = np.array(ci.ref_wfn.nconf, dtype=int)

    # Bitci eigenvector scratch file numbers
    ciunit   = 0
    ciunits = []
    
    # Loop over irreps
    for irrep in range(nirr):
    
        # Number of roots for the current irrep
        nroots = ci.n_states(irrep)

        # Call to the bitci reference space diagonalisation routine
        args = (irrep, nroots, confunits, nconf, ciunit)
        (nroots, ciunit) = libs.lib_func('ref_diag_mrci',args)

        # Bitci eigenvector scratch number
        ciunits.append(ciunit)
    
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
            args = (ciunits[irrep], nroots, ener[irrep, :nroots])
            (ener[irrep, :nroots]) = \
                    libs.lib_func('retrieve_energies', args)
        
    # Save the list of bitci eigenvector scratch numbers
    ci.ref_wfn.set_ciunits(ciunits)
    
    # Save the reference space state energies
    ci.ref_wfn.set_ener(np.transpose(ener))
    
    # Print the section header
    output.print_refdiag_summary(mol, ci)
    
    # Stop timing
    timing.stop('ref_diag')
    
    return
