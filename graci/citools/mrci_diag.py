"""
Module for performing the MRCI space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert

def diag(ci):
    """Diagonalisation of the MRCI Hamiltonian"""

    # Start timing
    timing.start('mrci_diag')

    # Get the information needed to run the MRCI diagonalisation
    # calculation
    ialg, tol, niter, blocksize, deflate = diag_vars(ci)

    # nirr is given by the length of the nstates vector in ci obj
    nirr = len(ci.nstates)

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(ci.mrci_wfn.conf_units, dtype=int)

    # Initialise the list to hold the eigenvector scratch file
    # numbers
    ciunit  = 0
    ciunits = []
    
    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci.n_states(irrep)

        # Block size for the current irrep
        nblock = blocksize[irrep]
        
        args = (irrep, nroots, ci_confunits, ciunit, ialg, tol, niter, 
                blocksize, deflate)
        ciunit = libs.lib_func('diag_mrci', args)

        # Bitci eigenvector scratch number
        ciunits.append(ciunit)

    # Retrieve the MRCI energies
    maxroots = max(ci.nstates)
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci.n_states(irrep) > 0:

            # Number of roots for the current irrep
            nroots = ci.n_states(irrep)

            args = (ciunits[irrep], nroots, ener[irrep,:nroots])
            (ener[irrep,:nroots]) = \
                    libs.lib_func('retrieve_energies', args)

    # Save the list of bitci eigenvector scratch numbers
    ci.mrci_wfn.set_ciunits(ciunits)

    # Retrieve the MRCI eigenvector scratch file names
    ciname = []
    name    = ''
    for irrep in range(nirr):
        args = (ci.mrci_wfn.ci_units[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        ciname.append(name)

    ci.mrci_wfn.set_ciname(ciname)

    # Save the MRCI state energies
    ci.mrci_wfn.set_ener(np.transpose(ener))
    ciunits = np.array(ci.mrci_wfn.ci_units, dtype=int)
    nstates = ci.nstates
    
    # Print the report of the MRCI states
    args = (ci_confunits, ciunits, nstates)
    libs.lib_func('print_mrci_states', args)

    # Stop timing
    timing.stop('mrci_diag')
    
    return

def diag_vars(ci):
    """ Returns the variables needed for the MRCI diagonalisation"""

    # Iterative diagonalisation algorithm
    if ci.diag_method == 'gendav':
        ialg = 1
    elif ci.diag_method == 'blockdav':
        ialg = 2
    else:
        sys.exit(
            '\n Unrecognised iterative diagonalisation algorithm:'
            +' '+ci.diag_method)

    # Maximum no. iterations
    niter = ci.diag_iter

    # Blocksizes
    if len(ci.diag_blocksize) == 0:
        # Blocksizes not given: asign sensible default values
        ci.diag_blocksize = [int(round(max(n+5, n*1.3))) 
                             for n in ci.nstates]
    
    blocksize = np.array(ci.diag_blocksize, dtype=int)

    # Deflation flag
    deflate = ci.diag_deflate

    # Residual norm convergence threshold
    # If deflation is being used, then this needs to be tight
    if ci.diag_deflate and ci.diag_tol > 1.e-6:
        print('\n Deflation is being used -> '+
                 'adjusting convergence threshold to 1e-6')
        ci.diag_tol = 1.e-6

    # If this is the first MRCI iteration, then use a loose
    # convergence threshold
    if ci.niter == 1 and not ci.diag_deflate:
        tol = 1.e-3
    else:
        tol = ci.diag_tol

    return ialg, tol, niter, blocksize, deflate

