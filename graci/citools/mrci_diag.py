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
    nirr, confscr, vecscr1, ialg, tol, niter, blocksize, deflate = \
            diag_vars(ci)

    # Initialise the list to hold the eigenvector scratch file
    # numbers
    vecscr  = []
    
    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci.n_states(irrep)

        # Block size for the current irrep
        nblock = blocksize[irrep]
        
        args = (irrep, nroots, confscr, vecscr1, ialg, tol, niter, 
                blocksize, deflate)
        (irrep, nroots, confscr, vecscr1, ialg, tol, niter,
                blocksize, deflate) = libs.lib_func('diag_mrci', args)

        # Bitci eigenvector scratch number
        vecscr.append(vecscr1)

    # Retrieve the MRCI energies
    maxroots = max(ci.nstates)
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci.n_states(irrep) > 0:

            # Number of roots for the current irrep
            nroots = ci.n_states(irrep)

            # Retrieve the energies
            ener1   = np.zeros(nroots, dtype=float)
            vecscr1 = vecscr[irrep]

            args = (vecscr1, nroots, ener1)
            (vecscr1, nroots, ener1) = \
                    libs.lib_func('retrieve_energies', args)
            ener[irrep,:nroots] = ener1

    # Save the list of bitci eigenvector scratch numbers
    ci.mrci_conf.set_vecscr(vecscr)

    # Retrieve the MRCI eigenvector scratch file names
    vecname = []
    name    = ''
    for irrep in range(nirr):
        #scrnum = convert.convert_ctypes(ci.mrci_conf.vecscr[i],
        #                                dtype='int32')
        vecscr1 = ci.mrci_conf.vecscr[irrep]
        args = (vecscr1, name)
        (vecscr1, name) = libs.lib_func('retrieve_filename', args)
        vecname.append(name)

    ci.mrci_conf.set_vecname(vecname)

    # Save the MRCI state energies
    ci.mrci_conf.set_ener(np.transpose(ener))
    vscr    = np.array(ci.mrci_conf.vecscr, dtype=int)
    nstates = ci.nstates
    
    # Print the report of the MRCI states
    args = (confscr, vscr, nstates)
    (confscr, vscr, nstates) = \
            libs.lib_func('print_mrci_states', args)

    # Stop timing
    timing.stop('mrci_diag')
    
    return

def diag_vars(ci):
    """ Returns the variables needed for the MRCI diagonalisation"""

    # nirr is given by the length of the nstates vector in ci obj
    nirr = len(ci.nstates)

    # Bitci MRCI configuration scratch file numbers
    confscr = np.array(ci.mrci_conf.confscr, dtype=int)

    # Bitci eigenvector scratch file numbers
    vecscr1 = 0

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

    return nirr, confscr, vecscr1, ialg, tol, niter, blocksize, deflate

