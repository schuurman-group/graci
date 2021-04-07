"""
Module for performing the MRCI space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.io.convert as convert

def diag(ci, lib_bitci):
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
    for i in range(nirr):

        # Number of roots for the current irrep
        nroots = ctypes.c_int32(ci.nstates[i])

        # Block size for the current irrep
        nblock = blocksize[i]
        
        # Irrep number
        irrep = ctypes.c_int32(i)

        # Call to the bitci MRCI diagonalisation routine
        lib_bitci.diag_mrci(ctypes.byref(irrep),
                            ctypes.byref(nroots),
                            confscr,
                            ctypes.byref(vecscr1),
                            ctypes.byref(ialg),
                            ctypes.byref(tol),
                            ctypes.byref(niter),
                            ctypes.byref(blocksize),
                            ctypes.byref(deflate))

        # Bitci eigenvector scratch number
        vecscr.append(vecscr1.value)

    # Retrieve the MRCI energies
    maxroots = max(ci.nstates)
    ener = np.zeros((nirr, maxroots), dtype=float)
    for i in range(nirr):
        if ci.nstates[i] > 0:

            # Energies
            ener1 = np.zeros(ci.nstates[i], dtype=float)
            ener1 = convert.convert_ctypes(ener1, dtype='double')
    
            # Number of roots for the current irrep
            nroots = ctypes.c_int32(ci.nstates[i])
    
            # Eigenvector Scratch file number
            vecscr1 = ctypes.c_int32(vecscr[i])
    
            # Retrieve the energies
            lib_bitci.retrieve_energies(ctypes.byref(vecscr1),
                                        ctypes.byref(nroots),
                                        ener1)
            ener[i][:ci.nstates[i]] = ener1[:]

    # Save the list of bitci eigenvector scratch numbers
    ci.mrci_conf.set_vecscr(vecscr)

    # Retrieve the MRCI eigenvector scratch file names
    name    = convert.convert_ctypes(' '*255, dtype='string')
    vecname = []
    for i in range(nirr):
        scrnum = convert.convert_ctypes(ci.mrci_conf.vecscr[i],
                                        dtype='int32')
        lib_bitci.retrieve_filename(ctypes.byref(scrnum), name)
        vecname.append(bytes.decode(name.value))
    ci.mrci_conf.set_vecname(vecname)

    # Save the MRCI state energies
    ci.mrci_conf.set_ener(np.transpose(ener))

    # Print the report of the MRCI states
    vecscr = convert.convert_ctypes(np.array(ci.mrci_conf.vecscr, dtype=int),
                                 dtype='int32')
    nstates = convert.convert_ctypes(ci.nstates, 
                                 dtype='int32')
    lib_bitci.print_mrci_states(confscr, vecscr, nstates)
    
    # Stop timing
    timing.stop('mrci_diag')
    
    return

def diag_vars(ci):
    """ Returns the variables needed for the MRCI diagonalisation"""

    # nirr is given by the length of the nstates vector in ci obj
    nirr = len(ci.nstates)

    # Bitci MRCI configuration scratch file numbers
    confscr = np.array(ci.mrci_conf.confscr, dtype=int)
    confscr = convert.convert_ctypes(confscr, dtype='int32')

    # Bitci eigenvector scratch file numbers
    vecscr1 = 0
    vecscr1 = ctypes.c_int32(vecscr1)

    # Iterative diagonalisation algorithm
    if ci.diag_method == 'gendav':
        ialg = convert.convert_ctypes(1, dtype='int32')
    elif ci.diag_method == 'blockdav':
        ialg = convert.convert_ctypes(2, dtype='int32')
    else:
        sys.exit(
            '\n Unrecognised iterative diagonalisation algorithm:'
            +' '+ci.diag_method)

    # Maximum no. iterations
    niter = convert.convert_ctypes(ci.diag_iter,
                                dtype='int32')
    
    # Blocksizes
    if len(ci.diag_blocksize) == 0:
        # Blocksizes not given: asign sensible default values
        ci.diag_blocksize = [int(round(max(n+5, n*1.3))) 
                             for n in ci.nstates]
    
    blocksize = convert.convert_ctypes(np.array(ci.diag_blocksize),
                                       dtype='int32')
    # Deflation flag
    deflate = convert.convert_ctypes(ci.diag_deflate, dtype='logical')

    # Residual norm convergence threshold
    # If deflation is being used, then this needs to be tight
    if ci.diag_deflate and ci.diag_tol > 1e-6:
        print('\n Deflation is being used -> '+
                 'adjusting convergence threshold to 1e-6')
        ci.diag_tol = 1e-6

    # If this is the first MRCI iteration, then use a loose
    # convergence threshold
    if ci.niter == 1 and not ci.diag_deflate:
        tol = convert.convert_ctypes(1e-3, dtype='double')
    else:
        tol = convert.convert_ctypes(ci.diag_tol, dtype='double')

    return nirr, confscr, vecscr1, ialg, tol, niter, blocksize, deflate
