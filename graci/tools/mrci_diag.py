"""
Module for performing the MRCI space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.methods.params as params
import graci.molecule.molecule as molecule
import graci.io.convert

def diag(mol, confsd, lib_bitci):
    """Diagonalisation of the MRCI Hamiltonian"""

    # Start timing
    timing.start('mrci_diag')

    # Get the information needed to run the MRCI diagonalisation
    # calculation
    nirrep, confscr, vecscr1, ialg, tol, niter, \
        blocksize, deflate = diag_vars(mol, confsd)

    # Initialise the list to hold the eigenvector scratch file
    # numbers
    vecscr  = []
    
    # Loop over irreps
    for i in range(nirrep):

        # Number of roots for the current irrep
        nroots = ctypes.c_int32(params.mol_param['nstates'][i])

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
    maxroots = max(params.mol_param['nstates'])
    ener = np.zeros((nirrep, maxroots), dtype=float)
    for i in range(nirrep):
        if params.mol_param['nstates'][i] > 0:

            # Energies
            ener1 = np.zeros(params.mol_param['nstates'][i], dtype=float)
            ener1 = convert.convert_ctypes(ener1, dtype='double')
    
            # Number of roots for the current irrep
            nroots = ctypes.c_int32(params.mol_param['nstates'][i])
    
            # Eigenvector Scratch file number
            vecscr1 = ctypes.c_int32(vecscr[i])
    
            # Retrieve the energies
            lib_bitci.retrieve_energies(ctypes.byref(vecscr1),
                                        ctypes.byref(nroots),
                                        ener1)
            ener[i][:params.mol_param['nstates'][i]] = ener1[:]

    # Save the list of bitci eigenvector scratch numbers
    confsd.set_vecscr(vecscr)
    
    # Save the MRCI state energies
    confsd.set_ener(np.transpose(ener))

    # Print the report of the MRCI states
    vecscr = convert.convert_ctypes(np.array(confsd.vecscr, dtype=int),
                                 dtype='int32')
    nstates = convert.convert_ctypes(params.mol_param['nstates'], 
                                 dtype='int32')
    lib_bitci.print_mrci_states(confscr, vecscr, nstates)
    
    # Stop timing
    timing.stop('mrci_diag')
    
    return

def diag_vars(mol, confsd):
    """ Returns the variables needed for the MRCI diagonalisation"""

    # Number of irreps
    if mol.comp_sym != 'c1':
        nirrep = molecule.nirrep[mol.sym_indx]
    else: 
        nirrep = 1

    # Bitci MRCI configuration scratch file numbers
    confscr = np.array(confsd.confscr, dtype=int)
    confscr = convert.convert_ctypes(confscr, dtype='int32')

    # Bitci eigenvector scratch file numbers
    vecscr1 = 0
    vecscr1 = ctypes.c_int32(vecscr1)

    # Iteratice diagonalisation algorithm
    if params.mrci_param['diag_algorithm'] == 'gendav':
        ialg = convert.convert_ctypes(1, dtype='int32')
    elif params.mrci_param['diag_algorithm'] == 'blockdav':
        ialg = convert.convert_ctypes(2, dtype='int32')
    else:
        sys.exit(
            '\n Unrecognised iterative diagonalisation algorithm:'
            +' '+var.d3_inp['diag_algorithm'])

    # Maximum no. iterations
    niter = convert.convert_ctypes(params.mrci_param['diag_iter'],
                                dtype='int32')
    
    # Blocksizes
    if len(params.mrci_param['diag_blocksize']) == 0:
        # Blocksizes not given: asign sensible default values
        params.mrci_param['diag_blocksize'] = \
            [int(round(max(n+5, n*1.3))) for n in params.mol_param['nstates']]
    
    blocksize = convert.convert_ctypes(np.array(params.mrci_param['diag_blocksize']),
                                    dtype='int32')
    
    # Deflation flag
    deflate = convert.convert_ctypes(params.mrci_param['diag_deflate'],
                                  dtype='logical')

    # Residual norm convergence threshold
    # If deflation is being used, then this needs to be tight
    if params.mrci_param['diag_deflate'] and params.mrci_param['diag_tol'] > 1e-6:
        print('\n Deflation is being used -> adjusting convergence threshold to 1e-6')
        params.mrci_param['diag_tol'] = 1e-6
    
    tol = convert.convert_ctypes(params.mrci_param['diag_tol'],
                              dtype='double')
    
    return nirrep, confscr, vecscr1, ialg, tol, niter, blocksize, deflate
