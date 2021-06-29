"""
Module for performing the MRCI space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert

def diag(ci_method):
    """Diagonalisation of the MRCI Hamiltonian"""

    # Start timing
    timing.start('mrci_diag')

    # nirr is given by the length of the nstates vector in ci obj
    nirr = ci_method.n_irrep()

    # the bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # Get the information needed to run the MRCI diagonalisation
    # calculation
    iguess, ialg, tol, niter, blocksize, deflate = diag_vars(ci_method)

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)

    # Bitci ref space eigenvector scratch file numbers
    ref_ciunits = np.array(ci_method.ref_wfn.ci_units, dtype=int)

    # Initialise the list to hold the eigenvector scratch file
    # numbers
    ciunit  = 0
    ciunits = []

    # Q-space energy corrections (only needed if pruning is being
    # used)
    equnits     = np.array(ci_method.mrci_wfn.eq_units, dtype=int)

    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci_method.n_state_sym(irrep)

        # Block size for the current irrep
        nblock = blocksize[irrep]

        args = (irrep, nroots, ci_confunits, ciunit, ialg, tol, niter, 
                blocksize, deflate, iguess, ref_ciunits)
        ciunit = libs.lib_func('diag_mrci', args)

        # Bitci eigenvector scratch number
        ciunits.append(ciunit)

    # Retrieve the MRCI energies
    maxroots = max(ci_method.n_state_sym())
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci_method.n_state_sym(irrep) > 0:

            # Number of roots for the current irrep
            nroots = ci_method.n_state_sym(irrep)

            args = (ciunits[irrep], nroots, ener[irrep,:nroots])
            (ener[irrep,:nroots]) = \
                    libs.lib_func('retrieve_energies', args)
            
    # Retrieve the MRCI eigenvector scratch file names
    ciname = []
    name    = ''
    for irrep in range(nirr):
        args = (ciunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        ciname.append(name)

    # Apply the Q-space energy corrections
    if ci_method.prune != 'off':
        for irrep in range(nirr):
            nstates  = ci_method.n_state_sym()
            qcorr    = np.zeros(nstates[irrep], dtype=float)
            maxovrlp = np.zeros(nstates[irrep], dtype=float)
            nextra = ci_method.nextra['prune'][irrep]
            
            args = (irrep,
                    ciunits[irrep],
                    ref_ciunits[irrep],
                    ci_confunits[irrep],
                    equnits[irrep],
                    nstates[irrep],
                    nextra,
                    qcorr,
                    maxovrlp)
            (qcorr, maxovrlp) = libs.lib_func('retrieve_qcorr', args)

            ener[irrep,:nstates[irrep]] += qcorr
            
    # Print the report of the MRCI states
    ciunits = np.array(ciunits, dtype=int)
    nstates = ci_method.n_state_sym()
    nextra  = np.array(ci_method.nextra['prune'], dtype=int)
    if ci_method.prune == 'off':
        args = (ci_confunits, ciunits, nstates)
        libs.lib_func('print_mrci_states', args)
    else:
        args = (ci_confunits, ciunits, ref_ciunits, equnits, nstates,
                nextra)
        libs.lib_func('print_pmrci_states', args)
        
    # Stop timing
    timing.stop('mrci_diag')
    
    return ciunits, ciname, ener 

def diag_vars(ci_method):
    """ Returns the variables needed for the MRCI diagonalisation"""

    # Guess vector generation
    if ci_method.diag_guess == 'subdiag':
        iguess = 1
    elif ci_method.diag_guess == 'enpt2':
        iguess = 2
    else:
        sys.exit(
            '\n Unrecognised guess vector type:'
            +' '+ci_method.diag_guess)

    # Iterative diagonalisation algorithm
    if ci_method.diag_method == 'gendav':
        ialg = 1
    elif ci_method.diag_method == 'blockdav':
        ialg = 2
    else:
        sys.exit(
            '\n Unrecognised iterative diagonalisation algorithm:'
            +' '+ci_method.diag_method)

    # Maximum no. iterations
    niter = ci_method.diag_iter

    # Blocksizes
    if len(ci_method.diag_blocksize) == 0:
        # Blocksizes not given: asign sensible default values
        ci_method.diag_blocksize = [int(round(max(n+5, n*1.3))) 
                                    for n in ci_method.nstates]
    
    blocksize = np.array(ci_method.diag_blocksize, dtype=int)

    # Deflation flag
    deflate = ci_method.diag_deflate

    # Residual norm convergence threshold
    # If deflation is being used, then this needs to be tight
    if ci_method.diag_deflate and ci_method.diag_tol > 1.e-6:
        print('\n Deflation is being used -> '+
                 'adjusting convergence threshold to 1e-6')
        ci_method.diag_tol = 1.e-6

    # If this is the first MRCI iteration, then use a loose
    # convergence threshold
    if ci_method.niter == 1 and not ci_method.diag_deflate:
        tol = 1.e-3
    else:
        tol = ci_method.diag_tol

    return iguess, ialg, tol, niter, blocksize, deflate

