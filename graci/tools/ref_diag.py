"""
Module for performing the reference space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

def diag(mol, nstates, conf0, lib_bitci):
    """Diagonalisation of the reference space Hamiltonian"""

    # length of nstates vector is the number of irreps
    nirr = len(nstates)

    # Construct the molecule object
    timing.start('ref_diag')

    # Print the section header
    output.print_refdiag_header()
    
    # Bitci reference configuration scratch file numbers
    confscr = np.array(conf0.confscr, dtype=int)
    confscr = convert.convert_ctypes(confscr, dtype='int32')

    # Numbers of configurations
    nconf = np.array(conf0.nconf, dtype=int)
    nconf = convert.convert_ctypes(nconf, dtype='int32')
    
    # Bitci eigenvector scratch file numbers
    vecscr1 = 0
    vecscr1 = ctypes.c_int32(vecscr1)
    vecscr  = []
    
    # Loop over irreps
    for i in range(nirr):
    
        # Number of roots for the current irrep
        nroots = ctypes.c_int32(nstates[i])
    
        # Irrep number
        irrep = ctypes.c_int32(i)
        
        # Call to the bitci reference space diagonalisation routine
        lib_bitci.ref_diag_mrci(ctypes.byref(irrep),
                                ctypes.byref(nroots),
                                confscr,
                                nconf,
                                ctypes.byref(vecscr1))
        
        # Bitci eigenvector scratch number
        vecscr.append(vecscr1.value)
    
        # If the number of reference space configurations for the
        # current irrep is less than the requested number of roots
        # then reset nstates accordingly
        if nroots.value < nstates[i]:
            nstates[i] = nroots.value
    
    # Retrieve the reference space energies
    maxroots = max(nstates)
    ener = np.zeros((nirr, maxroots), dtype=float)
    for i in range(nirr):
        if nstates[i] > 0:
    
            # Energies
            ener1 = np.zeros(nstates[i], dtype=float)
            ener1 = convert.convert_ctypes(ener1, dtype='double')
    
            # Number of roots for the current irrep
            nroots = ctypes.c_int32(nstates[i])
    
            # Eigenvector Scratch file number
            vecscr1 = ctypes.c_int32(vecscr[i])
    
            # Retrieve the energies
            lib_bitci.retrieve_energies(ctypes.byref(vecscr1),
                                        ctypes.byref(nroots),
                                        ener1)
            ener[i][:nstates[i]] = ener1[:]
            
    # Save the list of bitci eigenvector scratch numbers
    conf0.set_vecscr(vecscr)
    
    # Save the reference space state energies
    conf0.set_ener(np.transpose(ener))
    
    # Print the section header
    output.print_refdiag_summary(mol, nstates, conf0)
    
    # Stop timing
    timing.stop('ref_diag')
    
    return
