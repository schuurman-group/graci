"""
Module for performing the reference space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output

def diag(mol, ci, lib_bitci):
    """Diagonalisation of the reference space Hamiltonian"""

    # length of nstates vector is the number of irreps
    nirr = len(ci.nstates)

    # Construct the molecule object
    timing.start('ref_diag')

    # Print the section header
    output.print_refdiag_header()
    
    # Bitci reference configuration scratch file numbers
    confscr = np.array(ci.ref_conf.confscr, dtype=int)
    confscr = convert.convert_ctypes(confscr, dtype='int32')

    # Numbers of configurations
    nconf = np.array(ci.ref_conf.nconf, dtype=int)
    nconf = convert.convert_ctypes(nconf, dtype='int32')
    
    # Bitci eigenvector scratch file numbers
    vecscr1 = 0
    vecscr1 = ctypes.c_int32(vecscr1)
    vecscr  = []
    
    # Loop over irreps
    for i in range(nirr):
    
        # Number of roots for the current irrep
        nroots = ctypes.c_int32(ci.nstates[i])
    
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
        if nroots.value < ci.nstates[i]:
            ci.nstates[i] = nroots.value
    
    # Retrieve the reference space energies
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
    ci.ref_conf.set_vecscr(vecscr)
    
    # Save the reference space state energies
    ci.ref_conf.set_ener(np.transpose(ener))
    
    # Print the section header
    output.print_refdiag_summary(mol, ci)
    
    # Stop timing
    timing.stop('ref_diag')
    
    return
