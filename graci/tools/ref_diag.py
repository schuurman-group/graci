"""
Module for performing the reference space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.methods.params as params
import graci.molecule.molecule as molecule
import graci.io.convert
import graci.io.output

def diag(mol, conf0, lib_bitci):
    """Diagonalisation of the reference space Hamiltonian"""

    # Construct the molecule object
    timing.start('ref_diag')

    # Print the section header
    output.print_refdiag_header()
    
    # Number of irreps
    if mol.comp_sym != 'c1':
        nirrep = molecule.nirrep[mol.sym_indx]
    else:
        nirrep = 1
    
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
    for i in range(nirrep):
    
        # Number of roots for the current irrep
        nroots = ctypes.c_int32(params.mol_param['nstates'][i])
    
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
        if nroots.value < params.mol_param['nstates'][i]:
            params.mol_param['nstates'][i] = nroots.value
    
    # Retrieve the reference space energies
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
    conf0.set_vecscr(vecscr)
    
    # Save the reference space state energies
    conf0.set_ener(np.transpose(ener))
    
    # Print the section header
    output.print_refdiag_summary(mol, conf0)
    
    # Stop timing
    timing.stop('ref_diag')
    
    return
