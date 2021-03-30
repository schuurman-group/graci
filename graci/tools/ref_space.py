"""
Module for constructing the reference space configurations
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.methods.molecule as molecule
import graci.methods.wavefunction as wavefunction
import graci.io.output as output
import graci.io.convert as convert

def generate(mol, ci, lib_bitci):
    """generate the reference space object"""

    # Initialise the reference space wavefunction object
    conf0 = wavefunction.Wavefunction()

    # Optional automated determination of the RAS MO
    # spaces via the analysis of the DFT/CIS eigenvectors
    if ci.autoras:
        autoras(mol, ci, conf0, lib_bitci)
    
    # Generate the reference space configurations
    genconf(mol, ci, conf0, lib_bitci)
    
    return conf0


def genconf(mol, ci, conf0, lib_bitci):
    """generate the reference space configurations"""
        
    # Number of orbitals in RAS1, RAS2, and RAS3
    m1 = convert.convert_ctypes(ci.ras1.size, dtype='int32')
    m2 = convert.convert_ctypes(ci.ras2.size, dtype='int32')
    m3 = convert.convert_ctypes(ci.ras3.size, dtype='int32')

    # Number of electrons in RAS2
    n2 = 0
    for i in ci.ras2:
        n2 += int(mol.orb_occ[i-1])
    n2 = convert.convert_ctypes(n2, dtype='int32')

    # Maximum number of holes in RAS1
    nh1 = convert.convert_ctypes(ci.nhole1, dtype='int32')    
    
    # Maximum number of electrons in RAS3
    ne3 = convert.convert_ctypes(ci.nelec3, dtype='int32')
    
    # Array of RAS1 orbital indices
    iras1 = np.zeros(mol.nmo, dtype=int)
    k = -1
    for i in ci.ras1:
        k +=1
        iras1[k] = i
    iras1 = convert.convert_ctypes(iras1, dtype='int32')
    
    # Array of RAS2 orbital indices    
    iras2 = np.zeros(mol.nmo, dtype=int)
    k = -1
    for i in ci.ras2:
        k +=1
        iras2[k] = i
    iras2 = convert.convert_ctypes(iras2, dtype='int32')
    
    # Array of RAS3 orbital indices
    iras3 = np.zeros(mol.nmo, dtype=int)
    k = -1
    for i in ci.ras3:
        k +=1
        iras3[k] = i
    iras3 = convert.convert_ctypes(iras3, dtype='int32')

    # Number of irreps
    if mol.comp_sym != 'c1':
        nirrep = molecule.nirrep[mol.sym_indx]
    else:
        nirrep = 1

    # Number of reference space configurations per irrep
    nconf0 = np.zeros(nirrep, dtype=int)
    nconf0 = convert.convert_ctypes(nconf0, dtype='int32')

    # Reference space configurations scratch file number
    confscr = np.zeros(nirrep, dtype=int)
    confscr = convert.convert_ctypes(confscr, dtype='int32')
    
    # Construct the reference space configurations
    lib_bitci.generate_ref_confs(iras1,
                                iras2,
                                iras3,
                                ctypes.byref(nh1),
                                ctypes.byref(m1),
                                ctypes.byref(n2),
                                ctypes.byref(m2),
                                ctypes.byref(ne3),
                                ctypes.byref(m3),
                                nconf0,
                                confscr)
    
    # Convert the nconf0 and confscr ctypes array to lists
    # (note that slicing a ctypes array will automatically
    # produce a list)
    nconf0=nconf0[:]
    confscr=confscr[:]

    # Set the number of reference space configurations
    conf0.set_nconf(nconf0)
    
    # Set the reference space configuration scratch file number
    conf0.set_confscr(confscr)
    
    return


def autoras(mol, ci, conf0, lib_bitci):
    """determination of the RAS subspaces via preliminary
    DFT/CIS calculations"""

    # Print the section header
    output.print_autoras_header()
    
    # Number of irreps
    if mol.comp_sym != 'c1':
        nirrep = molecule.nirrep[mol.sym_indx]
    else:
        nirrep = 1

    # Number of extra roots
    n_extra = 2
    
    #
    # Diagonalisation of the DFT/CIS Hamiltonian
    #
    # Bitci eigenvector scratch file numbers
    vecscr1 = ctypes.c_int32(0)
    vecscr  = []

    # Loop over irreps
    for i in range(nirrep):

        # Number of roots for the current irrep
        nroots = ctypes.c_int32(ci.nstates[i] + n_extra)

        # Irrep number
        irrep = ctypes.c_int32(i)
        
        # Call the the bitci DFT/CIS routine
        lib_bitci.diag_dftcis(ctypes.byref(irrep),
                              ctypes.byref(nroots),
                              ctypes.byref(vecscr1))

        # Bitci eigenvector scratch number
        vecscr.append(vecscr1.value)

    #
    # RAS1 and RAS3 guess
    #
    # Dominant particle/hole indices
    iph  = np.zeros(mol.nmo, dtype=int)
    iph1 = convert.convert_ctypes(np.zeros(mol.nmo, dtype=int),
                               dtype='int32')
    # Loop over irreps
    for i in range(nirrep):

        # Number of roots for the current irrep
        nroots = ctypes.c_int32(ci.nstates[i] + n_extra)

        # Irrep number
        irrep = ctypes.c_int32(i)

        # Eigenpair scratch file number
        vecscr1 = ctypes.c_int32(vecscr[i])
        
        # Get the particle/hole MOs corresponding to the
        # dominant DFT/CIS CSFs
        lib_bitci.ras_guess_dftcis(ctypes.byref(irrep),
                                   ctypes.byref(nroots),
                                   ctypes.byref(vecscr1),
                                   iph1)

        # Update the array of dominant particle/hole indices
        iph = np.maximum(iph, np.array(iph1[:], dtype=int))

    # RAS1 & RAS3 MO indices
    ras1 = []
    ras3 = []
    for n in range(mol.nmo):
        if iph[n] == 1:
            if mol.orb_occ[n] == 0:
                # RAS3 MO
                ras3.append(n+1)
            else:
                # RAS1 MO
                ras1.append(n+1)
    ci.ras1 = np.array(ras1, dtype=int)
    ci.ras3 = np.array(ras3, dtype=int)

    # Output the selected RAS1 & RAS3 spaces
    print('\n Selected RAS MOs:', flush=True)
    print('\n RAS1',[n for n in ci.ras1], flush=True)
    print(' RAS3',[n for n in ci.ras3], flush=True)

    return
