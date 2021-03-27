"""
Module for computing DFT/MRCI energies
"""
import graci.methods.scf as scf
import graci.tools.init_libs as init_libs
import graci.tools.ref_space as ref_space
import graci.tools.ref_diag as ref_diag
import graci.tools.mrci_space as mrci_space
import graci.tools.mrci_diag as mrci_diag
import graci.tools.mrci_refine as mrci_refine

def energy(nroots):
""" compute the DFT/MRCI energy for nroots """

    # run the KS-DFT computation 
    scf.run(mol)

    # initialize int_pyscf
    lib_intpyscf = init_libs.init_intpyscf(mol)

    # initialize bitci
    lib_bitci = init_libs.init_bitci(mol)

    # generate the reference space configurations
    conf0 = ref_space.generate(mol, lib_bitci)

    # Perform the MRCI iterations, refining the reference space
    # as we go
    for i in range(var.d3_inp['refiter']):
        # reference space diagonalisation
        ref_diag.diag(mol, conf0, lib_bitci)

        # generate the MRCI configurations
        confsd = mrci_space.generate(mol, conf0, lib_bitci)

        # MRCI diagonalisation
        mrci_diag.diag(mol, confsd, lib_bitci)

        # refine the reference space
        min_norm = mrci_refine.refine_ref_space(mol, conf0, confsd, lib_bitci)

        # break if the reference space is converged
        if min_norm > 0.95 and i > 0:
            print('\n * Reference Space Converged *', flush=True)
            break


    return 


def density(states):
""" computes the density matrices for the states in 
    the array 'states'"""

    return

def slater_dets(state)
""" return the slater determinant list for state 'state'"""

    return

def csfs(state)
""" return the CSF list for state 'state'"""

    return

##########################################################################



