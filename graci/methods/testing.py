"""
Module for computing DFT/MRCI energies
"""
import graci.methods.scf as scf
import graci.methods.params as params
import graci.io.output as output

def energy(nroots):
""" compute the DFT/MRCI energy for nroots """

    # evaluate the KS DFT energy
    scf.run(mol)

    # if this is an interface file generation run,
    # then write it and finish here
    write_interface(mol)
    output.write_cleanup()
    sys.exit()

    return 

def density(states):
""" computes the density matrices for the states in 
    the array 'states'"""

    return

def slater_dets(state)
""" return the slater determinant list for state 'state'"""

    return

def csf(state)
""" return the CSF list for state 'state'"""

    return

################################################################

def write_interface(mol):
    """writes an ASCII bitci interface file"""

    # Open the output file
    outfile = open('ci.info', 'w')

    # No. electrons
    outfile.write('-- Nel\n')
    outfile.write(str(mol.nel)+'\n')

    # Multiplicity
    outfile.write('\n-- Multiplicity\n')
    outfile.write(str(mol.mult)+'\n')

    # Point group index
    isym  = mol.sym_indx + 1 if mol.sym_indx > 0 else 1
    outfile.write('\n-- Point_group\n')
    outfile.write(str(isym)+'\n')

    # No. MOs
    outfile.write('\n-- Nmo\n')
    outfile.write(str(mol.nmo)+'\n')

    # Density fitting
    outfile.write('\n-- Density_fitting\n')
    outfile.write(str(params.mol_param['use_df'])+'\n')

    # No. auxiliary basis functions
    outfile.write('\n-- Naux\n')
    outfile.write(str(mol.naux)+'\n')

    # MO irreps and energies
    outfile.write('\n-- MO_info\n')
    for i in range(len(mol.orb_sym)):
        outfile.write(str(mol.orb_sym[i])+'  '+str(mol.orb_ener[i])+'\n')

    # Nuclear repulsion energy
    outfile.write('\n-- Enuc\n')
    outfile.write(str(mol.enuc)+'\n')

    # Integrals threshold
    outfile.write('\n-- Integrals_threshold\n')
    outfile.write('1E-14'+'\n')

    # Hamiltonian index
    outfile.write('\n-- Hamiltonian\n')
    outfile.write(params.mrci_param['hamiltonian']+'\n')

    # No. roots per irrep
    outfile.write('\n-- Nroots\n')
    for i in range(len(params.mrci_param['nstates'])):
        outfile.write(str(params.mrci_param['nstates'][i])+'\n')

    # RAS spaces
    n2 = 0
    for i in params.mrci_param['ras2']:
        n2 += int(mol.orb_occ[i-1])
    outfile.write('\n-- RAS\n')
    outfile.write(str(len(params.mrci_param['ras1']))
                  +' '+str(params.mrci_param['nhole1'])+'\n')
    outfile.write(str(len(params.mrci_param['ras2']))+' '+str(n2)+'\n')
    outfile.write(str(len(params.mrci_param['ras3']))
                  +' '+str(params.mrci_param['nelec3'])+'\n')
    if (len(params.mrci_param['ras1']) > 0):
        outfile.write(str(params.mrci_param['ras1'])[1:-1]+'\n')
    else:
        outfile.write('-\n')
    if (len(params.mrci_param['ras2']) > 0):
        outfile.write(str(params.mrci_param['ras2'])[1:-1]+'\n')
    else:
        outfile.write('-\n')
    if (len(params.mrci_param['ras3']) > 0):
        outfile.write(str(params.mrci_param['ras3'])[1:-1]+'\n')
    else:
        outfile.write('-\n')

    return


