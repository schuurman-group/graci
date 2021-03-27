"""Module for performing file operations"""

import sys
import ctypes
import numpy as np
import graci.utils.constants as constants
import graci.utils.timing as timing
import graci.molecule.molecule as molecule
import graci.methods.params as params

#
def print_header(mol):
    """print the output log file header"""

    header =(" -----------------------------------------------------\n"+
             "                                                     \n"+
             " GRaCI                                               \n"+
             " General Reference Configuration Interaction         \n"+
             "                                                     \n"+
             " Authors: Simon Neville and Michael Schuurman        \n"+
             "                                                     \n"+
             " -----------------------------------------------------\n")
    
    inp_key =" Input Parameters \n"

    # Read input file. Small enough to gulp the whole thing
    with open(params.d3_inp['out_file'], 'w') as outfile:
        print(header+'\n\n', flush=True)
        
        print(inp_key+'\n', flush=True)
        
        print(' $geometry section\n -----------------\n', flush=True)
        ostr = ''
        for iatm in range(mol.n_atoms()):
            cstr = '   '.join(['{:10.4f}'.format(mol.geom[iatm,j])
                   for j in range(3)])
            ostr = ' ' + str(mol.atoms[iatm]) + cstr
            print(ostr+'\n', flush=True)
            
        print('\n $molecule section\n -----------------\n', flush=True)
        ostr = ''
        for kword in params.mol_param:
            ostr += ' '+kword.ljust(12,' ')+' = '+str(params.mol_param[kword])+'\n'
        print(ostr, flush=True)
        
        print('\n $scf section\n ------------\n', flush=True)
        ostr = ''
        for kword in params.scf_param:
            ostr += ' '+kword.ljust(12,' ')+' = '+str(params.scf_param[kword])+'\n'
        print(ostr, flush=True)
        
        print('\n $mrci section\n -------------\n', flush=True)
        ostr = ''
        for kword in params.mrci_param:
            ostr += ' '+kword.ljust(12,' ')+' = '+str(params.mrci_param[kword])+'\n'
        print(ostr, flush=True)
        
        print('\n', flush=True)
        print(' Full symmetry:     '+str(mol.full_sym)+'\n', flush=True)
        print(' Abelian sub-group: '+str(mol.comp_sym)+'\n', flush=True)
        print('\n', flush=True)
        
    return

#
def print_scf_header():
    """print the SCF header"""

    print(' SCF Computation with PySCF\n', flush=True)
    print(' -----------------------------\n\n', flush=True)
    
    return

#
def print_scf_summary(scf_energy, mol):
    """print summary of the SCF computation"""

    print(' DFT energy = {:16.10f}\n\n'.format(scf_energy), flush=True)
    print(' Orbital Energies and Occupations\n', flush=True)
    print(' --------------------------------\n', flush=True)
    print(' {:>5}  {:>9}  {:>10}  {:>8}\n'.
          format('Index','Symmetry','Energy','Occ'), flush=True)

    if mol.comp_sym != 'c1':
        orb_cnt = np.zeros(molecule.nirrep[mol.sym_indx], dtype=int)
    else:
        orb_cnt = np.zeros(1, dtype=int)

    for iorb in range(mol.nmo):
        sym_indx = mol.orb_sym[iorb]
        orb_cnt[sym_indx] += 1
        print(' {:5d}  {:>4d}({:>3})  {:10.5f}  {:8.4}\n'.
              format(iorb+1,
                     orb_cnt[sym_indx],
                     mol.orb_irrep[iorb],
                     mol.orb_ener[iorb],
                     mol.orb_occ[iorb]),
              flush=True)
    
    return

#
def print_refdiag_header():
    """print the reference space diagonalisation header"""
    print('\n Reference Space Diagonalisation\n', flush=True)
    print(' -------------------------------', flush=True)
    return

#
def print_refdiag_summary(mol, refdets):
    """print the summary of the reference space diagonalisation"""

    if mol.comp_sym != 'c1':
        nirrep = molecule.nirrep[mol.sym_indx]
    else:
        nirrep = 1

    mine = np.amin(refdets.ener)

    print('\n Reference state energies', flush=True)
    print(' -------------------------', flush=True)
    
    for i in range(nirrep):
        if params.d3_inp['nstates'][i] > 0:
            print('\n', flush=True)
            for n in range(params.d3_inp['nstates'][i]):
                print(' {:<3d} {:3} {:10.6f} {:10.6f}'
                      .format(n+1, mol.irreplbl[i],
                             refdets.ener[n][i],
                             (refdets.ener[n][i]-mine)*constants.au2ev),
                      flush=True)
    
    return

#
def print_mrcispace_header():
    """print the MRCI space generation header"""
    print('\n MRCI Configuration Generation\n', flush=True)
    print(' -------------------------------', flush=True)

#
def print_autoras_header():
    """print the automatic RAS space generation header"""
    print('\n Automatic RAS Space Generation\n', flush=True)
    print(' -------------------------------', flush=True)
    
#
def print_cleanup():
    """shutdown the timers and print timing information"""

    ostr = timing.print_timings()

    print(ostr, flush=True)
        
    return
