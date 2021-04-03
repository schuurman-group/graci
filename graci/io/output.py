"""Module for performing file operations"""

import sys
import ctypes
import numpy as np
from contextlib import contextmanager
import graci.utils.constants as constants
import graci.utils.timing as timing
import graci.methods.params as params
import graci.methods.molecule as molecule

#
file_names = {'input_file'   : '',
              'out_file'     : '',
              'pyscf_out'    : '',
              '1ei'          : '',
              '2ei'          : ''}

@contextmanager
def output_file(file_name, mode):
    """return the file handle if file_name string is not None, else
       return stdout"""

    if file_name is None:
        yield sys.stdout
    else:
        with open(file_name, mode) as out_file:
            yield out_file

def print_header(run_list):
    """print the output log file header"""
    global file_names

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
    with output_file(file_names['out_file'], 'w') as outfile:
        outfile.write(header+'\n\n')
        outfile.write(inp_key+'\n')

        for calc_obj in run_list:

            calc_type = params.method_name(calc_obj)

            # if molecule object, we need to pull out the 
            # geometry
            if calc_type == 'molecule':
                mol = calc_obj
                outfile.write(' $geometry section\n -----------------\n')
                for iatm in range(mol.n_atoms()):
                    cstr = '   '.join(['{:10.4f}'.format(mol.geom[iatm,j])
                           for j in range(3)])
                    ostr = ' ' + str(mol.atoms[iatm]) + cstr
                    outfile.write(ostr+'\n')
            
            # else just outfile.write the keyword input
            outfile.write('\n $'+calc_type+' section\n ------------\n')
            ostr = ''
            for kword in params.kwords[calc_type].keys():
                ostr += ' '+kword.ljust(12,' ')+\
                        ' = '+str(getattr(calc_obj,kword))+'\n'
            outfile.write(ostr)

        outfile.write('\n')
        outfile.write(' Full symmetry:     '+str(mol.full_sym)+'\n')
        outfile.write(' Abelian sub-group: '+str(mol.comp_sym)+'\n')
        outfile.write('\n')
        outfile.flush()
        
    return

#
def print_scf_header():
    """print the SCF header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write(' SCF Computation with PySCF\n')
        outfile.write(' -----------------------------\n\n')
        outfile.flush()

    return

#
def print_scf_summary(mol, scf):
    """print summary of the SCF computation"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write(' SCF energy = {:16.10f}\n\n'.format(scf.energy))
        outfile.write(' Orbital Energies and Occupations\n')
        outfile.write(' --------------------------------\n')
        outfile.write(' {:>5}  {:>9}  {:>10}  {:>8}\n'.
              format('Index','Symmetry','Energy','Occ'))

        orb_cnt = np.zeros(mol.n_irrep(), dtype=int)

        for iorb in range(scf.nmo):
            sym_indx = scf.orb_sym[iorb]
            orb_cnt[sym_indx] += 1
            outfile.write(' {:5d}  {:>4d}({:>3})  {:10.5f}  {:8.4}\n'.
                  format(iorb+1,
                         orb_cnt[sym_indx],
                         scf.orb_irrep[iorb],
                         scf.orb_ener[iorb],
                         scf.orb_occ[iorb]))
        outfile.flush()

    return

#
def print_refdiag_header():
    """print the reference space diagonalisation header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n Reference Space Diagonalisation\n')
        outfile.write(' -------------------------------')
        outfile.flush()

    return

#
def print_refdiag_summary(mol, ci):
    """print the summary of the reference space diagonalisation"""
    global file_names

    mine = np.amin(ci.ref_conf.ener)

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n Reference state energies')
        outfile.write('\n -------------------------')
    
        for i in range(len(ci.nstates)):
            if ci.nstates[i] > 0:
                outfile.write('\n')
                for n in range(ci.nstates[i]):
                    outfile.write('\n {:<3d} {:3} {:10.6f} {:10.6f}'
                          .format(n+1, mol.irreplbl[i],
                            ci.ref_conf.ener[n][i],
                            (ci.ref_conf.ener[n][i]-mine)*constants.au2ev))
        outfile.write('\n')
        outfile.flush()

    return

#
def print_mrcispace_header():
    """print the MRCI space generation header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n MRCI Configuration Generation\n')
        outfile.write(' -------------------------------')
        outfile.flush()

#
def print_autoras_header():
    """print the automatic RAS space generation header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n Automatic RAS Space Generation\n')
        outfile.write(' -------------------------------')
        outfile.flush()
    
#
def print_cleanup():
    """shutdown the timers and print timing information"""
    global file_names

    ostr = timing.print_timings()

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write(ostr)
        outfile.flush()
        
    return
