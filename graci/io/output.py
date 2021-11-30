"""Module for performing file operations"""

import sys
import ctypes
import os as os
import numpy as np
from contextlib import contextmanager
import graci.utils.constants as constants
import graci.utils.timing as timing
import graci.core.params as params

#
file_names = {'input_file'   : '',
              'out_file'     : '',
              'chkpt_file'   : '',
              'pyscf_out'    : ''}

orb_formats = {'molden', 'gamess'}

@contextmanager
def output_file(file_name, mode):
    """return the file handle if file_name string is not None, else
       return stdout"""

    if file_name is None:
        yield sys.stdout
    else:
        with open(file_name, mode) as out_file:
            yield out_file

#
def print_message(message_str):
    """Print a message string to the output file"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('ATTN: '+str(message_str)+'\n')

    return


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

            calc_name = type(calc_obj).__name__
            outfile.write('\n $'+calc_name+' section\n ------------\n')

            # if geometry object, we need to pull out the 
            # geometry
            if calc_name == 'geometry':
                gm  = calc_obj.geom()
                atm = calc_obj.atoms()

                for iatm in range(len(atm)):
                    cstr = '   '.join(['{:12.8f}'.format(gm[iatm,j]) 
                               for j in range(3)])
                    ostr = ' ' + str(atm[iatm]) + cstr
                    outfile.write(ostr+'\n')

            # outfile.write the keyword input
            ostr = '\n'
            for kword in params.kwords[calc_name].keys():
                ostr += ' '+kword.ljust(12,' ')+\
                        ' = '+str(getattr(calc_obj,kword))+'\n'
            outfile.write(ostr)

        outfile.write('\n Symmetry Information\n ---------------\n')
        for calc_obj in run_list:
            # pull out the molecule objets and print symmetry
            # information for each molecule
            if type(calc_obj).__name__ == 'Molecule':

                outfile.write(' $molecule '+str(calc_obj.label)+'\n')
                outfile.write(' Full symmetry:     '+
                        str(calc_obj.full_sym)+'\n')
                outfile.write(' Abelian sub-group: '+
                        str(calc_obj.comp_sym)+'\n')
                outfile.write('\n')
                outfile.flush()
        
    return

#
def print_scf_header(scf):
    """print the SCF header"""
    global file_names

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'SCF Computation with PySCF, label = '+str(scf.label)
        lpad = int(0.5*(max(0,LLEN-len(title))))
        pstr = str('*'.ljust(lpad)+title)
        pstr = pstr.ljust(LLEN-1)+'*'

        outfile.write('\n\n '+str('*'*LLEN))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str(pstr))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')                                   
        outfile.write(  '\n '+str('*'*LLEN))

        if scf.restart:
            outfile.write('\n\n **** RESTART ACTIVATED ****\n\n')
            outfile.write(' Extracting SCF result from checkpoint file:'
                          + str(file_names['chkpt_file'])+'\n\n')

        outfile.flush()
      
    return

#
def print_scf_summary(scf):
    """print summary of the SCF computation"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        if scf.mol.use_df:
            outfile.write('\n\n density fitting employed, basis: '+
                            str(scf.auxbasis)+'\n\n')

        outfile.write(' REF energy        = {:16.10f}\n'.format(scf.energy))
        outfile.write(' Nuclear Repulsion = {:16.10f}\n\n'.format(scf.mol.enuc))
        outfile.write(' Orbital Energies and Occupations\n')
        outfile.write(' --------------------------------\n')
        outfile.write(' {:>5}  {:>9}  {:>10}  {:>8}\n'.
              format('Index','Symmetry','Energy','Occ'))

        orb_cnt = np.zeros(scf.mol.n_irrep(), dtype=int)

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
def print_dftmrci_header(label):
    global file_names

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'DFT/MRCI computation, label = '+str(label)
        lpad = int(0.5*(max(0,LLEN-len(title))))
        pstr = str('*'.ljust(lpad)+title)
        pstr = pstr.ljust(LLEN-1)+'*'

        outfile.write('\n\n '+str('*'*LLEN))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str(pstr))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')                                   
        outfile.write(  '\n '+str('*'*LLEN))
        outfile.flush()

    return

#
def print_dftmrenpt2_header(label):
    global file_names

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'DFT/MR-ENPT2 computation, label = '+str(label)
        lpad = int(0.5*(max(0,LLEN-len(title))))
        pstr = str('*'.ljust(lpad)+title)
        pstr = pstr.ljust(LLEN-1)+'*'

        outfile.write('\n\n '+str('*'*LLEN))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str(pstr))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str('*'*LLEN))
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
def print_refdiag_summary(ci_method):
    """print the summary of the reference space diagonalisation"""
    global file_names

    mine = np.amin(ci_method.ref_ener)

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n Reference state energies')
        outfile.write('\n -------------------------')
    
        nirr = len(ci_method.nstates)
        for i in range(nirr):
            if ci_method.n_state_sym(i) > 0:
                outfile.write('\n')
                for n in range(ci_method.n_state_sym(i)):
                    outfile.write('\n {:<3d} {:3} {:10.6f} {:10.6f}'
                        .format(n+1, ci_method.scf.mol.irreplbl[i],
                        ci_method.ref_ener[i,n],
                        (ci_method.ref_ener[i,n]-mine)*constants.au2ev))
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
def print_dftmrci_states_header():
    """print the DFT/MRCI eigenstate report header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n DFT/MRCI Eigenstates\n')
        outfile.write(' -------------------------------')
        outfile.flush()

#
def print_dftmrenpt2_states_header():
    """print the DFT/MR-ENPT2 eigenstate report header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n DFT/MR-ENPT2 Eigenstates\n')
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

#
def print_moments_header():
    """prints out header for moments section"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        ostr =    '\n\n\n --------------------------------------------'
        ostr = ostr + '\n Dipole and Quadrupole Moments'
        ostr = ostr + '\n --------------------------------------------'
        outfile.write(ostr)

    return

#
def print_moments(st, irr, mu, q2):
    """prints out the dipole moment vector"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        st_str = str(st+1)+' ('+str(irr)+')'

        outfile.write('\n state: '+st_str)
        outfile.write(  '\n ----------')
        ostr = '\n dipole vector[au]: {:10.6f} {:10.6f} {:10.6f} '+\
                'total: {:10.6f} D'
        outfile.write(ostr.format(mu[0],mu[1],mu[2],
                        np.linalg.norm(mu)*constants.au2debye))

        ostr = '\n <'+st_str+' | r^2 | '+st_str+'> (a.u.): {:10.6f}\n'
        outfile.write(ostr.format(q2))

    return

#
def print_quad(st, irr, qtensor):
    """prints out the quadrupole moment tensor"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n\n quadrupole moment tensor -- \n')
        ostr = '\n {:10.6f}   {:10.6f}   {:10.6f}'
        for i in range(3):
            outfile.write(ostr.format(qtensor[i,0],qtensor[i,1],qtensor[i,2]))
        
    return

#
def print_transition_table(init_st, init_sym, final_st, final_sym, 
                                   exc_ener, f0l, f2l, f0v, f2v, f0xyz): 
    """print out the summary files for the transition moments"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        # print a table of oscillator strengths and transition 
        # dipole vectors

        outfile.write('\n\n\n ------------------------------------------')
        outfile.write(    '\n Transition Properties')
        outfile.write(    '\n ------------------------------------------')

        header  = '\n\n  Initial     Final    Exc Ener                  '
        header += '                        Oscillator Strength (V)\n'
        undr_str = '-' * (len(header))
        header += '  State       State      (eV)      f0(L)    f2(L)'
        header += '    f0(V)    f2(V)       x        y        z'

        fstr   = '\n {:3d}({:>3}) -> {:3d}({:>3}) {:7.2f}'+ \
                    ' {:9.4f}{:9.4f}{:9.4f}{:9.4f}  {:9.4f}{:9.4f}{:9.4f}'

        outfile.write(header)
        outfile.write('\n '+undr_str)
        outfile.write('\n')

        for i in range(len(final_st)):
            outfile.write(fstr.format(
                             init_st,
                             init_sym,
                             final_st[i],
                             final_sym[i],
                             exc_ener[i]*constants.au2ev, 
                             f0l[i], 
                             f2l[i], 
                             f0v[i], 
                             f2v[i],
                             *f0xyz[i][:]))

    return

