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
        outfile.flush()

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
    calc_types = [type(calc_obj).__name__ for calc_obj in run_list]

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

        if 'Molecule' in calc_types:
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

def print_rydano_header(rydano):
    """
    print header for Rydberg orbital generation
    """
    global file_names

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'Rydberg Orbital Generation, label = '+str(rydano.label)
        lpad = int(0.5*(max(0,LLEN-len(title))))
        pstr = str('*'.ljust(lpad)+title)
        pstr = pstr.ljust(LLEN-1)+'*'

        outfile.write('\n\n '+str('*'*LLEN))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str(pstr))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')                                   
        outfile.write(  '\n '+str('*'*LLEN))

        ostr = '\n\n -- Automatic Construction of an ANO basis for ' + \
               ' Rydberg states using KBJ\n' + \
               '    uncontracted primitives and a density '+ \
               ' constructed from the low-lying\n' + \
               '    virtual orbitals of the corresponding cation.\n'

        ostr+= '\n       See: K. Kaufmann, W. Baumeistert, M. Jungen'+\
               '\n            J. Phys. B: At. Mol. Opt. Phys., 22, '+\
               ' 2223-2240, (1989).\n\n'

        outfile.write(ostr)
        outfile.flush()

    return  

def print_rydano_summary(exps, occ, nos, con_str):
    """print the result of the Rydberg ANO basis construction

    Args:
        prims:    the KBJ primitives used in the expansion
        occ:      the 'occupations' of the NOs 
        nos:      the contracted ANOs
        contract: the requested contraction to be added to basis

    Returns:
        None
    """

    with output_file(file_names['out_file'], 'a+') as outfile:
        ostr  = '\n Summary of ANO basis set generation\n'
        ostr += ' '+'-'*35

        ostr += '\n\n Exponents\n ----------\n'
        for l in range(len(exps)):
            ostr += '\n  l = '+str(l)+' | '+ \
                    ''.join(['{:10.6f}'.format(exps[l][i]) \
                      for i in range(len(exps[l]))])        

        ostr += '\n\n Contractions\n ------------'
        for l in range(len(nos)):
            ncon = len(occ[l])
            ostr += '\n\n angular momentum l='+str(l)+'\n'
            ostr += ' occ:'+' '.join(['{:10.4f}'.format(occ[l][i]) 
                                               for i in range(ncon)])
            ostr += '\n '+'-'*(13*ncon+5)+'\n'
            for n in range(nos[l].shape[0]):
                ostr += '\n     '
                ostr += ' '.join(['{:10.4f}'.format(nos[l][n,i]) 
                                                 for i in range(ncon)])
          
        ostr += '\n\n Contractions to be applied: '+str(con_str)

        outfile.write(ostr)
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
            outfile.write(' density fitting employed, basis: '+
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
            outfile.write(' {:5d}  {:>4d}({:>3})  {:10.5f}  {:8.4f}\n'.
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
def print_dftmrci2_header(label):
    global file_names

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'DFT/MRCI(2) computation, label = '+str(label)
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
            if ci_method.n_states_sym(i) > 0:
                outfile.write('\n')
                for n in range(ci_method.n_states_sym(i)):
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
def print_dftmrci_states_header(prune):
    """print the DFT/MRCI eigenstate report header"""
    global file_names

    title = 'DFT/MRCI Eigenstates'
    if prune:
        title = 'p-'+title
    
    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n '+title+'\n')
        outfile.write(' -------------------------------')
        outfile.flush()

#
def print_dftmrci2_states_header():
    """print the DFT/MRCI(2) eigenstate report header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n DFT/MRCI(2) Eigenstates\n')
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
def print_promotion(ref, states, irr, pd, pa):
    """print attachment and detachment numbers for a set of states,
       relative to state 'ref'"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        ostr =      '\n\n --------------------------------------------'
        ostr = ostr + '\n Promotion Numbers'
        ostr = ostr + '\n --------------------------------------------'
        outfile.write(ostr)

        ostr = '\n\n {:8s} -> {:8s}  {:10s}  {:10s}'
        outfile.write(ostr.format('state'.rjust(8), 
                                  'state'.rjust(8),
                                  'P(detach)'.rjust(10), 
                                  'P(attach)'.rjust(10)))

        for ist in range(len(states)):
            st_str1 = str(states[ref]+1)+' ('+str(irr[ref])+')'
            st_str2 = str(states[ist]+1)+' ('+str(irr[ist])+')'
            ostr = '\n {:8s} -> {:8s}  {:10.4f}  {:10.4f}'
            outfile.write(ostr.format(st_str1.rjust(8),
                                      st_str2.rjust(8),
                                      pd[ist], pa[ist]))

    return

#
def print_moments(states, irr, momts):
    """prints out the dipole moment vector"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        ostr =    '\n\n\n --------------------------------------------'
        ostr = ostr + '\n Electronic Moments'
        ostr = ostr + '\n --------------------------------------------'
        outfile.write(ostr)

        ostr =      '\n\n Dipole Moments'
        ostr = ostr + '\n --------------'
        outfile.write(ostr)
        ostr =      '\n\n Dipole moment components are given in a.u.,'+ \
                         ' the total moment is in Debye'
        outfile.write(ostr)
        ostr = '\n\n {:8s} {:10s} {:10s} {:10s}    {:10s}'
        outfile.write(ostr.format('state'.rjust(8), 'x   '.rjust(10),
                                  'y   '.rjust(10), 'z   '.rjust(10),
                                  'total (D)'.rjust(10)))

        for ist in range(len(states)):
            st_str = str(states[ist]+1)+' ('+str(irr[ist])+')'
            mu_vec = momts.dipole(states[ist])
            mu     = np.linalg.norm(mu_vec)*constants.au2debye

            ostr = '\n {:8s} {:10.6f} {:10.6f} {:10.6f}    {:10.6f}'
            outfile.write(ostr.format(st_str.rjust(8), mu_vec[0], 
                                      mu_vec[1], mu_vec[2], mu))

        ostr =      '\n\n Quadrupole Moments'
        ostr = ostr + '\n ------------------'
        outfile.write(ostr)
        ostr =  '\n\n All quadrupole tensor elements are given in a.u.'
        ostr += '\n Qij is traceless: Qij = 3rirj - r^2 * delta[i,j]'
        outfile.write(ostr)
        ostr = '\n\n {:8s} {:9s} {:9s} {:9s} {:9s} {:9s} {:9s} {:9s}' 
        outfile.write(ostr.format('state'.rjust(8), 'xx  '.rjust(9), 
                                  'xy  '.rjust(9),  'xz  '.rjust(9), 
                                  'yy  '.rjust(9),  'yz  '.rjust(9), 
                                  'zz  '.rjust(9),  'r^2 '.rjust(9)))

        for ist in range(len(states)):
            st_str = str(states[ist]+1)+' ('+str(irr[ist])+')'
            q_tens = momts.quadrupole(states[ist])
            q2     = momts.second_moment(states[ist])

            ostr = '\n {:8s} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f} {:9.4f}'
            outfile.write(ostr.format(st_str.rjust(8), q_tens[0,0], 
                                      q_tens[0,1], q_tens[0,2],
                                      q_tens[1,1], q_tens[1,2],
                                      q_tens[2,2], q2))

    return

#
def print_transition_header(label):
    """ print out Transition section header"""

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'Transition, label = '+str(label)
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
def print_transition_table(init_st, init_sym, final_st, final_sym, 
                            exc_ener, f0l, f2l, f0v, f2v, f0xyz, promo): 
    """print out the summary files for the transition moments"""

    with output_file(file_names['out_file'], 'a+') as outfile:
        # print a table of oscillator strengths and transition 
        # dipole vectors

        fstr = '\n\n Transitions, initial state = {:3d}({:>3})'
        outfile.write(fstr.format(init_st, init_sym))
        outfile.write(    '\n ------------------------------------------')

        header  = '\n\n  Initial     Final    Exc Ener                  '
        header += '             Oscillator Strength (V)'
        if len(promo) == len(final_st):
            header+= '    Promotion Numbers\n'

        undr_str = '-' * (len(header))

        header += '  State       State      (eV)     f0(L)    f0(V)'
        header += '    f2(V)       x        y        z'
        if len(promo) == len(final_st):
            header += '      attach     detach'

        #f2(L) is actually mixed gauge and often gives nonsensical results:
        # removing it for now
        #header += '  State       State      (eV)      f0(L)    f2(L)'
        #header += '    f0(V)    f2(V)       x        y        z'

        fstr   = '\n {:3d}({:>3}) -> {:3d}({:>3}) {:7.2f}'+ \
                ' {:9.4f}{:9.4f}{:9.4f}  {:9.4f}{:9.4f}{:9.4f} {:8.4f}   {:8.4f}'

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
                             f0v[i], 
                             f2v[i],
                             *f0xyz[i][:],
                             *promo[i]))

    return


def print_spinorbit_header(label):
    """ print out Spinorbit section header"""

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'Spinorbit, label = '+str(label)
        lpad = int(0.5*(max(0,LLEN-len(title))))
        pstr = str('*'.ljust(lpad)+title)
        pstr = pstr.ljust(LLEN-1)+'*'

        outfile.write('\n\n '+str('*'*LLEN))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str(pstr))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str('*'*LLEN)+'\n')
        outfile.flush()

    return

def print_spinorbit_table(hsoc, hdim, stlbl, thrsh):
    """print out the summary files for the SOC matrix elements"""

    # max bra/ket label length
    llen = max([len(lbl[0]) for lbl in stlbl])

    # table header
    delim = ' '+'-'*(53+2*llen)
    print('\n'+delim, flush=True)
    fstr = '  {:<'+str(llen)+'}' \
        + ' {:>3}' \
        + ' {:>4}' \
        + '   ' \
        + ' {:<'+str(llen)+'}' \
        + ' {:>3}' \
        + ' {:>4}' \
        + ' {:>11}' \
        + ' {:>11}' \
    
    print(fstr.format('Bra',
                      'I',
                      'M',
                      'Ket',
                      'I',
                      'M',
                      'Re <SOC>',
                      'Im <SOC>'), flush=True)
    print(delim, flush=True)

    # matrix elements
    fstr = '  {:<'+str(llen)+'}' \
        + ' {:3d}' \
        + ' {:4.1f}' \
        + '   ' \
        + ' {:<'+str(llen)+'}' \
        + ' {:3d}' \
        + ' {:4.1f}' \
        + ' {:11.4f}' \
        + ' {:11.4f}' \
        + ' {:>4}'
    
    for i in range(hdim):
        for j in range(0, i):

            soc_cm = hsoc[i,j] * constants.au2cm

            if np.abs(soc_cm) > thrsh:
                print(fstr.format(stlbl[i][0],
                                  stlbl[i][2] + 1,
                                  stlbl[i][3],
                                  stlbl[j][0],
                                  stlbl[j][2] + 1,
                                  stlbl[j][3],
                                  np.real(soc_cm),
                                  np.imag(soc_cm),
                                  'cm-1'),
                      flush=True)
    
    # footer
    print(delim, flush=True)
    
    return


def print_hsoc_eig(eig, vec, hdim, stlbl):
    """print out the summary of the eigenpairs of H_SOC"""

    # max bra/ket label length
    llen = max([len(lbl[0]) for lbl in stlbl])

    thrsh = 1e-3

    delim = ' '+'-'*(31+llen)
    
    fstr_en = '  State {:3d}:' \
        + ' {:10.4f},' \
        + ' {:10.4f} eV'

    fstr_vec = '  ({:7.4f}, {:7.4f})' \
        + ' {:<'+str(llen)+'}' \
        + ' {:3d}' \
        + ' {:4.1f}'

    
    for i in range(hdim):

        print('\n'+delim, flush=True)
        
        print(fstr_en.format(i+1,
                             eig[i],
                             (eig[i] - eig[0]) * constants.au2ev),
              flush=True)

        print(delim, flush=True)

        indx = np.flip(np.argsort(np.abs(vec[:, i])))

        for j in range(hdim):
            if np.abs(vec[indx[j], i]) > thrsh:
                print(fstr_vec.format(np.real(vec[indx[j], i]),
                                      np.imag(vec[indx[j], i]),
                                      stlbl[indx[j]][0],
                                      stlbl[indx[j]][2] + 1,
                                      stlbl[indx[j]][3],),
                      flush=True)
                
        print(delim, flush=True)
                
    return


def print_overlap_header(label):
    """print out Overlap section header"""

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'Overlap, label = '+str(label)
        lpad = int(0.5*(max(0,LLEN-len(title))))
        pstr = str('*'.ljust(lpad)+title)
        pstr = pstr.ljust(LLEN-1)+'*'

        outfile.write('\n\n '+str('*'*LLEN))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str(pstr))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str('*'*LLEN)+'\n')
        outfile.flush()

    return


def print_overlaps(trans_list, overlaps, bra_label, ket_label,
                   irreplbl, bra_state_sym, ket_state_sym):
    """
    Prints the table of wave function overlaps
    """

    # max bra/ket label length
    llen = max(len(bra_label), len(ket_label))
    
    # table header
    delim = ' '+'-'*(36)
    print('\n'+delim, flush=True)

    fstr = '  {:'+str(10+llen)+'}'
    print(fstr.format('Bra Label: '+bra_label), flush=True)
    print(fstr.format('Ket Label: '+ket_label), flush=True)
    
    print(delim, flush=True)
    
    fstr = '  {:<12} {:<12} {:<12}'
    
    print(fstr.format('Bra State',
                      'Ket State',
                      'Overlap'),
          flush=True)
    
    print(delim, flush=True)

    # Overlaps
    fstr = '{:4d} {:<7} {:4d} {:<8} {:9.6f}'
    for indx in range(len(trans_list)):
        
        sij = overlaps[indx]
        
        bk_st = trans_list[indx]
            
        #[birr, bst]    = bra_obj.state_sym(bk_st[0])
        #[kirr, kst]    = ket_obj.state_sym(bk_st[1])
        
        [birr, bst] = bra_state_sym[bk_st[0]]
        [kirr, kst] = ket_state_sym[bk_st[1]]
                
        if birr != kirr:
            continue

        if np.abs(sij) < 1e-6:
            continue
        
        irr    = birr
        irrlbl = irreplbl[irr]
        
        print(fstr.format(bk_st[0]+1, '('+irrlbl+')',
                          bk_st[1]+1, '('+irrlbl+')',
                          sij),
              flush=True)

    # table footer
    print(delim, flush=True)
        
    return

def print_dyson_header(label):
    """print out Dyson section header"""

    LLEN = 76

    with output_file(file_names['out_file'], 'a+') as outfile:
        title = 'Dyson, label = '+str(label)
        lpad = int(0.5*(max(0,LLEN-len(title))))
        pstr = str('*'.ljust(lpad)+title)
        pstr = pstr.ljust(LLEN-1)+'*'

        outfile.write('\n\n '+str('*'*LLEN))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str(pstr))
        outfile.write(  '\n '+str('*'.ljust(LLEN-1))+'*')
        outfile.write(  '\n '+str('*'*LLEN)+'\n')
        outfile.flush()

    return

#
def print_dyson_table(init_st, init_sym, final_st, final_sym, 
                      exc_ener, sqnorm): 
    """
    prints out the summary of a Dyson orbital calculation for
    a single initial state
    """

    with output_file(file_names['out_file'], 'a+') as outfile:

        delim = ' '+50*'-'

        fstr = '\n Ionisation probabilities, initial state = {:3d}({:>3})'
        outfile.write('\n\n'+delim)
        outfile.write(fstr.format(init_st, init_sym))
        outfile.write('\n'+delim)

        header  = '\n  Initial     Final    Exc Ener    <psi_D|psi_D>'
        header += '\n  State       State      (eV)'

        fstr   = '{:3d}({:>3}) -> {:3d}({:>3}) {:7.2f}'+ \
                '    {:9.4f} \n'

        outfile.write(header)
        outfile.write('\n '+delim)
        outfile.write('\n')

        for i in range(len(final_st)):
            outfile.write(fstr.format(init_st,
                                      init_sym,
                                      final_st[i],
                                      final_sym[i],
                                      exc_ener[i]*constants.au2ev, 
                                      sqnorm[i]))

        outfile.write(delim)
            
    return
    
#
def print_param_header(target_data, ci_objs, ref_states, hparams):
    """
    print header for reparameterization run
    """

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n Parameterization Optimization Run\n')
        outfile.write(' -----------------------------------\n')

        try:
            n_omp = os.environ['OMP_NUM_THREADS']
        except KeyError:
            n_omp = 1
        outfile.write('\n Maximum number of parallel processes:  '+ 
                                 str(params.nproc))
        outfile.write('\n Number of threads requested:           '+ 
                                 str(n_omp))
        outfile.write('\n N_proc x N_thread:                     '+
                                 str(int(params.nproc)*int(n_omp)))

        outfile.write('\n\n Reference Data -------------\n\n')
        for molecule, states in target_data.items():
            outfile.write(' '+str(molecule)+': '+str(states)+'\n')

        outfile.write('\n Found Reference States -----\n\n')
        for molecule in ci_objs.keys():
            outfile.write(' '+str(molecule) + ': ' + 
                          str(ci_objs[molecule]) + ': ' + 
                          str(ref_states[molecule])+'\n')

        outfile.write('\n Initial Parameter Values -----\n')
        outfile.write(' '+str(hparams)+'\n\n')

        outfile.flush()

    return

#
def print_param_iter(cur_iter, params, error):
    """
    Print results of current parameterization iterations
    """

    with output_file(file_names['out_file'], 'a+') as outfile:

        nparam = len(params)
        args   = [cur_iter] + list(params) + [error]
        fstr   =  ' iteration {:>5d} |'
        fstr   += ' parameters: ' + ' '.join(['{:10.8f}']*nparam)
        fstr   += ' |error| = {:10.8f}\n'

        outfile.write(fstr.format(*args))
        outfile.flush()

    return

#
def print_param_results(p_final, res, target, init_ener, final_ener):
    """
    print result of a parameterization run
    """

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n\n Results\n')
        outfile.write(' -----------------------------------------')
        outs = [str(res['message']), res['fun'], res['nfev']]

        outfile.write('\n Status:            {:>50s}'.format(outs[0]))
        outfile.write('\n Norm of Error:     {:50.8f}'.format(outs[1]))
        outfile.write('\n # of Evaluations:  {:50d}'.format(outs[2]))

        outfile.write('\n\n Final Parameter Values')
        outfile.write('\n -----------------------------------------')
        fstr = '\n '+' '.join(['{:10.8f} ']*len(p_final))
        outfile.write(fstr.format(*list(p_final)))

        outfile.write('\n\n Reference Data')
        outfile.write('\n -----------------------------------------')

        tstr = '\n {:<20s} {:>10s} {:>10s}'+' '.join(['{:>10s}']*4)+'\n'
        fstr = '\n {:<20s} {:>10s} {:>10s}'+' '.join(['{:10.5f}']*4)
        estr = '\n {:<52s} {:10.5f} {:10.5f}'

        outfile.write(tstr.format('Molecule', 'istate', 'fstate', 
                                  'Reference', 'Initial', 'Final', 
                                  '|Change|'))

        au_ev = constants.au2ev

        rmsd = np.asarray([0.,0.], dtype=float)
        mae  = np.asarray([0.,0.], dtype=float)
        n    = 0
        for molecule, states in target.items():
            for trans, ener in target[molecule].items():
                init, final  = trans.strip().split()
                exc_i     = (init_ener[molecule][final] - 
                               init_ener[molecule][init]) * au_ev
                exc_f    = (final_ener[molecule][final] - 
                               final_ener[molecule][init]) * au_ev
                err_i    = abs(exc_i - ener)
                err_f    = abs(exc_f - ener)

                mae     += np.asarray([err_i, err_f], dtype=float)
                rmsd    += np.asarray([err_i**2, err_f**2], dtype=float)
                n       += 1
                outfile.write(fstr.format(molecule, init, final, ener, 
                                        exc_i, exc_f, err_f - err_i ))

        mae  = mae / n
        rmsd = np.sqrt(rmsd/n)

        outfile.write('\n'+'-'*(86))
        outfile.write(estr.format('RMSD', rmsd[0], rmsd[1]))
        outfile.write(estr.format('MAE', mae[0], mae[1]))
        outfile.write('\n')

        outfile.flush()

    return


#
def print_bdd_header():
    """print the block diagonalisation diabatisation header"""
    global file_names

    with output_file(file_names['out_file'], 'a+') as outfile:
        outfile.write('\n Block Diagonalisation Diabatisation\n')
        outfile.write(' -----------------------------------')
        outfile.flush()

    return

def print_diabpot(diabpot, nroots, nirr, irrlbl):
    """
    Prints the table of diabatic potential matrix elements for each irrep
    """

    delim = ' '+'-'*(36)

    fstr = '{:4d} {:4d}     {:10.6f}'

    # table header
    print('\n'+delim, flush=True)
    print('  Diabatic potential matrix elements', flush=True)
    
    # loop over irreps
    for irr in range(nirr):
    
        # sub-table header
        print(delim, flush=True)
        print('  '+irrlbl[irr]+' block', flush=True)
        print(delim, flush=True)

        # matrix elements for this irrep
        for i in range(nroots[irr]):
            for j in range(i,nroots[irr]):
                print(fstr.format(i+1, j+1, diabpot[irr][i,j]),
                      flush=True)
                
    # table footer
    print(delim, flush=True)

    return

def print_coords(crds, asym):
    """prints a set Cartesian coordinates and atom labels"""

    fstr = ' {:<}'+3*'  {:10.7f}'

    print('\n\n Cartesian Coordinates', flush=True)
    print(' ---------------------', flush=True)
    for i in range(crds.shape[0]):
        print(fstr.format(asym[i],
                          crds[i,0],
                          crds[i,1],
                          crds[i,2]),
              flush=True)
    
    return
