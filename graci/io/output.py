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
        else:
            outfile.write('\n\n')

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
    print('\n'+delim)
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
                      'Im <SOC>'))
    print(delim)

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
                                  'cm-1'))
    
    # footer
    print(delim)
    
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

        print('\n'+delim)
        
        print(fstr_en.format(
            i+1,
            eig[i],
            (eig[i] - eig[0]) * constants.au2ev))

        print(delim)

        indx = np.flip(np.argsort(np.abs(vec[:, i])))

        for j in range(hdim):
            if np.abs(vec[indx[j], i]) > thrsh:
                print(fstr_vec.format(
                    np.real(vec[indx[j], i]),
                    np.imag(vec[indx[j], i]),
                    stlbl[indx[j]][0],
                    stlbl[indx[j]][2] + 1,
                    stlbl[indx[j]][3],))

        print(delim)
                
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
    print('\n'+delim)

    fstr = '  {:'+str(10+llen)+'}'
    print(fstr.format('Bra Label: '+bra_label))
    print(fstr.format('Ket Label: '+ket_label))
    
    print(delim)
    
    fstr = '  {:<12} {:<12} {:<12}'
    
    print(fstr.format('Bra State',
                      'Ket State',
                      'Overlap'))
    
    print(delim)

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
                          sij))

    # table footer
    print(delim)
        
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
    print('\n'+delim)
    print('  Diabatic potential matrix elements')
    
    # loop over irreps
    for irr in range(nirr):
    
        # sub-table header
        print(delim)
        print('  '+irrlbl[irr]+' block')
        print(delim)

        # matrix elements for this irrep
        for i in range(nroots[irr]):
            for j in range(i,nroots[irr]):
                print(fstr.format(i+1, j+1, diabpot[irr][i,j]))
                
    # table footer
    print(delim)

    return
