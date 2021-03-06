"""
A set of standard structures to store basis set
of MO information
"""
import os
import numpy as np
import graci.io.moinfo as moinfo

gam_cart_ao_ordr = [['s'],
                   ['px','py','pz'],
                   ['dxx','dyy','dzz','dxy','dxz','dyz'],
                   ['fxxx','fyyy','fzzz','fxxy','fxxz',
                    'fxyy','fyyz','fxzz','fyzz','fxyz']]

gam_cart_ao_norm = [[1.],
                   [1.,1.,1.],
                   [1.,1.,1.,np.sqrt(3.),np.sqrt(3.),np.sqrt(3.)],
                   [1.,1.,1.,np.sqrt(5.),np.sqrt(5.),np.sqrt(5.),
                             np.sqrt(5.),np.sqrt(5.),np.sqrt(5.),np.sqrt(15.)]]

gam_sph_ao_ordr  =  [['s'],
                    ['px','py','pz'],
                    ['dxy','dyz','dz2','dxz','dx2-y2'],
                    ['f3x2y-y3','fxyz','fyz2','fz3','fxz2',
                                      'fx2z-y2z','fx3-xy2']]
gam_sph_ao_norm  = [[1.],
                    [1.,1.,1.],
                    [1.,1.,1.,1.,1.],
                    [1.,1.,1.,1.,1.,1.,1]]

def write_orbitals(file_name, mol, orbs,
                    occ=None, sym_lbl=None, ener=None, cart=None):
    """print the orbitals to gamess dat file format. Code assumes
       input orbs are in pyscf format. Default is to output orbitals in
       cartesian AOs"""

    gam_geom  = moinfo.create_geom(mol)
    gam_basis = moinfo.create_basis(mol, gam_geom)
    gam_mos   = moinfo.create_orbitals(mol, orbs, occ, ener, cart)

    # if 'cart' not specified, default is to use the AOs (cart or 
    # spherical) that were used in original calculation
    if cart is None:
        cart = mol.pymol().cart

    # construct the order array for the output and provide the AO
    # normalization factors
    if cart:
        out_ordr = gam_cart_ao_ordr
        out_nrm  = gam_cart_ao_norm
    else:
        out_ordr = gam_sph_ao_ordr
        out_nrm  = gam_sph_ao_norm

    # reorder to the gamess AO ordering format
    gam_mos.reorder(gam_basis, out_ordr)

    # right now, default is to _not_ scale, as superdyson accounts 
    # for normalization factor
    #scalevec = gam_basis.construct_scalevec(out_nrm, cart)
    #gam_mos.scale(scalevec)

    # print the MOs file
    with open(file_name, 'w') as orb_file:
        print_basis_set(orb_file, gam_basis)
        print_orbitals(orb_file, gam_mos)
        if gam_mos.occ is not None:
            print_occ(orb_file, gam_mos.occ)
        if gam_mos.ener is not None:
            print_ener(orb_file, gam_mos.ener)
    return


def print_atom(file_handle, atom):
    """Prints the atom in GAMESS file format."""
    ofmt   = ('{:2s}'+'{:>5s}'+
              ''.join('{:18.10f}' for i in range(3))+'\n')
    w_data = [atom.symbol, atom.number] + atom.coords.tolist()
    file_handle.write(ofmt.format(*w_data))

    return

#
def print_geom(file_handle, geom):
    """Prints geometry in GAMESS format."""
    for i in range(geom.natoms()):
        geoms.atoms[i].print_atom(file_handle)

    return 

#
def print_orbitals(file_handle, gam_mos):
    """Prints orbitals in gamess style VEC format."""
    # default is to print all the orbitals
    (n_ao, n_orb) = gam_mos.mo_vectors.shape

    # open file_name, append if file already exists
    file_handle.write(' $VEC\n')
    for i in range(n_orb):
        print_movec(file_handle, gam_mos.mo_vectors[:,i], i)
    file_handle.write(' $END\n')

    return

#
def print_movec(file_handle, movec, mo_i):
    """Prints an orbital vector."""
    n_aos  = len(movec)
    n_col  = 5 # this is set by GAMESS format
    n_row  = int(np.ceil(n_aos / n_col))
    mo_lab = (mo_i+1) % 100 

    for i in range(n_row):
        r_data = [mo_lab, i+1]
        n_end = min(n_aos,5*(i+1))
        mo_row = '{:>2d}{:>3d}' + (n_end - 5*i)*'{:15.8E}' + '\n'
        r_data += movec[5*i:n_end].tolist()
        file_handle.write(mo_row.format(*r_data))

    return

#
def print_occ(file_handle, occ):
    """Prints orbital occupations."""
    n_orb = len(occ) 

    n_col = 5
    n_row = int(np.ceil(n_orb / n_col))

    file_handle.write(' $OCCNO\n')
    for i in range(n_row):
        n_end = min(n_orb, 5*(i+1))
        oc_row = (n_end - 5*i)*'{:16.10f}' + '\n'
        r_data = occ[5*i:n_end]
        file_handle.write(oc_row.format(*r_data))
    file_handle.write(' $END\n')

    return

#
def print_ener(file_handle, ener):
    """Prints orbital energies."""
    n_orb = len(ener)

    n_col = 5
    n_row = int(np.ceil(n_orb / n_col))

    file_handle.write(' $ENER\n')
    for i in range(n_row):
        n_end = min(n_orb, 5*(i+1))
        e_row = (n_end - 5*i)*'{:16.10f}' + '\n'
        r_data = ener[5*i:n_end]
        file_handle.write(e_row.format(*r_data))
    file_handle.write(' $END\n')

    return

#
def print_basis_function(file_handle, bf):
    """Prints the basis function in gamess format."""
    ofmt1 = ('{:1s}'+'{:>6d}'+'\n')
    ofmt2 = ('{:>3d}'+'{:>15.7f}'+'{:>23.7f}'+'\n')

    w_data = [bf.ang_sym, bf.n_prim]
    file_handle.write(ofmt1.format(*w_data))
    for i in range(bf.n_prim):
        w_data = [i+1, bf.exps[i], bf.coefs[i]]
        file_handle.write(ofmt2.format(*w_data))

    return

def print_basis_set(file_handle, basis_set):
    """Print the data section of a GAMESS style input file."""

    # if file doesn't exist, create it. Overwrite if it does exist
    file_handle.write(' $DATA\n'+
                   'Comment line | basis='+str(basis_set.label)+'\n'+
                   'C1'+'\n')

    for i in range(basis_set.geom.natoms()):
        print_atom(file_handle, basis_set.geom.atoms[i])
        for j in range(basis_set.n_bf[i]):
            print_basis_function(file_handle, basis_set.basis_funcs[i][j])
        # each atomic record ends with an empty line
        file_handle.write('\n')

    file_handle.write(' $END\n')

    return

