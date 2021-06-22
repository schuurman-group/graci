"""
A set of standard structures to store basis set
of MO information
"""
import os
import numpy as np
import graci.io.moinfo as moinfo

ao_ordr     = [['s'],
               ['px','py','pz'],
               ['dxx','dyy','dzz','dxy','dxz','dyz'],
               ['fxxx','fyyy','fzzz','fxxy','fxxz',
                'fxyy','fyyz','fxzz','fyzz','fxyz']]

ao_norm     = [[1.],
               [1.,1.,1.],
               [1.,1.,1.,np.sqrt(3.),np.sqrt(3.),np.sqrt(3.)],
               [1.,1.,1.,np.sqrt(5.),np.sqrt(5.),np.sqrt(5.),
                         np.sqrt(5.),np.sqrt(5.),np.sqrt(5.),np.sqrt(15.)]]

def write_orbitals(file_name, mol, orbs, occ):
    """print the orbitals to gamess dat file format. Code assumes
       input orbs are in pyscf format"""

    # construt geometry
    gam_geom = moinfo.Geom()
    atms     = mol.atoms()
    carts    = mol.cart()
    for iatm in range(mol.n_atoms()):
        gam_atom = moinfo.Atom(atms[iatm], carts[iatm,:])
        gam_geom.add_atom(gam_atom)

    # pyscf mol object contains all basis set info
    pymol     = mol.pymol()
    gam_basis = moinfo.BasisSet(mol.basis, gam_geom)

    # load the basis function info into our standard
    # basis set objects
    for atm_id in range(len(atms)):
        atm_bas_ids = pymol.atom_shell_ids(atm_id)
        for bas_id in range(len(atm_bas_ids)):
            ang    = pymol.bas_angular(atm_bas_ids[bas_id])
            n_cont = pymol.bas_nctr(atm_bas_ids[bas_id])
            n_prim = pymol.bas_nprim(atm_bas_ids[bas_id])
            b_funcs = [moinfo.BasisFunction(ang) for m in range(n_cont)]
            expons  = pymol.bas_exp(atm_bas_ids[bas_id])
            coefs   = pymol.bas_ctr_coeff(atm_bas_ids[bas_id])

            for iprim in range(n_prim):
                expon = expons[iprim]
                for icont in range(n_cont):
                    b_funcs[icont].add_primitive(expon,
                                                 coefs[iprim, icont])
            
            for icont in range(n_cont):
                gam_basis.add_function(atm_id, b_funcs[icont])


    # put the orbitals into a GAMESS orbital object
    # we will necessarily convert orbitals to cartesians 
    # if need be
    print("pymol.cart: "+str(pymol.cart))
    if pymol.cart:
        orb_cart   = orbs
    else:
        # this generates (nmo, nao) transformation matrix
        sph2cart = np.linalg.pinv(pymol.cart2sph_coeff())
        orb_cart = np.matmul(sph2cart.T, orbs)
    (nao, nmo) = orb_cart.shape

    gam_mos = moinfo.Orbitals(nao, nmo)
    gam_mos.mo_vectors = orb_cart
    gam_mos.occ = occ

    # construct mapping array from pyscf to GAMESS
    pygam_map, scale_py, scale_gam = gam_basis.construct_map(ao_ordr, ao_norm)

    # remove the pyscf normalization factors
    gam_mos.scale(scale_py)

    # re-sort orbitals to GAMESS ordering
    gam_mos.sort_aos(pygam_map)

    # apply the GAMESS normalization factors
    gam_mos.scale(scale_gam)

    # print the MOs file
    with open(file_name, 'w') as orb_file:
        print_basis_set(orb_file, gam_basis)
        print_orbitals(orb_file, gam_mos)
        if gam_mos.occ is not None:
            print_occ(orb_file, gam_mos.occ)

    return


def print_atom(file_handle, atom):
    """Prints the atom in GAMESS file format."""
    ofmt   = ('{:2s}'+'{:>5s}'+
              ''.join('{:18.10f}' for i in range(3))+'\n')
    w_data = [atom.symbol, atom.number] + atom.coords.tolist()
    file_handle.write(ofmt.format(*w_data))

    return

def print_geom(file_handle, geom):
    """Prints geometry in GAMESS format."""
    for i in range(geom.natoms()):
        geoms.atoms[i].print_atom(file_handle)

    return 

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

