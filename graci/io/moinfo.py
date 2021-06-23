"""
A set of standard structures to store basis set
of MO information
"""

import sys as sys
import numpy as np
from pyscf.gto import mole

a_symbols   = ['H','HE',
               'LI','BE','B','C','N','O','F','NE',
               'NA','MG','AL','SI','P','S','CL','AR']

a_numbers   = [ '1.0', '2.0',
                '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0','10.0',
               '11.0','12.0','13.0','14.0','15.0','16.0','17.0','18.0']

ang_mom_sym = ['S','P','D','F','G','H','I']

nfunc_cart  = [1, 3, 6, 10]

nfunc_sph   = [1, 3, 5, 7]

# this is the standard ao ordering scheme used by libcint/pyscf
cart_ao_ordr = [['s'],
                ['px','py','pz'],
                ['dxx','dxy','dxz','dyy','dyz','dzz'],
                ['fxxx','fxxy','fxxz','fxyy','fxyz',
                 'fxzz','fyyy','fyyz','fyzz','fzzz']]

sph_ao_ordr  =  [['s'],
                 ['px','py','pz'],
                 ['dxy','dyz','dz2','dxz','dx2-y2'],
                 ['f3x2y-y3','fxyz','fyz2','fz3','fxz2',
                                   'fx2z-y2z','fx3-xy2']]

cart_ao_norm = [[1.],
                [1.,1.,1.],
                [1.,1.,1.,1.,1.,1],
                [1.,1.,1.,1.,1.,1.,1.,1.,1.,1.]]

sph_ao_norm = [[1.],
               [1.,1.,1.],
               [1.,1.,1.,1.,1.],
               [1.,1.,1.,1.,1.,1.,1]]

# generate geom, basis objects from GRaCI molecule object
def create_geom(mol):
    """create a Geometry object from a GRaCI molecule object"""

    new_geom = Geom()
    atms     = mol.atoms()
    carts    = mol.cart()
    for iatm in range(mol.n_atoms()):
        new_atom = Atom(atms[iatm], carts[iatm,:])
        new_geom.add_atom(new_atom)

    return new_geom

# generate a basis function object from GRaCI molecule object
def create_basis(mol, geom):
    """create a BasisSet object from a GRaCI molecule object"""

    new_basis = BasisSet(mol, geom)

    return new_basis

# create orbital object
def create_orbitals(mol, orbs, occ, cart):
    """create an Orbitals object from input orbitals"""

    pymol = mol.pymol()

    # if AOs already cartesian and we want cartesian orbitals,
    # do nothing
    if cart and pymol.cart:
        orb_out   = orbs
    elif cart and pymol.cart is False:
        # this generates (nmo, nao) transformation matrix
        s2c = np.linalg.pinv(pymol.cart2sph_coeff(normalized='sp'))
        orb_out = np.matmul(s2c.T, orbs)
    elif cart is False and pymol.cart is False:
        orb_out = orbs
    elif cart is False and pymol.cart:
        c2s = pymol.cart2sph_coeff(normalized='sp')
        orb_out = np.matmul(c2s.T, orbs)

    (nao, nmo) = orb_out.shape

    new_mos            = Orbitals(nao, nmo, cart_aos=cart)
    new_mos.mo_vectors = orb_out
    new_mos.occ        = occ

    return new_mos


class Atom:
    """Wrapper for all requisite data about an atom."""
    def __init__(self, label, coords=None):
        # atomic symbol
        self.symbol    = label

        # set atomic number based on symbol
        try:
            self.number = a_numbers[
                    a_symbols.index(self.symbol.strip().upper())]
        except:
            print('Atom: '+str(label)+' not found,'+
                  'setting atomic number to \'0\'\n')
            self.number = '0.0'

        # coordinates
        if coords is None:
            self.coords = np.zeros(3, dtype=float)
        else:
            self.coords = np.array(coords, dtype=float)

    def set_coords(self, coords):
        """Defines the coordinates of the atom."""
        self.coords = np.array(coords, dtype=float)

class Geom:
    """Class to hold geometry information."""
    def __init__(self, atoms=None):
        # if atoms is None, initialize empty atom list
        if atoms is None:
            self.atoms = np.array([], dtype=object)
        else:
            self.atoms = np.array(atoms)
        self.atom_map = np.arange(self.natoms())

    def add_atom(self, atom):
        """Adds an atom to the geometry."""
        nat = self.natoms()
        self.atoms = np.hstack((self.atoms, atom))
        self.atom_map = np.hstack((self.atom_map, nat))

    def natoms(self):
        """Gets the number of atoms."""
        return len(self.atoms)

    def reorder(self, other_geom):
        """Reorders a geometry to match another geometry and keeps
        the mapping to use on mos."""
        fmt = '{:10.4e}{:10.4e}{:10.4e}'
        scoords = [fmt.format(*atm.coords.tolist()) for atm in self.atoms]
        ocoords = [fmt.format(*atm.coords.tolist()) for atm in other_geom.atoms]
        self.atom_map = self._find_map(ocoords, scoords)
        self.atoms = self.atoms[self.atom_map]

    def _find_map(self, arr1, arr2):
        """Finds a map between two arrays."""
        o1 = np.argsort(arr1)
        o2 = np.argsort(arr2)
        q = np.empty_like(o1)
        q[o1] = np.arange(len(o1))
        return o2[q]

class Orbitals:
    """Wrapper to hold information about the orbitals."""
    def __init__(self, n_aos, n_mos, cart_aos=True):
        # number of AOs
        self.naos       = n_aos
        # number of MOs
        self.nmos       = n_mos
        # cartesian AOs or spherical harmonic AOs
        self.cart_aos   = cart_aos

        # matrix holding MOs
        self.mo_vectors = np.zeros((n_aos, n_mos), dtype=float)
        # vector holding occupations
        self.occ        = None

    def add(self, mo_vec):
        """Adds an orbital to the end of the list (i.e. at first column
        with norm==0)"""
        # only add orbital if the number of aos is correct
        if len(mo_vec) != self.naos:
            print("Cannot add orbital, wrong number of naos: "+
                   str(len(mo_vec))+"!="+str(self.naos))
            return None

        i = 0
        while(np.linalg.norm(self.mo_vectors[:,i]) > 1.e-10):
            i += 1
        self.mo_vectors[:,i] = mo_vec

    def insert(self, mo_vec, mo_i):
        """Inserts an orbital."""
        self.mo_vectors = np.insert(self.mo_vectors,mo_i,mo_vec,axis=1)

    def delete(self, mo_i):
        """Deletes an orbital."""
        self.mo_vectors = np.delete(self.mo_vectors,mo_i,axis=1)

    def scale(self, fac_vec):
        """Scales each MO by the vector fac_vec."""
        scale_fac       = np.array(fac_vec)
        old_mos         = self.mo_vectors.T
        new_mos         = np.array([mo * scale_fac for mo in old_mos]).T
        self.mo_vectors = new_mos

    def sort_aos(self, map_lst):
        """Re-sorts the MOs, ordering the AO indices via map_lst."""
        for i in range(self.nmos):
            vec_srt = self.mo_vectors[map_lst,i]
            self.mo_vectors[:,i] = vec_srt

    def sort_mos(self, map_lst):
        """Re-sorts the MO ordering"""
        vec_srt = self.mo_vectors[:,map_lst]
        self.mo_vectors = vec_srt

    def norm(self, mo_i):
        """Takes the norm of an orbital."""
        np.linalg.norm(self.mo_vectors[:,mo_i])

    def convert(self, basis_set, new_ao_ordr, new_ao_norm):
        """convert the orbitals to new order specified by ao_ordr and new
           norm specified by ao_norm. The basis set used is given by
           basis_set"""

        pyout_map, scale_py, scale_new = basis_set.construct_map(
                                                     new_ao_ordr, 
                                                     new_ao_norm, 
                                                     self.cart_aos)

        # remove the pyscf normalization factors
        self.scale(scale_py)

        # re-sort orbitals to new ordering
        self.sort_aos(pyout_map)

        # apply the new normalization factors
        self.scale(scale_new)

        return


class BasisFunction:
    """An object to hold information about a single basis function
    including angular momentum, exponents and contraction coefficients."""
    def __init__(self, ang_mom):
        # angular momentum of basis function
        self.ang_mom = ang_mom
        # text symbol of the angular momentum shell
        self.ang_sym = ang_mom_sym[self.ang_mom]
        # number of primitives in basis
        self.n_prim  = 0
        # list of exponents
        self.exps    = []
        # list of coefficients
        self.coefs   = []

    def add_primitive(self, expo, coef):
        """Adds a primitive to the basis function."""
        self.exps.extend([expo])
        self.coefs.extend([coef])
        self.n_prim += 1

class BasisSet:
    """Contains the basis set information for a single atom."""
    def __init__(self, mol, geom):
        # nuclear coordinates (and coordinates of centers of gaussians
        self.geom        = geom
        # string label of basis set name
        self.label       = mol.basis
        # pyscf mole object
        self.pymol       = mol.pymol()
        # total number of cartesian functions
        self.n_bf_cart   = 0
        # total number of spherical functions
        self.n_bf_sph    = 0
        # total number of contractions for atom i
        self.n_bf        = [0 for i in range(self.geom.natoms())]
        # list of basis function objects
        self.basis_funcs = [[] for i in range(self.geom.natoms())]

        # load the basis function info into our standard
        # basis set objects
        atms     = mol.atoms()

        for atm_id in range(len(atms)):
            atm_bas_ids = self.pymol.atom_shell_ids(atm_id)
            for bas_id in range(len(atm_bas_ids)):
                ang    = self.pymol.bas_angular(atm_bas_ids[bas_id])
                n_cont = self.pymol.bas_nctr(atm_bas_ids[bas_id])
                n_prim = self.pymol.bas_nprim(atm_bas_ids[bas_id])
                b_funcs = [BasisFunction(ang) for m in range(n_cont)]
                expons  = self.pymol.bas_exp(atm_bas_ids[bas_id])
                coefs   = self.pymol.bas_ctr_coeff(atm_bas_ids[bas_id])

                for iprim in range(n_prim):
                    expon = expons[iprim]
                    print("l="+str(ang)+" norm="+str(mole.gto_norm(ang,expon)))
                    for icont in range(n_cont):
                        b_funcs[icont].add_primitive(expon,
                                                 coefs[iprim, icont])

                for icont in range(n_cont):
                    self.add_function(atm_id, b_funcs[icont])

    def ang_mom_lst(self):
        """Returns an array containing the value of angular momentum of the
        corresponding at that index."""
        ang_mom_cart = []
        ang_mom_sph  = []

        for i in range(self.geom.natoms()):
            for j in range(self.n_bf[i]):
                ang_mom = self.basis_funcs[i][j].ang_mom
                ang_mom_cart.extend([ang_mom for k in range(nfunc_cart[ang_mom])])
                ang_mom_sph.extend([ang_mom for k in range(nfunc_sph[ang_mom])])

        return ang_mom_cart, ang_mom_sph

    def add_function(self, atom_i, bf):
        """Adds a basis function to the basis_set -- always keeps "like"
        angular momentum functions together."""
        ang_mom = bf.ang_mom
        if len(self.basis_funcs[atom_i])>0:
            bf_i = 0
            while(ang_mom >= self.basis_funcs[atom_i][bf_i].ang_mom):
                bf_i += 1
                if bf_i == len(self.basis_funcs[atom_i]):
                    break
        else:
            bf_i = len(self.basis_funcs[atom_i])

        self.basis_funcs[atom_i].insert(bf_i, bf)
        self.n_bf[atom_i] += 1
        self.n_bf_cart    += nfunc_cart[bf.ang_mom]
        self.n_bf_sph     += nfunc_sph[bf.ang_mom]

    def construct_map(self, new_ordr, new_norm, cart=True):
        """Constructs a map array to convert the nascent pyscf AO ordering to
        an output AO ordering specified in out_ordr. Also returns the normalization 
        factors for the nascent and corresponding output-ordered AO basis."""
        if cart:
            pyordr = cart_ao_ordr
            pynorm = cart_ao_norm
            nfunc  = nfunc_cart
        else:
            pyordr = sph_ao_ordr
            pynorm = sph_ao_ordr
            nfunc  = nfunc_sph

        map_arr       = []
        scale_pyscf   = [] 
        scale_new     = []
        ao_map        = [[pyordr[i].index(new_ordr[i][j])
                         for j in range(nfunc[i])]
                         for i in range(len(new_ordr))]

        nf_cnt = 0
        for i in self.geom.atom_map:
            for j in range(len(self.basis_funcs[i])):
                ang_mom = self.basis_funcs[i][j].ang_mom
                nf      = nfunc[ang_mom]
                map_bf  = [nf_cnt + ao_map[ang_mom][k] for k in range(nf)]
                map_arr.extend(map_bf)
                scale_pyscf.extend([pynorm[ang_mom][k] for k in range(nf)])
                scale_new.extend([1./new_norm[ang_mom][k] for k in range(nf)])
                nf_cnt += nf

        return map_arr, scale_pyscf, scale_new
    
