"""
A set of standard structures to store basis set
of MO information
"""
import os
import numpy as np
import graci.core.molecule as molecule

a_symbols   = [' H','HE',
               'LI','BE',' B',' C',' N',' O',' F','NE',
               'NA','MG','AL','SI',' P',' S','CL','AR']

a_numbers   = [ '1.0', '2.0',
                '3.0', '4.0', '5.0', '6.0', '7.0', '8.0', '9.0','10.0',
               '11.0','12.0','13.0','14.0','15.0','16.0','17.0','18.0']

ang_mom_sym = ['S','P','D','F','G','H','I']

nfunc_cart  = [1, 3, 6, 10]

nfunc_sph   = [1, 3, 5, 7]

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

# how to convert from spherical to cartesian basis functions
# (in molden ordering)
# s -> s | px -> px | py -> py | pz -> pz
# d ordering: d0, d1+, d1-, d2+, d2-
# dxx ->  -d0 / 2 + sqrt(3.) * d2+ / 2
# dyy -> -d0 / 2 - sqrt(3.) * d2+ / 2
# dzz ->  d0
# dxy ->  d2-
# dxz ->  d1+
# dyz ->  d1-
# f ordering: f0, f1+, f1-, f2+, f2-, f3+, f3-
# fxxx -> -f1+
# fyyy -> -f1-
# fzzz -> -f0
# fxyy -> -(3/(2*sqrt(5)))*f1+ + (sqrt(3)/2)*f3+
# fxxy -> -(3/(2*sqrt(5)))*f1- + (sqrt(3)/2)*f3-
# fxxz -> -(3/(2*sqrt(5)))*f0  + (sqrt(3)/2)*f2+
# fxzz -> -(3/(2*sqrt(5)))*f1+ - (sqrt(3)/2)*f3+
# fyzz -> -(3/(2*sqrt(5)))*f1- - (sqrt(3)/2)*f3-
# fyyz -> -(3/(2*sqrt(5)))*f0  - (sqrt(3)/2)*f2+
# fxyz -> f2-
a        = 1./2.
b        = np.sqrt(3.)
c        = np.sqrt(5.)
sph2cart = [
    [[[0], [1.]]],                           # conversion for s orbitals
    [[[0], [1.]], [[1], [1.]], [[2], [1.]]], # conversion for p orbitals
    [[[0, 3], [-a,  a*b]],                   # conversion for d orbitals
     [[0, 3], [-a, -a*b]],
     [[0], [1.]],
     [[4], [1.]],
     [[1], [1.]],
     [[2], [1.]]],
    [[[1, 5], [1., 0.]],                     # conversion for f orbitals
     [[2, 6], [1., 0.]],
     [[0], [1.]],
     [[1, 5], [-3.*a/c , a*b ]],
     [[2, 6], [-3.*a/c,  a*b ]],
     [[0, 3], [-3.*a/c,  a*b ]],
     [[1, 5], [-3.*a/c, -a*b ]],
     [[2, 6], [-3.*a/c, -a*b ]],
     [[0, 3], [-3.*a/c, -a*b ]],
     [[4], [1.]]]
            ]


def write_orbitals(file_name, mol, orbs, occ):
    """print the orbitals to gamess dat file format"""

    # construt geometry
    gam_geom = self.Geom()
    atms     = mol.atoms()
    cart     = mol.cart()
    for iatm in range(mol.n_atoms()):
        atom = self.Atom(atms[iatm], cart[iatm,:])
        geom.add_atom(atom)

    # pyscf mol object contains all basis set info
    pymol     = mol.pymol()
    gam_basis = self.BasisSet(mol.basis, geom)

    atm_cnt = 0
    for atom,bset in pymol._basis.items():
        # identify correct atom in geom
        atom.strip().lower() != atms[atm_cnt].strip().lower()
        print("atom error: "+atom.strip().lower()+"!=" + 
                              atms[icnt].strip().lower())
        for iset in range(len(bset)):
            ang     = bset[iset][0]
            # first element is the angulatr momentum
            n_prim  = len(bset[iset])-1
            # first element is the exponent
            n_cont  = len(bset[iset][1])-1
            b_funcs = [self.BasisFunction(ang) for m in range(n_cont)]
            for iprim in range(n_prim):
                expon = bset[iset][iprim+1][0]
                for icont in range(n_cont):
                    cf = bset[iset][iprim+1][icont+1]
                    b_funcs[icont].add_primitive(expon, cf)

            for icont in range(n_cont):
                gam_basis.add_function(atm_cnt, b_funcs[icont])

        atm_cnt += 1

    # put the orbitals into a GAMESS orbital object
    # we will necessarily convert orbitals to cartesians 
    # if need be
    if pymol.cart:
        orb_cart   = orbs
    else:
        orb_cart = sph2cart(gam_basis, orbs, s2c)
    (nao, nmo) = orb_cart.shape

    gam_orb = moinfo.Orbitals(nao, nmo)
    gam_orb.mo_vectors = orb_cart
    gam_orb.occ = occ

    # construct mapping array from COLUMBUS to GAMESS
    pygam_map, scale_py, scale_gam = basis.construct_map(ao_ordr, ao_norm)

    # remove the COLUMBUS normalization factors
    gam_orb.scale(scale_py)

    # re-sort orbitals to GAMESS ordering
    gam_orb.sort_aos(pygam_map)

    # apply the GAMESS normalization factors
    gam_orb.scale(scale_gam)

    # print the MOs file
    gam_basis.print_basis_set(file_name)
    gam_mos.print_orbitals(file_name)
    if gam_mos.occ is not None:
        gam_mos.print_occ(file_name)

    return


def sph2cart(basis, mo_sph, s2c):
    """Converts orbitals from spherical to cartesian AO basis. Assumes
    input orbitals are numpy array."""
    nao_sph, nmo = mo_sph.shape
    l_cart, l_sph = basis.ang_mom_lst()

    orb_trans = [[] for i in range(nmo)]
    iao_sph   = 0
    while(iao_sph < nao_sph):
        lval = l_sph[iao_sph]
        for imo in range(nmo):
            cart_orb = [sum([mo_sph[iao_sph+
                        s2c[lval][j][0][k],imo]*s2c[lval][j][1][k]
                        for k in range(len(s2c[lval][j][0]))])
                        for j in range(nfunc_cart[lval])]
            orb_trans[imo].extend(cart_orb)
        iao_sph += nfunc_sph[lval]

    if iao_sph != basis.n_bf_sph:
        raise ValueError('Error in sph2cart: {:d} '.format(iao_sph) +
                         '!= {:d}'.format(basis.n_bf_sph))

    return np.array(orb_trans).T


class Atom:
    """Wrapper for all requisite data about an atom."""
    def __init__(self, label, coords=None):
        # atomic symbol
        self.symbol    = label

        # set atomic number based on symbol
        try:
            self.number = a_numbers[a_symbols.index(self.symbol)]
        except:
            print('Atom: '+str(label)+' not found,'+
                  'setting atomic number to \'0\'\n')
            self.number = '0.0'

        # coordinates
        if coords is None:
            self.coords = [0. for i in range(3)]
        else:
            self.coords = coords

    def set_coords(self, coords):
        """Defines the coordinates of the atom."""
        self.coords = coords

    def print_atom(self, file_handle):
        """Prints the atom in GAMESS file format."""
        ofmt   = ('{:2s}'+'{:>5s}'+
                  ''.join('{:18.10f}' for i in range(3))+'\n')
        w_data = [self.symbol, self.number] + self.coords
        file_handle.write(ofmt.format(*w_data))


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
        scoords = [fmt.format(*atm.coords) for atm in self.atoms]
        ocoords = [fmt.format(*atm.coords) for atm in other_geom.atoms]
        self.atom_map = self._find_map(ocoords, scoords)
        self.atoms = self.atoms[self.atom_map]

    def print_geom(self, file_handle):
        """Prints geometry in GAMESS format."""
        for i in range(self.natoms()):
            self.atoms[i].print_atom(file_handle)

    def _find_map(self, arr1, arr2):
        """Finds a map between two arrays."""
        o1 = np.argsort(arr1)
        o2 = np.argsort(arr2)
        q = np.empty_like(o1)
        q[o1] = np.arange(len(o1))
        return o2[q]


class Orbitals:
    """Wrapper to hold information about the orbitals."""
    def __init__(self, n_aos, n_mos):
        # number of AOs
        self.naos       = n_aos
        # number of MOs
        self.nmos       = n_mos
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

    def print_orbitals(self, file_name, n_orb=None):
        """Prints orbitals in gamess style VEC format."""
        # default is to print all the orbitals
        if n_orb is None:
            n_orb = self.nmos

        # open file_name, append if file already exists
        with open(file_name, 'a') as mo_file:
            mo_file.write(' $VEC\n')
            for i in range(n_orb):
                self.print_movec(mo_file, i)
            mo_file.write(' $END\n')

    def print_occ(self, file_name, n_orb=None):
        """Prints orbital occupations."""
        if n_orb is None:
            n_orb = self.nmos

        n_col = 5
        n_row = int(np.ceil(n_orb / n_col))
        with open(file_name, 'a') as mo_file:
            mo_file.write(' $OCCNO\n')
            for i in range(n_row):
                n_end = min(n_orb, 5*(i+1))
                oc_row = (n_end - 5*i)*'{:16.10f}' + '\n'
                r_data = self.occ[5*i:n_end]
                mo_file.write(oc_row.format(*r_data))
            mo_file.write(' $END\n')

    def print_movec(self, file_handle, mo_i):
        """Prints an orbital vector."""
        n_col  = 5 # this is set by GAMESS format
        n_row  = int(np.ceil(self.naos / n_col))
        mo_lab = (mo_i+1) % 100

        for i in range(n_row):
            r_data = [mo_lab, i+1]
            n_end = min(self.naos,5*(i+1))
            mo_row = '{:>2d}{:>3d}' + (n_end - 5*i)*'{:15.8E}' + '\n'
            r_data += self.mo_vectors[5*i:n_end,mo_i].tolist()
            file_handle.write(mo_row.format(*r_data))


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

    def print_basis_function(self, file_handle):
        """Prints the basis function in gamess format."""
        ofmt1 = ('{:1s}'+'{:>6d}'+'\n')
        ofmt2 = ('{:>3d}'+'{:>15.7f}'+'{:>23.7f}'+'\n')

        w_data = [self.ang_sym, self.n_prim]
        file_handle.write(ofmt1.format(*w_data))
        for i in range(self.n_prim):
            w_data = [i+1, self.exps[i],self.coefs[i]]
            file_handle.write(ofmt2.format(*w_data))


class BasisSet:
    """Contains the basis set information for a single atom."""
    def __init__(self, label, geom):
        # string label of basis set name
        self.label       = label
        # nuclear coordinates (and coordinates of centers of gaussians
        self.geom        = geom
        # total number of cartesian functions
        self.n_bf_cart   = 0
        # total number of spherical functions
        self.n_bf_sph    = 0
        # total number of contractions for atom i
        self.n_bf        = [0 for i in range(self.geom.natoms())]
        # list of basis function objects
        self.basis_funcs = [[] for i in range(self.geom.natoms())]

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

    def construct_map(self, orig_ordr, orig_norm):
        """Constructs a map array to convert the nascent AO ordering to
        GAMESS AO ordering. Also returns the normalization factors
        for the nascent and corresponding GAMESS-ordered AO basis."""
        map_arr       = []
        scale_nascent = []
        scale_gam     = []
        ao_map        = [[orig_ordr[i].index(ao_ordr[i][j])
                         for j in range(nfunc_cart[i])]
                         for i in range(len(ao_ordr))]

        nf_cnt = 0
        for i in self.geom.atom_map:
            for j in range(len(self.basis_funcs[i])):
                ang_mom = self.basis_funcs[i][j].ang_mom
                nfunc   = nfunc_cart[ang_mom]
                map_bf = [nf_cnt + ao_map[ang_mom][k] for k in range(nfunc)]
                map_arr.extend(map_bf)
                scale_nascent.extend([orig_norm[ang_mom][k] for k in range(nfunc)])
                scale_gam.extend([1./ao_norm[ang_mom][k] for k in range(nfunc)])
                nf_cnt += nfunc

        return map_arr, scale_nascent, scale_gam

    def print_basis_set(self, file_name):
        """Print the data section of a GAMESS style input file."""

        # if file doesn't exist, create it. Overwrite if it does exist
        with open(file_name, 'x') as dat_file:
            dat_file.write(' $DATA\n'+
                           'Comment line | basis='+str(self.label)+'\n'+
                           'C1'+'\n')

            for i in range(self.geom.natoms()):
                self.geom.atoms[i].print_atom(dat_file)
                for j in range(self.n_bf[i]):
                    self.basis_funcs[i][j].print_basis_function(dat_file)
                # each atomic record ends with an empty line
                dat_file.write('\n')

            dat_file.write(' $END\n')
