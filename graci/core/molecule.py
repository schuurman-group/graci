"""
The Molecule object and its associated functions.
"""
import sys as sys
import numpy as np
import copy as copy
import graci.io.output as output
import graci.utils.basis as basis
import graci.utils.constants as constants
from pyscf.lib import logger
from pyscf import gto, df

atom_name = ['Ghost', 'X' ,'H' , 'He', 
             'Li', 'Be', 
             'B' , 'C' , 'N' , 'O' , 'F', 'Ne', 
             'Na', 'Mg', 
             'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 
             'K', 'Ca', 
             'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni','Cu', 'Zn', 
             'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 
             'Rb', 'Sr', 
             'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
             'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe', 
             'Cs', 'Ba', 
             'La', 'Ce','Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 
             'Ho', 'Er','Tm', 'Yb', 'Lu', 
             'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt','Au', 'Hg', 
             'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 
             'Fr', 'Ra',
             'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
             'Es', 'Fm', 'Md', 'No', 'Lr', 
             'Rf', 'Db', 'Sg', 'Bh', 'Hs','Mt', 'Ds', 'Rg', 'Cn', 'Nh', 
             'Fl', 'Mc', 'Lv', 'Ts', 'Og']

# Atomic masses taken from PySCF elements.py
atom_mass = [ 0., 0., 1.008, 4.002602, 6.94, 9.0121831, 10.81, 12.011,
              14.007, 15.999, 18.998403163, 20.1797, 22.98976928,
              24.305, 26.9815385, 28.085, 30.973761998, 32.06, 35.45,
              39.948, 39.0983, 40.078, 44.955908, 47.867, 50.9415,
              51.9961, 54.938044, 55.845, 58.933194, 58.6934, 63.546,
              65.38, 69.723, 72.630, 74.921595, 78.971, 79.904,
              83.798, 85.4678, 87.62, 88.90584, 91.224, 92.90637,
              95.95, 97.90721, 101.07, 102.90550, 106.42, 107.8682,
              112.414, 114.818, 118.710, 121.760, 127.60, 126.90447,
              131.293, 132.90545196, 137.327, 138.90547, 140.116,
              140.90766, 144.242, 144.91276, 150.36, 151.964, 157.25,
              158.92535, 162.500, 164.93033, 167.259, 168.93422,
              173.054, 174.9668, 178.49, 180.94788, 183.84, 186.207,
              190.23, 192.217, 195.084, 196.966569, 200.592, 204.38,
              207.2, 208.98040, 208.98243, 209.98715, 222.01758,
              223.01974, 226.02541, 227.02775, 232.0377, 231.03588,
              238.02891, 237.04817, 244.06421, 243.06138, 247.07035,
              247.07031, 251.07959, 252.0830, 257.09511, 258.09843,
              259.1010, 262.110, 267.122, 268.126, 271.134, 270.133,
              269.1338, 278.156, 281.165, 281.166, 285.177, 286.182,
              289.190, 289.194, 293.204, 293.208, 294.214]

# No. core electrons. Incomplete for now: only goes up to Kr
atom_ncore = [0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 10, 10, 10, 10, 10, 10,
              10, 10, 10, 10, 10, 16, 16, 16, 16, 16, 16, 16, 16, 16,
              16, 26, 26, 26, 26, 26]

point_grps = ['c1','ci','c2','cs','c2h','c2v','d2','d2h']
nirrep     = [1, 2, 2, 2, 4, 4, 4, 8]

class Molecule:
    """Class constructor for the Molecule object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.xyz_file = None
        self.units    = 'Angstrom'
        self.use_sym  = False
        self.sym_grp  = None
        self.basis    = dict() 
        self.ri_basis = None 
        self.ao_cart  = False
        self.use_df   = True
        self.label    = 'default'
        self.add_rydberg = None

        # the following are determined based on user
        # input
        self.multi_geom = False
        self.charge     = None 
        self.mult       = None
        self.asym       = []
        self.crds       = None 
        self.masses     = []
        self.full_sym   = ''
        self.comp_sym   = ''
        self.sym_indx   = -1
        self.irreplbl   = None
        self.enuc       = 0.
        self.mol_obj    = None
        self.nao        = None
        self.basis_obj  = None

    def copy(self):
        """create of deepcopy of self"""
        new = Molecule()

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            if key == 'mol_obj':
                setattr(new, key, None)
            else:
                setattr(new, key, copy.deepcopy(value))

        if self.mol_obj is not None:
            new.run()

        return new


    def run(self):
        """return a gto.Molecule object: 
           compute it if one doesn't exist, or force_make=True"""
        
        # if at this point we _still_ don't have a structure, pull
        # the plug
        if len(self.asym) == 0 or self.crds is None:
            print("No geometry found for Molecule section")
            sys.exit(1)

        atms = self.asym
        cart = self.crds
        mol_str = ';'.join([atms[i]+'   '+
             ' '.join([str(cart[i,j]) for j in range(3)])
                       for i in range(len(atms))])
        
        # set the atom masses
        atm_strip = [atm.strip() for atm in atms]
        if any([atm not in atom_name for atm in atm_strip]):
            sys.exit('atom in '+str(atm_strip)+" not recognized...")
        self.asym  = atm_strip
        self.masses = np.array([atom_mass[atom_name.index(atm_strip[i])]
                          for i in range(len(self.asym))], dtype=float)

        # make the basis set objects from the string alias basis names
        self.make_basis_obj()

        # if we're using symmetry, check if pt grp is specified
        if self.sym_grp is not None:
            pt_grp  = self.sym_grp
        else:
            pt_grp  = self.use_sym

        # set the pyscf spin to self.mult-1, if not None
        if self.mult is None:
            spin = None
        else:
            spin = self.mult - 1

        self.mol_obj = gto.M(
                     dump_input = False,
                     parse_arg  = False,
                     verbose    = logger.NOTE,
                     atom       = mol_str,
                     charge     = self.charge,
                     spin       = spin,
                     output     = None,
                     cart       = self.ao_cart,
                     basis      = self.basis_obj,
                     symmetry   = pt_grp,
                     unit       = self.units)

        # do a quick check on symmetry: we will convert Coov and 
        # Dooh to C2v and D2h, respectively`
        if self.mol_obj.topgroup.lower() == 'coov':
            self.mol_obj.symmetry_subgroup  = 'C2v'
        if self.mol_obj.topgroup.lower() == 'dooh':
            self.mol_obj.symmetry_subgroup = 'D2h'
        self.mol_obj.build()

        # the nuclear repulsion energy
        self.enuc     = self.mol_obj.get_enuc()
        # full point group symmetry of the molecule
        self.full_sym = self.mol_obj.topgroup.lower()
        # point group to be used for computation
        self.comp_sym = self.mol_obj.groupname.lower()
        # number of atomic orbitals
        self.nao      = self.mol_obj.nao_nr()

        if self.mol_obj.symmetry:
            self.sym_indx = point_grps.index(self.comp_sym)
            self.irreplbl = self.mol_obj.irrep_name
        else:
            self.sym_indx = -1
            self.irreplbl = ['A']

        return

    # 
    def add_atom(self, atom_obj):
        """add an atom to existing molecule object, which 
           includes coordinates and basis set

        Args:

        Returns:
            None
        """
    
        return

    #
    def make_basis_obj(self):
        """make basis set objects using alias names"""

        # check the basis set exists in pyscf library
        self.basis_obj    = {}
        if self.ri_basis is None:
            self.ri_basis = {}

        # first go through and make sure asym and basis atom
        # labels are consistent
        
        for iasym in self.asym:
            acase = [iasym.lower(), iasym.upper(), iasym.capitalize()]
            for atom, bname in self.basis.items():
                if atom in acase and atom != iasym:
                    self.basis[iasym] = self.basis[atom]
        for atom,bname in self.basis.items():
            if atom not in self.asym:
                del self.basis[atom]

        # now go through and make the basis object
        for atom, bname in self.basis.items():
            alias = bname.lower().replace('-','').replace('_','')

            # highest priority is the user-specified local directory 
            if alias in basis.local_basis_sets(local_dir=True, 
                                               return_alias=True):
                self.basis_obj[atom] = basis.load_basis(atom, alias,
                                                        local_dir=True)

            # ...next the GRACI source directory basis sets
            elif alias in basis.local_basis_sets(return_alias=True):
                self.basis_obj[atom] = basis.load_basis(atom, alias)

            # next is alias if found in PySCF library, use those
            elif alias in gto.basis.ALIAS.keys():
                self.basis_obj[atom] = gto.basis.load(bname, atom)

            else:
                # lastly, try loading as a basis format string using
                # pyscf basis loader
                try:
                    self.basis_obj[atom] = gto.basis.parse(bname)
                except:
                    sys.exit('Basis: ' + str(bname) +
                             ' for ' + atom + ' not found.')

            # if using density-fitting, set up auxiliary basis
            if self.use_df:

                # if not specified, make a default
                if atom not in self.ri_basis.keys():

                    if alias in df.addons.DEFAULT_AUXBASIS.keys():
                        aux_name = df.addons.DEFAULT_AUXBASIS[alias][0]
                        self.ri_basis[atom] = aux_name
                    else:
                        self.ri_basis[atom] = None

                # else, if specified, make sure it something that is
                # supported
                else:

                    # None is OK -- will result in even-tempered
                    if self.ri_basis[atom] is not None:
                        a_name  = self.ri_basis[atom].lower()
                        a_alias = a_name.replace('-','').replace('_','')

                        # is in pyscf directory
                        if a_alias not in gto.basis.ALIAS.keys():
                            sys.exit('Auxiliay basis: ' + str(a_alias) +
                                     ' for atom '+ atom + ' not found.')

        # if we don't have RI-basis for each atom, use even-tempered 
        # for all atoms -- else PySCF squawks
        if any(b is None for b in self.ri_basis.values()):
            self.ri_basis = None

        return

    #
    def read_xyz(self):
        """
        read the geometry from the xyz_file specified by 'xyz_file'
        """

        # parse contents of xyz file
        try:
            with open(self.xyz_file, 'r') as xyzfile:
                xyz_gm = xyzfile.readlines()
        except:
            output.print_message('xyz_file: '
                  +str(self.xyz_file)+' not found.')
            sys.exit()

        # use the number of atoms rather than number of lines in file
        natm       = int(xyz_gm[0].strip())
        self.asym  = []
        xyz        = []

        # do we have a multi-geometry xyz file?
        ngm = int(len(xyz_gm) / (natm+2))
        if ngm > 1:
            self.multi_geom = True

        # parse the geometry
        for i in range(2, natm+2):
            line = xyz_gm[i].strip().split()
            try:
                name_capitalized = line[0][0].upper()+line[0][1:]                
                atm_indx = atom_name.index(name_capitalized)
                self.asym.append(atom_name[atm_indx])
            except ValueError:
                sys.exit('atom '+str(line.strip()[0])+' not found.')

            try:
                xyz.append([float(line[j]) for j in range(1,4)])
            except:
                sys.exit('Cannot interpret input as a geometry')

        self.crds = np.array(xyz, dtype=float)

        return

    # 
    def set_geometry(self, atms, coords):
        """set the atoms and cartesian coordinates of the 
           Molecule object"""

        # if an xyz file is specified, this is default for reading
        # coordinates
        if self.xyz_file is not None and \
                              (len(atms)==0 and len(coords)==0):
            self.read_xyz()
            return

        # else use the values in atms and coords
        if all([atm.capitalize() in atom_name for atm in atms]):
            self.asym = [atm.capitalize() for atm in atms]
        else:
            print("ATOM list: "+str(atms)+
                  " contains unrecognized atom. Valid atoms are: "+
                  str(atom_name))
            sys.exit(1)

        try:
            self.crds = np.asarray(coords, dtype=float)
        except:
            print("coords need to be interpretable as a numpy array: "+
                    str(coords))
            sys.exit(1)

        if self.crds.shape != (len(self.asym), 3):
            print(" coordinate array needs to dimensions (dim,3): "+
                  str(self.crds.shape))
            sys.exit(1)         
            
        return

    # 
    def coords_updated(self):
        """return false if the pymol coordinates don't matach the 
           graci.molecule coordinates"""

        thrsh = 1.e-6
        pyc   = self.pymol().atom_coords(unit=self.units)

        if pyc.shape != self.crds.shape or \
                        np.linalg.norm(self.crds - pyc) > thrsh:
            return True
        else:
            return False

    # 
    def pymol(self):
        """return the pymol object"""
        return self.mol_obj

    # 
    def n_irrep(self):
        """return the number of irreps"""
        if self.comp_sym != 'c1':
            return nirrep[self.sym_indx]
        else:
            return 1

    # 
    def cart(self):
        """returns the cartesian geometry as an np array"""
        return self.crds

    #
    def atoms(self):
        """ returns a list of atoms"""
        return self.asym

    #
    def masses(self):
        """return the atomic masses"""
        return self.masses

    # 
    def n_atoms(self):
        """returns the number of atoms in the molecule"""
        return len(self.asym)

    #
    def nuc_charge_center(self):
        """
        determine the center of nuclear charge
        """
        coc_xyz = np.zeros(3, dtype=float)
        coords  = self.mol_obj.atom_coords()
        chg     = self.mol_obj.atom_charges()

        for ixyz in range(len(self.mol_obj.atom_coords())):
            coc_xyz += chg[ixyz]*coords[ixyz]

        if self.units == 'Angstrom':
            coc_xyz *= constants.bohr2ang
 
        return coc_xyz / np.sum(chg)

