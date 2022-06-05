"""
The Molecule object and its associated functions.
"""
import sys as sys
import numpy as np
import graci.io.output as output
from pyscf.lib import logger
from pyscf import gto

atom_name = ['X' ,'H' , 'He', 'Li', 'Be', 'B' , 'C' , 'N' , 'O' , 'F' ,
             'Ne','Na', 'Mg', 'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', 'K' ,
             'Ca','Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
             'Zn','Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y' ,
             'Zr','Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
             'Sn','Sb', 'Te', 'I' , 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr',
             'Nd','Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
             'Yb','Lu', 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au',
             'Hg','Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac',
             'Th','Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es',
             'Fm','Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
             'Ds','Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

# Atomic masses taken from PySCF elements.py
atom_mass = [
    0.,                 # X (ghost
    1.008,              # H [1.00784, 1.00811]
    4.002602,           # He
    6.94,               # Li [6.938, 6.997]
    9.0121831,          # Be
    10.81,              # B [10.806, 10.821]
    12.011,             # C [12.0096, 12.0116]
    14.007,             # N [14.00643, 14.00728]
    15.999,             # O [15.99903, 15.99977]
    18.998403163,       # F
    20.1797,            # Ne
    22.98976928,        # Na
    24.305,             # Mg [24.304, 24.307]
    26.9815385,         # Al
    28.085,             # Si [28.084, 28.086]
    30.973761998,       # P
    32.06,              # S [32.059, 32.076]
    35.45,              # Cl [35.446, 35.457]
    39.948,             # Ar
    39.0983,            # K
    40.078,             # Ca
    44.955908,          # Sc
    47.867,             # Ti
    50.9415,            # V
    51.9961,            # Cr
    54.938044,          # Mn
    55.845,             # Fe
    58.933194,          # Co
    58.6934,            # Ni
    63.546,             # Cu
    65.38,              # Zn
    69.723,             # Ga
    72.630,             # Ge
    74.921595,          # As
    78.971,             # Se
    79.904,             # Br [79.901, 79.907]
    83.798,             # Kr
    85.4678,            # Rb
    87.62,              # Sr
    88.90584,           # Y
    91.224,             # Zr
    92.90637,           # Nb
    95.95,              # Mo
    97.90721,           # 98Tc
    101.07,             # Ru
    102.90550,          # Rh
    106.42,             # Pd
    107.8682,           # Ag
    112.414,            # Cd
    114.818,            # In
    118.710,            # Sn
    121.760,            # Sb
    127.60,             # Te
    126.90447,          # I
    131.293,            # Xe
    132.90545196,       # Cs
    137.327,            # Ba
    138.90547,          # La
    140.116,            # Ce
    140.90766,          # Pr
    144.242,            # Nd
    144.91276,          # 145Pm
    150.36,             # Sm
    151.964,            # Eu
    157.25,             # Gd
    158.92535,          # Tb
    162.500,            # Dy
    164.93033,          # Ho
    167.259,            # Er
    168.93422,          # Tm
    173.054,            # Yb
    174.9668,           # Lu
    178.49,             # Hf
    180.94788,          # Ta
    183.84,             # W
    186.207,            # Re
    190.23,             # Os
    192.217,            # Ir
    195.084,            # Pt
    196.966569,         # Au
    200.592,            # Hg
    204.38,             # Tl [204.382, 204.385]
    207.2,              # Pb
    208.98040,          # Bi
    208.98243,          # Po
    209.98715,          # At
    222.01758,          # Rn
    223.01974,          # Fr
    226.02541,          # Ra
    227.02775,          # Ac
    232.0377,           # Th
    231.03588,          # Pa
    238.02891,          # U
    237.04817,          # Np
    244.06421,          # Pu
    243.06138,          # Am
    247.07035,          # Cm
    247.07031,          # Bk
    251.07959,          # Cf
    252.0830,           # Es
    257.09511,          # Fm
    258.09843,          # Md
    259.1010,           # No
    262.110,            # Lr
    267.122,            # Rf
    268.126,            # Db
    271.134,            # Sg
    270.133,            # Bh
    269.1338,           # Hs
    278.156,            # Mt
    281.165,            # Ds
    281.166,            # Rg
    285.177,            # Cn
    286.182,            # Nh
    289.190,            # Fl
    289.194,            # Mc
    293.204,            # Lv
    293.208,            # Ts
    294.214,            # Og
]

point_grps = ['c1','ci','c2','cs','c2h','c2v','d2 ','d2h']
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
        self.basis    = dict() 
        self.ri_basis = None 
        self.ao_cart  = False
        self.use_df   = False
        self.label    = 'Molecule'

        # the following are determined based on user
        # input
        self.charge   = 0
        self.mult     = 1
        self.asym     = []
        self.crds     = None 
        self.masses   = []
        self.full_sym = ''
        self.comp_sym = ''
        self.sym_indx = -1
        self.irreplbl = None
        self.enuc     = 0.
        self.mol_obj  = None
        self.nao      = None

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

        self.mol_obj = gto.M(
                     dump_input = False,
                     parse_arg  = False,
                     verbose    = logger.NOTE,
                     atom       = mol_str,
                     charge     = self.charge,
                     spin       = self.mult-1,
                     output     = None,
                     cart       = self.ao_cart,
                     basis      = self.basis,
                     symmetry   = self.use_sym,
                     unit       = self.units)

        # the nuclear repulsion energy
        self.enuc     = self.mol_obj.energy_nuc()
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
    def read_xyz(self):
        """read the xyz_file specified by 'xyz_file'"""

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

        for i in range(2, natm+2):
            line = xyz_gm[i].strip().split()
            try:
                atm_indx = atom_name.index(line[0].upper())
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
        if self.xyz_file is not None:
            self.read_xyz()
            return

        # else use the values in atms and coords
        if all([atm in atom_name for atm in atms]):
            self.asym = atms
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
 
