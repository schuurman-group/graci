"""
The Molecule object and its associated functions.
"""
import os as os
import numpy as np

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

class Geometry:
    """Class constructor for the Molecule object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.xyz_file  = None
        self.units     = 'angstrom'
        self.label     = 'Geometry'

        # the following are determined from a combination
        # of user input
        self.cartesian = None
        self.asym      = None
        self.masses    = None

    #
    def read_xyz():
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
        atoms      = []
        cartesian  = []

        for i in range(2, natm+2):
            line = xyz_gm[i]
            try:
                atm_indx = self.atom_name.index(line[0].upper())
                atoms.append(self.atom_name[atm_indx])
            except ValueError:
                sys.exit('atom '+str(line[0])+' not found.')
            try:
                cartesian.append([float(line[i] * conv)
                    for i in range(1,4)])
            except:
                sys.exit('Cannot interpret input as a geometry')

        self.set_atoms(atoms)
        self.set_geom(cartesian)

        return

    #
    def set_geom(self, geom):
        """Set the cartesian geometry"""
        self.cartesian = np.array(geom, dtype=float)
        return

    #
    def set_atoms(self, atoms):
        """Set the atomic labels"""
        atm_strip = [atm.strip() for atm in atoms]
        if any([atm not in atom_name for atm in atm_strip]):
            sys.exit('atom in '+str(atm_strip)+" not recognized...")
        self.asym  = atm_strip
        self.masses = np.array([atom_mass[atom_name.index(atm_strip[i])]
                          for i in range(len(self.asym))], dtype=float)
        return

    #
    def set_masses(self, masses):
        """set the atomic masses"""
        self.masses = np.array(masses, dtype=float)
        return

    # 
    def geom(self):
        """returns the cartesian geometry as an np array"""
        return self.cartesian

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

    
