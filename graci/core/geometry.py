"""
The Molecule object and its associated functions.
"""
import os as os
import numpy as np
import copy as copy

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
atom_mass = [ 0., 1.008, 4.002602, 6.94, 9.0121831, 10.81, 12.011,
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
atom_ncore = [0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 10, 10, 10, 10, 10, 10,
              10, 10, 10, 10, 10, 16, 16, 16, 16, 16, 16, 16, 16, 16,
              16, 26, 26, 26, 26, 26]

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

    def copy(self):
        """create of deepcopy of self"""
        new = Geometry()

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            setattr(new, key, copy.deepcopy(value))

        return new

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

    
