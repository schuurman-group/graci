"""
The Molecule object and its associated functions.
"""
import os as os
import numpy as np

atom_name  = ['X', 'H', 'D', 'T', 'He', 
              'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
              'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
atom_mass  = [0.000000000, 1.007825037, 2.015650074, 3.023475111, 4.00260325,
              7.0160045, 9.0121825, 11.0093053, 12.0, 14.003074008, 15.99491464,
              18.99840325, 19.9924391, 22.9897697, 23.9850450, 26.9815413,
              27.9769284, 30.9737634, 31.9720718, 34.968852729, 39.9623831]
atom_anum  = [0., 1., 1., 1., 2., 
              3., 4., 5., 6., 7., 8., 9., 10.,
              11., 12., 13., 14., 15., 16., 17., 18.]

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

    
