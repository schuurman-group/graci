"""
The Molecule object and its associated functions.
"""

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
point_grps = ['c1','ci','c2','cs','c2h','c2v','d2 ','d2h']
nirrep     = [1, 2, 2, 2, 4, 4, 4, 8]

class Molecule:
    """Class constructor for the Molecule object."""
    def __init__(self):
        self.geom      = None
        self.atoms     = None
        self.mult      = 1
        self.spin      = 0.
        self.charge    = 0.
        self.nel       = 0
        self.use_sym   = False
        self.full_sym  = ''
        self.comp_sym  = ''
        self.sym_indx  = -1
        self.irreplbl  = None
        self.hij       = None
        self.enuc      = 0.
        self.orbs      = None
        self.orb_occ   = None
        self.orb_ener  = None
        self.orb_irrep = []
        self.orb_sym   = []
        self.nmo       = 0
        self.naux      = 0
        self.units     = 'bohr'
        self.basis     = ''
        self.ri_basis  = ''
        
    #
    def set_geom(self, geom):
        """Set the cartesian geometry"""
        self.geom = np.array(geom)
        return

    # 
    def geom(self):
        """returns the cartesian geometry"""
        return self.geom

    #
    def set_atoms(self, atoms):
        """Set the atomic labels"""
        atm_strip = [atm.strip() for atm in atoms]
        if any([atm not in atom_name for atm in atm_strip]):
            sys.exit('atom in '+str(atm_strip)+" not recognized...")
        self.atoms = atm_strip
        return

    #
    def set_orbs(self, orbs):
        """set the molecular orbitals"""
        self.orbs = orbs
        return

    # 
    def set_occ(self, occ):
        """sets the natural orbital occupations"""
        self.occ = occ
        return

    #
    def set_ener(self, ener):
        """sets the orbital energies"""
        self.ener = ener
        return

    #
    def set_mult(self, mult):
        """sets the multiplicity"""
        self.mult = mult
        self.spin = 0.5*(mult - 1.)
        return

    #
    def set_charge(self, charge):
        """ses the molecular charge"""
        self.charge = charge
        return

    #
    def sort_orbs(self, method):
        """sorts the orbitals based on value of method"""
        return

    # 
    def n_atoms(self):
        """returns the number of atoms in the molecule"""
        return len(self.atoms)

    
