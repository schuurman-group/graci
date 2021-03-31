"""
The Molecule object and its associated functions.
"""
import numpy as np
import graci.io.output as output
from pyscf.lib import logger
from pyscf import gto

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
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.mult     = 1
        self.charge   = 0.
        self.use_sym  = False
        self.units    = 'bohr'
        self.basis    = ''
        self.ri_basis = ''
        self.use_df   = False
        self.name     = ''

        # the following are determined based on user
        # input
        self.nel      = 0
        self.spin     = 0
        self.geom     = None
        self.atoms    = None
        self.full_sym = ''
        self.comp_sym = ''
        self.sym_indx = -1
        self.irreplbl = None
        self.enuc     = 0.
        self.mol_obj  = None

    def pymol(self, force_make=False):
        """return a gto.Molecule object: 
           compute it if one doesn't exist, or force_make=True"""
  
        if self.mol_obj is None or force_make:
            atms = self.atoms
            cart = self.geom
            mol_str = ';'.join([atms[i]+'   '+
                 ' '.join([str(cart[i,j]) for j in range(3)])
                           for i in range(self.n_atoms())])

            self.mol_obj = gto.M(
                         dump_input = False,
                         parse_arg  = False,
                         verbose    = logger.NOTE,
                         atom       = mol_str,
                         charge     = self.charge,
                         spin       = self.spin,
                         output     = output.file_names['pyscf_out'],
                         basis      = self.basis,
                         symmetry   = self.use_sym,
                         unit       = self.units)

        return self.mol_obj

    # 
    def n_irrep(self):
        """return the number of irreps"""
        if self.comp_sym != 'c1':
            return nirrep[self.sym_indx]
        else:
            return 1

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
    def n_atoms(self):
        """returns the number of atoms in the molecule"""
        return len(self.atoms)

    
