"""
The Molecule object and its associated functions.
"""
import numpy as np
import graci.io.output as output
from pyscf.lib import logger
from pyscf import gto

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
        self.basis    = ''
        self.ri_basis = ''
        self.use_df   = False
        self.use_rrdf = False
        self.rrdf_fac = 3
        self.label    = 'default'

        # the following are determined based on user
        # input
        self.nel      = 0
        self.spin     = 0
        self.full_sym = ''
        self.comp_sym = ''
        self.sym_indx = -1
        self.irreplbl = None
        self.enuc     = 0.
        self.mol_obj  = None

    def name(self):
        """ return the name of the class object as a string"""
        return 'molecule'

    def set_pymol(self, gm_obj):
        """return a gto.Molecule object: 
           compute it if one doesn't exist, or force_make=True"""
  
        atms = gm_obj.atoms()
        cart = gm_obj.geom()
        mol_str = ';'.join([atms[i]+'   '+
             ' '.join([str(cart[i,j]) for j in range(3)])
                       for i in range(gm_obj.n_atoms())])

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
                     unit       = gm_obj.units)

        # the nuclear repulsion energy
        self.enuc     = self.mol_obj.energy_nuc()
        # the number of electrons
        self.nel      = int(sum(self.mol_obj.nelec))
        # full point group symmetry of the molecule
        self.full_sym = self.mol_obj.topgroup.lower()
        # point group to be used for computation
        self.comp_sym = self.mol_obj.groupname.lower()
        if self.mol_obj.symmetry:
            self.sym_indx = point_grps.index(self.comp_sym)
            self.irreplbl = self.mol_obj.irrep_name
        else:
            self.sym_indx = -1
            self.irreplbl = ['A']

        return

    # checks whether there is a valid pyscf molecule object set
    def pymol_exists(self):
        """return True if self.mol_obj is not 'None' """
        if self.mol_obj is None:
            return False
        else:
            return True

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

    
