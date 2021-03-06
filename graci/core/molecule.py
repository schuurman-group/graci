"""
The Molecule object and its associated functions.
"""
import sys as sys
import numpy as np
import graci.io.output as output
from pyscf.lib import logger
from pyscf import gto

atom_name  = ['X', 'H', 'D', 'T', 'He',
              'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
              'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
              'Fe']
atom_mass  = [0.000000000, 1.007825037, 2.015650074, 3.023475111, 4.00260325,
              7.0160045, 9.0121825, 11.0093053, 12.0, 14.003074008, 15.99491464,
              18.99840325, 19.9924391, 22.9897697, 23.9850450, 26.9815413,
              27.9769284, 30.9737634, 31.9720718, 34.968852729, 39.9623831,
              55.847]
atom_anum  = [0., 1., 1., 1., 2.,
              3., 4., 5., 6., 7., 8., 9., 10.,
              11., 12., 13., 14., 15., 16., 17., 18.,
              26.]
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
 
