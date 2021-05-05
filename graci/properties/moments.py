"""
module for compute moments of the dipole operator
for a given electronic state
"""
import numpy as np
import graci.utils.constants as constants

class Moments:
    """Moment class for determing permanent and transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # user defined quantities
        nstates           = None
        
        # these variables are set internally
        self.dipole_vec  = None
        self.quad_tensor = None
        self.second_momt = None
        self.label       = 'default'

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'moments'

    #
    def run(self, mol, scf, method):
        """return the dipole moments for states in 'states'"""

        # compute the dipole moment integrals dimensions = (3,nmo,nmo)
        mu_ao = mol.pymol().intor('int1e_r')
        # compute the quadrupole tensor, dimensions = (9, nmo, nmo)
        q_ao  = mol.pymol().intor('int1e_rr')
        # compute the second moment
        q2_ao = mol.pymol().intor('int1e_r2')

        # load the dipole vector and quadrupole arrays
        mu_shape = (len(method.nstates), max(method.nstates), 3) 
        q_shape  = (len(method.nstates), max(method.nstates), 3, 3)
        q2_shape = (len(method.nstates), max(method.nstates))

        self.nstates     = method.nstates
        self.dipole_vec  = np.zeros(mu_shape, dtype=float)
        self.quad_tensor = np.zeros(q_shape, dtype=float)
        self.second_momt = np.zeros(q2_shape, dtype=float)

        for irr in range(len(method.nstates)):
            for st in range(method.nstates[irr]):
                occ, no_ao = method.natural_orbs(irr, st, basis='ao')
                rdm_ao = np.matmul(np.matmul(no_ao, np.diag(occ)), no_ao.T)
                self.second_momt[irr, st] = np.sum(rdm_ao*q2_ao)
                for i in range(3):
                    self.dipole_vec[irr,st,i]  = np.sum(rdm_ao*mu_ao[i])
                    for j in range(3):
                        ij = 3*i + j
                        self.quad_tensor[irr,st,i,j] = np.sum(rdm_ao*q_ao[ij])
        return

    #
    def dipole(self, irrep, state):
        """return the dipole array"""

        if self.dipole_vec is None:
            print("Must call moments.run() to populate dipole vector")
            return None

        return self.dipole_vec[irrep, state, :]

    #
    def quadrupole(self, irrep, state):
        """return the quadrupole tensor"""

        if self.quad_tensor is None:
            print("Must call moments.run() to populate quad tensor")
            return None

        return self.quad_tensor[irrep, state, :, :]

    #
    def second_moment(self, irrep, state):
        """return the second moment"""

        if self.second_momt is None:
            print("Must call moments.run() to determine second moment")
            return None

        return self.second_momt[irrep, state]
