"""
module for compute moments of the dipole operator
for a given electronic state
"""
import numpy as np
import graci.utils.timing as timing

class Moments:
    """Moment class for determing permanent and transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self, mol, natocc, natorb_ao):
        # user defined quantities
        self.mol         = mol
        self.natocc      = natocc
        self.natorb      = natorb_ao
        
        # these variables are set internally
        self.dipole_vec  = None
        self.quad_tensor = None
        self.second_momt = None
        self.label       = 'Moments'

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'moments'

    #
    @timing.timed
    def run(self):
        """return the dipole moments for states in 'states'"""

        # compute the dipole moment integrals dimensions = (3,nmo,nmo)
        mu_ao = self.mol.pymol().intor('int1e_r')
        # compute the quadrupole tensor, dimensions = (9, nmo, nmo)
        q_ao  = self.mol.pymol().intor('int1e_rr')
        # compute the second moment
        q2_ao = self.mol.pymol().intor('int1e_r2')

        # load the dipole vector and quadrupole arrays
        (nst, nmo1, nmo2) = self.natorb.shape
        self.dipole_vec  = np.zeros((nst, 3), dtype=float)
        self.quad_tensor = np.zeros((nst, 3, 3), dtype=float)
        self.second_momt = np.zeros((nst,), dtype=float)

        for ist in range(nst):
            no_ao = self.natorb[ist, :, :]
            occ   = self.natocc[ist, :]
            rdm_ao = np.matmul(np.matmul(no_ao, np.diag(occ)), no_ao.T)
            self.second_momt[ist]  = np.sum(rdm_ao*q2_ao)
            self.dipole_vec[ist,:] = np.einsum('xij,ij->x',mu_ao, rdm_ao)
            
            qtens = np.einsum('xij,ij->x',q_ao.reshape((9,nmo1,nmo2)), 
                                                                  rdm_ao)
            self.quad_tensor[ist,:,:]  = qtens.reshape((3,3))

        return

    #
    def dipole(self, ist):
        """return the dipole array"""

        if self.dipole_vec is None:
            print("Must call moments.run() to populate dipole vector")
            return None

        return self.dipole_vec[ist, :]

    #
    def quadrupole(self, ist):
        """return the quadrupole tensor"""

        if self.quad_tensor is None:
            print("Must call moments.run() to populate quad tensor")
            return None

        return self.quad_tensor[ist, :, :]

    #
    def second_moment(self, ist):
        """return the second moment"""

        if self.second_momt is None:
            print("Must call moments.run() to determine second moment")
            return None

        return self.second_momt[ist]
