"""
module for compute moments of the dipole operator
for a given electronic state
"""
import copy as copy
import numpy as np
import graci.utils.timing as timing

class Moments:
    """Moment class for determing permanent and transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # these variables are set internally
        self.dipole_vec  = None
        self.quad_tensor = None
        self.second_momt = None
        self.label       = 'Moments'

    def copy(self):
        """create of deepcopy of self"""
        new = Moments() 

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))

        return new

    #
    @timing.timed
    def run(self, mol, natocc, natorb):
        """return the dipole moments for states in 'states'"""

        # compute the dipole moment integrals dimensions = (3,nao,nao)
        mu_ao = mol.pymol().intor('int1e_r')
        # compute the quadrupole tensor, dimensions = (9, nao, nao)
        q_ao  = mol.pymol().intor('int1e_rr')
        # compute the second moment
        q2_ao = mol.pymol().intor('int1e_r2')

        # load the dipole vector and quadrupole arrays
        (nst, nao, nmo)  = natorb.shape
        self.dipole_vec  = np.zeros((nst, 3), dtype=float)
        self.quad_tensor = np.zeros((nst, 3, 3), dtype=float)
        self.second_momt = np.zeros((nst,), dtype=float)

        for ist in range(nst):
            no_ao = natorb[ist, :, :]
            occ   = natocc[ist, :]
            rdm_ao = np.matmul(np.matmul(no_ao, np.diag(occ)), no_ao.T)

            self.dipole_vec[ist,:] = np.einsum('xij,ij->x',mu_ao, rdm_ao)

            self.second_momt[ist]  = np.sum(rdm_ao*q2_ao)

            qtens = np.einsum('xij,ij->x',q_ao.reshape((9,nao,nao)), 
                                                                  rdm_ao)

            # traceless quad tensor defined as: 
            # Q = 3*xi*xj - r^2 * delta_ij
            self.quad_tensor[ist,:,:]  = 3*qtens.reshape((3,3)) - \
                  self.second_momt[ist] * np.identity(3, dtype=float)

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
