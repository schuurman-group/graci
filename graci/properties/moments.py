"""
module for compute moments of the dipole operator
for a given electronic state
"""

class Moments:
    """Moment class for determing permanent and transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # user defined quantities
        states           = None
        
        # these variables are set internally
        self.dipole_vec  = None
        self.quad_tensor = None
        self.label       = 'default'

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'moments'

    #
    def run(self, mol, scf, wfn):
        """return the dipole moments for states in 'states'"""

        # compute the dipole moment integrals
        mu_ao = self.mol.pymol().intor('int1e_z')
        q_ao  = self.mol.pymol().intor('int1e_zz')

        # contract with the MOs
        mu_mo = np.matmul(np.matmul(scf.orbs.T, mu_ao), scf.orbs)
        q_mo  = np.matmul(np.matmul(scf.orbs.T, q_mo), scf.orbs)

        # what comes next depends on libSI
 

        # load the dipole vector and quadrupole arrays
        for st in states:
            self.dipole_vec[st]  = 0.
            self.quad_tensor[st] = 0.


        return

    #
    def dipole(self, state):
        """return the dipole array"""
        return self.dipole_vec[state]

    #
    def quadrupole(self, state):
        """return the quadrupole tensor"""
        return self.quad_tensor[state]

