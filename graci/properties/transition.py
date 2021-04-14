"""
module for compute moments of the dipole operator
for a given electronic state
"""


class Transition:
    """Moment class for determing permanent and transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # user defined quanties 
        self.init        = None
        self.final       = None

        # passed and global variables
        self.trans_list  = None
        self.trans_mom   = None
        self.osc_str     = None
        self.label       = 'default'

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'transition'

    #
    def run(self, mol, scf, bra, ket):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used."""

        # compute the dipole moment integrals
        mu_ao = self.mol.pymol().intor('int1e_z')
        q_ao  = self.mol.pymol().intor('int1e_zz')

        # contract with the MOs
        mu_mo = np.matmul(np.matmul(scf.orbs.T, mu_ao), scf.orbs)
        q_mo  = np.matmul(np.matmul(scf.orbs.T, q_mo), scf.orbs)

        # construct the trans_list array
        # currently store them as [initial state, final state]
        self.trans_list = []
        for ket in np.sort(self.init):
            for bra in np.sort(self.final):
                self.trans_list.append([ket, bra])

        # what happens depends on the libSI interface




        # compute the oscillator strengths from the transition dipoles
        self.osc_str = []
        cnt       = 0
        for exc_pair in self.trans_list:
            exc_ener = bra.energy(exc_pair[1]) - ket.energy(exc_pair[0])
            self.osc_str.extend(2.*exc_ener*(self.trans_mom[cnt]**2)/3.)
            cnt += 1

        return

    def trans_dipole(self, bra, ket):
        """return transition dipole matrix element"""

         # find the state pair in the transition list
         ind = self.tran_list.index([bra, ket])

        return self.trans_mom[ind]

    # 
    def osc_strength(self, bra, ket):
        """return the transition dipole moments for all pairs of
           bra_states and ket_states. If one or both are not specified,
           return all transition moments"""

        # find the state pair in the transition list
        ind = self.tran_list.index([bra, ket])

        return self.osc_str[ind]



