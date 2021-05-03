"""
Module for computing DFT/MRCI energies
"""


class Dftcis:
    """Class constructor for SCF object"""
    def __init__(self):
        self.label   = 'dftcis'

    def name(self):
        """ return the name of the class object as a string"""
        return 'dftcis'

    def run(mol):
        """runs a dftcis calculation"""

        return

    def density(mol):
        """return the density matrix for state i"""

        return

    def slater_dets(state):
        """return the slater determinant list for state 'state'"""

        return

    def csfs(state):
        """csf list for state 'state'"""

        return

    #
    class Wavefunction:
        """wavefunction object for DFT/MRCI method"""
        def __init__(self):
            # List of numbers of configurations per irrep
            self.nconf       = None
            # bitci configuration scratch file numbers
            self.confscr     = 0
            # List of bitci eigenvector scratch file numbers (one per irrep)
            self.vecscr      = None
            # State energies
            self.ener        = None
            # Hamiltonian integer label
            self.hamiltonian = None

        #
        def set_nconf(self, nconf):
            """Sets the numbers of configurations"""
            self.nconf = nconf
            return

        #
        def set_confscr(self, confscr):
            """Sets the bitci configuration scratch file numbers"""
            self.confscr = confscr
            return

        #
        def set_vecscr(self, vecscr):
            """Adds the list of bitci eigenvector scratch file numbers"""
            self.vecscr = vecscr
            return

        #
        def set_ener(self, ener):
            """Adds the array of state energies"""
            self.ener = ener
            return

        #
        def set_hamiltonian(self, lbl):
            """Set the Hamiltonian integer label"""
            self.hamiltonian = hamiltonians.index(lbl)
            return

