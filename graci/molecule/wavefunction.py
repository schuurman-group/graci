"""
The Wavefunction object and its associated functions.
"""

class Wavefunction:
    """Class constructor for the Wavefunction object"""
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
