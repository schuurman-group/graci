"""
Module for computing DFT/MRCI energies
"""
import graci.methods.scf as scf
import graci.tools.init_libs as init_libs
import graci.tools.ref_space as ref_space
import graci.tools.ref_diag as ref_diag
import graci.tools.mrci_space as mrci_space
import graci.tools.mrci_diag as mrci_diag
import graci.tools.mrci_refine as mrci_refine

class Cvsdftmrci:
    """Class constructor for DFT/MRCI object"""
    def __init__(self):
        # user defined quanties
        self.nstates        = []
        self.hamiltonian    = 'canonical'
        self.de_select      = 0.8
        self.ras1           = []
        self.ras2           = []
        self.ras3           = []
        self.nhole1         = 0
        self.nelec3         = 0
        self.autoras        = False
        self.ciorder        = 2
        self.refiter        = 5
        self.asci           = 'off'
        self.diag_method    = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False
        self.name           = ''

        # class variables
        self.asci_thresh    = {'tight'  : 1e-5,
                               'normal' : 1e-4,
                               'loose'  : 1e-3}
        self.ref_conf       = None
        self.mrci_conf      = None

    def run(self, mol, scf):
        """ compute the DFT/MRCI energy for nroots """

        # run the KS-DFT computation 
        scf.run(mol)

        # initialize int_pyscf
        lib_intpyscf = init_libs.init_intpyscf(mol, scf)

        # initialize bitci
        lib_bitci = init_libs.init_bitci(mol, scf, self.hamiltonian)

        # generate the reference space configurations
        self.ref_conf = self.Wavefunction()
        ref_space.generate(scf, self, lib_bitci)

        # Perform the MRCI iterations, refining the reference space
        # as we go
        for i in range(self.refiter):
            # reference space diagonalisation
            ref_diag.diag(mol, self, lib_bitci)

            # generate the MRCI configurations
            self.mrci_conf = self.Wavefunction()
            mrci_space.generate(self, lib_bitci)

            # MRCI diagonalisation
            mrci_diag.diag(self, lib_bitci)

            # refine the reference space
            min_norm = mrci_refine.refine_ref_space(self, lib_bitci)

            # break if the reference space is converged
            if min_norm > 0.95 and i > 0:
                print('\n * Reference Space Converged *', flush=True)
                break

        return 

    def density(self, state):
        """ computes the density matrices for the states in 
        the array 'states'"""

        return

    def slater_dets(self, state):
        """ return the slater determinant list for state 'state'"""

        return

    def csfs(self, state):
        """ return the CSF list for state 'state'"""

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

