"""
Module for computing DFT/MRCI energies
"""
import sys as sys
import graci.io.convert as convert
import ctypes as ctypes

import graci.core.loadlibs as loadlibs
import graci.citools.ref_space as ref_space
import graci.citools.ref_diag as ref_diag
import graci.citools.mrci_space as mrci_space
import graci.citools.mrci_diag as mrci_diag
import graci.citools.mrci_refine as mrci_refine

# MRCI and DFT/MRCI Hamiltonian labels
hamiltonians   = ['canonical',
                  'grimme_standard',
                  'grimme_short',
                  'lyskov_standard',
                  'lyskov_short',
                  'heil17_standard',
                  'heil18_short']

class Dftmrci:
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
        self.icvs           = []
        self.ciorder        = 2
        self.refiter        = 5
        self.asci           = 'off'
        self.diag_method    = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False
        self.label          = 'default'

        # class variables
        self.asci_thresh    = {'tight'  : 1e-4,
                               'normal' : 1e-3,
                               'loose'  : 3e-3}
        self.niter          = 0
        self.ref_conf       = None
        self.mrci_conf      = None

    def name(self):
        """ return the name of the class object as a string"""
        return 'dftmrci'

    def run(self, mol, scf):
        """ compute the DFT/MRCI energy for nroots """

        # run the KS-DFT computation 
        scf.run(mol)

        # initialize int_pyscf
        lib_intpyscf = loadlibs.init_intpyscf(mol, scf)

        # initialize bitci
        lib_bitci = loadlibs.init_bitci(mol, scf, self)

        # generate the reference space configurations
        self.ref_conf = self.Wavefunction()
        ref_space.generate(scf, self, lib_bitci)
        
        # Perform the MRCI iterations, refining the reference space
        # as we go
        for i in range(self.refiter):
            # Update the MRCI iteration number
            self.niter += 1
            
            # reference space diagonalisation
            ref_diag.diag(mol, self, lib_bitci)

            # generate the MRCI configurations
            self.mrci_conf = self.Wavefunction()
            mrci_space.generate(scf, self, lib_bitci)

            # MRCI diagonalisation
            mrci_diag.diag(self, lib_bitci)

            # refine the reference space
            min_norm = mrci_refine.refine_ref_space(self, lib_bitci)

            # break if the reference space is converged
            if min_norm > 0.9025 and i > 0:
                print('\n * Reference Space Converged *', flush=True)
                break
        
        return 

    def nstates(self):
        """number of states to compute"""

        return self.nstates

    def energy(self, state):
        """return the energy of state 'state'"""

        return self.mrci_conf.ener[state]

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
            # bitci configuration scratch file names
            self.confname    = 0
            # List of bitci eigenvector scratch file names (one per irrep)
            self.vecname     = None
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
        def set_confname(self, confname):
            """Sets the bitci configuration scratch file names"""
            self.confname = confname
            return

        #
        def set_vecname(self, vecname):
            """Adds the list of bitci eigenvector scratch file names"""
            self.vecname = vecname
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

