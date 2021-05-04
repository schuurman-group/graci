"""
Module for computing DFT/MRCI energies
"""
import sys as sys
import graci.io.convert as convert
import ctypes as ctypes
import graci.core.libs as libs
import graci.citools.ref_space as ref_space
import graci.citools.ref_diag as ref_diag
import graci.citools.mrci_space as mrci_space
import graci.citools.mrci_diag as mrci_diag
import graci.citools.mrci_refine as mrci_refine
import graci.citools.mrci_1rdm as mrci_1rdm

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
        self.refiter        = 3
        self.prune          = 'off'
        self.diag_method    = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False
        self.label          = 'dftmrci'

        # class variables
        self.prune_thresh   = {'tight'  : 0.9950,
                               'normal' : 0.9925,
                               'loose'  : 0.9900}
        self.niter          = 0
        self.ref_wfn        = None
        self.mrci_wfn       = None

    def name(self):
        """ return the name of the class object as a string"""
        return 'dftmrci'

    def run(self, mol, scf):
        """ compute the DFT/MRCI eigenpairs for all irreps """

        # run the KS-DFT computation 
        scf.run(mol)

        # initialize bitci
        libs.init_bitci(mol, scf, self)
        
        # generate the reference space configurations
        self.ref_wfn = self.Wavefunction()
        ref_space.generate(scf, self)

        # Perform the MRCI iterations, refining the reference space
        # as we go
        for i in range(self.refiter):
            # Update the MRCI iteration number
            self.niter += 1
            
            # reference space diagonalisation
            ref_diag.diag(mol, self)
            
            # generate the MRCI configurations
            self.mrci_wfn = self.Wavefunction()
            mrci_space.generate(scf, self)
        
            # MRCI diagonalisation
            mrci_diag.diag(self)
        
            # refine the reference space
            min_norm = mrci_refine.refine_ref_space(self)
        
            # break if the reference space is converged
            if min_norm > 0.9025 and i > 0:
                print('\n * Reference Space Converged *', flush=True)
                break

        # Compute the 1-RDMs for all states
        mrci_1rdm.rdm(self, scf)
            
        # Finalise the bitCI library
        libs.finalise_bitci()

        return 
    
    #
    def n_states(self, irrep):
        """number of states to compute"""

        if irrep < len(self.nstates):
            return self.nstates[irrep]
        else:
            print("irrep > nirrep, irrep = "+str(irrep))
            sys.exit()

    #
    def energy(self, state):
        """return the energy of state 'state'"""

        return self.mrci_wfn.ener[state]

    #
    def density(self, states, irr, mol, scf):
        """return the density matrices for the states in 
        the array 'states' for the irrep 'irr'"""
        
        return 

    #
    def slater_dets(self, state):
        """ return the slater determinant list for state 'state'"""

        return

    #
    def csfs(self, state):
        """ return the CSF list for state 'state'"""

        return

    #
    class Wavefunction:
        """wavefunction object for DFT/MRCI method"""
        def __init__(self):
            # List of numbers of configurations per irrep
            self.nconf       = None
            # bitci configuration scratch file names
            self.conf_name   = None
            # bitci configuration scratch file numbers
            self.conf_units  = None
            # List of bitci eigenvector scratch file names (one per irrep)
            self.ci_name     = None
            # List of bitci eigenvector scratch file numbers (one per irrep)
            self.ci_units    = None
            # State energies
            self.ener        = None
            # Hamiltonian integer label
            self.hamiltonian = None
            # Density matrices (list of numpy arrays, one per irrep,
            # shape: (nirr,nmo,nmo,nstates))
            self.dmat        = None
            
        #
        def set_nconf(self, nconf):
            """Sets the numbers of configurations"""
            self.nconf = nconf
            return

        #
        def set_confunits(self, confunits):
            """Sets the bitci configuration scratch file numbers"""
            self.conf_units = confunits
            return

        #
        def set_ciunits(self, ciunits):
            """Adds the list of bitci eigenvector scratch file numbers"""
            self.ci_units = ciunits
            return

        #
        def set_confname(self, confname):
            """Sets the bitci configuration scratch file names"""
            self.conf_name = confname
            return

        #
        def set_ciname(self, ciname):
            """Adds the list of bitci eigenvector scratch file names"""
            self.ci_name = ciname
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

        #
        def set_dmat(self, dmat):
            """Adds the list of 1-RDMs"""
            self.dmat = dmat
            return
