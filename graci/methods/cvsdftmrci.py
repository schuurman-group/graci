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
        # user defined variables
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
        self.diag_algorithm = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False
        self.name           = ''

        # class variables
        self.asci_thresh    = {'tight'  : 1e-5,
                               'normal' : 1e-4,
                               'loose'  : 1e-3}

    def run(self, mol, scf):
        """ compute the DFT/MRCI energy for nroots """

        # run the KS-DFT computation 
        scf.run(mol)

        # initialize int_pyscf
        lib_intpyscf = init_libs.init_intpyscf(mol, scf)

        # initialize bitci
        lib_bitci = init_libs.init_bitci(mol, scf, self.hamiltonian)

        # number of irreps -- should check this is constant with 
        # length self.nstates?
        nirr = mol.n_irrep()
  
        # generate the reference space configurations
        conf0 = ref_space.generate(scf, self, lib_bitci)

        # Perform the MRCI iterations, refining the reference space
        # as we go
        for i in range(self.refiter):
            # reference space diagonalisation
            ref_diag.diag(ci.nstates, conf0, lib_bitci)

            # generate the MRCI configurations
            confsd = mrci_space.generate(self, conf0, lib_bitci)

            # MRCI diagonalisation
            mrci_diag.diag(self, confsd, lib_bitci)

            # refine the reference space
            min_norm = mrci_refine.refine_ref_space(self, conf0, confsd, lib_bitci)

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




