"""
Module for computing DFT/MRCI energies
"""
import sys as sys
import graci.io.convert as convert
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.citools.ref_space as ref_space
import graci.citools.ref_diag as ref_diag
import graci.citools.mrci_space as mrci_space
import graci.citools.mrci_diag as mrci_diag
import graci.citools.mrci_refine as mrci_refine
import graci.citools.mrci_1rdm as mrci_1rdm
import graci.io.output as output
import graci.properties.moments as moments

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
        self.prune_thresh   = {'tight'  : 1e-4,
                               'normal' : 1e-3,
                               'loose'  : 3e-3}
        self.niter          = 0
        self.ref_wfn        = None
        self.mrci_wfn       = None
        self.print_quad     = False

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
        self.ref_wfn = self.Wavefunction(orbs=scf.orbs)
        ref_space.generate(scf, self)

        # Perform the MRCI iterations, refining the reference space
        # as we go
        for i in range(self.refiter):
            # Update the MRCI iteration number
            self.niter += 1
            
            # reference space diagonalisation
            ref_diag.diag(mol, self)
            
            # generate the MRCI configurations
            self.mrci_wfn = self.Wavefunction(orbs=scf.orbs)
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

        # by default, print out the natural orbitals for
        # each state
        for irr in range(len(self.nstates)):
            for st in range(self.nstates[irr]):
                occ,orb = self.natural_orbs(irr, st)
                output.print_nos_molden(mol, irr, st, occ, orb)

        # we'll also compute 1-electron properties by
        # default.
        momts = moments.Moments()
        momts.run(mol, scf, self)
        for irr in range(len(self.nstates)):
            output.print_moments_header(mol.irreplbl[irr])
            for st in range(self.nstates[irr]):
                output.print_moments(mol.irreplbl[irr], st, 
                        momts.dipole(irr,st), 
                        momts.second_moment(irr,st))
                if self.print_quad:
                    output.print_quad(mol.irreplbl[irr], st, 
                            momts.quadrupole(irr,st))

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
    def rdm1(self, irr, state):
        """return the density matrix for the state in 
        the array 'states' for the irrep 'irr'"""

        if self.mrci_wfn.dmat is not None:
            return self.mrci_wfn.dmat[irr][ :, :, state]
        else:
            print("rdm1 called in method=dftmrci, "+
                  "but density does not exist")
            return self.mrci_wfn.dmat

    #
    def natural_orbs(self, irr, state, basis='ao'):
        """print the natural orbitals for irrep 
           'irr' and state 'state'"""

        # if natural orbitals don't exist, create them --
        # we might as well create all of them, unless a truly
        # crazy of number of states have been computed
        if (self.mrci_wfn.natorb_mo is None and 
                         self.mrci_wfn.natorb_ao is None):
            # if the 1RDMs don't exist, then we have a problem
            if self.mrci_wfn.dmat is None:
                print("can't return natural orbitals: "+
                      "1RDMs don't exist")
                return None
            
            # this is hacky
            nd = len(self.mrci_wfn.dmat)
            (nmo1, nmo2, n_max) = self.mrci_wfn.dmat[0].shape
            nirr = len(self.nstates)

            # sanity check some things
            if irr >= nd or state >= n_max:
                print("irrep/state = "+str(irr)+","+str(state)+
                      " requested, but dmat dimensions are "+
                      str(nd)+","+str(n_max))
                return None
            if nirr != nd:
                print("diagreement on number of irreps in dftmrci")
                return None

            self.mrci_wfn.natocc     = np.zeros((nirr,nmo1,n_max),
                                             dtype=float)
            self.mrci_wfn.natorb_mo  = np.zeros((nirr,nmo1,nmo2,n_max),
                                             dtype=float)
            self.mrci_wfn.natorb_ao  = np.zeros((nirr,nmo1,nmo2,n_max),
                                             dtype=float)

            # by default, compute the natural orbitals for all states
            # requested
            for irrep in range(nirr):
                for st in range(self.nstates[irrep]):
                    rdm        = self.rdm1(irrep, st)
                    occ,nos_mo = np.linalg.eigh(rdm)
                    nos_ao     = np.matmul(self.mrci_wfn.orbs, nos_mo)
                    
                    sort_ordr = self.order_orbs(occ)
                    occ       = occ[sort_ordr]
                    nos_mo    = nos_mo[:,sort_ordr]
                    nos_ao    = nos_ao[:,sort_ordr]

                    self.mrci_wfn.natocc[irrep,:,st]      = occ
                    self.mrci_wfn.natorb_mo[irrep,:,:,st] = nos_mo
                    self.mrci_wfn.natorb_ao[irrep,:,:,st] = nos_ao


        # always return the contents of the natorb/natocc arrays in wfn
        occ = self.mrci_wfn.natocc[irr, :, state]

        if basis == 'ao':
            nos = self.mrci_wfn.natorb_ao[irr, :, :, state]
        else:
            nos = self.mrci_wfn.natorb_mo[irr, :, :, state]

        return occ, nos

    #
    def slater_dets(self, state):
        """ return the slater determinant list for state 'state'"""

        return

    #
    def csfs(self, state):
        """ return the CSF list for state 'state'"""

        return


#########################################################################
    #
    def order_orbs(self, occ):
        """sort the occupation vector 'occ'. This is a dedicated 
           function for orbitals because I can envision more elaborate 
           sorting criteria than just occupation number (i.e. core orbitals 
           first, etc.)"""

        # this sorts ascending
        ordr = np.argsort(occ)

        # flip so largest occupations come first
        return np.flip(ordr)

    #
    class Wavefunction:
        """wavefunction object for DFT/MRCI method"""
        def __init__(self, orbs=None):
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
            # one particle functions used in the expansion
            self.orbs        = orbs
            # Density matrices (list of numpy arrays, one per irrep,
            # shape: (nirr,nmo,nmo,nstates))
            self.dmat        = None
            # natural orbital occupations
            self.natocc      = None
            # natural orbitals, MO basis (same format as the density matrices)
            self.natorb_mo   = None
            # natural orbitals, AO basis
            self.natorb_ao   = None

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
        def set_orbs(self, orbs):
            """set the one particle functions"""
            self.orbs = orbs
            return

        #
        def set_dmat(self, dmat):
            """Adds the list of 1-RDMs"""
            self.dmat = dmat
            return
