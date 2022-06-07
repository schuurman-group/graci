"""
Module for computing DFT/MRCI(2) states and energies
"""
import sys as sys
import numpy as np
import copy as copy
import graci.utils.timing as timing
import graci.methods.cimethod as cimethod
import graci.core.bitciwfn as bitciwfn
import graci.io.output as output
import graci.bitcitools.bitci_init as bitci_init
import graci.bitcitools.ref_space as ref_space
import graci.bitcitools.ref_diag as ref_diag
import graci.bitcitools.ref_prune as ref_prune
import graci.bitcitools.mrci_space as mrci_space
import graci.bitcitools.gvvpt2 as gvvpt2
import graci.bitcitools.gvvpt2_refine as gvvpt2_refine
import graci.bitcitools.mrci_1rdm as mrci_1rdm
import graci.bitcitools.mrci_wf as mrci_wf

# MRCI and DFT/MRCI Hamiltonian labels
hamiltonians   = ['canonical',
                  'grimme_standard',
                  'grimme_short',
                  'lyskov_standard',
                  'lyskov_short',
                  'heil17_standard',
                  'heil17_short',
                  'heil18_standard',
                  'heil18_short']

class Dftmrci2(cimethod.Cimethod):
    """Class constructor for DFT/MRCI(2) objects"""
    def __init__(self):        
         # parent attributes
        super().__init__()

        # user defined quanties
        self.truncate        = True
        self.truncate_thresh = 0.9
        self.shift           = 0.005
        self.hamiltonian     = 'heil17_standard'
        self.hparam          = None
        self.ras1            = []
        self.ras2            = []
        self.ras3            = []
        self.nhole1          = 0
        self.nelec3          = 0
        self.autoras         = False
        self.icvs            = []
        self.refiter         = 3
        self.ref_prune       = True
        self.nbuffer         = []
        self.label           = 'Dftmrci2'

        # class variables
        # No. extra ref space roots needed
        self.nextra         = {}
        # Max no. iterations of the ref space refinement 
        self.niter          = 0
        # ref space bitci wave function object
        self.ref_wfn        = None
        # MRCI bitci wave function object
        self.mrci_wfn       = None
        # reference space energies
        self.ref_ener       = None
        # Q-space norm file numbers
        self.qunits         = None
        # Damped strong perturbers file numbers
        self.dspunits       = None
        # dictionary of bitci wfns
        self.bitciwfns      = {}

# Required functions #############################################################

    @timing.timed
    def run(self, scf):
        """ compute the DFT/MRCI(2) eigenpairs for all irreps """

        # set the scf object
        self.set_scf(scf)

        if self.scf.mol is None or self.scf is None:
            sys.exit('ERROR: mol and scf objects not set in dftmrci')

        # write the output logfile header for this run
        output.print_dftmrci2_header(self.label)

        # initialize bitci
        bitci_init.init(self)

        # create the reference space and mrci wfn objects
        self.ref_wfn  = bitciwfn.Bitciwfn()
        self.mrci_wfn = bitciwfn.Bitciwfn()

        # set the number of extra roots to be calculated
        if len(self.nbuffer) == 0:
            # default: 10 roots per irrep
            self.nextra = {'pt2' : [10 for n in range(self.n_irrep())],
                           'max' : [10 for n in range(self.n_irrep())]}
        else:
            # user-specified values
            self.nextra = {'pt2' : np.ndarray.tolist(self.nbuffer),
                           'max' : np.ndarray.tolist(self.nbuffer)}
        
        # generate the initial reference space configurations
        n_ref_conf, ref_conf_units = ref_space.generate(self)
        # set the number of configurations and the scratch file numbers
        self.ref_wfn.set_nconf(n_ref_conf)
        self.ref_wfn.set_confunits(ref_conf_units)

        # Perform the DFT/MRCI(2) iterations, refining the reference
        # space as we go
        self.niter = 0
        for self.niter in range(self.refiter):

            # reference space diagonalisation
            ref_ci_units, ref_ener = ref_diag.diag(self)
            # set the ci files and reference energies
            self.ref_wfn.set_ciunits(ref_ci_units)
            self.ref_ener = ref_ener
            output.print_refdiag_summary(self)

            # optional removal of deadwood from the
            # guess reference space
            if self.ref_prune and self.niter == 0:
                # remove the deadwood
                n_ref_conf = ref_prune.prune(self)
                # set the new no. ref confs
                self.ref_wfn.set_nconf(n_ref_conf)
                # re-diagonalise
                ref_ci_units, ref_ener = ref_diag.diag(self)
                # set the ci files and reference energies
                self.ref_wfn.set_ciunits(ref_ci_units)
                self.ref_ener = ref_ener
                output.print_refdiag_summary(self)

            # generate the MRCI configurations
            n_mrci_conf, mrci_conf_units, mrci_conf_files, dummy = \
                mrci_space.generate(self)

            # set the number of mrci config, and the mrci unit numbers
            # and unit names
            self.mrci_wfn.set_nconf(n_mrci_conf)
            self.mrci_wfn.set_confunits(mrci_conf_units)
            self.mrci_wfn.set_confname(mrci_conf_files)

            # DFT/MRCI(2) calculation
            mrci_ci_units, mrci_ci_files, mrci_ener_sym, q_units, dsp_units, n_conf_new = \
                gvvpt2.diag_heff(self)

            # set the new number of mrci confs if wave function
            # truncation is being used
            if self.truncate:
                self.mrci_wfn.set_nconf(n_conf_new)            
                        
            # set the wfn unit numbers, file names, energies,
            # Q-space info, and damped strong perturber unit numbers
            self.mrci_wfn.set_ciunits(mrci_ci_units)
            self.mrci_wfn.set_ciname(mrci_ci_files)
            self.energies_sym = mrci_ener_sym
            self.qunits       = q_units
            self.dspunits     = dsp_units
            # generate the energies sorted by value, and their
            # corresponding states
            self.order_energies()
                
            # refine the reference space
            min_norm, n_ref_conf, ref_conf_units = \
                    gvvpt2_refine.refine_ref_space(self)
            self.ref_wfn.set_nconf(n_ref_conf)
            self.ref_wfn.set_confunits(ref_conf_units)

            # break if the reference space is converged
            if min_norm > 0.9025 and self.niter > 0:
                print('\n * Reference Space Converged *', flush=True)
                break

        # save the determinant expansions of the wave functions
        # if requested
        if self.save_wf:
            mrci_wf.extract_wf(self)

        # construct density matrices
        dmat_sym = mrci_1rdm.rdm(self)

        # store them in adiabatic energy order
        n_tot = self.n_states()
        (nmo1, nmo2, n_dum) = dmat_sym[0].shape  

        self.dmats = np.zeros((n_tot, nmo1, nmo2), dtype=float)
        for istate in range(n_tot):
            irr, st = self.state_sym(istate)
            self.dmats[istate, :, :] = dmat_sym[irr][:, :, st]

        # Finalise the bitCI library
        bitci_init.finalize()

        # Finalize the bitCI library
        bitci_init.finalize()

        # build the natural orbitals in AO basis by default
        self.build_nos()

        # only print if user-requested
        if self.print_orbitals:
            self.print_nos()

        # determine promotion numbers if ref_state != -1
        if self.ref_state != -1:
            self.print_promotion(self.ref_state)

        # print the moments
        self.print_moments()

        return 
        
    #
    def bitci_ref(self):
        """returns the bitci reference space wfn"""
        return self.ref_wfn

    #
    def bitci_mrci(self):
        """returns the bitci mrci wfn"""
        return self.mrci_wfn

#########################################################################
