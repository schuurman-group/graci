"""
Module for computing DFT/MRCI(2) states and energies
"""
import sys as sys
import numpy as np
import copy as copy
import graci.utils.timing as timing
import graci.methods.cimethod as cimethod
import graci.core.params as params
import graci.core.bitciwfn as bitciwfn
import graci.io.output as output
import graci.interfaces.bitci.bitci_init as bitci_init
import graci.interfaces.bitci.ref_space as ref_space
import graci.interfaces.bitci.ref_diag as ref_diag
import graci.interfaces.bitci.ref_prune as ref_prune
import graci.interfaces.bitci.mrci_space as mrci_space
import graci.interfaces.bitci.gvvpt2 as gvvpt2
import graci.interfaces.bitci.gvvpt2_diab as gvvpt2_diab
import graci.interfaces.bitci.gvvpt2_refine as gvvpt2_refine
import graci.interfaces.bitci.gvvpt2_truncate as gvvpt2_truncate
import graci.interfaces.bitci.mrci_1rdm as mrci_1rdm
import graci.interfaces.bitci.mrci_wf as mrci_wf
import graci.interfaces.overlap.bdd as bdd

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
        self.regularizer     = 'isa'
        self.regfac          = None
        self.hamiltonian     = 'heil17_standard'
        self.ras1            = []
        self.ras2            = []
        self.ras3            = []
        self.nhole1          = 2
        self.nelec3          = 2
        self.autoras         = False
        self.icvs            = []
        self.refiter         = 3
        self.ref_prune       = True
        self.diabatic        = False
        self.adt_type        = 'bdd'
        self.nbuffer         = []
        self.refsel          = 'dynamic'

        # class variables
        # MO energy cutoff: MOs above this value excluded
        self.mo_cutoff       = 1.
        # allowed regularizers
        self.allowed_regularizer = ['isa', 'sigmap']
        # allowed ref conf selection algorithms
        self.allowed_refsel      = ['static', 'dynamic']
        # default regularization factors
        self.default_regfac      = {'isa' : 0.005, 'sigmap' : 0.025}
        # no. extra ref space roots needed
        self.nextra              = {}
        # max no. iterations of the ref space refinement 
        self.niter               = 0
        # ref space bitci wave function object
        self.ref_wfn             = None
        # MRCI bitci wave function object
        self.mrci_wfn            = None
        # reference space energies
        self.ref_ener            = None
        # Q-space norm file numbers
        self.qunits              = None
        # damped strong perturbers file numbers
        self.dspunits            = None
        # A-vector file numbers
        self.Aunits              = None
        # allowed ADT types
        self.allowed_adt_type    = ['bdd', 'qdpt']
        # dictionary of bitci wfns
        self.bitciwfns           = {}
        
        
# Required functions #############################################################
    def copy(self):
        """create of deepcopy of self"""
        new = Dftmrci2()

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))

        return new

    @timing.timed
    def run(self, scf, guess):
        """ compute the DFT/MRCI(2) eigenpairs for all irreps """

        # set the scf object
        self.set_scf(scf)

        # sanity check on the input variables
        self.sanity_check(scf, guess)

        # if a guess CI object has been passed, compute the
        # MO overlaps
        if guess is not None:
            self.smo = self.scf.mo_overlaps(guess.scf)[:guess.nmo,:self.nmo]

        # set the regularization factor to something sensible if it has
        # not been given
        if self.regfac is None:
            self.regfac = self.default_regfac[self.regularizer]
            
        # write the output log file header for this run
        if self.verbose:
            output.print_dftmrci2_header(self.label)
            
        # write the Cartesian coordinate to the log file
        if self.verbose:
            output.print_coords(self.scf.mol.crds,
                                self.scf.mol.asym)
        
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
        if guess is not None:
            n_ref_conf, ref_conf_units, ref_conf_files = \
                ref_space.propagate(self, guess)
        else:
            n_ref_conf, ref_conf_units, ref_conf_files = \
                ref_space.generate(self)
        
        # set the number of configurations and scratch file numbers
        # and names
        self.ref_wfn.set_nconf(n_ref_conf)
        self.ref_wfn.set_confunits(ref_conf_units)
        self.ref_wfn.set_confname(ref_conf_files)
        
        # Perform the DFT/MRCI(2) iterations, refining the reference
        # space as we go
        self.niter = 0
        for self.niter in range(self.refiter):

            # reference space diagonalisation
            if self.diabatic:
                ref_ci_units, ref_ener = ref_diag.diag_follow(self, guess)
            else:
                ref_ci_units, ref_ener = ref_diag.diag(self)
            
            # set the ci files and reference energies
            self.ref_wfn.set_ciunits(ref_ci_units)
            self.ref_ener = ref_ener

            if self.verbose:
                output.print_refdiag_summary(self)

            # optional removal of deadwood from the
            # guess reference space
            if self.ref_prune and self.niter == 0 \
               and guess is None:
                # remove the deadwood
                n_ref_conf = ref_prune.prune(self)
                # set the new no. ref confs
                self.ref_wfn.set_nconf(n_ref_conf)
                # re-diagonalise
                ref_ci_units, ref_ener = ref_diag.diag(self)
                # set the ci files and reference energies
                self.ref_wfn.set_ciunits(ref_ci_units)
                self.ref_ener = ref_ener
                if self.verbose:
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
            if self.diabatic:
                mrci_ci_units, mrci_ci_files, mrci_ener_sym, \
                    q_units, dsp_units, A_units \
                    = gvvpt2.diag_heff_follow(self, guess)
            else:
                mrci_ci_units, mrci_ci_files, mrci_ener_sym, \
                    q_units, dsp_units = gvvpt2.diag_heff(self)
                
            # set the wfn unit numbers, file names, energies,
            # Q-space info, and damped strong perturber unit numbers
            self.mrci_wfn.set_ciunits(mrci_ci_units)
            self.mrci_wfn.set_ciname(mrci_ci_files)
            self.energies_sym = mrci_ener_sym
            self.qunits       = q_units
            self.dspunits     = dsp_units
            if self.diabatic:
                self.Aunits   = A_units
                
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
                if self.verbose:
                    print('\n * Reference Space Converged *', flush=True)
                break

        # save the determinant expansions of the wave functions
        # if requested
        if self.save_wf or self.diabatic:
            mrci_wf.extract_wf(self)

        # diabatisation
        if self.diabatic:

            if self.adt_type == 'bdd':
                # block diagonalisation diabatisation
                adt_matrices = bdd.adt(guess, self)
                self.adt     = adt_matrices
                self.diabatize()

            elif self.adt_type == 'qdpt':
                # QDPT diabatisation
                diabpots, ciunits, cinames, \
                    confunits, confnames, nconfs = \
                        gvvpt2_diab.diabpot(guess, self)
                self.diabpot = diabpots
                self.mrci_wfn.set_ciunits(ciunits, rep='diabatic')
                self.mrci_wfn.set_ciname(cinames, rep='diabatic')
                self.mrci_wfn.set_confunits(confunits, rep='diabatic')
                self.mrci_wfn.set_confname(confnames, rep='diabatic')
                self.mrci_wfn.set_nconf(nconfs, rep='diabatic')
                mrci_wf.extract_wf(self, rep='diabatic')
                
            # output the diabatic potentials
            if self.verbose:
                nroots = [self.n_states_sym(irr)
                          for irr in range(self.n_irrep())]
                output.print_diabpot(self.diabpot, nroots,
                                     self.n_irrep(),
                                     self.scf.mol.irreplbl)

        # removal of deadwood configurations
        # (this must be called _after_ the CSF-to-det
        #  transformation)
        if self.truncate:
            for rep in ['adiabatic', 'diabatic']:
                if self.mrci_wfn.nconf[rep] is not None:
                    n_conf_new = gvvpt2_truncate.truncate_wf(self, rep)
                    self.mrci_wfn.set_nconf(n_conf_new, rep)

        # construct density matrices
        self.get_dmats()
        
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

        # Finalize the bitCI library
        bitci_init.finalize()
        
        return 

    #
    def sanity_check(self, scf, guess):
        """sanity check on user specified variables"""

        # mol and scf objects
        if self.scf.mol is None or self.scf is None:
            sys.exit('ERROR: mol and scf objects not set in dftmrci2')

        # diabatic guess object
        if self.diabatic and guess is None:
            sys.exit('ERROR: diabatisation requires a guess dftmrci2 object')

        # ref conf selection algorithm
        if self.refsel not in self.allowed_refsel:
            sys.exit('\n ERROR: unrecognised value of refsel: '
                     +self.refsel + '\n allowed values: '
                     +str(self.allowed_refsel))

        # GVVPT2 regularizer
        if self.regularizer not in self.allowed_regularizer:
            sys.exit('\n ERROR: unrecognised value of regularizer: '
                     +self.regularizer + '\n allowed values: '
                     +str(self.allowed_regularizer))

        # ADT type
        if self.diabatic and self.adt_type not in self.allowed_adt_type:
            sys.exit('\n ERROR: unrecognised adt_type: '
                     +self.adt_type + '\n allowed values: '
                     +str(self.allowed_adt_type))
            
        return

    #
    def get_dmats(self):
        """Constructs the adiabatic and, optionally, diabatic
           density matrices"""

        # diabatic WFs only currently available for QDPT diabatisation
        if self.diabatic and self.adt_type == 'qdpt':
            do_diabatic = True
        else:
            do_diabatic = False
        
        # calculate the 1-RDMs
        dmat_sym = {'adiabatic' : None, 'diabatic' : None}
        dmat_sym['adiabatic'] = mrci_1rdm.rdm(self)
        if do_diabatic:
            dmat_sym['diabatic'] = mrci_1rdm.rdm(self, rep='diabatic')        
        
        # store them in adiabatic energy order
        n_tot = self.n_states()
        for rep in dmat_sym.keys():
            if dmat_sym[rep] is not None:
                (nmo1, nmo2, n_dum) = dmat_sym[rep][0].shape  
                self.dmats[rep] = np.zeros((n_tot, nmo1, nmo2), dtype=float)
                for istate in range(n_tot):
                    irr, st = self.state_sym(istate)
                    self.dmats[rep][istate, :, :] = dmat_sym[rep][irr][:, :, st]
                    
        return
        
    #
    def bitci_ref(self):
        """returns the bitci reference space wfn"""
        return self.ref_wfn

    #
    def bitci_mrci(self):
        """returns the bitci mrci wfn"""
        return self.mrci_wfn
    
