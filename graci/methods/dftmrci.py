"""
Module for computing DFT/MRCI energies
"""
import copy as copy
import sys as sys
import numpy as np
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
import graci.interfaces.bitci.mrci_diag as mrci_diag
import graci.interfaces.bitci.mrci_refine as mrci_refine
import graci.interfaces.bitci.mrci_1rdm as mrci_1rdm
import graci.interfaces.bitci.mrci_wf as mrci_wf

class Dftmrci(cimethod.Cimethod):
    """Class constructor for DFT/MRCI object"""
    def __init__(self):
        # parent attributes
        super().__init__()

        # user defined quanties
        self.hamiltonian    = 'heil17_standard'
        self.ras1           = []
        self.ras2           = []
        self.ras3           = []
        self.nhole1         = 0
        self.nelec3         = 0
        self.autoras        = False
        self.icvs           = []
        self.refiter        = 3
        self.ref_prune      = True
        self.prune          = False
        self.prune_thresh   = 0.9
        self.prune_qcorr    = True
        self.prune_extra    = 10
        self.diag_guess     = 'subdiag'
        self.diag_method    = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False
        self.label          = 'Dftmrci'

        # No. extra ref space roots needed
        # for various tasks
        self.nextra         = {}
        # Max no. iterations of the ref space refinement 
        self.niter          = 0
        # ref space bitci wave function object
        self.ref_wfn        = None
        # MRCI bitci wave function object
        self.mrci_wfn       = None
        # reference space energies
        self.ref_ener       = None
        # dictionary of bitci wfns
        self.bitciwfns      = {}


# Required functions #############################################################
    def copy(self):
        """create of deepcopy of self"""
        new = Dftmrci()

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
        """ compute the DFT/MRCI eigenpairs for all irreps """
        
        # set the scf object 
        self.set_scf(scf)

        if self.scf.mol is None or self.scf is None:
            sys.exit('ERROR: mol and scf objects not set in dftmrci')
            
        # write the output logfile header for this run
        if self.verbose:
            output.print_dftmrci_header(self.label)

        # if a guess CI object has been passed, compute the
        # MO overlaps
        if guess is not None:
            self.smo = self.scf.mo_overlaps(guess.scf)
       
        # initialize bitci
        bitci_init.init(self)
        
        # create the reference space and mrci wfn objects
        self.ref_wfn  = bitciwfn.Bitciwfn()
        self.mrci_wfn = bitciwfn.Bitciwfn()

        # determine the no. extra ref space roots needed
        self.nextra = ref_diag.n_extra(self)

        # generate the initial reference space configurations
        if guess is not None:
            n_ref_conf, ref_conf_units, ref_conf_files = \
                ref_space.propagate(self, guess)
        else:
            n_ref_conf, ref_conf_units, ref_conf_files = \
                ref_space.generate(self)
        
        # set the number of configurations and the scratch file numbers
        # and names
        self.ref_wfn.set_nconf(n_ref_conf)
        self.ref_wfn.set_confunits(ref_conf_units)
        self.ref_wfn.set_confname(ref_conf_files)
        
        # Perform the MRCI iterations, refining the reference space
        # as we go
        self.niter = 0
        for self.niter in range(self.refiter):
            
            # reference space diagonalisation
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
            n_mrci_conf, mrci_conf_units, mrci_conf_files, \
                eq_units = mrci_space.generate(self)
            # set the number of mrci config, the mrci unit numbers and
            # unit names, the Q-space energy correction unit numbers,
            # and the damped strong perturber unit numbers
            self.mrci_wfn.set_nconf(n_mrci_conf)
            self.mrci_wfn.set_confunits(mrci_conf_units)
            self.mrci_wfn.set_confname(mrci_conf_files)
            self.mrci_wfn.set_equnits(eq_units)
            
            # MRCI diagonalisation
            mrci_ci_units, mrci_ci_files, mrci_ener_sym = \
                    mrci_diag.diag(self)
            # set the mrci wfn unit numbers, file names, and mrci 
            # energies
            self.mrci_wfn.set_ciunits(mrci_ci_units)
            self.mrci_wfn.set_ciname(mrci_ci_files)
            self.energies_sym = mrci_ener_sym
            # generate the energies sorted by value, and their
            # corresponding states
            self.order_energies()
            
            # refine the reference space
            min_norm, n_ref_conf, ref_conf_units = \
                mrci_refine.refine_ref_space(self)
            self.ref_wfn.set_nconf(n_ref_conf)
            self.ref_wfn.set_confunits(ref_conf_units)

            # break if the reference space is converged
            if min_norm > 0.9025 and self.niter > 0:
                if self.verbose:
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

        # Finalize the bitCI library
        bitci_init.finalize()

        # build the natural orbitals in AO basis by default
        self.build_nos()

        # only print if user-requested
        if self.print_orbitals:
            self.print_nos()

        # determine promotion numbers if ref_state != -1
        if self.ref_state != -1:
            self.print_promotion(self.ref_state-1)

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

