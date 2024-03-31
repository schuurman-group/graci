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
    def __init__(self, ci_obj=None):
        # parent attributes
        super().__init__()

        # user defined quanties
        self.order          = 2
        self.hamiltonian    = 'qe8'
        self.ras1           = []
        self.ras2           = []
        self.ras3           = []
        self.nhole1         = 2
        self.nelec3         = 2
        self.autoras        = False
        self.icvs           = []
        self.refiter        = 3
        self.ref_prune      = True
        self.prune          = False
        self.prune_thresh   = 0.90
        self.prune_qcorr    = True
        self.prune_extra    = 10
        self.regfac         = 0.005
        self.diag_guess     = 'subdiag'
        self.diag_method    = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False

        # MO energy cutoff: MOs above this value excluded
        self.mo_cutoff      = 1.
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
                
        if isinstance(ci_obj, cimethod.Cimethod):
            for name,obj in ci_obj.__dict__.items():
                if hasattr(self, name):
                    # pass a copy of mutable objects
                    if isinstance(obj, (dict, list, np.ndarray)):
                        setattr(self, name, copy.deepcopy(obj))
                    # pass a copy of graci object
                    elif obj.__class__.__name__ in params.valid_objs:
                        setattr(self, name, obj.copy())
                    else:
                        setattr(self, name, obj)

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
    def run(self, scf, guess=None, mo_ints=None):
        """ compute the DFT/MRCI eigenpairs for all irreps """

        # set the scf object
        # currently the guess object is ignored, since propagation of
        # reference space in absense of root-following is not a good
        # idea. Root-following for davidson might be implemented in 
        # the future, though, so for the sake of API consistentcy with
        # DFT/MRCI(2) class, the argument can be passed, but is currently
        # ignored 
        guess      = None
        scf_energy = self.set_scf(scf, ci_guess=guess, mo_ints=mo_ints)

        if scf_energy is None:
            return None

        # set the Hamiltonian
        self.set_hamiltonian()

        if self.scf.mol is None or self.scf is None:
            sys.exit('ERROR: mol and scf objects not set in dftmrci')

        # write the output logfile header for this run
        if self.verbose:
            output.print_dftmrci_header(self.label)

        # write the Cartesian coordinate to the log file
        if self.verbose:
            output.print_coords(self.scf.mol.crds,
                                self.scf.mol.asym)

        # write the Hamiltonian information to the log file
        if self.verbose:
            output.print_hamiltonian(self.hamiltonian)
        
        # if a guess CI object has been passed, compute the
        # MO overlaps
        if guess is not None:
            self.smo = self.scf.mo_overlaps(guess.scf)[:guess.nmo,:self.nmo]

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
            # guess reference space:
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
            mrci_ci_units, mrci_ci_files, mrci_ener_sym, \
                avii_units, avii_files= mrci_diag.diag(self)

            # set the mrci wfn unit numbers, file names, mrci 
            # energies, and spin-couping averaged Hii unit numbers
            self.mrci_wfn.set_ciunits(mrci_ci_units)
            self.mrci_wfn.set_ciname(mrci_ci_files)
            self.energies_sym = mrci_ener_sym
            self.mrci_wfn.set_aviiunits(avii_units)
            self.mrci_wfn.set_aviiname(avii_files)

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

        # Finalize the bitCI library
        bitci_init.finalize()

        # store the density matrices in adiabatic energy order
        n_tot = self.n_states()
        indx  = self.irreps_nonzero()[0]
        (nmo1, nmo2, n_dum) = dmat_sym[indx].shape  

        self.dmats['adiabatic'] = np.zeros((n_tot, nmo1, nmo2), dtype=float)
        for istate in range(n_tot):
            irr, st = self.state_sym(istate)
            self.dmats['adiabatic'][istate, :, :] = dmat_sym[irr][:, :, st]

        # build the natural orbitals in AO basis by default
        self.build_nos()

        # only print if user-requested
        if self.print_orbitals:
            self.print_nos()

        # determine promotion numbers if ref_state != -1
        if self.ref_state > 0:
            self.print_promotion(self.ref_state-1)

        # print the moments
        self.print_moments()

        return True
    
    #
    def bitci_ref(self):
        """returns the bitci reference space wfn"""
        return self.ref_wfn

    #
    def bitci_mrci(self):
        """returns the bitci mrci wfn"""
        return self.mrci_wfn

