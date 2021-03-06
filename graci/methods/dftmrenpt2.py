"""
Module for computing DFT/MR-ENPT2 energies
"""
import os as os
import sys as sys
import numpy as np
import copy as copy
import graci.utils.timing as timing
import graci.methods.cimethod as cimethod
import graci.core.bitciwfn as bitciwfn
import graci.bitcitools.bitci_init as bitci_init
import graci.bitcitools.ref_space as ref_space
import graci.bitcitools.ref_diag as ref_diag
import graci.bitcitools.ref_prune as ref_prune
import graci.bitcitools.mrci_space as mrci_space
import graci.bitcitools.mrenpt2 as mrenpt2
import graci.bitcitools.mrenpt2_refine as mrenpt2_refine
import graci.bitcitools.mrci_1rdm as mrci_1rdm
import graci.bitcitools.mrci_wf as mrci_wf
import graci.io.output as output
import graci.properties.moments as moments

class Dftmrenpt2(cimethod.Cimethod):
    """Class constructor for DFT/MR-ENPT2 objects"""
    def __init__(self):        
         # parent attributes
        super().__init__()

        # user defined quanties
        self.truncate        = True
        self.truncate_thresh = 0.9
        self.shift           = 0.005
        self.multistate      = False
        self.hamiltonian     = 'heil17_standard'
        self.ras1            = []
        self.ras2            = []
        self.ras3            = []
        self.nhole1          = 0
        self.nelec3          = 0
        self.autoras         = False
        self.icvs            = []
        self.refiter         = 3
        self.ref_prune       = True
        self.save_wf         = False
        self.label           = 'Dftmrenpt2'

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
        """ compute the DFT/MR-ENPT2 eigenpairs for all irreps """

        # set the scf object
        self.set_scf(scf)

        if self.scf.mol is None or self.scf is None:
            sys.exit('ERROR: mol and scf objects not set in dftmrci')

        # write the output logfile header for this run
        output.print_dftmrenpt2_header(self.label)

        # initialize bitci
        bitci_init.init(self)

        # create the reference space and mrci wfn objects
        self.ref_wfn  = bitciwfn.Bitciwfn()
        self.mrci_wfn = bitciwfn.Bitciwfn()

        # set the number of extra roots to be calculated
        # (hard-wired for now)
        self.nextra = {'enpt2' : [10 for n in range(self.n_irrep())],
                       'max'   : [10 for n in range(self.n_irrep())]}
        
        # generate the initial reference space configurations
        n_ref_conf, ref_conf_units = ref_space.generate(self)
        # set the number of configurations and the scratch file numbers
        self.ref_wfn.set_nconf(n_ref_conf)
        self.ref_wfn.set_confunits(ref_conf_units)

        # Perform the MR-ENPT2 iterations, refining the reference
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

            # MR-ENPT2 calculation
            mrci_ci_units, mrci_ci_files, mrci_ener_sym, q_units, dsp_units, n_conf_new = \
                mrenpt2.corrections(self)

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
                    mrenpt2_refine.refine_ref_space(self)
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

        # build the natural orbitals
        self.build_nos()

        # Finalise the bitCI library
        bitci_init.finalize()

        # print orbitals if requested
        if self.print_orbitals:
            self.export_orbitals(orb_format='molden')

        # build the list and symmetries for subsequent printing
        states = [i for i in range(n_tot)]
        syms   = [self.scf.mol.irreplbl[self.state_sym(i)[0]]
                  for i in range(n_tot)]
            
        # also compute attachment and detachment numbers
        # (relative to ground state)
        ndo, ndo_wt = self.build_ndos(0, basis='mo')
        pd, pa      = self.promotion_numbers(ndo, ndo_wt)
        output.print_promotion(0, states, syms, pd, pa)

        # we'll also compute 1-electron properties by
        # default.
        momts = moments.Moments(self.scf.mol, self.natocc, self.natorb_ao)
        momts.run()
        states = [i for i in range(n_tot)]
        syms   = [self.scf.mol.irreplbl[self.state_sym(i)[0]]
                  for i in range(n_tot)]
        output.print_moments(states, syms, momts)

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
