"""
Module for computing DFT/MRCI energies
"""
import os as os
import sys as sys
import importlib
import graci.io.convert as convert
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.core.libs as libs
import graci.core.bitciwfn as bitciwfn
import graci.citools.ref_space as ref_space
import graci.citools.ref_diag as ref_diag
import graci.citools.ref_prune as ref_prune
import graci.citools.mrci_space as mrci_space
import graci.citools.mrci_diag as mrci_diag
import graci.citools.mrci_refine as mrci_refine
import graci.citools.mrci_1rdm as mrci_1rdm
import graci.citools.mrci_wf as mrci_wf
import graci.io.output as output
import graci.properties.moments as moments

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

class Dftmrci:
    """Class constructor for DFT/MRCI object"""
    def __init__(self):
        # user defined quanties
        self.nstates        = []
        self.hamiltonian    = 'canonical'
        self.ras1           = []
        self.ras2           = []
        self.ras3           = []
        self.nhole1         = 0
        self.nelec3         = 0
        self.autoras        = False
        self.icvs           = []
        self.refiter        = 3
        self.ref_prune      = True
        self.prune          = 'off'
        self.prune_thresh   = 1.
        self.prune_qcorr    = True
        self.prune_extra    = 10
        self.diag_guess     = 'subdiag'
        self.diag_method    = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False
        self.print_orbitals = False
        self.save_wf        = False
        self.label          = 'Dftmrci'

        # class variables
        # KS SCF object
        self.scf            = None 
        # Pruning variables
        self.pmrci          = False
        self.prune_dict     = {'tight'  : 0.9900,
                               'normal' : 0.9500,
                               'loose'  : 0.9000}
        # No. extra ref space roots needed
        # for various tasks
        self.nextra         = {}
        # Max no. iterations of the ref space refinement 
        self.niter          = 0
        # ref space bitci wave function object
        self.ref_wfn        = None
        # MRCI bitci wave function object
        self.mrci_wfn       = None
        # print quadrupoles
        self.print_quad     = False
        # reference space energies
        self.ref_ener       = None
        # ci energies, by adiabatic state
        self.mrci_ener      = None
        # ci energies by irrep and root index
        self.mrci_ener_sym  = None
        # irrep,state pairs sorted by energy
        self.sym_sorted     = None
        # Density matrices (list of numpy arrays, one per irrep,
        # shape: (nirr,nmo,nmo,nstates))
        self.dmat           = None
        # natural orbital occupations (shape = (nirr, nmo, nstates)
        self.natocc         = None
        # natural orbitals, MO basis (same shape as the density matrices)
        self.natorb_mo      = None
        # natural orbitals, AO basis (same shape as the density matrices)
        self.natorb_ao      = None
        # symmetry of the natural orbitals
        self.natorb_sym     = None
        # dictionary of bitci wfns
        self.bitciwfns      = {}


# Required functions #############################################################

    def set_scf(self, scf):
        """set the scf object for the dftmrci class object"""
        self.scf = scf
        return

    def scf_exists(self):
        """return true if scf object is not None"""
        try:
            return type(self.scf).__name__ is 'Scf'
        except:
            return False

    @timing.timed
    def run(self):
        """ compute the DFT/MRCI eigenpairs for all irreps """

        if self.scf.mol is None or self.scf is None:
            sys.exit('ERROR: mol and scf objects not set in dftmrci')

        # write the output logfile header for this run
        output.print_dftmrci_header(self.label)

        # initialize bitci
        libs.init_bitci(self)
        
        # create the reference space and mrci wfn objects
        self.ref_wfn  = bitciwfn.Bitciwfn()
        self.mrci_wfn = bitciwfn.Bitciwfn()

        # determine the no. extra ref space roots needed
        self.nextra = ref_diag.n_extra(self)

        # set the pruning variables
        self.set_prune_vars()
        
        # generate the initial reference space configurations
        n_ref_conf, ref_conf_units = ref_space.generate(self)
        # set the number of configurations and the scratch file numbers
        self.ref_wfn.set_nconf(n_ref_conf)
        self.ref_wfn.set_confunits(ref_conf_units)

        # Perform the MRCI iterations, refining the reference space
        # as we go
        self.niter = 0
        for self.niter in range(self.refiter):
            
            # reference space diagonalisation
            ref_ci_units, ref_ener = ref_diag.diag(self)
            # set the ci files and reference energies
            self.ref_wfn.set_ciunits(ref_ci_units)
            self.ref_ener = ref_ener
            output.print_refdiag_summary(self)


            sys.exit('Remove this sys.exit()')
            
            
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
            n_mrci_conf, mrci_conf_units, mrci_conf_files, eq_units = \
                    mrci_space.generate(self)
            # set the number of mrci config, the mrci unit numbers and
            # unit names, and the Q-space energy correction unit numbers
            self.mrci_wfn.set_nconf(n_mrci_conf)
            self.mrci_wfn.set_confunits(mrci_conf_units)
            self.mrci_wfn.set_confname(mrci_conf_files)
            self.mrci_wfn.set_equnits(eq_units)
            
            # MRCI diagonalisation
            mrci_ci_units, mrci_ci_files, mrci_ener_sym = \
                    mrci_diag.diag(self)
            # set the mrci wfn unit numbers, file names and mrci 
            # energies
            self.mrci_wfn.set_ciunits(mrci_ci_units)
            self.mrci_wfn.set_ciname(mrci_ci_files)
            self.mrci_ener_sym = mrci_ener_sym
            # generate the energies sorted by value, and their
            # corresponding states
            self.order_ener()

            # refine the reference space
            min_norm, n_ref_conf, ref_conf_units = \
                    mrci_refine.refine_ref_space(self)
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
        n_tot = self.n_state()
        (nmo1, nmo2, n_dum) = dmat_sym[0].shape  

        self.dmat = np.zeros((n_tot, nmo1, nmo2), dtype=float)
        for istate in range(n_tot):
            irr, st = self.state_sym(istate)
            self.dmat[istate, :, :] = dmat_sym[irr][:, :, st]

        # build the natural orbitals
        self.build_nos()

        # Finalise the bitCI library
        libs.finalise_bitci()

        # print orbitals if requested
        if self.print_orbitals:
            self.export_orbitals(orb_format='molden')

        # we'll also compute 1-electron properties by
        # default.
        momts = moments.Moments(self.scf.mol, self.natocc, self.natorb_ao)
        momts.run()
 
        output.print_moments_header()
        for ist in range(n_tot):
            [irr, st_sym] = self.state_sym(ist)
            output.print_moments(
                        ist,
                        self.scf.mol.irreplbl[irr],
                        momts.dipole(ist),
                        momts.second_moment(ist))
  
            if self.print_quad:
                output.print_quad(
                        ist,
                        self.scf.mol.irreplbl[irr],
                        momts.quadrupole(ist))

        return

    #
    def set_prune_vars(self):
        """Handles the setting of the MRCI pruning logical flag
        and threshold variables"""

        # Pruning not requested
        if self.prune == 'off' and self.prune_thresh == 1.:
            self.pmrci        = False
            self.prune_thresh = 1.
            return

        # Pruning reqested: return True and the requested threshold
        # Note that the 'prune_thresh' keyword takes precedent over
        # the 'prune' keyword
        if self.prune_thresh != 1.:
            self.pmrci = True
            return
        else:
            self.pmrci        = True
            self.prune_thresh = self.prune_dict[self.prune]
            return
    
    # 
    def n_irrep(self):
        """return the number of irreps"""
        return len(self.nstates)
 
    #
    def n_state(self):
        """total number of states"""

        return sum(self.n_state_sym())

    #
    def n_state_sym(self, irrep=None):
        """number of states to compute"""
 
        if irrep is None:
            return self.nstates

        if irrep < len(self.nstates):
            return self.nstates[irrep]
        else:
            print("irrep > nirrep, irrep = "+str(irrep))
            sys.exit()

    #
    def state_index(self, irrep, state):
        """return adiabatic state label for a given irrep and  root"""

        if [irrep, state] in self.sym_sorted:
            return self.sym_sorted.index([irrep, state])
        else:
            return None

    #
    def state_sym(self, n):
        """return the irrep and state index of the root corresponding
           to the energy of the nth root"""

        if n < len(self.sym_sorted):
            return self.sym_sorted[n]
        else:
            return None

    #
    def energy(self, istate):
        """return the energy of state 'state'"""
        return self.mrci_ener[istate]

    #
    def energy_sym(self, irrep, state):
        """return the energy of state by irrep, st index"""
        return self.mrci_ener_sym[irrep, state]

    #
    def rdm1(self, istate):
        """return the density matrix for the state istate"""

        if self.dmat is not None:
            return self.dmat[istate, :, :]
        else:
            print("rdm1 called in method=dftmrci, "+
                  "but density does not exist")
            return self.dmat

    #
    def rdm1_sym(self, irrep, state):
        """return the density matrix for the state in 
        the array 'states' for the irrep 'irr'"""

        istate = self.state_index(irrep, state)

        if self.dmat is not None:
            return self.dmat[istate, :, :]
        else:
            print("rdm1 called in method=dftmrci, "+
                  "but density does not exist")
            return self.dmat

    #
    @timing.timed
    def build_nos(self):
        """print the natural orbitals for irrep 
           'irr' and state 'state'"""

        # if the 1RDMs don't exist, then we have a problem
        if self.dmat is None:
            print("can't return natural orbitals: "+
                  "1RDMs don't exist")
            return None
            
        # check that istate is less than the total number of states:
        (n_tot, nmo1, nmo2) = self.dmat.shape
        self.natocc     = np.zeros((n_tot, nmo2), dtype=float)
        self.natorb_mo  = np.zeros((n_tot, nmo1, nmo2), dtype=float)
        self.natorb_ao  = np.zeros((n_tot, nmo1, nmo2), dtype=float)
        self.natorb_sym = np.zeros((n_tot, nmo2), dtype=int)

        # by default, compute the natural orbitals for all states
        # requested
        for ist in range(n_tot):
            rdm        = self.rdm1(ist)
            occ,nos_mo = np.linalg.eigh(rdm)
            nos_ao     = np.matmul(self.scf.orbs, nos_mo)

            sort_ordr = self.order_orbs(occ)
            occ       = occ[sort_ordr]
            nos_mo    = nos_mo[:,sort_ordr]
            nos_ao    = nos_ao[:,sort_ordr]
            syms      = self.get_orb_sym(nos_ao)

            self.natocc[ist, :]        = occ
            self.natorb_mo[ist, :, :]  = nos_mo
            self.natorb_ao[ist, :, :]  = nos_ao
            self.natorb_sym[ist, :]    = syms

        return
  
    # 
    def natural_orb(self, istate, basis='ao'):
        """return natural orbitals and occupations for state 'istate' """

        if basis == 'ao' and self.natorb_ao is None:
            return None

        if basis == 'mo' and self.natorb_mo is None:
            return None

        (nst, nmo) = self.natocc.shape
        if istate >= nst:
            return None

        # always return the contents of the natorb/natocc arrays in wfn
        occ = self.natocc[istate, :]
        if basis == 'ao':
            nos = self.natorb_ao[istate, :, :]
        else:
            nos = self.natorb_mo[istate, :, :]

        return occ, nos

    #
    def natural_sym(self, istate):
        """return the symmetry of the natural orbitals"""

        # if natural orbitals haven't been generated, generate them now
        if self.natorb_sym is None:
            occ, nos = self.natural_orb(istate)

        return self.natorb_sym[istate, :]

    #
    def slater_dets(self, state):
        """ return the slater determinant list for state 'state'"""

        return

    #
    def csfs(self, state):
        """ return the CSF list for state 'state'"""

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
    #
    def order_ener(self):
        """orders the mrci_ener_sym array to create mrci_ener, a 1D array
           of energies in ascending order, and a corresponding sym_sorted
           array of [irrep, st] pairs for each energy in mrci_sorted."""

        if self.mrci_ener_sym is None:
            sys.exit('cannot sort dftmrci.mrci_ener')

        nirr    = self.n_irrep()
        istate  = np.zeros((nirr), dtype=int)
        n_tot   = self.n_state()

        self.mrci_ener  = np.zeros((n_tot), dtype=float)
        self.sym_sorted = []
        # mrci_ener_sym is an irrep x maxroots array, with trailing values 
        # of 'zero'. 
        mrci_vals         = np.pad(self.mrci_ener_sym, ((0,0),(0,1)), 
                                   'constant', 
                                   constant_values=((0,0),(0,0)))
       
        n_srt = 0
        while n_srt < n_tot:
            eners = np.array([mrci_vals[irr, istate[irr]] 
                              for irr in range(nirr)])
            iirr  = np.argsort(eners)[0]

            self.mrci_ener[n_srt] = mrci_vals[iirr, istate[iirr]]
            self.sym_sorted.append([iirr, istate[iirr]])
            # increment the number of elements in the sorted array,
            # and the irrep state index of irrep just chosen
            n_srt         += 1
            istate[iirr] += 1

        return

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
    def get_orb_sym(self, orbs):
        """determine the symmetry of each orb in orbs. Right now we just use
           the symmetry of the scf/ks orbitals to determine that. We can 
           get more elaborate later if need be"""

        nmo      = len(orbs[0,:])
        sym_indx = np.zeros((nmo), dtype=int)

        # take the overlap of the orb with each KS orbital. We will
        # assign the symmetry of this orbital to be the same as the 
        # symmetry of the KS orbital with maximum overlap
        for orb in range(nmo):
            olap = [abs(np.dot(orbs[:,orb], self.scf.orbs[:,i])) 
                    for i in range(nmo)]
            maxo = olap.index(max(olap))
            sym_indx[orb] = self.scf.orb_sym[maxo]
 
        # make sure the symmetry orbital counts match up
        oval, ocnts = np.unique(sym_indx, return_counts=True)
        val,  cnts  = np.unique(self.scf.orb_sym, return_counts=True)

        #if not (ocnts==cnts).all():
        #    sys.exit('ERROR: cannot assign symmetry of natural orbs')

        return sym_indx

    # 
    def export_orbitals(self, orb_format='molden', orb_dir=True):
        """export natural orbitals for all states to file"""

        for ist in range(self.n_state()):
            self.export_orbitals_state(ist, orb_format=orb_format, 
                                                   orb_dir=orb_dir)
        return

    #
    def export_orbitals_state(self, state, orb_format='molden', 
                                         orb_dir=True, cart=None):
        """export orbitals of a single state to various file formats"""

        if state >= self.n_state():
            return

        [irr, st_sym] = self.state_sym(state)
        occ, orb      = self.natural_orb(state)
        syms          = self.natural_sym(state)
        sym_lbl       = [self.scf.mol.irreplbl[syms[i]]
                        for i in range(len(syms))]
        fname = 'nos.'+str(state+1)+ \
                '_'+str(self.scf.mol.irreplbl[irr].lower())+ \
                '_'+str(orb_format).lower()+"_"+str(self.label)

        if orb_dir:
            fname = 'orbs/'+fname
            # if dir_name directory doesn't exist, create it
            if not os.path.exists('orbs'):
                os.mkdir('orbs')

        # if a file called fname exists, delete it
        if os.path.isfile(fname):
            os.remove(fname)

        # import the appropriate library for the file_format
        if orb_format in output.orb_formats:
            orbtype = importlib.import_module('graci.io.'+orb_format)
        else:
            print('orbital format type=' + orb_format +
                                        ' not found. exiting...')
            sys.exit(1)

        orbtype.write_orbitals(fname, self.scf.mol, orb, 
                               occ=occ, sym_lbl=sym_lbl, cart=None)

        return

