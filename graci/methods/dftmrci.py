"""
Module for computing DFT/MRCI energies
"""
import sys as sys
import graci.io.convert as convert
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.core.bitciwfn as bitciwfn
import graci.citools.ref_space as ref_space
import graci.citools.ref_diag as ref_diag
import graci.citools.ref_prune as ref_prune
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
        self.prune_extra    = 10
        self.diag_method    = 'gendav'
        self.diag_tol       = 0.0001
        self.diag_iter      = 50
        self.diag_blocksize = []
        self.diag_deflate   = False
        self.label          = 'dftmrci'

        # class variables
        # KS SCF object
        self.scf            = None 
        # molecule object
        self.mol            = None
        # Pruning thresholds
        self.prune_thresh   = {'tight'  : 0.9900,
                               'normal' : 0.9500,
                               'loose'  : 0.9000}
        
        self.niter          = 0
        self.ref_wfn        = None
        self.mrci_wfn       = None
        self.print_quad     = False 
        # reference space energies
        self.ref_ener       = None
        # ci energies, by irrep and root index
        self.mrci_ener      = None
        # ci energies, sorted by value
        self.mrci_sorted    = None
        # corresponding states of sorted energy list
        self.state_sorted   = None
        # Density matrices (list of numpy arrays, one per irrep,
        # shape: (nirr,nmo,nmo,nstates))
        self.dmat           = None
        # natural orbital occupations (shape = (nirr, nmo, nstates)
        self.natocc         = None
        # natural orbitals, MO basis (same shape as the density matrices)
        self.natorb_mo      = None
        # natural orbitals, AO basis (same shape as the density matrices)
        self.natorb_ao      = None
        # dictionary of bitci wfns
        self.bitciwfns      = {}


# Required functions #############################################################

    def set_mol(self, mol):
        """set the mol object for the dftmrci class object"""
        self.mol = mol
        return

    def set_scf(self, scf):
        """set the scf object for the dftmrci class object"""
        self.scf = scf
        return

    def name(self):
        """ return the name of the class object as a string"""
        return 'dftmrci'

    def run(self):
        """ compute the DFT/MRCI eigenpairs for all irreps """

        if self.mol is None or self.scf is None:
            sys.exit('ERROR: mol and scf objects not set in dftmrci')

        # run the KS-DFT computation 
        self.scf.set_mol(self.mol)
        self.scf.run()

        # initialize bitci
        libs.init_bitci(self)
        
        # create the reference space and mrci wfn objects
        self.ref_wfn  = bitciwfn.Bitciwfn()
        self.mrci_wfn = bitciwfn.Bitciwfn()

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

            # if this is the first iteration, then remove any
            # deadwood from the reference space
            if self.niter == 0:
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
            mrci_ci_units, mrci_ci_files, mrci_ener = \
                    mrci_diag.diag(self)
            # set the mrci wfn unit numbers, file names and mrci 
            # energies
            self.mrci_wfn.set_ciunits(mrci_ci_units)
            self.mrci_wfn.set_ciname(mrci_ci_files)
            self.mrci_ener = mrci_ener
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

        # Compute the 1-RDMs for all states
        self.dmat = mrci_1rdm.rdm(self)

        # Finalise the bitCI library
        libs.finalise_bitci()

        # by default, print out the natural orbitals for
        # each state
        nirr = len(self.nstates)

        for irr in range(nirr):
            for st in range(self.nstates[irr]):
                occ, orb = self.natural_orbs(irr, st)
                syms     = self.natural_sym(irr, st)
                sym_lbl  = [self.mol.irreplbl[syms[i]] 
                            for i in range(len(syms))]
                fname = 'nos.'+str(self.state_index(irr,st)+1)+ \
                        '_'+str(self.mol.irreplbl[irr].lower())
                output.print_nos_molden(fname, self.mol, orb, occ, 
                                                       sym=sym_lbl)

        # we'll also compute 1-electron properties by
        # default.
        momts = moments.Moments(self.mol, self.scf, self)
        momts.run()
        for irr in range(nirr):
            output.print_moments_header(self.mol.irreplbl[irr])
            for st in range(self.nstates[irr]):
                output.print_moments(
                        self.mol.irreplbl[irr], 
                        st, 
                        momts.dipole(irr,st), 
                        momts.second_moment(irr,st))
                if self.print_quad:
                    output.print_quad(
                            self.mol.irreplbl[irr], 
                            st, 
                            momts.quadrupole(irr,st))
        return 

    # 
    def n_irrep(self):
        """return the number of irreps"""
        return len(self.nstates)
 
    #
    def n_states(self, irrep=None):
        """number of states to compute"""
 
        if irrep is None:
            return self.nstates

        if irrep < len(self.nstates):
            return self.nstates[irrep]
        else:
            print("irrep > nirrep, irrep = "+str(irrep))
            sys.exit()

    #
    def energy(self, irrep, state):
        """return the energy of state 'state'"""

        return self.mrci_ener[irrep, state]

    #
    def state_index(self, irrep, state):
        """return adiabatic state label for a given irrep and  root"""

        if [irrep, state] in self.state_sorted:
            return self.state_sorted.index([irrep, state])
        else:
            return None

    #
    def energy_n(self, n):
        """return the energy of the nth root"""

        return self.mrci_sorted[n]

    #
    def state_n(self, n):
        """return the irrep and state index of the root corresponding
           to the energy of the nth root"""

        return self.state_sorted[n]

    #
    def rdm1(self, irrep, state):
        """return the density matrix for the state in 
        the array 'states' for the irrep 'irr'"""

        if self.dmat is not None:
            return self.dmat[irrep][ :, :, state]
        else:
            print("rdm1 called in method=dftmrci, "+
                  "but density does not exist")
            return self.dmat

    #
    def natural_orbs(self, irrep, state, basis='ao'):
        """print the natural orbitals for irrep 
           'irr' and state 'state'"""

        # if natural orbitals don't exist, create them --
        # we might as well create all of them, unless a truly
        # crazy of number of states have been computed
        if (self.natorb_mo is None and self.natorb_ao is None):
            # if the 1RDMs don't exist, then we have a problem
            if self.dmat is None:
                print("can't return natural orbitals: "+
                      "1RDMs don't exist")
                return None
            
            # this is hacky
            nd = len(self.dmat)
            (nmo1, nmo2, n_max) = self.dmat[0].shape
            nirr = len(self.nstates)

            # sanity check some things
            if irrep >= nd or state >= n_max:
                print("irrep/state = "+str(irr)+","+str(state)+
                      " requested, but dmat dimensions are "+
                      str(nd)+","+str(n_max))
                return None
            if nirr != nd:
                print("diagreement on number of irreps in dftmrci")
                return None

            self.natocc     = np.zeros((nirr, n_max, nmo1),
                                             dtype=float)
            self.natorb_mo  = np.zeros((nirr, n_max, nmo1, nmo2),
                                             dtype=float)
            self.natorb_ao  = np.zeros((nirr, n_max, nmo1, nmo2),
                                             dtype=float)
            self.natorb_sym = np.zeros((nirr, n_max, nmo1), 
                                             dtype=int)

            # by default, compute the natural orbitals for all states
            # requested
            for irr in range(nirr):
                for st in range(self.nstates[irr]):
                    rdm        = self.rdm1(irr, st)
                    occ,nos_mo = np.linalg.eigh(rdm)
                    nos_ao     = np.matmul(self.scf.orbs, nos_mo)
                    
                    sort_ordr = self.order_orbs(occ)
                    occ       = occ[sort_ordr]
                    nos_mo    = nos_mo[:,sort_ordr]
                    nos_ao    = nos_ao[:,sort_ordr]
                    sym       = self.get_orb_sym(nos_ao)

                    self.natocc[irr, st, :]       = occ
                    self.natorb_mo[irr, st, :, :] = nos_mo
                    self.natorb_ao[irr, st, :, :] = nos_ao
                    self.natorb_sym[irr, st, :]   = sym 

        # always return the contents of the natorb/natocc arrays in wfn
        occ = self.natocc[irrep, state, :]
        if basis == 'ao':
            nos = self.natorb_ao[irrep, state, :, :]
        else:
            nos = self.natorb_mo[irrep, state, :, :]

        return occ, nos
   
    #
    def natural_sym(self, irrep, state):
        """return the symmetry of the natural orbitals"""

        # if natural orbitals haven't been generated, generate them now
        if self.natorb_sym is None:
            occ, nos = self.natural_orbs(irrep, state)

        return self.natorb_sym[irrep, state, :]

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
        """orders the mrci_ener array to create mrci_sorted, a 1D array
           of energies in ascending order, and a corresponding state_sorted
           array of [irrep, st] pairs for each energy in mrci_sorted."""

        if self.mrci_ener is None:
            sys.exit('cannot sort dftmrci.mrci_ener')

        nirr   = self.mol.n_irrep()
        istate = np.zeros((nirr),dtype=int)
        ntot   = sum(self.n_states())

        self.mrci_sorted  = np.zeros((ntot), dtype=float)
        self.state_sorted = []
        # mrci_ener is an irrep x maxroots array, with trailing values 
        # of 'zero'. 
        mrci_vals         = np.pad(self.mrci_ener, ((0,0),(0,1)), 
                                   'constant', 
                                   constant_values=((0,0),(0,0)))
       
        nsrt = 0
        while nsrt < ntot:
            eners = np.array([mrci_vals[irr, istate[irr]] 
                              for irr in range(nirr)])
            iirr  = np.argsort(eners)[0]

            self.mrci_sorted[nsrt]    = mrci_vals[iirr, istate[iirr]]
            self.state_sorted.append([iirr, istate[iirr]])
            # increment the number of elements in the sorted array,
            # and the irrep state index of irrep just chosen
            nsrt         += 1
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
 
        # max sure the symmetry orbital counts match up
        oval, ocnts = np.unique(sym_indx, return_counts=True)
        val,  cnts  = np.unique(self.scf.orb_sym, return_counts=True)

        #if not (ocnts==cnts).all():
        #    sys.exit('ERROR: cannot assign symmetry of natural orbs')

        return sym_indx

