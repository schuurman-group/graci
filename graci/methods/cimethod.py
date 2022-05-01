"""
Module for computing CI and CI+PT2 energies
"""
import os as os
import sys as sys
import importlib
import numpy as np
import copy as copy
import graci.utils.timing as timing
import graci.io.output as output

# Allowed MRCI and DFT/MRCI Hamiltonian labels
# (not sure if this is the correct home for this, but we
#  will put it here for now)
hamiltonians   = ['canonical',
                  'grimme_standard',
                  'grimme_short',
                  'lyskov_standard',
                  'lyskov_short',
                  'heil17_standard',
                  'heil17_short',
                  'heil18_standard',
                  'heil18_short',
                  'cvs_standard']

class Cimethod:
    """Class constructor for CI objects"""
    def __init__(self):
        # user defined quanties
        self.charge         = 0
        self.mult           = None
        self.nstates        = []
        self.print_orbitals = False
        self.save_wf        = False
        self.ddci           = True

        # class variables
        # total number of electrons
        self.nel            = 0
        # SCF object
        self.scf            = None 
        # reference occupation 
        self.ref_occ        = None
        # Energies, by adiabatic state
        self.energies       = None
        # ci energies by irrep and root index
        self.energies_sym   = None
        # irrep,state pairs sorted by energy
        self.sym_sorted     = None
        # shape: (nirr,nmo,nmo,nstates))
        self.dmats          = None
        # natural orbital occupations (shape = (nirr, nmo, nstates)
        self.natocc         = None
        # natural orbitals, MO basis (same shape as the density matrices)
        self.natorb_mo      = None
        # natural orbitals, AO basis (same shape as the density matrices)
        self.natorb_ao      = None
        # symmetry of the natural orbitals
        self.natorb_sym     = None


# Required functions #############################################################

    def set_scf(self, scf):
        """set the scf object for the dftmrci class object"""

        self.scf = scf
        self.nel = scf.nel

        # Default spin multiplicity: inherited from the scf object
        if self.mult is None:
            self.mult = scf.mult

        # Default charge: inherited from the scf object
        if self.charge is None or self.charge == scf.charge:
            self.charge  = scf.charge
            self.ref_occ = copy.copy(scf.orb_occ)

        # if the number of electrons have changed update ref_occ
        if self.charge != scf.charge:
            self.nel = scf.nel + (scf.charge - self.charge)
            # set the occupation vector to yield the multiplicity
            # with a maximum number of closed shells
            nopen = self.mult - 1
            nclsd = int(0.5 * (self.nel - nopen))
            if 2*nclsd + nopen != self.nel:
                sys.exit('Inconsistent charge='+str(self.charge)+ 
                                      ' / multp='+str(self.mult))
            self.ref_occ = np.zeros(self.scf.nmo, dtype=float)
            self.ref_occ[:nclsd]            = 2.
            self.ref_occ[nclsd:nclsd+nopen] = 1.

        return

    def scf_exists(self):
        """return true if scf object is not None"""
        try:
            return type(self.scf).__name__ is 'Scf'
        except:
            return False

    # 
    def n_irrep(self):
        """return the number of irreps"""
        return len(self.nstates)
 
    #
    def n_states(self):
        """total number of states"""

        return sum(self.n_states_sym())

    #
    def n_states_sym(self, irrep=None):
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
        return self.energies[istate]

    #
    def energy_sym(self, irrep, state):
        """return the energy of state by irrep, st index"""
        return self.energies_sym[irrep, state]

    #
    def rdm(self, istate):
        """return the density matrix for the state istate"""

        if self.dmats is not None and istate < self.n_states():
            return self.dmats[istate, :, :]
        else:
            print("rdm called but density does not exist")
            return None

    #
    def rdm_sym(self, irrep, state):
        """return the density matrix for the state in 
        the array 'states' for the irrep 'irr'"""

        istate = self.state_index(irrep, state)

        if self.dmats is not None:
            return self.dmats[istate, :, :]
        else:
            print("rdm_sym called but density does not exist")
            return self.dmats

    #
    @timing.timed
    def build_nos(self):
        """print the natural orbitals for irrep irr and state state"""

        # if the 1RDMs don't exist, then we have a problem
        if self.n_states() > 0 and self.rdm(0) is None:
            print("can't return natural orbitals: 1RDMs don't exist")
            return None
            
        # check that istate is less than the total number of states:
        n_tot           = self.n_states()
        (nmo1, nmo2)    = self.rdm(0).shape

        self.natocc     = np.zeros((n_tot, nmo2), dtype=float)
        self.natorb_mo  = np.zeros((n_tot, nmo1, nmo2), dtype=float)
        self.natorb_ao  = np.zeros((n_tot, nmo1, nmo2), dtype=float)
        self.natorb_sym = np.zeros((n_tot, nmo2), dtype=int)

        # by default, compute the natural orbitals for all states
        # requested
        for ist in range(n_tot):
            rdmi       = self.rdm(ist)
            occ,nos_mo = np.linalg.eigh(rdmi)
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
    @timing.timed
    def build_ndos(self, ref_state, basis='mo'):
        """ build natural difference orbitals"""

        # first check if bra and ket objects define an
        # rdm method. If not, return an empty dictionary

        nmo    = self.scf.nmo
        nao    = self.scf.mol.nao
        mos    = self.scf.orbs

        ndo    = np.zeros((nao, nmo, self.n_states()), dtype=float)
        ndo_wt = np.zeros((nmo, self.n_states()), dtype=float)

        for ist in range(self.n_states()):

            if ist == ref_state:
                continue

            # delta is the different 1RDM between ist and
            # reference state (likely the ground state)
            delta = self.rdm(ist) - self.rdm(ref_state) 

            # form different natural orbitals and transform
            wt, ndo_mo = np.linalg.eigh(delta)

            # sort NDO wts by increasing magnitude (hole 
            # orbitals to start, then particle
            ordr           = np.argsort(wt)
            ndo_wt[:, ist] = wt[ordr].copy()
            if basis == 'ao':
                ndo_ao         = mos @ ndo_mo
                ndo[:, :, ist] = ndo_ao[:,ordr].copy()
            else:
                ndo[:, :, ist] = ndo_mo[:,ordr].copy()            

        return ndo, ndo_wt

    # 
    def promotion_numbers(self, ndos, ndo_wts):
        """compute the detachment and attachment promotion numbers
           for a given set of NDOs"""

        (nao, nmo, nst) = ndos.shape
        p_detach        = np.zeros(nst, dtype=float)
        p_attach        = np.zeros(nst, dtype=float)        

        for ist in range(nst):

            wt_mat = np.diag([wt if wt < 0. else 0. for
                                                wt in ndo_wts[:, ist]])
            d_mat  = ndos[:,:,ist] @ wt_mat @ ndos[:,:,ist].transpose()
            p_detach[ist] = np.trace(d_mat)

            wt_mat = np.diag([wt if wt > 0. else 0. for
                                                wt in ndo_wts[:, ist]])
            d_mat  = ndos[:,:,ist] @ wt_mat @ ndos[:,:,ist].transpose()
            p_attach[ist] = np.trace(d_mat)

        return p_detach, p_attach

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

#########################################################################
    #
    def order_energies(self):
        """orders the self.energies_sym array to create self.energies, a 1D array
           of energies in ascending order, and a corresponding sym_sorted
           array of [irrep, st] pairs for each energy in mrci_sorted."""

        if self.energies_sym is None:
            sys.exit('error in order_ener: energies_sym is None')

        nirr    = self.n_irrep()
        istate  = np.zeros((nirr), dtype=int)
        n_tot   = self.n_states()

        self.energies   = np.zeros((n_tot), dtype=float)
        self.sym_sorted = []
        # mrci_ener_sym is an irrep x maxroots array, with trailing values 
        # of 'zero'. 
        ener_vals         = np.pad(self.energies_sym, ((0,0),(0,1)), 
                                   'constant', 
                                   constant_values=((0,0),(0,0)))
       
        n_srt = 0
        while n_srt < n_tot:
            eners = np.array([ener_vals[irr, istate[irr]] 
                              for irr in range(nirr)])
            iirr  = np.argsort(eners)[0]

            self.energies[n_srt] = ener_vals[iirr, istate[iirr]]
            self.sym_sorted.append([iirr, istate[iirr]])
            # increment the number of elements in the sorted array,
            # and the irrep state index of irrep just chosen
            n_srt        += 1
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

        for ist in range(self.n_states()):
            self.export_orbitals_state(ist, orb_format=orb_format, 
                                                   orb_dir=orb_dir)
        return

    #
    def export_orbitals_state(self, state, orb_format='molden', 
                                         orb_dir=True, cart=None):
        """export orbitals of a single state to various file formats"""

        if state >= self.n_states():
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
