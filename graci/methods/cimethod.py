"""
Module for computing CI and CI+PT2 energies
"""
import os as os
import sys as sys
import numpy as np
import copy as copy
import graci.core.params as params
import graci.core.orbitals as orbitals
import graci.core.ao2mo as ao2mo
import graci.utils.timing as timing
import graci.io.output as output
import graci.properties.moments as moments
from pyscf import gto

class Cimethod:
    """Class constructor for CI objects"""
    def __init__(self):
        # user defined quanties
        self.charge         = 0
        self.mult           = None
        self.nstates        = None 
        self.print_orbitals = False
        self.save_wf        = False
        self.wf_thresh      = 1.
        self.ref_state      = -1
        self.ddci           = True
        self.guess_label    = None
        self.verbose        = True
        self.hparam         = None
        self.scf_label      = 'default'
        self.label          = 'default'
        self.precision      = 'single'

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
        self.dmats          = {'adiabatic' : None, 'diabatic' : None}
        # natural orbital occupations (shape = (nirr, nmo, nstates)
        self.natocc         = {'adiabatic' : None, 'diabatic' : None}
        # natural orbitals, AO basis (same shape as the density matrices)
        self.natorb_ao      = {'adiabatic' : None, 'diabatic' : None}
        # List of determinant bit strings (one per irrep)
        # for various representations (adiabatic, diabatic,...)
        self.det_strings     = {'adiabatic' : None, 'diabatic' : None}
        # Lists of eigenvectors in the determinant basis (one per irrep)
        # for various representations (adiabatic, diabatic,...)
        self.vec_det         = {'adiabatic' : None, 'diabatic' : None}
        # R_n-1 - R_n MO overlaps
        self.smo            = None
        # ADT matrices (one per irrep)
        self.adt            = None
        # Diabatic potential matrices (one per irrep)
        self.diabpot        = None
        # mos 
        self.mos            = None
        # number of MOs in CSF expansions
        self.nmo            = None
        # mo energies
        self.emo            = None
        # mo symmetries
        self.mosym          = None

# Required functions #############################################################

    def set_scf(self, scf):
        """set the scf object for the dftmrci class object"""

        # check that the coordinates in the graci.molecule
        # object agree with pymol. If not: we need re-run SCF
        if scf.mol.coords_updated():
            scf.run(scf.mol, None)
            self.update_eri()

        self.scf = scf
        self.nel = scf.nel

        # Default spin multiplicity: inherited from the scf object
        if self.mult is None:
            self.mult = scf.mult

        # Default charge: inherited from the scf object
        if self.charge is None or self.charge == scf.charge:
            self.charge  = scf.charge
            self.ref_occ = copy.copy(scf.orb_occ)[:self.nmo]

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
            self.ref_occ = np.zeros(self.nmo, dtype=float)
            self.ref_occ[:nclsd]            = 2.
            self.ref_occ[nclsd:nclsd+nopen] = 1.
        
        return

    def set_geometry(self, atms, crds, update_scf=False):
        """calls set_geometry of the corresponding molecule object"""
        self.scf.mol.set_geometry(atms, crds)
        if update_scf:
            self.set_scf(scf)

    def scf_exists(self):
        """return true if scf object is not None"""
        try:
            return type(self.scf).__name__ == 'Scf'
        except:
            return False

    def set_hamiltonian(self):
        """sets the default Hamiltonian based on the XC functional
        if the user has not specified one"""

        defaults = {'hf'               : 'abinitio',
                    'hyb_gga_xc_qtp17' : 'qe8',
                    'bhandhlyp'        : 'r2017'}
        
        if self.hamiltonian is None:
            try:
                self.hamiltonian = defaults[self.scf.xc.lower()]
            except:
                sys.exit('\n\n Error: could not determine a default'
                         + ' Hamiltonian for XC functional '
                         + self.scf.xc)
        
        return
        
    # 
    def n_irrep(self):
        """return the number of irreps"""
        return len(self.nstates)

    #
    def irreps_nonzero(self):
        """returns the indices of the irreps with a
        non-zero no. of roots"""        
        return np.where(self.nstates > 0)[0]
        
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
    def rdm(self, istate, rep='adiabatic'):
        """return the density matrix for the state istate"""

        if self.dmats[rep] is not None and istate < self.n_states():
            return self.dmats[rep][istate, :, :]
        else:
            print("rdm called but density does not exist")
            return None

    #
    def rdm_sym(self, irrep, state, rep='adiabatic'):
        """return the density matrix for the state in 
        the array 'states' for the irrep 'irr'"""

        istate = self.state_index(irrep, state)

        if self.dmats[rep] is not None:
            return self.dmats[rep][istate, :, :]
        else:
            print("rdm_sym called but density does not exist")
            return self.dmats[rep]

    # 
    def natorbs(self, istate, basis='ao', rep='adiabatic'):
        """return natural orbitals and occupations for state 'istate' """

        if basis == 'ao' and self.natorb_ao[rep] is None:
            return None

        if basis == 'mo':
            return None

        (nst, nmo) = self.natocc[rep].shape
        if istate >= nst:
            return None

        # always return the contents of the natorb/natocc arrays in wfn
        occ = self.natocc[rep][istate, :]
        nos = self.natorb_ao[rep][istate, :, :]

        return occ, nos

    #
    def update_hparam(self, hparams):
        """update the values of the Hamiltonian parameters"""

        self.hparam = hparams
        return

    #
    def update_eri(self, eri_mo=None):
        """update the MO integral information used by the CI 
           object. Particularly: the dimension of the MO basis

           Arguments:
               ao2mo: the Ao2mo class object used to generate 
                      the MO integrals

           Returns:
               None
        """
        if eri_mo == None:
            eri_mo         = ao2mo.Ao2mo()
            eri_mo.emo_cut = self.mo_cutoff
            eri_mo.run(self.scf, self.precision)

        self.nmo   = eri_mo.nmo
        self.emo   = eri_mo.emo
        self.mosym = eri_mo.mosym
        self.mos   = eri_mo.orbs
        #print(self.label+' update_eri, nmo, emo, mosym, mos='+str(self.nmo)+','+str(self.emo)+','+str(self.mosym)+','+str(self.mos.shape),flush=True)

        return

#########################################################################
    #
    def build_nos(self, rep='adiabatic'):
        """calls routines in orbitals module to construct natural
           orbitals and occupation vectors for each state.

           Arguments:
            None

           Returns:
            None
        """
        n_tot = self.n_states()
        if n_tot > 0:
            (dim1, nmo) = self.rdm(0).shape
            nao         = self.scf.mol.nao
        else:
            return

        self.natorb_ao[rep] = np.zeros((n_tot, nao, nmo), dtype=float)
        self.natocc[rep]    = np.zeros((n_tot, nmo), dtype=float)
        for i in range(n_tot):
            occ, nos = orbitals.build_nos(self.rdm(i), basis='ao',
                                                 mos=self.mos)
            self.natorb_ao[rep][i,:,:] = nos
            self.natocc[rep][i,:]      = occ

        return

    # 
    def print_nos(self, rep='adiabatic'):
        """Calls routines in orbitals module to print natural orbitals
           to file. Default file format is molden

           Arugments:
               None

           Returns:
               None
        """

        n_tot  = self.n_states()
        syms   = [self.scf.mol.irreplbl[self.state_sym(i)[0]]
                  for i in range(n_tot)]

        for i in range(n_tot):
            fname = 'nos_'+str(i+1)+'_'+syms[i]+'_molden'
            orbitals.export_orbitals(fname, self.scf.mol, 
                   self.natorb_ao[rep][i,:,:],
                   orb_occ=self.natocc[rep][i,:], 
                   orb_dir=str(type(self).__name__)+'.'+str(self.label),
                   fmt='molden')

        return

    #
    def print_promotion(self, refstate):
        """Print the promotion numbers relative to state 'refstate'

           Arguments:
              refstate: integer adiabatic state index

           Returns:
              None
        """

        # build the list and symmetries for subsequent printing
        n_tot  = self.n_states()
        states = [i for i in range(n_tot)]
        syms   = [self.scf.mol.irreplbl[self.state_sym(i)[0]]
                  for i in range(n_tot)]
        pd     = np.zeros(n_tot, dtype=float)
        pa     = np.zeros(n_tot, dtype=float)

        rdm_ref = self.rdm(refstate)
        for i in range(n_tot):

            if i == refstate:
                continue

            rdm_i = self.rdm(i)
            # also compute attachment and detachment numbers
            # (relative to ground state)
            wt, ndo  = orbitals.build_ndos(rdm_i, rdm_ref)
            pd[i], pa[i] = orbitals.promotion_numbers(wt, ndo)

        if self.verbose:
            output.print_promotion(refstate, states, syms, pd, pa)
        
        return

    #
    def print_moments(self, rep='adiabatic'):
        """ print the moments using the natural orbitals and occupation
            vector

            Arguments:
               None

            Returns:
               None
        """

        # we'll also compute 1-electron properties by
        # default.
        momts = moments.Moments()
        momts.run(self.scf.mol, self.natocc[rep], self.natorb_ao[rep])

        n_tot  = self.n_states()
        states = [i for i in range(n_tot)]
        syms   = [self.scf.mol.irreplbl[self.state_sym(i)[0]]
                  for i in range(n_tot)]

        if self.verbose:
            output.print_moments(states, syms, momts)

        return

    #
    def order_energies(self):
        """orders the self.energies_sym array to create self.energies,
        a 1D array of energies in ascending order, and a corresponding
        sym_sorted array of [irrep, st] pairs for each energy in 
        mrci_sorted."""

        if self.energies_sym is None:
            sys.exit('error in order_ener: energies_sym is None')

        nirr    = self.n_irrep()
        istate  = np.zeros((nirr), dtype=int)
        n_tot   = self.n_states()

        self.energies   = np.zeros((n_tot), dtype=float)
        self.sym_sorted = []
        # mrci_ener_sym is an irrep x maxroots array, with trailing values 
        # of 'zero'. 
        ener_vals = np.pad(self.energies_sym, ((0,0),(0,1)), 
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
    def diabatize(self):
        """
        Constructs the diabatic potential matrix using
        a prior-calculated ADT matrix
        """

        # number of irreps
        nirr = self.n_irrep()

        # initialise the list of diabatic potentials
        self.diabpot = [None for i in range(nirr)]

        # loop over irreps
        for irrep in self.irreps_nonzero():

            # no. states
            nstates = self.n_states_sym(irrep)

            # adiabatic potential matrix
            vmat = np.zeros((nstates, nstates), dtype=float)
            np.fill_diagonal(vmat, self.energies_sym[irrep][:nstates])

            # diabatic potential matrix
            wmat = np.matmul(self.adt[irrep].T,
                             np.matmul(vmat, self.adt[irrep]))

            self.diabpot[irrep] = wmat

        return
