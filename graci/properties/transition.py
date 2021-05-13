"""
module for compute moments of the dipole operator
for a given electronic state
"""
import numpy as np
import graci.core.libs as libs
import graci.citools.mrci_1tdm as mrci_1tdm
import graci.io.output as output

class Transition:
    """Transition class for determing transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # user defined quanties 
        self.init_states      = None
        self.final_states     = None
        self.istate_array     = None
        self.fstate_array     = None
        self.all_final_states = False
        self.init_method       = None
        self.final_method       = None

        # global variables
        self.tdms         = None
        self.trans_list   = None
        self.trans_mom    = None
        self.osc_str      = None
        self.label        = 'transition'

    #-----------------------------------------------------------

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'transition'

    #
    def set_final_method(self, final_method):
        """set the method used to compute the ket state(s)"""
        self.final_method = final_method
        return

    #
    def set_init_method(self, init_method):
        """set the method used to compute the bra state(s)"""
        self.init_method = init_method
        return

    #
    def run(self):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used."""

        mol_bra = self.final_method.mol
        mol_ket = self.final_method.mol

        scf_bra = self.final_method.scf
        scf_ket = self.final_method.scf
 
        # sanity check that orbitals and geometry are the same
        if np.any(scf_bra.orbs != scf_ket.orbs):
            sys.exit('transition moments require same bra/ket orbs')

        if mol_bra.pymol().atom != mol_ket.pymol().atom:
            sys.exit('transition moments require same geometry/basis')

        # initialize the bitsi library
        libs.init_bitsi(self)

        # compute the dipole moment integrals
        mu_ao = mol_ket.pymol().intor('int1e_r')

        # contract with the MOs
        mu_mo = [np.matmul(np.matmul(scf_ket.orbs.T, mu_ao[i]), 
                                    scf_ket.orbs) for i in range(3)]

        # construct the trans_list array
        # currently store them as [initial state, final state]
        self.build_pair_list()

        # grab the transition density matrices
        self.tdms = mrci_1tdm(self)

        # compute transition dipole moments and oscillator strengths
        self.build_trans_moments(mu_mo)

        # output oscillator strenths and transition dipole matrix 
        # elements
        output.print_transition(self)

        return

    #
    def tdm1(self, b_irr, b_st, k_irr, k_st):
        """return the tdm for the bra and ket state indicated"""

        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st,k_st])
        except:
            return None

        return self.tdms[b_irr][k_irr][:,:,ind]
 
    #
    def trans_dipole(self, b_irr, b_st, k_irr, k_st):
        """return transition dipole matrix element"""

        # return the transition dipole vector given the irreps and
        # state indices
        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st, k_st])
        except:
            None

        return self.trans_dipole[bra_irrep][ket_irrep][ind][:]

    # 
    def osc_strength(self,  b_irr, b_st, k_irr, k_st):
        """return the transition dipole moments for all pairs of
           bra_states and ket_states. If one or both are not specified,
           return all transition moments"""

        # return the transition dipole vector given the irreps and
        # state indices
        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st, k_st])
        except:
            None

        return self.osc_str[bra_irrep][ket_irrep][ind][:]

#--------------------------------------------------------------------------

    def build_trans_list(self):
        """built an array of initial states, ordered by irrep"""

        nirr = self.final_method.mol.n_irrep()
        if nirr != self.final_method.mol.n_irrep():
            sys.exit('initial and final states are different symmetry')

        # determine what initial states to include
        # -------------------------------------------------

        ket_list = [[] for irr in range(nirr)]
        ket_avail = self.final_method.n_states()

        # iterate through the initial states, adding roots
        # to various irreps as needed
        if self.init_states is not None:
            # iterate through the final states, adding roots
            # to various irreps as nee
            for state in self.init_states:
                istate = self.final_method.state_n(state)
                ket_list[istate[0]].append(istate[1])

        # use the istate array (format = irrep.state) to select
        # individual states
        elif self.istate_array:
            for state in self.fstate_array:
                irrep = int(state)
                st    = int(10*state) - 10*irrep
                if st <= st_avail[irrep]:
                    ket_list[irrep].append(st)

        else:
            sys.exit('no initial states selected in Transition')

        # determine what final states to include
        # -------------------------------------------------

        bra_list = [[] for irr in range(nirr)]
        bra_avail = self.final_method.n_states()

        # if all_final_states is true, that is overrides anything else
        if self.all_final_states:
            bra_list = [[st for st in range(bra_avail[irr])]
                            for irr in range(nirr)]

        # if the number of final roots is given, pick the lowest n roots
        # by energy
        elif self.final_states is not None:
            # iterate through the final states, adding roots
            # to various irreps as nee
            for state in self.final_states:
                fstate = self.final_method.state_n(state)
                bra_list[fstate[0]].append(fstate[1])

        # use the fstate array (format = irrep.state) to select
        # individual states
        elif self.fstate_array is not None:
            for state in self.fstate_array:
                irrep = int(state)
                st    = int(10*state) - 10*irrep
                if st <= bra_avail[irrep]:
                    bra_list[irrep].append(st)

        else:
            sys.exit('no final states selected in Transition')

        # create the pair list
        self.trans_list = [[[] for irr_bra in range(nirr)]
                               for irr_ket in range(nirr)]
        for b_irr in range(nirr):
            for k_irr in range(nirr):
                for b_st in bra_list[b_irr]:
                    for k_st in ket_list[k_irr]:
                        self.trans_list[b_irr][k_irr].append([b_st, k_st])

        return

    #
    def build_trans_moments(self, mu_mo):
        """builds the transition moment vector based on the based tdms"""

        # loop over the state pairs and construct transition 
        # dipole moments
        nirr_bra = self.bra_method.n_irrep()
        nirr_ket = self.ket_method.n_irrep()

        self.trans_dipole = [[[] for irr_bra in range(nirr_bra)]
                                 for irr_ket in range(nirr_ket)]

        self.osc_str      = [[[] for irr_bra in range(nirr_bra)]
                                 for irr_ket in range(nirr_ket)]

        for b_irr in range(nirr_bra):
            for k_irr in range(nirr_ket):
                for bra_ket in range(len(self.trans_list[b_irr][k_irr])):

                    states  = self.trans_list[b_irr][k_irr][bra_ket]
                    b_ener  = self.final_method.energy(b_irr, states[0])
                    k_ener  = self.init_method.energy(k_irr, states[1])

                    tdm_mat = self.tdms[b_irr][k_irr][:,:,bra_ket]
                    td_vec  = np.array([np.sum(tdm_mat*mu_mo[x])
                              for x in range(3)], dtype=float)
                    osc_vec = (2./3.)*(b_ener-k_ener)*td_vec**2

                    self.trans_dipole[b_irr][k_irr].append([td_vec])
                    self.osc_str[b_irr][k_irr].append([osc_vec])

        return
