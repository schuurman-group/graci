"""
module for compute moments of the dipole operator
for a given electronic state
"""
import numpy as np
import graci.core.libs as libs

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
        st_pairs = self.build_pair_list()

        # what happens depends on the libSI interface
        nirr = mol_bra.n_irrep()
        print("\n")
        for bra_irr in range(nirr):
            for ket_irr in range(nirr):

                pair_list = st_pairs[bra_irr][ket_irr]
                print("pair_list="+str(pair_list))
#                libs.lib_func('transition_density_mrci', )


        # compute the oscillator strengths from the transition dipoles
        #self.osc_str = []
        #cnt       = 0
        #for exc_pair in self.trans_list:
        #    exc_ener = bra.energy(exc_pair[1]) - ket.energy(exc_pair[0])
        #    self.osc_str.extend(2.*exc_ener*(self.trans_mom[cnt]**2)/3.)
        #    cnt += 1

        return

    def build_pair_list(self):
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
        st_pairs = [[[] for irr_bra in range(nirr)] 
                        for irr_ket in range(nirr)]
        for b_irr in range(nirr):
            for k_irr in range(nirr):
                for b_st in bra_list[b_irr]:
                    for k_st in ket_list[k_irr]:
                        st_pairs[b_irr][k_irr].append([b_st, k_st])

        return st_pairs

    #
    def tdm1(self, bra_irrep, bra_st, ket_irrep, ket_st):
        """return the tdm for the bra and ket state indicated"""

        return
 
    #
    def trans_dipole(self, bra, ket):
        """return transition dipole matrix element"""

        # find the state pair in the transition list
        ind = self.tran_list.index([bra, ket])

        return self.trans_mom[ind]

    # 
    def osc_strength(self, bra, ket):
        """return the transition dipole moments for all pairs of
           bra_states and ket_states. If one or both are not specified,
           return all transition moments"""

        # find the state pair in the transition list
        ind = self.tran_list.index([bra, ket])

        return self.osc_str[ind]



