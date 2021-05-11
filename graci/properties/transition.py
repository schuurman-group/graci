"""
module for compute moments of the dipole operator
for a given electronic state
"""
import graci.core.libs as libs

class Transition:
    """Transition class for determing transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # user defined quanties 
        self.init_states      = None
        self.final_states     = None
        self.nfinal_states    = None
        self.all_final_states = False

        # passed and global variables
        self.init_method  = None
        self.final_method = None
        self.trans_list   = None
        self.trans_mom    = None
        self.osc_str      = None
        self.label        = 'default'

    #-----------------------------------------------------------

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'transition'

    #
    def set_bra_method(self, bra_method):
        """set the method used to compute the ket state(s)"""
        self.bra_method = bra_method
        return

    #
    def set_ket_method(self, ket_method):
        """set the method used to compute the bra state(s)"""
        self.ket_method = ket_method
        return

    #
    def run(self):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used."""

        print("init")
        # initialize the bitsi library
        libs.init_bitsi(self)

        # compute the dipole moment integrals
        mu_ao = self.mol.pymol().intor('int1e_r')

        # contract with the MOs
        mu_mo = np.matmul(np.matmul(scf.orbs.T, mu_ao), scf.orbs)
        q_mo  = np.matmul(np.matmul(scf.orbs.T, q_mo), scf.orbs)

        # construct the trans_list array
        # currently store them as [initial state, final state]
        st_pairs = self.build_pair_list()

        # what happens depends on the libSI interface
#        for ket_irrep in init_irreps:
#            for bra_irrep in final_irreps:
#                libs.lib_func('transition_density_mrci', )


        # compute the oscillator strengths from the transition dipoles
        self.osc_str = []
        cnt       = 0
        for exc_pair in self.trans_list:
            exc_ener = bra.energy(exc_pair[1]) - ket.energy(exc_pair[0])
            self.osc_str.extend(2.*exc_ener*(self.trans_mom[cnt]**2)/3.)
            cnt += 1

        return

    def build_pair_list():
        """built an array of initial states, ordered by irrep"""

        print("init_states="+str(self.init_states))

        nirr = self.init_method.mol.n_irrep()
        if nirr != self.final_method.mol.n_irrep():
            sys.exit('initial and final states are different irrep')

        ket_list = [[] for irr in range(nirr)]
        ket_avail = self.init_method.n_states()

        for state in self.init_states:
            irrep = state[0]
            st    = state[1]
            if st <= st_avail[irrep]:
                ket_list[irrep].append(st)

        bra_list = [[] for irr in range(nirr)]
        bra_avail = self.final_method.n_states()

        # if all_final_states is true, that is overrides anything else
        if self.all_final_states:
            bra_list = [[st for st in range(bra_avail[irr])] 
                            for irr in range(nirr)]
        # if the number of final roots is given, pick the lowest n roots
        # by energy
        elif self.nfinal_states is not None: 
            # iterate through the final states, adding roots
            # to various irreps as needed
            nfinal = 0
            while (nfinal < self.nfinal_states and 
                    nfinal < sum(self.bra_method.n_states())):
                istate = self.bra_method.state_n(nfinal)
                bra_list[istate[0]].append(istate[1])
                nfinal += 1

        else:
            if self.final_states is None:
                sys.exit('no final states selected in Transition')
            for state in self.final_states:
                irrep = state[0]
                st    = state[1]
                if st <= bra_avail[irrep]:
                    bra_list[irrep].append(st)

        # create the pair list
        st_pairs = [[[] for irr_bra in range(nirr)] 
                        for irr_ket in range(nirr)]
        for b_irrep in range(nirr):
            for k_irrep in range(nirr):
                for b_state in bra_list[b_irrep]:
                    for k_state in ket_list[k_irrep]:
                        st_pairs[b_irrep][k_irrep].append(
                                                    [b_state, k_state])
                st_pairs[k_irrep][b_irrep] = st_pairs[b_irrep][k_irrep]

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



