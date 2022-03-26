"""
Parent state interaction class
"""
import numpy as np
import graci.utils.timing as timing

class Interaction:
    """Parent state interaction class to be inherited by any module
       computing matrix elements <psi_i|O|psi_j> for some operator O"""
    def __init__(self):

        # all state interaction sub-classes must have a 
        # bra object and a ket object

        # user set variables
        #-----------------------------------------
        # labels for those objects (set by user)
        self.bra_label    = None
        self.ket_label    = None
        self.bra_states   = []
        self.ket_states   = []

        # Internal class variables (should not be accessed directly)
        #-----------------------------------------------------------
        self.bra_obj      = None
        self.ket_obj      = None
        # list of objects from which to pull states
        self.objs         = {} 
        # state list -- list of groups of states
        self.states       = {}
        # symmetries of the states in the state list
        self.symmetries   = {}
        # spin info for the states in group
        self.spins        = {}

    class SpinInfo:
        """class to hold spin formation for a set of states"""
        def __init__(self, method_obj):
            # total spin
            self.mult  = method_obj.mult
            self.S     = (self.mult - 1.)/2.
            self.M     = np.array([-self.S + i for 
                                   i in range(self.mult)], dtype=float)

    def same_braket(self):
        """return true if the bra and ket objects are the same"""
        class_bra = type(self.bra_obj).__name__ 
        class_ket = type(self.ket_obj).__name__ 

        lbl_bra  = self.bra_obj.label
        lbl_ket  = self.ket_obj.label

        return class_bra == class_ket and lbl_bra == lbl_ket 

    #
    def add_states(self, lbl, obj, state_list):
        """append a list of states and corresponding object"""

        # add object and state list
        self.objs[lbl]   = obj
        self.states[lbl] = state_list

        # also store state symmetries
        self.symmetries[lbl] = []
        for state in state_list:
            # state_sym returns an [irrep, st_indx] pair
            self.symmetries[lbl].append(obj.state_sym(state)[0])

        return

    # 
    def add_spins(self, lbl, spins):
        """add a spin object to group 'lbl'"""
        self.spins[lbl] = spins
        return

    #
    def get_states(self, lbl):
        """return the states and symmetries corresponding to grp lbl"""

        state_list = None
        state_sym  = None

        if lbl in list(self.states.keys()):
            state_list = self.states[lbl]    

        if lbl in list(self.symmetries.keys()):
            state_sym = self.symmetries[lbl]

        return state_list, state_sym

    #
    def get_spins(self, lbl):
        """return the spin object associated with 'lbl'"""
        spin = None

        if lbl in list(self.spins.keys()):
            spin = self.spins[lbl]

        return spin

    #
    def get_bkobj(self, lbl):
        """return the method object corresponding to lbl"""

        obj = None

        if lbl in list(self.objs.keys()):
            obj = self.objs[lbl]

        return obj

    #
    @timing.timed 
    def build_pair_list(self, bra, ket, 
                              pairs = 'full', sym_blk = False):
        """
           builds an array of initial and final states, ordered by
           adiabatic energy. The symmetries of each state will be
           stored in a separatearray.

           If sym_blk is True, state pairs are grouped bra/ket irrep,
           If sym_blk is False, output is a list of state pairs
    
          The string list_type can be one of the following:
   
           'full'         <-> build the full matrix of bra and ket
                              state pairs, i.e. (bra, ket) and (ket, bra)
    
           'lower'        <-> build the lower triangle of initial and final
                              state pairs i.e. (bra, ket) only 
    
           'nodiag'       <-> build the lower triangle of initial and final
                              state pairs minus the on-diagonal elements
                              -- this is only active when bra == ket label
        """

        # if either bra or ket group not in stored group names, return None
        if (bra not in list(self.objs.keys()) or 
              ket not in list(self.objs.keys())):
            return None

        # if bra and ket objects are same,
        braket_same = self.same_braket()

        # initialize the pair list. If symblk, pair_list is an nirr x nirr
        # state pairs. Else, just a list of state pairs
        if sym_blk:
            nb = self.objs[bra].n_irrep()
            nk = self.objs[ket].n_irrep()
            pair_list = [[[] for b in range(nb)] for k in range(nk)]
        else:
            pair_list = []

        for k in range(len(self.states[ket])):
            kst = self.states[ket][k]
            ksm = self.objs[ket].state_sym(kst)
    
            for b in range(len(self.states[bra])):
                bst = self.states[bra][b]
                bsm = self.objs[bra].state_sym(bst)
   
                # exclude bra == ket pairs: 'lower_i>j'
                if braket_same and bst == kst and pairs == 'nodiag':
                    continue

                # only include unique pairs if 'lower' or 'nodiag'
                if b < k and pairs != 'full':
                    continue
    
                # add the state pair
                if sym_blk:
                    pair_list[bsm[0]][ksm[0]].append([bsm[1], ksm[1]])
                else:
                    pair_list.append([bst, kst])
    
        return pair_list

    #
    def set_ket(self, ket_method):
        """set the method used to compute the ket state(s)"""

        if ket_method.label != self.ket_label:
            print('WARNING: set_ket called, but ' + 
                   str(ket_method.label)+' != ' + str(self.ket_label))

        self.ket_obj = ket_method
        return

    # 
    def ket_exists(self):
        """return True is self.ket_obj is not None"""
        return self.ket_obj is not None

    #
    def set_bra(self, bra_method):
        """set the method used to compute the bra state(s)"""

        if bra_method.label != self.bra_label:
            print('WARNING: set_bra called, but ' + 
                   str(bra_method.label)+' != ' + str(self.bra_label))

        self.bra_obj = bra_method
        return

    # 
    def bra_exists(self):
        """return True is self.ket_obj is not None"""
        return self.bra_obj is not None
