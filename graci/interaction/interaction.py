"""
Parent state interaction class
"""

import graci.utils.timing as timing
import sys as sys

class Interaction:
    """Parent state interaction class to be inherited by any module
       computing matrix elements <psi_i|O|psi_j> for some operator O"""
    def __init__(self):
        # Class name: this is to be overridden by the Child
        # class name
        self.class_name     = 'interaction'

        # user defined quanties common to all state interaction
        # calculations
        self.init_states      = None
        self.final_states     = None
        self.init_states_sym  = None
        self.final_states_sym = None
        self.all_final_states = False
        self.init_label       = None
        self.final_label      = None

        # global variables common to all state interaction
        # calculations
        #
        # method object for the bra states
        self.bra_obj        = None
        # method object for the ket states
        self.ket_obj        = None
        # comprehensive list of lists of all transition pairs
        # for all interaction matrix blocks
        self.trans_list     = []
        # list of lists of transition pairs, grouped by irrep
        # (format for bitsi) for all interaction matrix blocks
        self.trans_list_sym = []
        # convenient to just keep track of initial states for each
        # interaction matrix block
        self.ket_list       = []
        # simple boolean that is true if bra and ket objects
        # are identical. If this is true, default behavior is 
        # to only consider unique bra/ket pairs
        self.braket_iden    = False

    #
    def name(self):
        """ return the name of the class object as a string"""
        return self.class_name

    #
    def set_ket(self, ket_method):
        """set the method used to compute the ket state(s)"""
        try:
            self.init_label = ket_method.label
        except:
            self.init_label = None
            return
        self.ket_obj = ket_method
        return
        
    # 
    def ket_exists(self):
        """return True is self.ket_obj is not None"""

        return self.ket_obj is not None

    #
    def set_bra(self, bra_method):
        """set the method used to compute the bra state(s)"""

        try:
            self.final_label = bra_method.label
        except:
            self.final_label = None
            return

        self.bra_obj = bra_method
        return

    #
    def bra_exists(self):
        """return True if self.bra_obj is not None"""

        return self.bra_obj is not None

    #
    def trans_index(self, b_st, k_st, iblock):
        """return the index in the trans_list corresponding to
           [b_st, k_st] for the iblock'th block of the interaction
           matrix"""

        try:
            indx = self.trans_list[iblock].index([b_st, k_st])
        except:
            indx = None

        return indx

    # 
    def trans_sym_index(self, b_irr, k_irr, b_st, k_st, iblock):
        """return the index in the trans_sym_list corresponding to
           b_st,k_st for the iblock'th block of the interaction
           matrix"""
        
        try:
            ind = self.trans_list_sym[iblock][b_irr][k_irr].index([b_st, k_st])
        except:
            return None

        return ind

    ##
    #@timing.timed
    #def build_trans_list(self, list_type='full'):
    #    """
    #    builds an array of initial and final states, ordered by
    #       adiabatic energy. The symmetries of each state will be
    #       stored in a separatearray.
    #
    #       The string list_type can be one of the following:
    #
    #       'full'      <-> build the full matrix of initial and final
    #                       state pairs
    #
    #       'lower'     <-> build the lower triangle of initial and final
    #                       state pairs
    #
    #       'lower_i>j' <-> build the lower triangle of initial and final
    #                       state pairs minus the on-diagonal elements
    #       """
    #
    #    # check on list_type
    #    allowed_types = ['full', 'lower', 'lower_i>j']
    #    if list_type not in allowed_types:
    #        sys.exit('\n Illegal list_type in build_trans_list.'
    #                 +' Allowed types: \n\n '+str(allowed_types))
    #
    #    # determine what initial states to include
    #    # -------------------------------------------------
    #        
    #    state_list = []
    #
    #    # iterate through the initial states, adding roots
    #    # to various irreps as needed -- these states are 
    #    # ordered by adiabatic energy
    #    ninit_total = self.ket_obj.n_state()
    #    if self.init_states is not None:
    #        # iterate through the initial states, adding roots
    #        # to various irreps as needed
    #        for istate in self.init_states:
    #            if istate < ninit_total:
    #                state_list.append(istate)
    #
    #    # use the istate array (format = irrep.state) to select
    #    # individual states - input is form irrep.state
    #    if self.init_states_sym is not None:
    #        for state in self.init_states_sym:
    #            irrep   = int(state)
    #            st      = int(10*state) - 10*irrep
    #            istate  = self.ket_obj.state_index(irrep,st)
    #            if istate is not None:
    #                state_list.append(istate)
    #
    #    # take the union of the two approaches
    #    self.ket_list = sorted(list(set(state_list)))
    #    ket_sym  = [self.ket_obj.state_sym(self.ket_list[i])[0] 
    #                                for i in range(len(self.ket_list))]
    #
    #    # determine what final states to include
    #    # -------------------------------------------------
    #
    #    state_list  = []
    #
    #    # iterate through the final states, adding roots
    #    # to various irreps as needed -- these states are 
    #    # ordered by adiabatic energy
    #    nfinal_total = self.bra_obj.n_state()
    #    if self.final_states is not None:
    #        # iterate through the final states, adding roots
    #        # to various irreps as nee
    #        for fstate in self.final_states:
    #            if fstate < nfinal_total:
    #                state_list.append(fstate)
    #
    #    # use the istate array (format = irrep.state) to select
    #    # individual states
    #    if self.final_states_sym is not None:
    #        for state in self.final_states_sym:
    #            irrep   = int(state)
    #            st      = int(10*state) - 10*irrep
    #            fstate  = self.bra_obj.state_index(irrep,st)
    #            if fstate is not None:
    #                state_list.append(fstate)
    #
    #    # take the union of the two approaches
    #    bra_list = sorted(list(set(state_list)))
    #    bra_sym  = [self.bra_obj.state_sym(bra_list[i])[0]
    #                           for i in range(len(bra_list))]
    #
    #    # create the pair list -- we build the list this way 
    #    # in order to facilitate call to the tdm code in bitci
    #    # which computes tmds for pairs of states in particular
    #    # irreps
    #    nirr_bra = self.bra_obj.n_irrep()
    #    nirr_ket = self.ket_obj.n_irrep()
    #
    #    self.trans_list_sym = [[[] for irr_bra in range(nirr_bra)]
    #                               for irr_ket in range(nirr_ket)]
    #
    #    # trans list is a list of pairs, with first dimension 
    #    # equal to the number of initial states
    #    self.trans_list     = []
    #
    #    for k_st in self.ket_list:
    #        ksym = self.ket_obj.state_sym(k_st)
    #
    #        for b_st in bra_list:
    #            bsym = self.bra_obj.state_sym(b_st)
    #
    #            # exclude bra == ket pairs: 'lower_i>j'
    #            if list_type == 'lower_i>j' and k_st == b_st:
    #                continue
    #
    #            # only include unique pairs: 'lower' and 'lower_i>j'
    #            if list_type != 'full' \
    #               and self.trans_index(k_st, b_st) is not None:
    #                continue
    #
    #            # add the state pair
    #            self.trans_list_sym[bsym[0]][ksym[0]].append(
    #                                       [bsym[1], ksym[1]])
    #            self.trans_list.append([b_st, k_st])
    #
    #    return

    #
    @timing.timed
    def build_trans_list_new(self, bstring='bra', kstring='ket',
                             list_type='full'):
        """
        builds an array of initial and final states, ordered by
           adiabatic energy. The symmetries of each state will be
           stored in a separatearray.

           The string list_type can be one of the following:

           'full'      <-> build the full matrix of initial and final
                           state pairs

           'lower'     <-> build the lower triangle of initial and final
                           state pairs

           'lower_i>j' <-> build the lower triangle of initial and final
                           state pairs minus the on-diagonal elements
           """

        # check on list_type
        allowed_types = ['full', 'lower', 'lower_i>j']
        if list_type not in allowed_types:
            sys.exit('\n Illegal list_type in build_trans_list.'
                     +' Allowed types: \n\n '+str(allowed_types))

        # set up the bra and ket variables
        allowed_strings = ['bra', 'ket']
        if bstring not in allowed_strings or kstring not in allowed_strings:
            sys.exit('\n Illegal bstring/kstring in build_trans_list.'
                     +' Allowed strings: \n\n '+str(allowed_strings))
            
        istring = 'init'  if kstring == 'ket' else 'final'
        fstring = 'final' if bstring == 'bra' else 'init'

        bra_obj          = getattr(self, bstring+'_obj')
        ket_obj          = getattr(self, kstring+'_obj')
        init_states      = getattr(self, istring+'_states')
        init_states_sym  = getattr(self, istring+'_states_sym')
        final_states     = getattr(self, fstring+'_states')
        final_states_sym = getattr(self, fstring+'_states_sym')

        # determine what initial states to include
        # -------------------------------------------------
            
        state_list = []

        # iterate through the initial states, adding roots
        # to various irreps as needed -- these states are 
        # ordered by adiabatic energy
        ninit_total = ket_obj.n_state()
        if init_states is not None:
            # iterate through the initial states, adding roots
            # to various irreps as needed
            for istate in init_states:
                if istate < ninit_total:
                    state_list.append(istate)

        # use the istate array (format = irrep.state) to select
        # individual states - input is form irrep.state
        if init_states_sym is not None:
            for state in init_states_sym:
                irrep   = int(state)
                st      = int(10*state) - 10*irrep
                istate  = ket_obj.state_index(irrep,st)
                if istate is not None:
                    state_list.append(istate)

        # take the union of the two approaches
        ket_list = sorted(list(set(state_list)))
        ket_sym  = [ket_obj.state_sym(ket_list[i])[0] 
                    for i in range(len(ket_list))]

        # determine what final states to include
        # -------------------------------------------------

        state_list  = []

        # iterate through the final states, adding roots
        # to various irreps as needed -- these states are 
        # ordered by adiabatic energy
        nfinal_total = bra_obj.n_state()
        if final_states is not None:
            # iterate through the final states, adding roots
            # to various irreps as nee
            for fstate in final_states:
                if fstate < nfinal_total:
                    state_list.append(fstate)

        # use the istate array (format = irrep.state) to select
        # individual states
        if final_states_sym is not None:
            for state in final_states_sym:
                irrep   = int(state)
                st      = int(10*state) - 10*irrep
                fstate  = bra_obj.state_index(irrep,st)
                if fstate is not None:
                    state_list.append(fstate)

        # take the union of the two approaches
        bra_list = sorted(list(set(state_list)))
        bra_sym  = [bra_obj.state_sym(bra_list[i])[0]
                    for i in range(len(bra_list))]

        # create the pair list -- we build the list this way 
        # in order to facilitate call to the tdm code in bitci
        # which computes tmds for pairs of states in particular
        # irreps
        nirr_bra = bra_obj.n_irrep()
        nirr_ket = ket_obj.n_irrep()

        trans_list_sym = [[[] for irr_bra in range(nirr_bra)]
                          for irr_ket in range(nirr_ket)]
        
        # trans list is a list of pairs, with first dimension 
        # equal to the number of initial states
        trans_list     = []

        for k_st in ket_list:
            ksym = ket_obj.state_sym(k_st)

            for b_st in bra_list:
                bsym = bra_obj.state_sym(b_st)
 
                # exclude bra == ket pairs: 'lower_i>j'
                if list_type == 'lower_i>j' and k_st == b_st:
                    continue

                # only include unique pairs: 'lower' and 'lower_i>j'
                if list_type != 'full':
                    try:
                        indx = self.trans_list.index([b_st, k_st])
                    except:
                        indx = None
                    if indx is not None:
                        continue

                # add the state pair
                trans_list_sym[bsym[0]][ksym[0]].append([bsym[1], ksym[1]])
                trans_list.append([b_st, k_st])
                
        return ket_list, trans_list, trans_list_sym
