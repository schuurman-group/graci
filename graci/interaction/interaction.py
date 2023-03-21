"""
Parent state interaction class
"""
import sys as sys
import numpy as np
import math as math
import graci.core.params as params
import graci.utils.timing as timing

class Cigroup:
    """ 
     Class to hold information on a group of states: the 
     ci method and the number of CI states, etc.
    """
    def __init__(self, ci_objs, ci_states=None, spin_states=False,
                                                       lbl='default'):
        self.label      = lbl
        self.grp_objs   = []
        self.grp_states = None
        self.grp_vecs   = None
        self.grp_energy = None
        self.grp_spin_states = spin_states
        self.ci_lbls    = []
        self.ci_objs    = {}
        self.ci_mult    = {}
        self.ci_ener    = {}
        self.ci_S       = {}
        self.ci_M       = {}
        self.ci_st      = {}
        self.ci_st_irr  = {}
        self.ci_st_sym  = {}
        self.n_ci       = 0
        self._nmult     = {}    

        # if this is a ci method, create a new group associated
        # with the CI states
        if ci_states is None:
            ci_states = [None] * len(ci_obj)

        for ci in range(len(ci_objs)):

            self.grp_objs.append(ci_objs[ci])

            # if this is a cimethod object, just add the ci object
            # to this group
            if type(ci_objs[ci]).__name__ in params.ci_objs:
                self._add_ciobj(ci_objs[ci], ci_states[ci])

            # if this is a postci object, should already extend
            # interaction and have groups set up. Add all of them
            # here
            elif type(ci_objs[ci]).__name__ in params.postci_objs:
                grp_names = list(ci_objs[ci].groups.keys())
                if len(grp_names) > 1:
                    msg = 'Cannot construct Cigroup from complex postci_obj'
                    sys.exit(msg)
                grp = grp_names[0]

                # add the ci objects from the passed group as they 
                # are      
                for lbl in ci_objs[ci].get_ci_lbls(grp):
                    self._add_ciobj(ci_objs[ci].get_ci_obj(grp, lbl),
                                    ci_objs[ci].get_ci_states(grp, lbl))

                # only add those group vectors requested by ci_states
                n_g_states = len(ci_objs[ci].get_energy(grp))
                if ci_states[ci] is None:
                    g_states = np.array(range(n_g_states), dtype=int)
                else:
                    if max(ci_states[ci]) > n_g_states:
                        sys.exit('Too many states requested from ' + 
                                  str(ci_objs[ci].label))
                    g_states = ci_states[ci]
                  
                self.grp_vecs = ci_objs[ci].get_vectors(grp)[:,g_states]
                self.grp_states = ci_objs[ci].get_states(grp)[g_states]
                self.grp_energy = ci_objs[ci].get_energy(grp)[g_states]
                self.grp_spin_states = \
                            ci_objs[ci].groups[grp].grp_spin_states

        # if the passed object does not define the group states/vectors,
        # this will default to: 
        # grp_states = CI states, grp_energies = CI energies 
        if self.grp_spin_states:
            self._nmult = {lbl: len(self.ci_M[lbl]) 
                                               for lbl in self.ci_lbls}
        else:
            self._nmult = {lbl: 1 for lbl in self.ci_lbls}

        if self.grp_vecs is None or self.grp_states is None:
            # ntot = total number of states, including M_s values if this 
            # group is holding spin states
            n_tot = sum([len(self.ci_st[lbl])*self._nmult[lbl]
                                          for lbl in self.ci_lbls])
            self.grp_vecs   = np.identity(n_tot, dtype=float)

            # THIS FEELS LIKE A HACK...
            # if there is a single CI object, and it's not a 
            # spin-state object, the 'grp_states' will be the 
            # same as the CI state indices.
            if len(self.ci_lbls)==1 and not self.grp_spin_states:
                self.grp_states = self.ci_st[self.ci_lbls[0]]
            else:
                self.grp_states = np.array(range(n_tot), dtype=int)    

        if self.grp_energy is None:
            # set up all energies, including degenerate m_s components
            # if group is holding spin states
            ener = []
            for lbl in self.ci_lbls:
                for ien in range(len(self.ci_ener[lbl])):
                    ener.extend([self.ci_ener[lbl][ien]]*self._nmult[lbl])
            ener.sort()
            self.grp_energy = np.array(ener, dtype=float)
 
    #
    def _add_ciobj(self, ciobj, states):
        """ 
        Add a cimethod object to the group list
        """
        lbl = ciobj.label
        while lbl in self.ci_lbls:
            lbl += 'x'

        if states is None:
            st = np.array(range(ciobj.n_states()), dtype=int)
        else:
            st = np.array(states, dtype=int)

        self.ci_lbls.append(lbl)
        self.ci_objs[lbl]   = ciobj
        self.ci_mult[lbl]   = ciobj.mult
        self.ci_st[lbl]     = st
        self.ci_ener[lbl]   = np.array([ciobj.energy(st[i]) 
                                    for i in range(len(st))])
        self.ci_st_irr[lbl] = np.array([ciobj.state_sym(st[i])[0] 
                                    for i in range(len(st))], dtype=int)
        self.ci_st_sym[lbl] = np.array([ciobj.state_sym(st[i])[1] 
                                    for i in range(len(st))], dtype=int)
        S                   = (ciobj.mult - 1.)/2.
        self.ci_S[lbl]      = S
        self.ci_M[lbl]      = np.linspace(-S, S, num=round(2*S+1))

        self.n_ci          += 1
        return

    def get_grp_indx(self, ci_lbl, ci_st, ci_M=None):
        """
        return the index in the grp_vec that the ci_lbl/ci_st pair 
        correspond to
        """
        if ci_lbl not in self.ci_lbls:
            return None
        else:
            (sind, ) = np.where(self.ci_st[ci_lbl] == ci_st)

            if len(sind) == 0:
                return None

            ciind    = self.ci_lbls.index(ci_lbl)
            cilen    = [len(self.ci_st[self.ci_lbls[i]]) 
                                  for i in range(ciind)]
            nci      = [self._nmult[self.ci_lbls[i]]
                                  for i in range(ciind+1)]
            off      = sum([cilen[i]*nci[i]  for i in range(ciind)])
           
            # if M is not specified, return all valid Ms indices
            if ci_M is None:
                return [off + sind[0]*nci[ciind] + i 
                                    for i in range(nci[ciind])]
            else:
                (mind,) = np.where(self.ci_M[ci_lbl] == ci_M)
                return off + sind[0]*nci[ciind] + mind[0]

    def get_ci_indx(self, grp_ind):
        """
        return the ci lbl and ci state index corresponding to the 
        group vec index grp_ind
        """

        # if we don't have any group vectors, or, the index is out
        # bounds, return None
        if self.grp_vecs is None or grp_ind > self.grp_vecs.shape[0]:
            return None
        else:
            i     = 0
            cind  = grp_ind
            while i < len(self.ci_lbls):
                ci_len = len(self.ci_st[self.ci_lbls[i]]) * \
                             self._nmult[self.ci_lbls[i]]

                if cind - ci_len >= 0:
                    cind -= ci_len
                    i    += 1
                else:
                    break
            
            lbl  = self.ci_lbls[i]
            sind = int(math.floor(cind / self._nmult[lbl]))
            if self.grp_spin_states:
                mval = self.ci_M[lbl][cind - sind*self._nmult[lbl]]
            else:
                mval = None

            return lbl, self.ci_st[lbl][sind], mval          

#
class Interaction:
    """Parent state interaction class to be inherited by any module
       computing matrix elements <psi_i|O|psi_j> for some operator O"""
    def __init__(self):

        # all state interaction sub-classes must have a 
        # bra object and a ket object

        # user set variables that are common to
        # _all_ child classes
        #-----------------------------------------
        self.representation = 'adiabatic'
        
        # Internal class variables (should not be accessed directly)
        #-----------------------------------------------------------
        # state list -- list of groups of states
        self.states         = {}
        # the objects that hold the ci states
        self.groups         = {}
        self.verbose        = True
        self.label          = 'default'
        self.allowed_reps   = ['adiabatic']

    #
    def add_group(self, grp_name, ci_objs, 
                                 states=None, spin=False, lbl='default'):
        """append a list of states and corresponding object"""

        self.groups[grp_name] = Cigroup(ci_objs, 
                                        ci_states=states, 
                                        spin_states=spin, 
                                        lbl=lbl)

        return

    #
    def check_representation(self):
        """checks that the user-specified representation is
        supported"""

        if self.representation not in self.allowed_reps:
            sys.exit('\n ERROR: unsupported representation: '
                     +self.representation
                     +'\n Allowed values: '+str(self.allowed_reps))

        return

    def same_group(self, grp_name1, grp_name2):
        """return true if the bra and ket objects are the same"""

        nci1 = self.groups[grp_name1].n_ci
        nci2 = self.groups[grp_name2].n_ci

        ci_lbl1 = self.get_ci_lbls(grp_name1)
        ci_lbl1.sort()
        ci_lbl2 = self.get_ci_lbls(grp_name2)
        ci_lbl2.sort()

        if nci1 != nci2 or ci_lbl1 != ci_lbl2:
            return False

        is_same = True
        for ci_lbl in ci_lbl1:

            obj1 = self.get_ci_obj(grp_name1, ci_lbl)
            obj2 = self.get_ci_obj(grp_name2, ci_lbl)

            is_same = is_same and self.same_ci_obj(obj1, obj2)

        return is_same 

    #
    def n_groups(self):
        """return the number of entries in the state group 
           dictionary"""

        return len(self.groups.keys())

    #
    def get_groups(self):
        """return a list of the group labels"""
        return list(self.groups.keys())
   
    #
    def get_group_objs(self, grp_name):
        """
        return the group object: could be cimethod or interaction 
        object
        """
        return self.groups[grp_name].grp_objs

    #
    def is_spin_states(self, grp_name):
        """ 
        return true if the group contains spin-states
        """
        return self.groups[grp_name].grp_spin_states

    #
    def get_states(self, grp_name):
        """return the states and symmetries corresponding to grp lbl"""
        return self.groups[grp_name].grp_states

    #
    def get_energy(self, grp_name, state=None):
        """ return group energy """
        if state is None:
            return self.groups[grp_name].grp_energy
        else:
            (ind,) = np.where(self.groups[grp_name].grp_states==state)
            return self.groups[grp_name].grp_energy[ind[0]]

    #
    def get_vectors(self, grp_name):
        """ return group vectors """
        return self.groups[grp_name].grp_vecs

    #
    def get_vector(self, grp_name, g_state):
        """return the gropu vector"""
        (ind,) = np.where(self.groups[grp_name].grp_states == g_state)
        return self.groups[grp_name].grp_vecs[:, ind[0]]

    #
    def get_syms(self, grp_name):
        """ return the symmetries of the grp states """

        # currently can't handle symmetry of SO objects: return
        # C1 symmetry labels
        if self.is_spin_states(grp_name):
            syms = [0]*len(self.get_energy(grp_name))
            return np.array(syms, dtype=int)

        # if group is a single set of spin-free states, return
        # the symmetries of the ci states
        elif len(self.get_ci_lbls(grp_name)) == 1:
            lbls = self.get_ci_lbls(grp_name)
            return self.groups[grp_name].ci_st_irr[lbls[0]]

        # else, we don't have a way to assign symmetry of group
        # states, return None
        else:
            return None

    # 
    def get_sym_lbl(self, grp_name):
        """ return the symmetry lables of the grp states """

        # currently can't handle symmetry of SO objects: return
        # C1 symmetry labels
        if self.is_spin_states(grp_name):
            return ['A']*len(self.get_states(grp_name))

        # if group is a single set of spin-free states, return
        # the symmetries of the ci states
        elif len(self.get_ci_lbls(grp_name)) == 1:
            lbls    = self.get_ci_lbls(grp_name)
            irr_lbl = self.get_ci_obj(grp_name,lbls[0]).scf.mol.irreplbl
            st_sym  = self.get_syms(grp_name)
            return [irr_lbl[st_sym[i]] for i in range(len(st_sym))]
            
        # else, we don't have a way to assign symmetry of group
        # states, return None
        else:
            return None

    #
    def set_group_states(self, grp_name, grp_vectors, grp_energy): 
        """ set the group vectors -- these may correspond to states
            that are linear combination of the CI states.
        """

        grp_states = np.array(range(grp_vectors.shape[1]), dtype=int)

        if grp_vectors.shape[1] != grp_energy.shape[0]:
            sys.exit('set_grp_states: # of states != # of energies')

        self.groups[grp_name].grp_vecs   = grp_vectors
        self.groups[grp_name].grp_states = grp_states
        self.groups[grp_name].grp_energy = grp_energy
        return

    #       
    def get_group_index(self, grp_name, lbl, state, M=None):
        """ return the index in the grp_vector corresponding 
            to the ci_obj, ci_st, ci_m value
        """
        return self.groups[grp_name].get_grp_indx(lbl, state, M)
   
    #
    def same_ci_obj(self, obj1, obj2):
        """
        check if same objects
        """
        class1 = type(obj1).__name__
        class2 = type(obj2).__name__

        lbl1  = obj1.label
        lbl2  = obj2.label

        return class1 == class2 and lbl1 == lbl2

    #
    def get_ci_index(self, grp_name, grp_index):
        """
            return the ci_lbl, ci_state, and ci_M value for a
            given index of the group vector
        """
        return self.groups[grp_name].get_ci_indx(grp_index)
 
    #
    def n_ci_objs(self, grp_name):
        """
            return the number of ci object comprising the 
            group with named 'grp_name'
        """

        return self.groups[grp_name].n_ci

    #
    def get_ci_lbls(self, grp_name):
        """
        return a list of ci_obj labels from specified group
        """
        if grp_name in self.get_groups():
            return self.groups[grp_name].ci_lbls
        else:
            return None

    #
    def get_ci_obj(self, grp_name, ci_lbl):
        """
        return the ci object from grp_name with label ci_lbl
        """
        if ci_lbl in self.get_ci_lbls(grp_name):
            return self.groups[grp_name].ci_objs[ci_lbl]
        else:
            return None

    #
    def get_ci_states(self, grp_name, ci_name):
        """
        get the states included in group 'grp_name' in ci object
        'ci_name'
        """
        return self.groups[grp_name].ci_st[ci_name]
    
    #
    def get_ci_energy(self, grp_name, ci_name, st):
        """
        get the ci energy of state st from grp grp_name 
        in ci object ci_name   
        """
        (ind,) = np.where(self.get_ci_states(grp_name, ci_name) == st)
        return self.groups[grp_name].ci_ener[ci_name][ind[0]]       
 
    #
    def get_ci_sym(self, grp_name, ci_name, st):
        """
        get the irrep and state index within the irrep
        of adiabatic state st
        """
        (ind,)  = np.where(self.get_ci_states(grp_name, ci_name) == st)
        irr  = self.groups[grp_name].ci_st_irr[ci_name][ind[0]]
        sind = self.groups[grp_name].ci_st_sym[ci_name][ind[0]]
        return [irr, sind]

    #
    def get_ci_spins(self, grp_name, ci_name):
        """return the symmetries of the states in group 'lbl'"""
        mult = self.groups[grp_name].ci_mult[ci_name]
        S    = self.groups[grp_name].ci_S[ci_name]
        M    = self.groups[grp_name].ci_M[ci_name]

        return mult, S, M 

    @timing.timed
    def grp_pair_list(self, bra_grp, ket_grp, pairs='full'):
        """
           builds an array of initial and final states, ordered by
           adiabatic energy. The symmetries of each state will be
           stored in a separate array.

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
        if (bra_grp not in list(self.groups.keys()) or
              ket_grp not in list(self.groups.keys())):
            return None

        braket_same = self.same_group(bra_grp, ket_grp)
        pair_list   = []

        for k in range(len(self.groups[ket_grp].grp_states)):
            kst = self.groups[ket_grp].grp_states[k]

            for b in range(len(self.groups[bra_grp].grp_states)):
                bst = self.groups[bra_grp].grp_states[b]

                # exclude bra == ket pairs: 'lower_i>j'
                if braket_same and bst == kst and pairs == 'nodiag':
                    continue

                # only include unique pairs if 'lower' or 'nodiag'
                if b < k and pairs != 'full':
                    continue

                pair_list.append([bst, kst])

        return pair_list

    #
    @timing.timed 
    def ci_pair_list(self, bra_grp, bra_ci, ket_grp, ket_ci, 
                                         pairs='full', sym_blk=False):
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
        if (bra_grp not in list(self.groups.keys()) or 
              ket_grp not in list(self.groups.keys())):
            return None

        # if bra and ket objects are same,
        braket_same = self.same_ci_obj(
                                self.groups[bra_grp].ci_objs[bra_ci], 
                                self.groups[ket_grp].ci_objs[ket_ci])

        # initialize the pair list. If symblk, pair_list is an nirr x nirr
        # state pairs. Else, just a list of state pairs
        if sym_blk:
            nb = self.groups[bra_grp].ci_objs[bra_ci].n_irrep()
            nk = self.groups[ket_grp].ci_objs[ket_ci].n_irrep()
            pair_list = [[[] for b in range(nb)] for k in range(nk)]
        else:
            pair_list = []

        for k in range(len(self.groups[ket_grp].ci_st[ket_ci])):
            kst = self.groups[ket_grp].ci_st[ket_ci][k]
    
            for b in range(len(self.groups[bra_grp].ci_st[bra_ci])):
                bst = self.groups[bra_grp].ci_st[bra_ci][b]
   
                # exclude bra == ket pairs: 'lower_i>j'
                if braket_same and bst == kst and pairs == 'nodiag':
                    continue

                # only include unique pairs if 'lower' or 'nodiag'
                if b < k and pairs != 'full':
                    continue
    
                # add the state pair
                if sym_blk:
                    ksym = self.groups[ket_grp].ci_st_irr[ket_ci][k]
                    kind = self.groups[ket_grp].ci_st_sym[ket_ci][k]
                    bsym = self.groups[bra_grp].ci_st_irr[bra_ci][b]
                    bind = self.groups[bra_grp].ci_st_sym[bra_ci][b]                       
                    
                    pair_list[bsym][ksym].append([bind, kind])
                else:
                    pair_list.append([bst, kst])
    
        return pair_list

    @timing.timed
    def build_grp_tensor(self, bra_grp, bra_ci, ket_grp, ket_ci,
                                            ci_pair, ci_tens, grp_tens):
        """
        Rotate a tensor in the basis of ci states into one in the
        basis of group states.
        """

        if ci_tens.shape[1:] != grp_tens.shape[1:]:
            sys.exit('ERROR tensor shape mis-match. ci_tens.shape = ' + 
                      str(ci_tens.shape[1:]) + '/ grp_tens.shape = ' + 
                      str(grp_tens.shape[1:]))

        # tensor is assumed to be in basis of ci states. Compute the 
        # contribution of each CI pair to each GROUP pair
        for i in range(len(ci_pair)):
            cpair   = ci_pair[i]
            bra_ind = self.get_group_index(bra_grp, bra_ci, cpair[0])
            ket_ind = self.get_group_index(ket_grp, ket_ci, cpair[1])

            for j in range(len(self.trans_list)):
                gpair = self.trans_list[j]
                bvec  = self.get_vector(bra_grp, gpair[0])
                kvec  = self.get_vector(ket_grp, gpair[1])

                b_cf = [bvec[b] for b in bra_ind]
                k_cf = [kvec[k] for k in ket_ind]

                cf = sum( [b*k for k in k_cf for b in b_cf] )

                if abs(cf) > 1.e-12:
                    grp_tens[j, ...] += cf * ci_tens[i, ...]

        return


