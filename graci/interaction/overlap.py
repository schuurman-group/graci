"""
Overlap class
"""

import sys as sys
import numpy as np
import copy
import graci.utils.timing as timing
import graci.interaction.interaction as interaction
import graci.interfaces.bitci.bitwf_init as bitwf_init
import graci.interfaces.bitci.wf_overlap as wf_overlap
import graci.io.output as output

class Overlap(interaction.Interaction):
    """
    Calculation of wave functions overlaps using the bitwf and
    overlap libraries
    """
    def __init__(self):
        super().__init__()

        # user defined quantities
        self.bra_label    = None
        self.ket_label    = None
        self.bra_states   = None
        self.ket_states   = None
        self.norm_thresh  = 0.999
        self.representation = 'adiabatic'

        # ----------------------------------------------------------
        # internal class variables -- should not be accessed
        # directly
        self.bra_wfunit    = None
        self.ket_wfunit    = None
        self.overlaps      = None
        self.allowed_reps  = ['adiabatic']        

    def copy(self):
        """create of deepcopy of self"""
        new = Overlap()

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))

        return new

    #
    @timing.timed
    def run(self, bra, ket):
        """
        computes the overlap matrix elements between the
        bra and ket states
        """
        # section header
        if self.verbose:
            output.print_overlap_header(self.label)

        # sanity check on the representation
        self.check_representation()

        # set the bra/ket objects and add the state groups associated
        # with each 
        self.bra_obj = bra
        self.ket_obj = ket
        self.add_group('bra', [bra], states = [self.bra_states])
        self.add_group('ket', [ket], states = [self.ket_states])

        # check on the requested calculation
        self.check_bra_ket()

        # if bra and ket are the same object, only compute the unique
        # overlaps 
        ltype = 'full'
        if self.same_group('bra', 'ket'):
            ltype = 'lower'

        # initialize the list of overlaps to compute
        self.trans_list = self.grp_pair_list('bra', 'ket', pairs=ltype)
        ntrans = len(self.trans_list)
        self.overlaps = np.zeros((ntrans), dtype=complex)      

        # loop over all pairs of spin free states in both the bra
        # and ket objects. In the future, may be prudent to check
        # init and final states and just which states are necessary
        for k_lbl in self.get_ci_lbls('ket'):
            ket_ci            = self.get_ci_obj('ket', k_lbl)
            k_mult, k_s, k_m  = self.get_ci_spins('ket', k_lbl)

            for b_lbl in self.get_ci_lbls('bra'):
                bra_ci            = self.get_ci_obj('bra', b_lbl)
                b_mult, b_s, b_m  = self.get_ci_spins('bra', b_lbl)

                # if the same object, only build lower diagonal
                pair_type = 'full'
                if self.same_obj(ket_ci, bra_ci):
                    pair_type = 'lower'

                # initialise the bitwf library
                bitwf_init.init(bra_ci, ket_ci, 'overlap', self.verbose)

                # this state pair list for the current pair of CI 
                # objects, stored by adiabatic label
                ci_tran = self.ci_pair_list('bra', b_lbl, 
                                            'ket', k_lbl,
                                             pairs=pair_type)

                # the state pairs, ordered by interacting irreps, is how
                # we call bitwf
                ci_tran_sym = self.ci_pair_list('bra', b_lbl,
                                                'ket', k_lbl,
                                                 pairs=pair_type,
                                                 sym_blk=True)

                # overlaps over the spin-free CI states
                ci_ovrlp =  self.build_ci_overlaps(b_lbl, k_lbl,
                                                  ci_tran, ci_tran_sym)

                # rotate overlaps to final state basis
                self.build_grp_tensor('bra', b_lbl,
                                      'ket', k_lbl,
                                       ci_tran, ci_ovrlp,
                                       self.overlaps)

                # finalize the bitwf library
                bitwf_init.finalize()


        # output the overlaps
        bra_state_sym = [self.bra_obj.state_sym(n)
                         for n in range(self.bra_obj.n_states())]
        ket_state_sym = [self.ket_obj.state_sym(n)
                         for n in range(self.ket_obj.n_states())]

        if self.verbose:
            output.print_overlaps(self.trans_list, self.overlaps,
                                  self.ket_label, self.bra_label,
                                  self.bra_obj.scf.mol.irreplbl,
                                  bra_state_sym, ket_state_sym)

    #
    def get_overlap(bra_st=None, ket_st=None, spin_state=False):
        """
        Return requested overlaps.
 
        If spin_state is False, real overlaps are assume and
        only the real component is returned. If spin_state
        is True, the complex overlap is returned.
        """
        if spin_state:
            olap = self.overlaps
        else:
            olap = self.overlaps.real


        if bra_st is None and ket_st is None:
            return olap
        else:
            if bra_st is None:
                chk_st  = ket_st
                chk_ind = 1
            else:
                chk_st  = bra_st
                chk_ind = 0
            lst = [self.trans_list.index(pair)
                  for pair in self.trans_list if pair[chk_ind]==chk_st]
            return olap[lst]


                
#----------------------------------------------------------------------
# "Private" class methods
#

    #
    def check_bra_ket(self):
        """
        Sanity check on the requested calculation
        """

        # make sure that we have equal numbers of bra and ket
        # electrons
        for b_lbl in self.get_ci_lbls('bra'):
            for k_lbl in self.get_ci_lbls('ket'):

                # make sure that we have equal numbers of bra and ket
                # electrons
                if (self.get_ci_obj('bra', b_lbl).nel != 
                        self.get_ci_obj('ket',k_lbl).nel):
                    sys.exit('Error: unequal bra and ket N_el')
        
                # bitwf currently only supports equal bra and ket point groups
                if (self.get_ci_obj('bra', b_lbl).scf.mol.sym_indx != 
                        self.get_ci_obj('ket',k_lbl).scf.mol.sym_indx):
                    sys.exit('Error: unequal bra and ket point groups')                

        return

    #
    def build_ci_overlaps(self, b_lbl, k_lbl, ci_trans, ci_trans_sym):
        """
        Grab the wave function overlaps from bitwf and then reshape the
        list of these into a more usable format
        """

        # get ci objects
        bra = self.get_ci_obj('bra', b_lbl)
        ket = self.get_ci_obj('ket', k_lbl)

        # extract the determinant representation of the wave functions
        self.bra_wfunit, self.ket_wfunit = wf_overlap.extract(bra, ket)

        # compute the overlaps
        overlap_list = wf_overlap.overlap(bra, ket,
                                          self.bra_wfunit,
                                          self.ket_wfunit,
                                          ci_trans_sym,
                                          self.norm_thresh)
        # make the overlap list
        npairs   = len(ci_trans)
        overlaps = np.zeros((npairs), dtype=float)

        for indx in range(npairs):
            bk_st          = ci_trans[indx]
            [bir, bst]     = self.get_ci_sym('bra', b_lbl, bk_st[0])
            [kir, kst]     = self.get_ci_sym('ket', k_lbl, bk_st[1])
            sym_indx       = ci_trans_sym[bir][kir].index([bst, kst])
            overlaps[indx] = overlap_list[birr][kirr][sym_indx]
                
        return overlaps

