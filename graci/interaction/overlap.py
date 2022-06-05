"""
Overlap class
"""

import sys as sys
import numpy as np
import graci.utils.timing as timing
import graci.interaction.interaction as interaction
import graci.bitcitools.bitwf_init as bitwf_init
import graci.bitcitools.wf_overlap as wf_overlap
import graci.io.output as output

class Overlap(interaction.Interaction):
    """
    Calculation of overlaps, Dyson orbitals and ADT matrices using
    the bitwf and overlap libraries
    """
    def __init__(self):
        super().__init__()

        # user defined quantities
        self.calc         = None
        self.label        = 'Overlap'
        self.init_label   = None
        self.final_label  = None
        self.init_states  = None
        self.final_states = None
        self.norm_thresh  = 0.995

        # ----------------------------------------------------------
        # internal class variables -- should not be accessed
        # directly
        self.allowed_calcs = ['overlap', 'dyson', 'adt']
        self.bra_obj       = None
        self.ket_obj       = None
        self.bra_wfunit    = None
        self.ket_wfunit    = None
        self.overlaps      = None
        
    #
    @timing.timed
    def run(self, obj_list):
        """
        computes the overlaps/Dyson orbitals/ADT matrix elements
        between the bra and ket states
        """

        # set the bra and ket objects
        self.bra_obj, self.ket_obj = self.set_objs(obj_list)

        # section header
        output.print_overlap_header(self.label)
        
        # check on the calculation type
        self.check_calc()
        
        # initialise the bitwf library
        bitwf_init.init(self.bra_obj, self.ket_obj, self.calc)

        # if bra and ket are the same object, only compute the unique
        # overlaps
        list_type = 'full'
        if self.same_obj(self.bra_obj, self.ket_obj):
            list_type = 'nodiag'

        # construct the list of state pairs between which we
        # will compute overlaps
        self.add_group('ket', self.ket_obj, self.init_states)
        self.add_group('bra', self.bra_obj, self.final_states)
        self.trans_list = self.build_pair_list('bra',
                                               'ket',
                                               pairs=list_type)
        self.trans_list_sym = self.build_pair_list('bra',
                                                   'ket',
                                                   pairs=list_type,
                                                   sym_blk=True)
        
        # compute the wave function overlaps
        self.overlaps = self.build_overlaps(self.bra_obj, self.ket_obj,
                                            self.trans_list,
                                            self.trans_list_sym)

        # output the overlaps
        output.print_overlaps(self.trans_list, self.overlaps,
                              self.bra_obj, self.ket_obj)
                
        # finalize the bitwf library
        bitwf_init.finalize()
         
#----------------------------------------------------------------------
# "Private" class methods
#

    #
    def set_objs(self, obj_list):
        """
        Determines the bra and ket CI objects based on their labels
        """

        bra_obj = None
        ket_obj = None

        if obj_list[0].label == self.init_label:
            bra_obj = obj_list[0]
            ket_obj = obj_list[1]
        else:
            ket_obj = obj_list[0]
            bra_obj = obj_list[1]
            
        return bra_obj, ket_obj

    #
    def check_calc(self):
        """
        Checks whether the requested calculation type is supported
        """

        # check the calculation type
        if self.calc not in self.allowed_calcs:
            print('\n Unsupported calc type in overlap: '+str(self.calc)
                  +'\n Allowed keywords: \n'+str(self.allowed_calcs))
            sys.exit()

        # bitwf currently only supports equal bra and ket point groups
        if self.bra_obj.scf.mol.sym_indx != self.ket_obj.scf.mol.sym_indx:
            sys.exit('Error: unequal bra and ket point groups')
        
        return

    #
    def build_overlaps(self, bra, ket, trans_list, trans_list_sym):
        """
        grab the wave function overlaps from bitwf and then reshape the
        list of these into a more usable format
        """

        # extract the determinant representation of the wave functions
        self.bra_wfunit, self.ket_wfunit = wf_overlap.extract(bra, ket)

        # compute the overlaps
        overlap_list = wf_overlap.overlap(bra, ket,
                                          self.bra_wfunit,
                                          self.ket_wfunit,
                                          trans_list_sym,
                                          self.norm_thresh)

        # make the overlap list
        npairs   = len(trans_list)
        overlaps = np.zeros((npairs), dtype=float)

        for indx in range(len(trans_list)):
            bk_st          = trans_list[indx]
            [birr, bst]    = bra.state_sym(bk_st[0])
            [kirr, kst]    = ket.state_sym(bk_st[1])
            sym_indx       = trans_list_sym[birr][kirr].index([bst,kst])
            overlaps[indx] = overlap_list[birr][kirr][sym_indx]
        
        return overlaps
