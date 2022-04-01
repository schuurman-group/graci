"""
Overlap class
"""

import sys as sys
import graci.utils.timing as timing
import graci.interaction.interaction as interaction
import graci.bitcitools.bitwf_init as bitwf_init
import graci.io.output as output

class Overlap(interaction.Interaction):
    """
    Calculation of overlaps, Dyson orbitals and ADT matrices using
    the bitwf library
    """
    def __init__(self):
        super().__init__()

        # user defined quantities
        self.calc       = None
        self.label      = 'Overlap'

        self.bra_label  = None
        self.ket_label  = None

        self.bra_states = None
        sekf.ket_states = None

    #
    @timing.timed
    def run(self):
        """
        computes the overlaps/Dyson orbitals/ADT matrix elements
        between the bra and ket states
        """

        # section header
        output.print_overlap_header(self.label)
        
        # check on the calculation type
        allowed_calcs = ['overlap', 'dyson', 'adt']
        if self.calc not in allowed_calcs:
            print('\n Unsupported calc type in overlap: '+str(self.calc)
                  +'\n Allowed keywords: \n'+str(allowed_calcs))
            sys.exit()

        # initialise the bitwf library
        bitwf_init.init(self.bra_obj, self.ket_obj, self.calc)
       
        # construct the list state pairs between which we
        # will compute overlaps
        self.add_group('bra', self.bra_obj, self.bra_states)
        self.add_group('ket', self.ket_obj, self.ket_states)

        # if bra/ket are the same, only consider unique pairs
        if self.same_obj(self.bra_obj, self.ket_obj):
            pair_blks = 'nodiag'
        else:
            pair_blks = 'full'

        s_list     = self.build_pair_list('bra', 'ket', 
                                           pairs=pair_blks)
        s_list_sym = self.build_pair_list('bra', 'ket', 
                                           pairs=pair_blks, 
                                           sym_blk=True)


        sys.exit()
