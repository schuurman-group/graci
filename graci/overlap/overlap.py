"""
Overlap class
"""

import sys as sys
import graci.utils.timing as timing
import graci.bitcitools.bitwf_init as bitwf_init
import graci.io.output as output

class Overlap:
    """
    Calculation of overlaps, Dyson orbitals and ADT matrices using
    the bitwf library
    """
    def __init__(self):

        # user defined quantities
        self.bra_states = None
        self.ket_states = None
        self.bra_label  = None
        self.ket_label  = None
        self.calc       = None
        self.label      = 'Overlap'

        # global variables
        self.class_name = 'overlap'
        # method object for the bra states
        self.bra_obj = None
        # method object for the ket states
        self.ket_obj = None
        # comprehensive list of all matrix elements
        # to be computed
        self.trans_list = []
        # list of lists of matrix elements, grouped by irrep
        # (format for bitwf)
        self.trans_list_sym = []
        # convenient to just keep track of initial states for each
        # matrix element
        self.ket_list       = []
        # simple boolean that is true if bra and ket objects
        # are identical. If this is true, default behavior is 
        # to only consider unique bra/ket pairs
        self.braket_iden    = False


    #
    def name(self):
        """
        returns the name of the class as a string
        """
        return self.class_name

    
    #
    def set_ket(self, ket_method):
        """set the method used to compute the ket state(s)"""
        try:
            self.ket_label = ket_method.label
        except:
            self.ket_label = None
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
            self.bra_label = bra_method.label
        except:
            self.bra_label = None
            return

        self.bra_obj = bra_method
        return

    #
    def bra_exists(self):
        """return True if self.bra_obj is not None"""

        return self.bra_obj is not None
    
    
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
        
        sys.exit()
