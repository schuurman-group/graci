"""
module to compute spin-orbit coupling matrix 
elements
"""

import graci.utils.timing as timing
import graci.core.libs as libs
import graci.bitcitools.bitsi_init as bitsi_init
import graci.io.output as output
import graci.utils.constants as constants
import sys as sys

class Spinorbit:
    """Spin orbit coupling class. For now will only accept calculations
       for which the 'Molecule and 'Scf' objects are the same"""
    def __init__(self):
        # user defined quanties 
        self.init_states      = None
        self.final_states     = None
        self.init_states_sym  = None
        self.final_states_sym = None
        self.all_final_states = False
        self.init_label       = None
        self.final_label      = None
        self.label      = 'spinorbit'

        # global variables
        # method object for the bra states
        self.bra_obj        = None
        # method object for the ket states
        self.ket_obj        = None
        # comprehensive list of all transition pairs
        self.trans_list     = None
        # list of transition pairs, grouped by irrep
        # (format for bitci)
        self.trans_list_sym = None
        # convenient to just keep track of initial states
        # for printing purposes
        self.ket_list       = None
        # simple boolean that is true if bra and ket objects
        # are identical. If this is true, default behavior is 
        # to only consider unique bra/ket pairs
        self.braket_iden    = False

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'spinorbit'

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
    @timing.timed
    def run(self):
        """returns the SOC matrix elements between the bra and ket
           states. If b_state and k_state are None, assume 
           SOC matrix elements from all states in method object should
           be used."""

        # Check on spin multiplicities:
        # (1) Delta S = -1, 0 or +1 must hold
        # (2) In the case of Delta S = 0, only S > 0 must hold
        S_bra = (self.bra_obj.mult-1)/2
        S_ket = (self.ket_obj.mult-1)/2
        if S_bra - S_ket not in [-1., 0., 1.]:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination' \
                     +'in spinorbit')
        if S_bra == S_ket and S_bra == 0.:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination ' \
                     +'in spinorbit')
            
        # first check to see if bra and ket are identical
        if (type(self.bra_obj).__name__ == type(self.ket_obj).__name__ 
            and self.bra_obj.label == self.ket_obj.label):
            self.braket_iden = True


        sys.exit('\n here...')
