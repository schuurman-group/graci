"""
module for compute spin-orbit coupling matrix 
elements
"""


class Spinorbit:
    """spin orbit coupling class"""
    def __init__(self):
        self.bra_states = None
        self.ket_states = None
        self.bra_wfn    = None
        self.ket_wfn    = None
        self.label      = 'default'


    def name(self):
        """ return the name of the class object as a string"""
        return 'spinorbit'


    
