"""
module for compute moments of the dipole operator
for a given electronic state
"""


class Moments:
    """Moment class for determing permanent and transition moments"""
    def __init__(self):
        self.bra_states = None
        self.ket_states = None
        self.bra_wfn    = None
        self.ket_wfn    = None
        self.label      = 'default'

    def name(self):
        """ return the name of the class object as a string"""
        return 'moments'


    def dipole(self):
        """return the dipole moments for states in 'states'"""


        return



    def tr_dipole(self):
        """return the transition dipole moments between the bra and
           ket states"""


        return


    def quadrupole(self):
        """return second moments for states in 'states' array"""


        return
