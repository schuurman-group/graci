"""
The Parameterize object and its associated functions.
"""
import numpy as np
import graci.core.driver as driver
import graci.io.parse


class Parameterize:
    """Class constructor for the Parameterize object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.algorithm = None
        self.label     = 'parameterize'


    def run(self):
        """re-parameterize the Hamiltonian"""

        return




