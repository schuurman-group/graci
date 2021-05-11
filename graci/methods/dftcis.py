"""
Module for computing DFT/MRCI energies
"""


class Dftcis:
    """Class constructor for SCF object"""
    def __init__(self):
        self.label   = 'dftcis'
        self.mol     = None
        self.scf     = None


#################################################################

    def set_mol(self, mol):
        """set the mol object for the dftcis method"""
        self.mol = mol
        return

    def set_scf(sel, scf):
        """set the scf object for the dftcis method"""
        self.scf = scf
        return

    def name(self):
        """ return the name of the class object as a string"""
        return 'dftcis'

    def run():
        """runs a dftcis calculation"""

        return

    def rdm1():
        """return the density matrix for state i"""

        return

    def slater_dets(irrep, state):
        """return the slater determinant list for state 'state'"""

        return

    def csfs(irrep, state):
        """csf list for state 'state'"""

        return

