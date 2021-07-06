"""
A set of standard structures to store basis set
of MO information
"""
import os
import numpy as np
from pyscf.tools import molden

def write_orbitals(file_name, mol, orbs, 
                    occ=None, sym_lbl=None, ener=None, cart=True):
    """print the orbitals to molden dat file format. Code assumes
       input orbs are in pyscf format. Default is to output orbitals in
       cartesian AOs"""

    molden.from_mo(mol.pymol(), file_name, orbs, 
            symm=sym_lbl, occ=occ, spin='Alpha', ene=ener, 
            ignore_h=True)

    return

