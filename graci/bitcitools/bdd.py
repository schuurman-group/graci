"""
Module for the calculation of an adiabatic-to-diabatic transformation
matrix using the propagative block diagonalisation diabatisation
scheme
"""

import sys as sys
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.core.molecule as molecule

@timing.timed
def bdd(ref_obj, disp_obj):
    """
    Calculation of the P-BDD ADT matrix
    Here, ref_obj corresponds to the previous geometry, R_n-1,
    and disp_obj to the current geometru, R_n
    """

    # (1) Get the ref_obj ADT matrix. If it doesn't exist, assume
    #     that it corresponds to R0 and set it to the unit matrix
    #
    # (2) Compute the ref-disp WF overlaps using norm_thrsh = 0.999
    #
    # (3) Compute the BDD ADT matrix
    
    print('\n Here')
    sys.exit()
    
    return
