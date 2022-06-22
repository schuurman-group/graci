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
import graci.io.output as output
import graci.overlaptools.overlap as overlap
import scipy as sp

@timing.timed
def adt(ref_obj, disp_obj):
    """
    Calculation of the P-BDD ADT matrices for all irreps
    Here, ref_obj corresponds to the previous geometry, R_n-1,
    and disp_obj to the current geometru, R_n
    """

    # section header
    output.print_bdd_header()
    
    # number of irreps
    nirr = ref_obj.n_irrep()

    # initialise the list of disp ADT matrices
    disp_obj.adt = []
    
    # loop over irreps
    for irr in range(nirr):

        # check on the number of roots
        nref  = ref_obj.vec_det[irr].shape[1]
        ndisp = disp_obj.vec_det[irr].shape[1]
        if nref != ndisp:
            sys.exit('\n ERROR: nroots_ref != nroots_disp in bdd.bdd')
            
        # get the ref_obj ADT matrix. If it doesn't exist, assume
        # that it corresponds to R0 and set it to the unit matrix
        nstates = nref
        if ref_obj.adt is None:
            adt_ref = np.eye(nstates)
        else:
            adt_ref = ref_obj.adt[irr]
        
        # wave function truncation threshold
        norm_thresh = 0.999

        # pairs of states for which we want overlaps
        npairs = nstates**2
        pairs  = np.array([[i,j] for i in range(nstates)
                           for j in range(nstates)], dtype=int)

        # get the overlaps
        Sij = overlap.overlap(ref_obj, disp_obj, disp_obj.smo,
                              pairs, irr, norm_thresh,
                              disp_obj.verbose)

        # reshape the wave function overlaps into a more
        # useful form
        Sij = np.reshape(Sij, (nstates, nstates))

        # transform using the ref ADT matrix
        Sij = np.matmul(adt_ref.T, Sij)
        
        # compute the disp BDD ADT matrix
        S_inv    = np.linalg.inv(Sij)
        SST      = np.matmul(Sij, Sij.T)
        SST_sqrt = sp.linalg.sqrtm(SST)
        adt      = np.matmul(S_inv, SST_sqrt)

        # save the disp ADT matrix
        disp_obj.adt.append(adt)
        
    return
