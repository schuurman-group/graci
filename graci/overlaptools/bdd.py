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
import scipy as sp

@timing.timed
def adt(ref_obj, disp_obj):
    """
    Calculation of the P-BDD ADT matrices for all irreps
    Here, ref_obj corresponds to the previous geometry, R_n-1,
    and disp_obj to the current geometru, R_n
    """
    
    # number of irreps
    nirr = ref_obj.n_irrep()

    # initialise the list of disp ADT matrices
    disp_obj.adt = []
    
    # loop over irreps
    for irr in range(nirr):
    
        # number of states
        nstates = ref_obj.n_states_sym(irr)
    
        # get the ref_obj ADT matrix. If it doesn't exist, assume
        # that it corresponds to R0 and set it to the unit matrix
        if ref_obj.adt is None:
            adt_ref = np.eye(nstates)
        else:
            adt_ref = ref_obj.adt[irr]
        
        # number of ref and disp MOs
        nmo_ref  = ref_obj.scf.nmo
        nmo_disp = disp_obj.scf.nmo

        # dimensions
        n_int_ref   = ref_obj.det_strings[irr].shape[0]
        n_int_disp  = disp_obj.det_strings[irr].shape[0]
        ndet_ref    = ref_obj.det_strings[irr].shape[2]
        ndet_disp   = disp_obj.det_strings[irr].shape[2]
        nroots_ref  = ref_obj.vec_det[irr].shape[1]
        nroots_disp = disp_obj.vec_det[irr].shape[1]

        # check on the number of roots
        if nroots_ref != nroots_disp:
            sys.exit('\n ERROR: nroots_ref != nroots_disp in bdd.bdd')
            
        # determinant bit strings
        det_ref  = np.reshape(ref_obj.det_strings[irr],
                              (n_int_ref*2*ndet_ref), order='F')
        det_disp = np.reshape(disp_obj.det_strings[irr],
                              (n_int_disp*2*ndet_disp), order='F')

        # eigenvectors in the determinant basis
        vec_ref  = np.reshape(ref_obj.vec_det[irr],
                              (ndet_ref*nroots_ref), order='F')
        vec_disp = np.reshape(disp_obj.vec_det[irr],
                              (ndet_disp*nroots_disp), order='F')

        # ref-diso MO overlaps
        smo = np.reshape(disp_obj.smo, (nmo_ref * nmo_disp), order='F')

        # norm-based wave function truncation threshold
        norm_thresh = 0.999

        # frozen core MO information
        ncore_el    = np.array([molecule.atom_ncore[n] for n in
                                [molecule.atom_name.index(m)
                                 for m in ref_obj.scf.mol.asym]])
        ncore       = int(np.sum(ncore_el/2))
        icore       = np.arange(0, ncore, 1) + 1
        delete_core = True

        # ref-disp wave function overlaps
        npairs = nroots_ref * nroots_disp
        Sij    = np.zeros((npairs), dtype=float)
        pairs  = np.array([[i,j] for i in range(nroots_ref)
                           for j in range(nroots_disp)],
                          dtype=int)
        pairs  = np.reshape(1+pairs, (npairs*2), order='F')

        # compute the ref-disp wave function overlaps
        args = (nmo_ref, nmo_disp, n_int_ref, n_int_disp, ndet_ref,
                ndet_disp, nroots_ref, nroots_disp, det_ref, det_disp,
                vec_ref, vec_disp, smo, norm_thresh, ncore, icore,
                delete_core, npairs, Sij, pairs)
        Sij  = libs.lib_func('overlap_c', args)

        # reshape the wave function overlaps into a more
        # useful form
        Sij = np.reshape(Sij, (nroots_ref, nroots_disp))

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
