"""
Module for the calculation of wave function overlaps using
the overlap library
"""
import sys as sys
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.core.molecule as molecule

@timing.timed
def overlap(bra_obj, ket_obj, smo, pairs, irrep, norm_thresh):
    """
    Computes the overlaps between the pairs of bra and ket
    wave functions specified in the pairs list for a single
    irrep

    pairs is an numpy array of dimension (npairs,2):

    pairs[n,1] <-> n'th bra state index
    pairs[n,2] <-> n'th ket state index
    """

    # number of bra-ket wave function pairs
    npairs = pairs.shape[0]

    # number of MOs
    nmo_bra = bra_obj.scf.nmo
    nmo_ket = ket_obj.scf.nmo

    # number of dets, roots, etc.
    n_int_bra  = bra_obj.det_strings[irrep].shape[0]
    n_int_ket  = ket_obj.det_strings[irrep].shape[0]
    ndet_bra   = bra_obj.det_strings[irrep].shape[2]
    ndet_ket   = ket_obj.det_strings[irrep].shape[2]
    nroots_bra = bra_obj.vec_det[irrep].shape[1]
    nroots_ket = ket_obj.vec_det[irrep].shape[1]

    # determinant bit strings
    det_bra = np.reshape(bra_obj.det_strings[irrep],
                         (n_int_bra*2*ndet_bra), order='F')

    det_ket = np.reshape(ket_obj.det_strings[irrep],
                          (n_int_ket*2*ndet_ket), order='F')

    # eigenvectors in the determinant basis
    vec_bra = np.reshape(bra_obj.vec_det[irrep],
                         (ndet_bra*nroots_bra), order='F')
    
    vec_ket = np.reshape(ket_obj.vec_det[irrep],
                          (ndet_ket*nroots_ket), order='F')

    # reshape the bra-ket MO overlap array
    smo1 = np.reshape(smo, (nmo_bra * nmo_ket), order='F')

    # frozen core MO information
    # for now we will just assume that the core is frozen
    ncore_el    = np.array([molecule.atom_ncore[n] for n in
                            [molecule.atom_name.index(m)
                             for m in bra_obj.scf.mol.asym]])
    ncore       = int(np.sum(ncore_el/2))
    icore       = np.arange(0, ncore, 1) + 1
    delete_core = True

    # bra-ket wave function overlaps
    Sij = np.zeros((npairs), dtype=float)

    # reshape the pairs list
    # N.B. the +1 is because the overlap library uses Fortran
    # indexing
    pairs1 = np.reshape(1+pairs, (npairs*2), order='F')

    # compute the bra-ket wave function overlaps
    args = (nmo_bra, nmo_ket, n_int_bra, n_int_ket,
            ndet_bra, ndet_ket, nroots_bra, nroots_ket,
            det_bra, det_ket, vec_bra, vec_ket, smo1,
            norm_thresh, ncore, icore, delete_core, npairs,
            Sij, pairs1)
    Sij  = libs.lib_func('overlap_c', args)
    
    return Sij
