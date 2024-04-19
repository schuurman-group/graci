"""
Module for the calculation of Dyson orbitals
"""

import sys as sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.core.molecule as molecule

@timing.timed
def dyson(bra, ket, mo_basis, n_basis, bra_wfunit,
          ket_wfunit, norm_thresh, det_thresh,
          trans_list):
    """
    Calculation of the Dyson orbitals for all pairs of states
    in trans_list using their expansion in terms of Slater
    determinants
    """

    # number of irreps
    nirr_ket = ket.n_irrep()
    nirr_bra = bra.n_irrep()

    # Dyson orbitals for all pairs of irreps
    dysorb = [[[] for i in range(nirr_bra)] for j in range(nirr_ket)]

    # Loop over pairs of bra and ket irreps
    for ket_irr in ket.irreps_nonzero():
        for bra_irr in bra.irreps_nonzero():

            # pairs of states for this bra irrep and ket irrep
            # bitwf uses Fortran indexing for these, hence the +1
            npairs = len(trans_list[bra_irr][ket_irr]) 

            if npairs == 0:
                continue

            dyson_pairs = 1 + np.reshape(
                np.array(trans_list[bra_irr][ket_irr], dtype=int),
                (2*npairs), order='F')

            # total number of bra and ket roots for this irrep
            bra_tot = bra.n_states_sym(bra_irr)
            ket_tot = ket.n_states_sym(ket_irr)
            
            # total number of bra and ket roots for this irrep
            bra_tot = bra.n_states_sym(bra_irr)
            ket_tot = ket.n_states_sym(ket_irr)

            # Dyson orbital array
            dysij = np.zeros((n_basis*npairs), dtype=np.float64)

            # bitwf wave function file numbers
            bra_unit = bra_wfunit[bra_irr]
            ket_unit = ket_wfunit[ket_irr]

            # compute the Dyson orbitals for all requested pairs
            # of states
            args = (bra_irr, ket_irr, bra_tot, ket_tot, npairs,
                    dyson_pairs, bra_unit, ket_unit, norm_thresh,
                    det_thresh, n_basis, dysij)

            dysij = libs.lib_func('detdyson', args)

            # Add the Dyson orbitals to the list
            dysorb[bra_irr][ket_irr] = np.reshape(dysij, (n_basis, npairs),
                                                  order='F')

    return dysorb
