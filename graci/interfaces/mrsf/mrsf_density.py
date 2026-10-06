"""
Module for the calculation of MRSF-TDDFT one-electron reduced density
matrices and transition density matrices (spin-summed, MO basis)
"""

import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing

@timing.timed
def rdm(ci_method):
    """
    State 1-RDMs for all computed states, one (nmo, nmo, nstates)
    array per irrep
    """

    nmo  = ci_method.nmo
    xdim = int(ci_method.xdim)
    mult = int(ci_method.mult)
    occ  = np.array(ci_method.occ_ref, dtype=float)

    dmat_sym = []
    for irrep in range(ci_method.n_irrep()):
        nst = ci_method.n_states_sym(irrep)
        if nst == 0:
            dmat_sym.append(np.zeros((nmo, nmo, 0), dtype=float))
            continue
        xvec = np.reshape(np.asarray(ci_method.amps['adiabatic'][irrep]),
                          (xdim*nst), order='F')
        dmat = np.zeros(nmo*nmo*nst, dtype=float)
        args = (nmo, occ, mult, xdim, nst, xvec, dmat)
        dmat = libs.lib_func('mrsf_density', args)
        dmat_sym.append(np.reshape(dmat, (nmo, nmo, nst), order='F'))

    return dmat_sym

@timing.timed
def tdm(bra, ket, trans_list, rep='adiabatic'):
    """
    1-TDMs <bra|E_pq|ket> for all pairs of states in trans_list, which
    is a nested list over [bra_irrep][ket_irrep] of [bra_st, ket_st]
    pairs (0-based, within-irrep indices), as built by
    Interaction.ci_pair_list(..., sym_blk=True)
    """

    nirr_bra = bra.n_irrep()
    nirr_ket = ket.n_irrep()
    nmo      = bra.nmo
    xdim     = int(bra.xdim)
    mult     = int(bra.mult)
    occ      = np.array(bra.occ_ref, dtype=float)

    if ket.mult != bra.mult:
        raise ValueError('MRSF transition densities require equal '
                         'multiplicities')

    rho = [[[] for i in range(nirr_ket)] for j in range(nirr_bra)]

    for bra_irr in bra.irreps_nonzero():
        for ket_irr in ket.irreps_nonzero():

            npairs = len(trans_list[bra_irr][ket_irr])
            if npairs == 0:
                continue

            # 1-based (bra, ket) pairs, Fortran order (2, npairs)
            pairs = 1 + np.reshape(np.array(trans_list[bra_irr][ket_irr],
                                            dtype=int), (2*npairs), order='C')

            xb = np.asarray(bra.amps[rep][bra_irr])
            xk = np.asarray(ket.amps[rep][ket_irr])
            nb = xb.shape[1]
            nk = xk.shape[1]
            xb = np.reshape(xb, (xdim*nb), order='F')
            xk = np.reshape(xk, (xdim*nk), order='F')

            rhoij = np.zeros(nmo*nmo*npairs, dtype=float)
            args  = (nmo, occ, mult, xdim, npairs, nb, nk, pairs, xb, xk,
                     rhoij)
            rhoij = libs.lib_func('mrsf_tdm', args)

            rho[bra_irr][ket_irr] = np.reshape(rhoij, (nmo, nmo, npairs),
                                               order='F')

    return rho
