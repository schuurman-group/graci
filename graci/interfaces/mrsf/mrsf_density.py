"""
Module for the calculation of MRSF-TDDFT one-electron reduced density
matrices and transition density matrices (spin-summed, MO basis), for
the standard and the extended (EMRSF) response spaces. The amplitudes,
occupations and output arrays are passed to the library by address
(F-contiguous float64, no element-wise conversion).
"""

import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing


def _ncol(ci_method):
    """number of (nvirb-long) columns of the amplitude vectors: nocca for
    standard MRSF, nocca + nC for the extended method (objects from older
    checkpoints have no ncol attribute)"""
    ncol = getattr(ci_method, 'ncol', None)
    return int(ci_method.nocca) if ncol is None else int(ncol)


def _famps(x):
    """amplitudes as an F-contiguous float64 (xdim, nstates) array"""
    return np.asfortranarray(np.asarray(x, dtype=np.float64))


def _stack(ci_method, rep='adiabatic'):
    """amplitudes of all irreps as one F-contiguous (xdim, nstates) array
    and the column offset of each irrep"""
    blocks, offs, n = [], {}, 0
    for irr in range(ci_method.n_irrep()):
        nst = ci_method.n_states_sym(irr)
        offs[irr] = n
        if nst > 0:
            blocks.append(np.asarray(ci_method.amps[rep][irr], dtype=np.float64)[:, :nst])
            n += nst
    if n == 0:
        return np.zeros((int(ci_method.xdim), 0), order='F'), offs
    return np.asfortranarray(np.concatenate(blocks, axis=1)), offs


@timing.timed
def rdm(ci_method):
    """
    State 1-RDMs for all computed states, one (nmo, nmo, nstates)
    array per irrep (views of one library call over all states)
    """

    nmo  = int(ci_method.nmo)
    xdim = int(ci_method.xdim)
    mult = int(ci_method.mult)
    ncol = _ncol(ci_method)
    occ  = np.asfortranarray(ci_method.occ_ref, dtype=np.float64)

    xvec, offs = _stack(ci_method)
    ntot = xvec.shape[1]
    dmat = np.zeros((nmo, nmo, ntot), dtype=np.float64, order='F')
    if ntot > 0:
        libs.lib_func('mrsf_density', (nmo, occ, mult, ncol, xdim, ntot, xvec, dmat))

    return [dmat[:, :, offs[irr]:offs[irr]+ci_method.n_states_sym(irr)]
            for irr in range(ci_method.n_irrep())]


@timing.timed
def tdm(bra, ket, trans_list, rep='adiabatic'):
    """
    1-TDMs <bra|E_pq|ket> for all pairs of states in trans_list, which
    is a nested list over [bra_irrep][ket_irrep] of [bra_st, ket_st]
    pairs (0-based, within-irrep indices), as built by
    Interaction.ci_pair_list(..., sym_blk=True). All irrep blocks are
    evaluated in one library call (one OpenMP region over all pairs);
    the returned blocks are views of the common output array.
    """

    nirr_bra = bra.n_irrep()
    nirr_ket = ket.n_irrep()
    nmo      = int(bra.nmo)
    xdim     = int(bra.xdim)
    mult     = int(bra.mult)
    ncol     = _ncol(bra)
    occ      = np.asfortranarray(bra.occ_ref, dtype=np.float64)

    if ket.mult != bra.mult:
        raise ValueError('MRSF transition densities require equal '
                         'multiplicities')
    if bool(getattr(bra, 'extended', False)) != bool(getattr(ket, 'extended', False)):
        raise ValueError('MRSF transition densities require two extended '
                         'or two standard MRSF-TDDFT objects')
    if int(ket.xdim) != xdim or _ncol(ket) != ncol:
        raise ValueError('MRSF transition densities require equal response '
                         'spaces for the bra and ket objects')

    xb, boff = _stack(bra, rep)
    xk, koff = _stack(ket, rep)

    # global 1-based (bra, ket) pairs in the order of the irrep blocks
    pairs, blocks, n = [], {}, 0
    for bra_irr in range(nirr_bra):
        for ket_irr in range(nirr_ket):
            blk = trans_list[bra_irr][ket_irr]
            if len(blk) == 0:
                continue
            blocks[(bra_irr, ket_irr)] = (n, n + len(blk))
            for bst, kst in blk:
                pairs.append([boff[bra_irr] + bst + 1, koff[ket_irr] + kst + 1])
            n += len(blk)

    rho = [[[] for i in range(nirr_ket)] for j in range(nirr_bra)]
    if n == 0:
        return rho

    pairs = np.reshape(np.array(pairs, dtype=int), (2*n), order='C')
    rhoij = np.zeros((nmo, nmo, n), dtype=np.float64, order='F')
    args  = (nmo, occ, mult, ncol, xdim, n, xb.shape[1], xk.shape[1], pairs,
             xb, xk, rhoij)
    libs.lib_func('mrsf_tdm', args)

    for (bra_irr, ket_irr), (i0, i1) in blocks.items():
        rho[bra_irr][ket_irr] = rhoij[:, :, i0:i1]

    return rho
