"""
Overlaps between MRSF-TDDFT states of two geometries through the libmrsf
overlap module (determinant factorisation; exact two-index determinants
or the truncated Leibniz formula of Lee, Kim, Lee, Choi, JCTC 15, 882
(2019)). PySCF supplies the AO overlap matrix between the two molecules;
everything else is evaluated in the library.
"""
import sys
import numpy as np
from pyscf import gto
import graci.core.libs as libs

methods = {'exact': 0, 'tlf0': 1, 'tlf1': 2, 'tlf2': 3}


def state_overlap(Sao, Cb, Ck, occ_b, occ_k, mult, xb, xk, method='exact', align=True):
    """overlap matrix (nb, nk) between the bra states with compressed
    amplitudes xb (xdim, nb) in the MOs Cb (nao_b, nmo) and the ket states
    xk (xdim, nk) in Ck (nao_k, nmo); Sao (nao_b, nao_k) is the AO overlap
    matrix between the two geometries; occ_b/occ_k the reference
    occupations; align: align the ket MOs to the bra MOs within each
    orbital class (orthogonal Procrustes rotations from the MO overlap
    blocks, ket amplitudes transformed accordingly). Returns (S, ierr): ierr
    0 ok, 1 orbital-class mismatch, 2 explicit minors used (near-singular
    core/hole block), 3 bad method, 4 SVD failure in the alignment."""
    if method not in methods:
        sys.exit('MRSF overlaps: unknown method ' + str(method) + ' (exact|tlf0|tlf1|tlf2)')
    Cb = np.asfortranarray(Cb, dtype=np.float64)
    Ck = np.asfortranarray(Ck, dtype=np.float64)
    Sao = np.asfortranarray(Sao, dtype=np.float64)
    xb = np.asfortranarray(np.reshape(xb, (xb.shape[0], -1)), dtype=np.float64)
    xk = np.asfortranarray(np.reshape(xk, (xk.shape[0], -1)), dtype=np.float64)
    nao_b, nmo = Cb.shape
    nao_k = Ck.shape[0]
    nb, nk = xb.shape[1], xk.shape[1]
    S = np.zeros((nb, nk), dtype=np.float64, order='F')
    args = (nao_b, nao_k, nmo, np.asfortranarray(occ_b, dtype=np.float64),
            np.asfortranarray(occ_k, dtype=np.float64), int(mult), methods[method], int(bool(align)), nb, nk,
            Cb, Ck, Sao, xb, xk, S, 0)
    S, ierr = libs.lib_func('mrsf_state_overlap', args)
    return np.array(S), int(ierr)


def overlap(bra, ket, trans_list, method='exact', align=True, rep='adiabatic'):
    """
    <bra_I|ket_J> for all pairs of states in trans_list, a nested list over
    [bra_irrep][ket_irrep] of [bra_st, ket_st] pairs (0-based, within-irrep
    indices) as built by Interaction.ci_pair_list(..., sym_blk=True);
    returns the same nested layout of overlap values (the layout expected
    by Overlap.build_ci_overlaps). bra and ket are Mrsftddft objects of
    two geometries (same basis, same reference occupation pattern).
    """
    nirr_bra = bra.n_irrep()
    nirr_ket = ket.n_irrep()
    out = [[[] for i in range(nirr_ket)] for j in range(nirr_bra)]
    if bra.mult != ket.mult:
        # exact by the (+ / -) mirror structure of the two manifolds
        for bi in range(nirr_bra):
            for ki in range(nirr_ket):
                out[bi][ki] = np.zeros(len(trans_list[bi][ki]))
        return out, 0
    Sao = gto.intor_cross('int1e_ovlp', bra.scf.mol.mol_obj, ket.scf.mol.mol_obj)
    Cb, Ck = np.asarray(bra.mos), np.asarray(ket.mos)
    # all requested states of both objects in one library call
    xb_blocks, xk_blocks, boff, koff = [], [], {}, {}
    for bi in bra.irreps_nonzero():
        x = np.asarray(bra.amps[rep][bi])
        if x.shape[1] == 0:
            continue
        boff[bi] = sum(b.shape[1] for b in xb_blocks)
        xb_blocks.append(x)
    for ki in ket.irreps_nonzero():
        x = np.asarray(ket.amps[rep][ki])
        if x.shape[1] == 0:
            continue
        koff[ki] = sum(b.shape[1] for b in xk_blocks)
        xk_blocks.append(x)
    S, ierr = state_overlap(Sao, Cb, Ck, bra.occ_ref, ket.occ_ref, bra.mult,
                            np.concatenate(xb_blocks, axis=1), np.concatenate(xk_blocks, axis=1),
                            method, align)
    if ierr == 1:
        sys.exit('MRSF overlaps: the bra and ket objects have different orbital-class structures')
    if ierr == 3:
        sys.exit('MRSF overlaps: unknown method ' + str(method))
    if ierr == 4:
        sys.exit('MRSF overlaps: SVD failure in the MO alignment')
    for bi in range(nirr_bra):
        for ki in range(nirr_ket):
            pairs = trans_list[bi][ki]
            vals = np.zeros(len(pairs))
            for n, (bst, kst) in enumerate(pairs):
                vals[n] = S[boff[bi] + bst, koff[ki] + kst]
            out[bi][ki] = vals
    return out, ierr
