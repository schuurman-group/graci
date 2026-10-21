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


# non-fatal library return codes: 2 explicit minors used for a near-singular
# core/hole/G block, 5 SOMO pair rotated or swapped between the geometries
# (extended objects: the alignment of the SOMO pair is restricted to sign
# flips), 7 = both
warning_codes = {2: 'near-singular core/hole MO overlap block; explicit minors used',
                 5: 'the SOMO pair is rotated or swapped between the two geometries; '
                    'the TLF overlaps of extended objects are unreliable (exact ones are not affected)',
                 7: 'near-singular MO overlap block (explicit minors used) and rotated SOMO pair '
                    '(TLF overlaps of extended objects unreliable)'}


def state_overlap(Sao, Cb, Ck, occ_b, occ_k, mult, xb, xk, method='exact', align=True, wpib=None, wpik=None):
    """overlap matrix (nb, nk) between the bra states with compressed
    amplitudes xb (xdim, nb) in the MOs Cb (nao_b, nmo) and the ket states
    xk (xdim, nk) in Ck (nao_k, nmo); Sao (nao_b, nao_k) is the AO overlap
    matrix between the two geometries; occ_b/occ_k the reference
    occupations; align: align the ket MOs to the bra MOs within each
    orbital class (orthogonal Procrustes rotations from the MO overlap
    blocks, ket amplitudes transformed accordingly). Each side is standard
    (xdim = nvirb*nocca) or extended (xdim = nvirb*(nocca+nC), the CV
    amplitudes in the last nC columns); a standard side paired with an
    extended one counts as an extended state without CV part. wpib/wpik:
    parent weights (w_G, w_D, w_S) of the extension space Pi of the bra and
    ket objects (ext_space = pi; exact method only; each side its own
    weights), None for the published parent G. Returns (S, ierr): ierr 0
    ok, 1 orbital-class mismatch, 2 explicit minors used (near-singular
    core/hole/G block), 3 bad method, 4 SVD failure in the alignment, 5 SOMO
    pair rotated/swapped (extended objects, warning), 6 TLF closed forms of
    the extended blocks need the orbital classes in index order, 7 = 2 and
    5, 8 invalid amplitude dimension, 9 TLF requested for the parent Pi,
    10 singular hole block for the parent Pi with too many CV slots for the
    explicit fallback."""
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
    occ_b = np.asfortranarray(occ_b, dtype=np.float64)
    occ_k = np.asfortranarray(occ_k, dtype=np.float64)
    # column counts of the two sides (standard: nocca, extended: nocca + nC)
    nC = int(np.sum(np.abs(occ_b - 2.0) < 1.0e-6))
    nocca, nvirb = nC + 2, nmo - nC
    ncol = []
    for side, x in (('bra', xb), ('ket', xk)):
        n = x.shape[0] // nvirb
        if x.shape[0] != n*nvirb or n not in (nocca, nocca + nC):
            sys.exit('MRSF overlaps: the ' + side + ' amplitude vectors have ' + str(x.shape[0])
                     + ' elements; expected ' + str(nvirb*nocca) + ' (standard) or '
                     + str(nvirb*(nocca + nC)) + ' (extended)')
        ncol.append(n)
    S = np.zeros((nb, nk), dtype=np.float64, order='F')
    pimode = 0 if (wpib is None and wpik is None) else 1
    wb = np.asfortranarray(wpib if wpib is not None else (1.0, 0.0, 0.0), dtype=np.float64)
    wk = np.asfortranarray(wpik if wpik is not None else (1.0, 0.0, 0.0), dtype=np.float64)
    args = (nao_b, nao_k, nmo, occ_b, occ_k, int(mult), methods[method], int(bool(align)), nb, nk,
            ncol[0], ncol[1], pimode, wb, wk, Cb, Ck, Sao, xb, xk, S, 0)
    S, ierr = libs.lib_func('mrsf_state_overlap', args)
    return np.array(S), int(ierr)


def align_vectors(Sao, Cb, Ck, occ_b, occ_k, mult, ext, mode, xk):
    """warm-start transformation of the vectors xk (xdim, nk) of a previous
    geometry (MOs Ck) to the present MOs Cb: the ket MOs are aligned to the
    bra MOs within the orbital classes (orthogonal Procrustes rotations of
    the MO overlap blocks; the SOMO pair only sign-aligned when ext) and
    the vectors are transformed accordingly. mode 0: compressed amplitude
    vectors (standard or extended layout); mode 1: Z-vector-type matrices
    in the (nvirb, nocca) layout (no reference sign). Sao (nao_b, nao_k) is
    the AO overlap matrix between the two geometries. Returns
    (aligned vectors, dmin (3) = smallest aligned diagonal MO overlap of the
    core, SOMO and virtual classes, ierr: 0 ok, 1 class mismatch, 4 SVD
    failure, 5 SOMO pair rotated/swapped (ext, warning), 8 bad dimension)."""
    Cb = np.asfortranarray(Cb, dtype=np.float64)
    Ck = np.asfortranarray(Ck, dtype=np.float64)
    Sao = np.asfortranarray(Sao, dtype=np.float64)
    xk = np.asfortranarray(np.reshape(xk, (np.shape(xk)[0], -1)), dtype=np.float64).copy(order='F')
    nao_b, nmo = Cb.shape
    nao_k = Ck.shape[0]
    nk = xk.shape[1]
    occ_b = np.asfortranarray(occ_b, dtype=np.float64)
    occ_k = np.asfortranarray(occ_k, dtype=np.float64)
    nC = int(np.sum(np.abs(occ_k - 2.0) < 1.0e-6))
    nocca, nvirb = nC + 2, nmo - nC
    ncol = xk.shape[0] // nvirb
    if xk.shape[0] != ncol*nvirb or ncol not in (nocca, nocca + nC):
        sys.exit('MRSF warm start: the vectors have ' + str(xk.shape[0]) + ' elements; expected '
                 + str(nvirb*nocca) + ' (standard) or ' + str(nvirb*(nocca + nC)) + ' (extended)')
    dmin = np.zeros(3, dtype=np.float64)
    args = (nao_b, nao_k, nmo, occ_b, occ_k, int(mult), int(bool(ext)), int(mode), nk, ncol,
            Cb, Ck, Sao, xk, dmin, 0)
    xk, dmin, ierr = libs.lib_func('mrsf_align_vectors', args)
    return np.array(xk), np.array(dmin), int(ierr)


def overlap(bra, ket, trans_list, method='exact', align=True, rep='adiabatic'):
    """
    <bra_I|ket_J> for all pairs of states in trans_list, a nested list over
    [bra_irrep][ket_irrep] of [bra_st, ket_st] pairs (0-based, within-irrep
    indices) as built by Interaction.ci_pair_list(..., sym_blk=True);
    returns the same nested layout of overlap values (the layout expected
    by Overlap.build_ci_overlaps). bra and ket are Mrsftddft objects of
    two geometries (same basis, same reference occupation pattern),
    standard or extended (a standard object paired with an extended one
    is treated as an extended state without CV amplitudes).
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
    # parent Pi of the extension space: each side with its own weights
    pi = [getattr(o, 'extended', False) and str(getattr(o, 'ext_space', 'published')).lower() == 'pi'
          for o in (bra, ket)]
    wpib = wpik = None
    if any(pi):
        from graci.interfaces.mrsf.mrsf_density import _wpi
        wpib, wpik = _wpi(bra), _wpi(ket)
    S, ierr = state_overlap(Sao, Cb, Ck, bra.occ_ref, ket.occ_ref, bra.mult,
                            np.concatenate(xb_blocks, axis=1), np.concatenate(xk_blocks, axis=1),
                            method, align, wpib, wpik)
    if ierr == 1:
        sys.exit('MRSF overlaps: the bra and ket objects have different orbital-class structures')
    if ierr == 3:
        sys.exit('MRSF overlaps: unknown method ' + str(method))
    if ierr == 4:
        sys.exit('MRSF overlaps: SVD failure in the MO alignment')
    if ierr == 6:
        sys.exit('MRSF overlaps: the TLF variants for extended objects require the orbital '
                 'classes in index order (core, SOMOs, virtuals); use mrsf_method = exact')
    if ierr == 8:
        sys.exit('MRSF overlaps: inconsistent amplitude dimensions')
    if ierr == 9:
        sys.exit('MRSF overlaps: only mrsf_method = exact is available for ext_space = pi')
    if ierr == 10:
        sys.exit('MRSF overlaps: the occupied MO spaces of the two geometries are (nearly) '
                 'orthogonal (singular hole block); the Pi overlaps are not available for this '
                 'pair of geometries')
    for bi in range(nirr_bra):
        for ki in range(nirr_ket):
            pairs = trans_list[bi][ki]
            vals = np.zeros(len(pairs))
            for n, (bst, kst) in enumerate(pairs):
                vals[n] = S[boff[bi] + bst, koff[ki] + kst]
            out[bi][ki] = vals
    return out, ierr
