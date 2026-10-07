"""
AO-side contraction of the density-fitted two-electron gradient term for
one or more MRSF states, in auxiliary-shell blocks with hole-width factor
pairs (the naux x nao^2 three-index density is never formed).

For each state the B-space (Cholesky aux basis) families are
  Ghh(nocca,nocca,naux), Fhp(nocca,nvirb,naux), Yf(nvirb,nocca,naux), Xt(nvirb,nocca), g(naux,naux)
and the symmetrised three-index density in the raw aux basis is
  Gamma^P = sum_f (L_f^P R_f^T + R_f L_f^P^T),
  f=1: L_1^P = 1/2 (C_H Ghh^P + C_P F^P^T), R_1 = C_H ;  f=2: L_2^P = C_P Yf^P, R_2 = C_P Xt
with Ghh^P etc. the L^-T rotated families.  The contribution to the gradient is
  dE/dR_A^x = -4 sum_{mu in A,nu,P} (nabla_x mu nu|P) Gamma^P_{mu nu}
              -2 sum_{P in A,mu nu} (mu nu|nabla_x P) Gamma^P_{mu nu}
              + sum_{P in A,Q} (nabla_x P|Q) (gamma + gamma^T)_{PQ},  gamma = L^-T g L^-1.
The derivative integrals are evaluated once per auxiliary block for all states.
"""
import numpy as np
import scipy.linalg
from pyscf.ao2mo.outcore import balance_partition
from pyscf.df.grad.rhf import _int3c_wrapper
import graci.core.libs as libs


def rotate_aux(fam, L):
    """fam(..., naux) F-ordered -> fam L^-1 over the last index."""
    shape = fam.shape
    n = int(np.prod(shape[:-1]))
    A = np.asfortranarray(fam).reshape(n, shape[-1], order='F')
    X = scipy.linalg.solve_triangular(L, A.T, lower=True, trans='T').T
    return np.asfortranarray(X).reshape(shape, order='F')


def df_grad_factorised(mol, auxmol, L, C_H, C_P, states, max_memory=2000):
    """states: list of dicts with keys Xt, Ghh, Fhp, Yf, g.
    Returns an array (nstates, natm, 3)."""
    nao, naux = mol.nao, auxmol.nao_nr()
    nbas = mol.nbas
    aux_loc = auxmol.ao_loc
    nst = len(states)
    nocca = C_H.shape[1]
    Linv = scipy.linalg.solve_triangular(L, np.eye(naux), lower=True)

    # fixed right factors (F-order for the C transformation driver)
    R1 = np.asfortranarray(C_H)
    fams = []
    for st in states:
        GhhP = rotate_aux(st['Ghh'], L)
        FP = rotate_aux(st['Fhp'], L)
        YfP = rotate_aux(st['Yf'], L)
        gam = Linv.T @ st['g'] @ Linv
        # L1(mu,j,P) = 1/2 [C_H Ghh^P + C_P F^P^T], L2(mu,j,P) = C_P Yf^P  (dgemm on F-ordered views)
        L1 = 0.5*(C_H @ GhhP.reshape(nocca, -1, order='F')).reshape(nao, nocca, naux, order='F')
        FPt = np.asfortranarray(FP.transpose(1, 0, 2))                        # (a, j, P)
        L1 = L1 + 0.5*(C_P @ FPt.reshape(C_P.shape[1], -1, order='F')).reshape(nao, nocca, naux, order='F')
        L2 = (C_P @ YfP.reshape(C_P.shape[1], -1, order='F')).reshape(nao, nocca, naux, order='F')
        fams.append({'L': [np.asfortranarray(L1), np.asfortranarray(L2)],
                     'R': [R1, np.asfortranarray(C_P @ st['Xt'])],
                     'active': [bool(np.any(L1)), bool(np.any(L2))],
                     'gam2': np.asfortranarray(gam + gam.T)})
        del GhhP, FP, YfP

    # auxiliary blocks: the library holds ip1, its transpose and the
    # unpacked ip2 block (3 x 3 nao^2 nP words) plus the half-transforms
    nP_blk = int(max(8, min(naux, max_memory*1e6*0.4/(8.0*9*nao*nao))))
    ao_ranges = balance_partition(aux_loc, nP_blk)
    get_ip1 = _int3c_wrapper(mol, auxmol, 'int3c2e_ip1', 's1')
    get_ip2 = _int3c_wrapper(mol, auxmol, 'int3c2e_ip2', 's2ij')
    # active (state, family) entries contracted together per block
    flist = [(k, f) for k, fam in enumerate(fams) for f in range(2) if fam['active'][f]]
    nf = len(flist)
    Rf = np.empty((nao, nocca, nf), order='F')
    for i, (k, f) in enumerate(flist):
        Rf[:, :, i] = fams[k]['R'][f]
    tf = np.zeros((nao, 3, nf), order='F')
    u = np.zeros((nst, 3, naux))
    for shl0, shl1, nL in ao_ranges:
        p0, p1 = aux_loc[shl0], aux_loc[shl1]
        nP = p1 - p0
        # libcint buffers are (nao, nao, nP, 3) / (npair, nP, 3) Fortran-ordered: no copy
        ip1 = np.asfortranarray(get_ip1((0, nbas, 0, nbas, shl0, shl1)).transpose(1, 2, 3, 0))
        ip2 = np.asfortranarray(get_ip2((0, nbas, 0, nbas, shl0, shl1)).transpose(1, 2, 0))
        Lf = np.empty((nao, nocca, nP, nf), order='F')
        for i, (k, f) in enumerate(flist):
            Lf[:, :, :, i] = fams[k]['L'][f][:, :, p0:p1]
        ublk = np.zeros((nP, 3, nf), order='F')
        libs.lib_func('mrsf_aograd_block', (nao, nocca, nP, nf, ip1, ip2, Lf, Rf, tf, ublk))
        for i, (k, f) in enumerate(flist):
            u[k, :, p0:p1] += ublk[:, :, i].T
    t = np.zeros((nst, 3, nao))
    for i, (k, f) in enumerate(flist):
        t[k] += tf[:, :, i].T
    ip2c = auxmol.intor('int2c2e_ip1', comp=3)
    aoslices = mol.aoslice_by_atom()[:, 2]
    auxslices = auxmol.aoslice_by_atom()[:, 2]
    de = np.zeros((nst, mol.natm, 3))
    for k, fam in enumerate(fams):
        w = np.einsum('xPQ,PQ->xP', ip2c, fam['gam2'])
        de[k] += -4.0*np.add.reduceat(t[k], aoslices, axis=1).T
        de[k] += (-2.0*np.add.reduceat(u[k], auxslices, axis=1) + np.add.reduceat(w, auxslices, axis=1)).T
    return de
