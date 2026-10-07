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
import ctypes
import numpy as np
import scipy.linalg
from pyscf import lib
from pyscf.ao2mo import _ao2mo
from pyscf.ao2mo.outcore import balance_partition
from pyscf.df.grad.rhf import _int3c_wrapper


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
        L1 = 0.5*(np.einsum('mi,ijP->mjP', C_H, GhhP) + np.einsum('ma,jaP->mjP', C_P, FP))
        L2 = np.einsum('ma,ajP->mjP', C_P, YfP)
        fams.append({'L': [np.asfortranarray(L1), np.asfortranarray(L2)],
                     'R': [R1, np.asfortranarray(C_P @ st['Xt'])],
                     'gam2': np.asfortranarray(gam + gam.T)})
        del GhhP, FP, YfP

    blksize = int(min(max(max_memory*.5e6/8/(nao**2*3), 20), naux, 240))
    ao_ranges = balance_partition(aux_loc, blksize)
    get_ip1 = _int3c_wrapper(mol, auxmol, 'int3c2e_ip1', 's1')
    get_ip2 = _int3c_wrapper(mol, auxmol, 'int3c2e_ip2', 's1')
    fmmm = _ao2mo.libao2mo.AO2MOmmm_bra_nr_s1
    fdrv = _ao2mo.libao2mo.AO2MOnr_e2_drv
    ftrans = _ao2mo.libao2mo.AO2MOtranse2_nr_s1
    null = lib.c_null_ptr()

    t = np.zeros((nst, 3, nao))
    u = np.zeros((nst, 3, naux))
    for shl0, shl1, nL in ao_ranges:
        p0, p1 = aux_loc[shl0], aux_loc[shl1]
        nP = p1 - p0
        ip1 = np.ascontiguousarray(get_ip1((0, nbas, 0, nbas, shl0, shl1)).transpose(0, 3, 2, 1))
        ip2 = np.ascontiguousarray(get_ip2((0, nbas, 0, nbas, shl0, shl1)).transpose(0, 3, 2, 1))
        M1 = np.empty((3, nP, nocca, nao))
        N1 = np.empty((3, nP, nocca, nao))
        for k, fam in enumerate(fams):
            for f in range(2):
                Rf, LfP = fam['R'][f], fam['L'][f][:, :, p0:p1]
                LfC = np.ascontiguousarray(LfP.transpose(2, 1, 0))             # (P, o, mu)
                fdrv(ftrans, fmmm, M1.ctypes.data_as(ctypes.c_void_p), ip1.ctypes.data_as(ctypes.c_void_p),
                     Rf.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(3*nP), ctypes.c_int(nao),
                     (ctypes.c_int*4)(0, nocca, 0, nao), null, ctypes.c_int(0))
                t[k] += np.einsum('xpom,pom->xm', M1, LfC, optimize=True)
                LfK = np.ascontiguousarray(LfP.transpose(0, 2, 1)).reshape(nao*nP, nocca, order='F')
                for x in range(3):
                    M2 = ip1[x].reshape(nP*nao, nao).T @ LfK
                    t[k, x] += np.einsum('mo,mo->m', M2, Rf)
                fdrv(ftrans, fmmm, N1.ctypes.data_as(ctypes.c_void_p), ip2.ctypes.data_as(ctypes.c_void_p),
                     Rf.ctypes.data_as(ctypes.c_void_p), ctypes.c_int(3*nP), ctypes.c_int(nao),
                     (ctypes.c_int*4)(0, nocca, 0, nao), null, ctypes.c_int(0))
                u[k, :, p0:p1] += 2.0*np.einsum('xpom,pom->xp', N1, LfC, optimize=True)
    ip2c = auxmol.intor('int2c2e_ip1', comp=3)
    aoslices = mol.aoslice_by_atom()[:, 2]
    auxslices = auxmol.aoslice_by_atom()[:, 2]
    de = np.zeros((nst, mol.natm, 3))
    for k, fam in enumerate(fams):
        w = np.einsum('xPQ,PQ->xP', ip2c, fam['gam2'])
        de[k] += -4.0*np.add.reduceat(t[k], aoslices, axis=1).T
        de[k] += (-2.0*np.add.reduceat(u[k], auxslices, axis=1) + np.add.reduceat(w, auxslices, axis=1)).T
    return de
