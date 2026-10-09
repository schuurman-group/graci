"""
Orbital Hessian of the ROKS/ROHF triplet reference of MRSF-TDDFT: the
exact second derivative of the reference energy with respect to the
orbital rotations C->O, C->V and O->V (exponential parametrisation,
one parameter per pair). It is the operator of the Z-vector equations of
the MRSF gradient and of the internal stability check of the reference.

The two-electron part (density fitting, MO basis) and the Fock couplings
are evaluated by the Z-vector operator of libmrsf, the XC kernel of the
trial densities by the library grid module, with the AO values cached
when they fit the memory budget.
"""

import numpy as np
import graci.interfaces.mrsf.mrsf_grad as mrsf_grad
import graci.interfaces.mrsf.mrsf_xc as mrsf_xc


class RoksHessian:
    """
    The library must hold the reference (mrsf_init.init) and its DF
    integrals (mrsf_init.init_ints) in double precision before the
    Hessian is constructed. Rotation pairs (Rp[k], Rq[k]) (particle p,
    hole q, MO indices) are ordered as in the library: for q in C: p in
    O; for q in C: p in V; for q in O: p in V.
    """

    def __init__(self, mf, C, occ, fa, fb, x_ref, eri_file, dims, mem_budget):
        """
        mf         : PySCF object of the reference (functional and grids)
        C, occ     : reference MOs (nao, nmo) and occupations (nmo)
        fa, fb     : alpha and beta Fock matrices in the MO basis
        x_ref      : fraction of HF exchange of the reference functional
        eri_file   : the DF integral file of the reference MOs
        dims       : (nocca, nvirb, naux)
        mem_budget : GB for the AO values cached on the grid
        """
        ierr = mrsf_grad.grad_init(eri_file)
        if ierr != 0:
            raise RuntimeError('mrsf_grad_init failed with error code %d' % ierr)

        C = np.asarray(C)
        occ = np.asarray(occ, dtype=float)
        self.nao = C.shape[0]
        self.nocca, self.nvirb = int(dims[0]), int(dims[1])
        self.x_ref = x_ref
        xc = getattr(mf, 'xc', None)
        self.is_dft = xc is not None and xc.lower() != 'hf'

        # orbital classes and the rotation space (same order as the library)
        Cidx = np.where(np.abs(occ - 2.0) < 1e-6)[0]
        Oidx = np.where(np.abs(occ - 1.0) < 1e-6)[0]
        Vidx = np.where(np.abs(occ) < 1e-6)[0]
        R = [(p, q) for q in Cidx for p in Oidx] \
            + [(p, q) for q in Cidx for p in Vidx] \
            + [(p, q) for q in Oidx for p in Vidx]
        self.Rp = np.array([p for p, q in R], dtype=int)
        self.Rq = np.array([q for p, q in R], dtype=int)
        self.lz = len(R)

        # Fock couplings and diagonal
        self.hdiag = np.asarray(mrsf_grad.zvec_setup(self.nao, self.lz, fa, fb))

        # XC kernel on the grid, with the libxc derivatives of the
        # reference cached once
        self.xcgrid = None
        if self.is_dft:
            occ_a = (occ > 0.5).astype(float)
            occ_b = (occ > 1.5).astype(float)
            H = np.concatenate([Cidx, Oidx])
            P = np.concatenate([Oidx, Vidx])
            self.xcgrid = mrsf_xc.XCGrid(mf, C, occ_a, occ_b, C[:, H], C[:, P],
                                         dims, mem_budget)

        # number of Hessian-vector products
        self.nprod = 0

    def apply(self, Z):
        """H Z for a batch of rotation vectors Z (lz, nvec)"""
        zF = mrsf_grad.farr(Z)
        nvec = zF.shape[1]
        if self.is_dft:
            Lf = mrsf_grad.zvec_factors(zF, self.nao, self.nocca)
            VHP, _ = self.xcgrid.kernel(Lf, Lf, np.ones((2, nvec), dtype=np.int32),
                                        True, False)
        else:
            VHP = mrsf_grad.fzeros(self.nocca, self.nvirb, 2, nvec)
        self.nprod += nvec
        return mrsf_grad.zvec_hessian(self.x_ref, zF, VHP, self.is_dft, self.lz)

    def diag(self):
        return self.hdiag

    def free(self):
        """release the grid cache and the library data of the Hessian"""
        if self.xcgrid is not None:
            self.xcgrid.free()
            self.xcgrid = None
        mrsf_grad.zvec_free()
        mrsf_grad.grad_free()
