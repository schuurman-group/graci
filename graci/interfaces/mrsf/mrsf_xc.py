"""
Exchange-correlation kernel contractions of the MRSF gradient through the
libmrsf grid module: PySCF supplies the quadrature grid, the libxc
functional derivatives and the AO values on the grid once; the kernel
action on factor-pair densities (D = L R^T + R L^T with hole-width
factors) is evaluated inside the library, projected directly onto the
hole-particle and hole-hole MO blocks.
"""
import time
import numpy as np
import graci.core.libs as libs
import graci.interfaces.mrsf.mrsf_grad as mrsf_grad


class XCGrid:

    def __init__(self, mf, C, occ_a, occ_b, C_H, C_P, dims, mem_budget_gb=1.0, block_bytes=1.0e6):
        """mf: PySCF (RO)KS object with built grids; C: MO coefficients;
        occ_a/occ_b: alpha/beta occupation vectors; C_H, C_P: hole and
        particle MO blocks; dims = (nocca, nvirb, naux). The AO values are
        cached in the library when they fit mem_budget_gb, otherwise they
        are re-evaluated block by block at every kernel call."""
        self.mf, self.mol, self.ni = mf, mf.mol, mf._numint
        self.time_kernel, self.ncalls_kernel, self.time_probe = 0.0, 0, 0.0
        self.dims = dims
        nocca, nvirb, naux = dims
        self.nao = C.shape[0]
        xctype = self.ni._xc_type(mf.xc)
        if xctype == 'LDA':
            self.ncomp, self.ao_deriv = 1, 0
        elif xctype == 'GGA':
            self.ncomp, self.ao_deriv = 4, 1
        else:
            raise NotImplementedError('MRSF gradients: XC type ' + xctype)
        grids = mf.grids
        self.coords = grids.coords
        ngrid = grids.weights.size
        self.ngrid = ngrid
        rho, vxc, fxc = self.ni.cache_xc_kernel(self.mol, grids, mf.xc, (C, C), (occ_a, occ_b), spin=1)
        nv = self.ncomp
        # PySCF vxc[s,x,g], fxc[s,x,t,y,g] (C order) = Fortran (g,x,s), (g,y,t,x,s)
        vxc = np.ascontiguousarray(np.asarray(vxc).reshape(2, nv, ngrid))
        fxc = np.ascontiguousarray(np.asarray(fxc).reshape(2, nv, 2, nv, ngrid))
        # grid blocks: ~block_bytes of AO values each
        nb = int(block_bytes/(8.0*self.nao*self.ncomp))
        nb = max(64, min(8192, (nb//64)*64))
        self.gbeg = np.arange(0, ngrid, nb)
        self.gend = np.minimum(self.gbeg + nb, ngrid)
        self.nblocks = len(self.gbeg)
        cache_bytes = 8.0*self.ncomp*ngrid*(self.nao + nocca)
        self.cache = bool(cache_bytes <= mem_budget_gb*1.0e9)
        libs.lib_func('mrsf_xc_init', (self.nao, int(ngrid), self.ncomp, nv, self.nblocks,
                                       np.asarray(self.gbeg + 1, dtype=np.int32),
                                       np.asarray(self.gend, dtype=np.int32),
                                       np.asfortranarray(grids.weights), fxc.T, vxc.T,
                                       np.asfortranarray(C_H), np.asfortranarray(C_P), self.cache))
        if self.cache:
            for ib in range(self.nblocks):
                libs.lib_func('mrsf_xc_add_block', (ib + 1, self.eval_block(ib)))

    @staticmethod
    def _ao_fortran(ao):
        """PySCF AO values (ncomp, nb, nao), natively stored as (ncomp, nao, nb),
        as the F-contiguous (nb, nao, ncomp) array the library expects (no copy)"""
        if ao.ndim == 2:
            ao = ao[None]
        return np.asarray(ao.transpose(0, 2, 1), order='C').T

    def eval_block(self, ib, deriv=None):
        """AO values on block ib, with derivatives up to deriv (default: those
        of the kernel), as an F-contiguous (nb, nao, ncomp) array"""
        if deriv is None:
            deriv = self.ao_deriv
        ao = self.ni.eval_ao(self.mol, self.coords[self.gbeg[ib]:self.gend[ib]], deriv=deriv)
        return self._ao_fortran(ao)

    def kernel(self, Lfac, Rfac, rch, want_hp=True, want_hh=False):
        """kernel potentials of the densities D^s_v = L R^T + R L^T,
        Lfac, Rfac: (nao, k, 2, nvec); rch (2, nvec): 1 where R is the
        hole block C_H (cached on the grid; then k = nocca).
        Returns VHP (nocca, nvirb, 2, nvec), VHH (nocca, nocca, 2, nvec)
        (zero when not requested)."""
        t0 = time.time()
        nocca, nvirb, naux = self.dims
        Lfac = np.asfortranarray(Lfac, dtype=np.float64)
        Rfac = np.asfortranarray(Rfac, dtype=np.float64)
        nao, k, two, nvec = Lfac.shape
        rch = np.asarray(np.reshape(rch, (2, nvec)), dtype=np.int32).flatten(order='F')
        libs.lib_func('mrsf_xc_begin', (nvec, k, Lfac, Rfac, rch))
        if self.cache:
            libs.lib_func('mrsf_xc_cached', ())
        else:
            for ib in range(self.nblocks):
                libs.lib_func('mrsf_xc_block', (ib + 1, self.eval_block(ib)))
        VHP = mrsf_grad.fzeros(nocca, nvirb, 2, nvec)
        VHH = mrsf_grad.fzeros(nocca, nocca, 2, nvec)
        libs.lib_func('mrsf_xc_end', (bool(want_hp), bool(want_hh), VHP, VHH))
        self.time_kernel += time.time() - t0
        self.ncalls_kernel += 1
        return VHP, VHH

    def probe(self, Lfac, Rfac, rch):
        """XC probe terms of the gradient on the fixed grid for the state
        densities P^s_v = L R^T + R L^T, Lfac/Rfac (nao, k, 2, nst), rch as
        in kernel(): returns t (nao, 3, nst+1) with entry 0 the reference
        (vxc, D) term and entry v+1 the (vxc, P_v) + (f_xc[P_v], D) term of
        state v; the gradient contributions are -2 sum_{mu in A} t_x(mu).
        The AO values with the derivatives needed are streamed from PySCF."""
        t0 = time.time()
        Lfac = np.asfortranarray(Lfac, dtype=np.float64)
        Rfac = np.asfortranarray(Rfac, dtype=np.float64)
        nao, k, two, nst = Lfac.shape
        rch = np.asarray(np.reshape(rch, (2, nst)), dtype=np.int32).flatten(order='F')
        libs.lib_func('mrsf_xc_probe_begin', (nst, k, Lfac, Rfac, rch))
        for ib in range(self.nblocks):
            aoF = self.eval_block(ib, deriv=self.ao_deriv + 1)
            libs.lib_func('mrsf_xc_probe_block', (ib + 1, int(aoF.shape[2]), aoF))
        t = mrsf_grad.fzeros(nao, 3, nst + 1)
        libs.lib_func('mrsf_xc_probe_end', (t,))
        self.time_probe += time.time() - t0
        return t

    def print_timings(self):
        """library-side timers of the grid kernel (diagnostics)"""
        libs.lib_func('mrsf_xc_timings', ())

    def free(self):
        libs.lib_func('mrsf_xc_free', ())
