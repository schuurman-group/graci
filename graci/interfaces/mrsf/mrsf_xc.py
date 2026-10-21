"""
Exchange-correlation kernel contractions of the MRSF gradient through the
libmrsf grid module: PySCF supplies the quadrature grid, the libxc
functional derivatives and the AO values on the grid once; the kernel
action on factor-pair densities (D = L R^T + R L^T with hole-width
factors) is evaluated inside the library, projected directly onto the
hole-particle and hole-hole MO blocks.
"""
import time
import sys
import numpy as np
import graci.core.libs as libs
import graci.interfaces.mrsf.mrsf_grad as mrsf_grad


class XCGrid:

    def __init__(self, mf, C, occ_a, occ_b, C_H, C_P, dims, mem_budget_gb=1.0, block_bytes=1.0e6, lazy=False):
        """mf: PySCF (RO)KS object with built grids; C: MO coefficients;
        occ_a/occ_b: alpha/beta occupation vectors; C_H, C_P: hole and
        particle MO blocks; dims = (nocca, nvirb, naux). The AO values are
        cached in the library when they fit mem_budget_gb, otherwise they
        are re-evaluated block by block at every kernel call. The kernel of
        the reference density is installed at once unless lazy (then on the
        first call of ensure_ref_kernel)."""
        self.mf, self.mol, self.ni = mf, mf.mol, mf._numint
        self.time_kernel, self.ncalls_kernel, self.time_probe = 0.0, 0, 0.0
        self.time_init_kernel, self.time_init_cache = 0.0, 0.0
        self.C = np.asarray(C)
        self.occ_a, self.occ_b = np.asarray(occ_a, float), np.asarray(occ_b, float)
        self.ref_kernel = False
        self.g_kernel, self.g_gidx, self.g_with_kxc, self.g_arrays = False, None, False, None
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
        nv = self.ncomp
        # grid blocks: ~block_bytes of AO values each
        nb = int(block_bytes/(8.0*self.nao*self.ncomp))
        nb = max(64, min(8192, (nb//64)*64))
        self.gbeg = np.arange(0, ngrid, nb)
        self.gend = np.minimum(self.gbeg + nb, ngrid)
        self.nblocks = len(self.gbeg)
        cache_bytes = 8.0*self.ncomp*ngrid*(self.nao + nocca)
        self.cache = bool(cache_bytes <= mem_budget_gb*1.0e9)
        libs.lib_func('mrsf_xc_init_grid', (self.nao, int(ngrid), self.ncomp, nv, self.nblocks,
                                            np.asarray(self.gbeg + 1, dtype=np.int32),
                                            np.asarray(self.gend, dtype=np.int32),
                                            np.asfortranarray(grids.weights),
                                            np.asfortranarray(C_H), np.asfortranarray(C_P), self.cache))
        if self.cache:
            t0 = time.time()
            for ib in range(self.nblocks):
                libs.lib_func('mrsf_xc_add_block', (ib + 1, self.eval_block(ib)))
            self.time_init_cache = time.time() - t0
        if not lazy:
            self.ensure_ref_kernel()

    def ensure_ref_kernel(self):
        """install the spin-resolved kernel of the reference density (needed
        by kernel() and probe()); computed once"""
        if self.ref_kernel:
            return
        t0 = time.time()
        nv, ngrid = self.ncomp, self.ngrid
        rho, vxc, fxc = self.ni.cache_xc_kernel(self.mol, self.mf.grids, self.mf.xc, (self.C, self.C),
                                                (self.occ_a, self.occ_b), spin=1)
        # PySCF vxc[s,x,g], fxc[s,x,t,y,g] (C order) = Fortran (g,x,s), (g,y,t,x,s)
        vxc = np.ascontiguousarray(np.asarray(vxc).reshape(2, nv, ngrid))
        fxc = np.ascontiguousarray(np.asarray(fxc).reshape(2, nv, 2, nv, ngrid))
        libs.lib_func('mrsf_xc_set_ref_kernel', (int(ngrid), nv, vxc.T, fxc.T))
        self.ref_kernel = True
        self.time_init_kernel += time.time() - t0

    def reset_counters(self):
        self.time_kernel, self.ncalls_kernel, self.time_probe = 0.0, 0, 0.0

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

    def probe(self, Lfac, Rfac, rch, gfac=None, g_kxc=True):
        """XC probe terms of the gradient on the fixed grid for the state
        densities P^s_v = L R^T + R L^T, Lfac/Rfac (nao, k, 2, nst), rch as
        in kernel(): returns t (nao, 3, nst+1) with entry 0 the reference
        (vxc, D) term and entry v+1 the (vxc, P_v) + (f_xc[P_v], D) term of
        state v; the gradient contributions are -2 sum_{mu in A} t_x(mu).
        gfac = (LM, RM, LS) registers, for the extended method, the
        G-density terms per state (see gprobe terms in the library) and the
        result is then (t, tG) with tG (nao, 3, nst). The AO values with the
        derivatives needed are streamed from PySCF."""
        t0 = time.time()
        Lfac = np.asfortranarray(Lfac, dtype=np.float64)
        Rfac = np.asfortranarray(Rfac, dtype=np.float64)
        nao, k, two, nst = Lfac.shape
        rch = np.asarray(np.reshape(rch, (2, nst)), dtype=np.int32).flatten(order='F')
        libs.lib_func('mrsf_xc_probe_begin', (nst, k, Lfac, Rfac, rch))
        if gfac is not None:
            LM, RM, LS = [np.asfortranarray(a, dtype=np.float64) for a in gfac]
            kM, kS = LM.shape[1], LS.shape[1]
            libs.lib_func('mrsf_xc_gprobe_set', (nst, kM, LM, RM, kS, LS, bool(g_kxc)))
        for ib in range(self.nblocks):
            aoF = self.eval_block(ib, deriv=self.ao_deriv + 1)
            libs.lib_func('mrsf_xc_probe_block', (ib + 1, int(aoF.shape[2]), aoF))
        tG = None
        if gfac is not None:
            tG = mrsf_grad.fzeros(nao, 3, nst)
            libs.lib_func('mrsf_xc_gprobe_get', (tG,))
        t = mrsf_grad.fzeros(nao, 3, nst + 1)
        libs.lib_func('mrsf_xc_probe_end', (t,))
        self.time_probe += time.time() - t0
        if gfac is not None:
            return t, tG
        return t

    def probe_q(self, Bexp, nb, idx, S, with_g=False, with_kxc=False, nsuper=32, with_w=False):
        """XC probe terms of the gradient for densities given as quadratic
        forms P = B S B^T in the basis [hole MOs (cached) | Bexp (nao, nexp)]:
        nb (4, nst) term lengths, idx (maxnb, 4, nst) 1-based basis indices
        (1..nocca the hole MOs, nocca+1.. the columns of Bexp), S (maxnb,
        maxnb, 4, nst); terms 0, 1 = the alpha/beta relaxed densities of
        each state, 2, 3 = the extended method's M and D_s (when with_g).
        The AO values with the derivatives needed are streamed from PySCF
        in superblocks of nsuper cache blocks and processed in parallel over
        the blocks. Returns t (nao, 3, nst+1) [and tG (nao, 3, nst)] as probe()."""
        t0 = time.time()
        nao = self.nao
        nb = np.asarray(nb, dtype=np.int32); idx = np.asarray(idx, dtype=np.int32)
        nkind = nb.shape[0]
        nst = nb.shape[1]; maxnb = idx.shape[0]
        Bexp = np.asfortranarray(np.asarray(Bexp, dtype=np.float64).reshape(nao, -1))
        nexp = Bexp.shape[1]
        if nexp == 0:
            Bexp = mrsf_grad.fzeros(nao, 1)
        S = np.asfortranarray(np.asarray(S, dtype=np.float64).reshape(maxnb, maxnb, nkind, nst))
        if nkind == 5:
            libs.lib_func('mrsf_xc_probe_begin_q2', (nst, nexp, Bexp, maxnb, nb.flatten(order='F'),
                                                     idx.flatten(order='F'), S, bool(with_g), bool(with_kxc),
                                                     bool(with_w)))
        else:
            if with_w:
                raise ValueError('probe_q: the E_W terms need the five-kind term arrays')
            libs.lib_func('mrsf_xc_probe_begin_q', (nst, nexp, Bexp, maxnb, nb.flatten(order='F'),
                                                    idx.flatten(order='F'), S, bool(with_g), bool(with_kxc)))
        ib = 0
        while ib < self.nblocks:
            ns = min(int(nsuper), self.nblocks - ib)
            g0, g1 = int(self.gbeg[ib]), int(self.gend[ib + ns - 1])
            ao = self.ni.eval_ao(self.mol, self.coords[g0:g1], deriv=self.ao_deriv + 1)
            aoF = self._ao_fortran(ao)
            libs.lib_func('mrsf_xc_probe_super', (ib + 1, ns, g1 - g0, int(aoF.shape[2]), aoF))
            ib += ns
        tG = None
        if with_g:
            tG = mrsf_grad.fzeros(nao, 3, nst)
            libs.lib_func('mrsf_xc_gprobe_get', (tG,))
        t = mrsf_grad.fzeros(nao, 3, nst + 1)
        libs.lib_func('mrsf_xc_probe_end', (t,))
        self.time_probe += time.time() - t0
        if with_g:
            return t, tG
        return t

    def set_g_kernel(self, C, occG, with_kxc):
        """closed-shell kernel set at the density of the configuration G
        (occG: occupations 2 on the G orbitals); see ensure_g_kernel"""
        self.ensure_g_kernel(np.where(np.asarray(occG) > 0.5)[0], with_kxc)

    def ensure_g_kernel(self, gidx, with_kxc, g_default=None):
        """closed-shell kernel set at the density of the configuration G of
        the extended method (orbitals gidx doubly occupied), in PySCF's spin-0
        convention (derivatives of the total density; trial densities
        spin-summed): v_xc, f_xc and, when with_kxc, k_xc. Computed once per
        (gidx, with_kxc) from the cached hole MOs when gidx is the library's
        default G (the first nocca - 1 hole MOs, g_default), else through
        PySCF. Installing the set with k_xc covers later requests without."""
        gidx = tuple(int(g) for g in gidx)
        if self.g_kernel and self.g_gidx == gidx and (self.g_with_kxc or not with_kxc):
            return
        t0 = time.time()
        xc, ni, mol, grids = self.mf.xc, self.ni, self.mol, self.mf.grids
        nv, ngrid = self.ncomp, self.ngrid
        xctype = ni._xc_type(xc)
        rho = None
        if self.cache and g_default is not None and tuple(int(g) for g in g_default) == gidx:
            rho_f = mrsf_grad.fzeros(ngrid, nv)
            libs.lib_func('mrsf_xc_g_density', (rho_f,))
            rho = np.ascontiguousarray(rho_f.T)
            if nv == 1:
                rho = rho.reshape(ngrid)
        if rho is None:
            occG = np.zeros(self.C.shape[1]); occG[list(gidx)] = 2.0
            rho = ni.cache_xc_kernel(mol, grids, xc, self.C, occG, spin=0)[0]
        out = ni.eval_xc_eff(xc, rho, deriv=3 if with_kxc else 2, xctype=xctype)
        vxc = np.ascontiguousarray(np.asarray(out[1]).reshape(nv, ngrid))
        fxc = np.ascontiguousarray(np.asarray(out[2]).reshape(nv, nv, ngrid))
        if with_kxc:
            kxc = np.ascontiguousarray(np.asarray(out[3]).reshape(nv, nv, nv, ngrid))
        else:
            kxc = np.zeros((nv, nv, nv, 1))
        libs.lib_func('mrsf_xc_g_set', (int(ngrid), nv, vxc.T, fxc.T, np.asfortranarray(kxc.T), bool(with_kxc)))
        self.g_kernel, self.g_gidx, self.g_with_kxc = True, gidx, bool(with_kxc)
        self.g_arrays = (rho, vxc, fxc, kxc if with_kxc else None)
        self.time_init_kernel += time.time() - t0

    def ensure_kernel_occ(self, occ_total, with_kxc=False):
        """closed-shell kernel set (spin-0 convention, see ensure_g_kernel) at
        the density of the (fractional) total occupations occ_total of the
        MOs C, through PySCF's density evaluation; keyed on the occupations"""
        key = ('occ',) + tuple(float(round(o, 10)) for o in np.asarray(occ_total, float))
        if self.g_kernel and self.g_gidx == key and (self.g_with_kxc or not with_kxc):
            return
        t0 = time.time()
        xc, ni, mol, grids = self.mf.xc, self.ni, self.mol, self.mf.grids
        nv, ngrid = self.ncomp, self.ngrid
        xctype = ni._xc_type(xc)
        rho = ni.cache_xc_kernel(mol, grids, xc, self.C, np.asarray(occ_total, float), spin=0)[0]
        out = ni.eval_xc_eff(xc, rho, deriv=3 if with_kxc else 2, xctype=xctype)
        vxc = np.ascontiguousarray(np.asarray(out[1]).reshape(nv, ngrid))
        fxc = np.ascontiguousarray(np.asarray(out[2]).reshape(nv, nv, ngrid))
        if with_kxc:
            kxc = np.ascontiguousarray(np.asarray(out[3]).reshape(nv, nv, nv, ngrid))
        else:
            kxc = np.zeros((nv, nv, nv, 1))
        libs.lib_func('mrsf_xc_g_set', (int(ngrid), nv, vxc.T, fxc.T, np.asfortranarray(kxc.T), bool(with_kxc)))
        self.g_kernel, self.g_gidx, self.g_with_kxc = True, key, bool(with_kxc)
        self.g_arrays = (rho, vxc, fxc, kxc if with_kxc else None)
        self.time_init_kernel += time.time() - t0

    def fxc_potential_ao(self, dm):
        """AO matrix of the linear XC response f_xc[rho0] rho_t for the total
        trial density matrix dm at the density of the installed closed-shell
        kernel set (PySCF spin-0 convention: this is the first-order change
        of the closed-shell potential for the total density change rho_t)"""
        rho0, vxc, fxc = self.g_arrays[0], self.g_arrays[1], self.g_arrays[2]
        ni, mol, grids, xc = self.ni, self.mol, self.mf.grids, self.mf.xc
        dm = 0.5*(dm + dm.T)
        v = ni.nr_rks_fxc(mol, grids, xc, None, dm, 0, 1, rho0, vxc, fxc)
        return np.asarray(v)

    def install_cs_kernel(self, fac=2.0):
        """closed-shell kernel array of the CV block of the extended method
        from the installed G set: w (f_aa + f_ab) = fac * w f^(0) (fac = 2)"""
        libs.lib_func('mrsf_xc_cs_from_g', (float(fac),))

    def vxc_g_ao(self):
        """AO matrix of the XC potential of the G density from the cached AO
        values and the installed G set (PySCF nr_rks convention)"""
        V = mrsf_grad.fzeros(self.nao, self.nao)
        libs.lib_func('mrsf_xc_g_vxc_ao', (V,))
        return V

    def potential_g(self, LM, RM, LS):
        """G-density kernel potentials of the extended method: returns
        (MM, MK, V1) with MM, MK (nocca, nao) the hole-MO projections of
        v[f rho(M)] and v[k rho(D_s)^2] (V(H, all) = MM C) and V1 (nao, nao)
        the AO matrix of v[f rho(D_s)]; M = LM RM^T + RM LM^T, D_s =
        LS C_C^T + C_C LS^T (LS has one column per doubly occupied MO)"""
        nocca = self.dims[0]
        LM = np.asfortranarray(LM, dtype=np.float64); RM = np.asfortranarray(RM, dtype=np.float64)
        LS = np.asfortranarray(LS, dtype=np.float64)
        kM, kS = LM.shape[1], LS.shape[1]
        MM = mrsf_grad.fzeros(nocca, self.nao); MK = mrsf_grad.fzeros(nocca, self.nao)
        V1 = mrsf_grad.fzeros(self.nao, self.nao)
        libs.lib_func('mrsf_xc_g_potential', (kM, LM, RM, kS, LS, MM, MK, V1))
        return MM, MK, V1

    def potential_g_q(self, Bexp, termM, termS, want_k):
        """G-density kernel potentials of the extended method for M and D_s
        given as quadratic forms (nb, idx, S) in the basis [hole MOs | Bexp]:
        returns MM, MK (nocca, nao) (projections of v[f rho(M)] and
        v[k rho(D_s)^2] on the hole MOs), VE = V1 Bexp (nao, nexp) and
        VH = V1 C_H (nao, nocca) with V1 = v[f rho(D_s)] (AO matrix)"""
        nocca = self.dims[0]
        nao = self.nao
        Bexp = np.asfortranarray(np.asarray(Bexp, dtype=np.float64).reshape(nao, -1))
        nexp = Bexp.shape[1]
        if nexp == 0:
            Bexp = mrsf_grad.fzeros(nao, 1)
        nbM, idxM, SM = termM; nbS, idxS, SS = termS
        maxnb = max(int(nbM), int(nbS), 1)
        iM = np.zeros(maxnb, dtype=np.int32); iM[:nbM] = idxM
        iS = np.zeros(maxnb, dtype=np.int32); iS[:nbS] = idxS
        SMf = mrsf_grad.fzeros(maxnb, maxnb); SMf[:nbM, :nbM] = SM
        SSf = mrsf_grad.fzeros(maxnb, maxnb); SSf[:nbS, :nbS] = SS
        MM = mrsf_grad.fzeros(nocca, nao); MK = mrsf_grad.fzeros(nocca, nao)
        VE = mrsf_grad.fzeros(nao, max(nexp, 1)); VH = mrsf_grad.fzeros(nao, nocca)
        libs.lib_func('mrsf_xc_g_potential_q', (nexp, Bexp, maxnb, int(nbM), iM, SMf, int(nbS), iS, SSf,
                                                bool(want_k), MM, MK, VE, VH))
        return MM, MK, VE[:, :nexp], VH

    def set_hole_weights(self, lam):
        """base-density weights n_i/2 of the hole MOs (the probe's D_base =
        2 C_H diag(lam) C_H^T): the published pattern (1 on the first
        nocca-1, 0 on O2) keeps the original code path"""
        lam = np.ascontiguousarray(np.asarray(lam, dtype=np.float64).reshape(-1))
        libs.lib_func('mrsf_xc_set_hole_weights', (int(lam.shape[0]), lam))

    def clear_hole_weights(self):
        libs.lib_func('mrsf_xc_clear_hole_weights', ())

    def set_dgamma(self, dg):
        """SOMO-pair density Dgamma (2x2, hole MOs nocca-1, nocca) of the
        linearised XC term of the Pi parent"""
        libs.lib_func('mrsf_xc_set_dgamma', (np.asfortranarray(np.asarray(dg, dtype=np.float64).reshape(2, 2)),))

    def clear_dgamma(self):
        libs.lib_func('mrsf_xc_clear_dgamma', ())

    def potential_g_q2(self, Bexp, termM, termS, termT, want_k, want_t):
        """potential_g_q plus the T^cv terms of the Pi parent: returns
        MM, MK, MT, MK2 (nocca, nao) (projections of v[f rho(M)],
        v[k rho(D_s)^2], v[f rho(T)] and v[k rho(T) rho(Dgamma)] on the
        hole MOs), VE = V1 Bexp and VH = V1 C_H with V1 = v[f rho(D_s)]"""
        nocca = self.dims[0]
        nao = self.nao
        Bexp = np.asfortranarray(np.asarray(Bexp, dtype=np.float64).reshape(nao, -1))
        nexp = Bexp.shape[1]
        if nexp == 0:
            Bexp = mrsf_grad.fzeros(nao, 1)
        nbM, idxM, SM = termM; nbS, idxS, SS = termS; nbT, idxT, ST = termT
        maxnb = max(int(nbM), int(nbS), int(nbT), 2)
        iM = np.zeros(maxnb, dtype=np.int32); iM[:nbM] = idxM
        iS = np.zeros(maxnb, dtype=np.int32); iS[:nbS] = idxS
        iT = np.zeros(maxnb, dtype=np.int32); iT[:nbT] = idxT
        SMf = mrsf_grad.fzeros(maxnb, maxnb); SMf[:nbM, :nbM] = SM
        SSf = mrsf_grad.fzeros(maxnb, maxnb); SSf[:nbS, :nbS] = SS
        STf = mrsf_grad.fzeros(maxnb, maxnb); STf[:nbT, :nbT] = ST
        MM = mrsf_grad.fzeros(nocca, nao); MK = mrsf_grad.fzeros(nocca, nao)
        MT = mrsf_grad.fzeros(nocca, nao); MK2 = mrsf_grad.fzeros(nocca, nao)
        VE = mrsf_grad.fzeros(nao, max(nexp, 1)); VH = mrsf_grad.fzeros(nao, nocca)
        libs.lib_func('mrsf_xc_g_potential_q2', (nexp, Bexp, maxnb, int(nbM), iM, SMf, int(nbS), iS, SSf,
                                                 int(nbT), iT, STf, bool(want_k), bool(want_t),
                                                 MM, MK, MT, MK2, VE, VH))
        return MM, MK, MT, MK2, VE[:, :nexp], VH

    def print_timings(self):
        """library-side timers of the grid kernel (diagnostics)"""
        libs.lib_func('mrsf_xc_timings', ())

    def free(self):
        libs.lib_func('mrsf_xc_free', ())


# ----------------------------------------------------------------------
# per-SCF grid session: one library grid initialisation (AO cache, hole
# MOs; kernels installed lazily), the closed-shell Fock matrix of the
# configuration G and the G-density kernel arrays, shared by the extended
# energy sections, the stability check and the gradients of one
# SCF object. The grid module is a process singleton, so one session is
# live at a time; a new geometry or new orbitals (stamp) replace it.
# ----------------------------------------------------------------------
_live = {'session': None}


def _stamp(ci):
    from graci.core.ao2mo import mo_stamp
    return mo_stamp(ci.scf.mol.crds, ci.mos)


def get_session(ci, mem_budget, block_bytes=1.0e6):
    """the grid session of the SCF object of the Mrsftddft object ci
    (created on first use with the given memory budget for the AO cache)"""
    key = (ci.scf.label, _stamp(ci))            # geometry + orbitals: the same SCF data, whatever the object
    sess = _live['session']
    if sess is not None and sess.key == key:
        return sess
    if sess is not None:
        sess.free()
    sess = GridSession(ci, key, mem_budget, block_bytes)
    _live['session'] = sess
    return sess


def free_session():
    sess = _live['session']
    if sess is not None:
        sess.free()
    _live['session'] = None


class GridSession:

    def __init__(self, ci, key, mem_budget, block_bytes):
        scf = ci.scf
        self.scf, self.key, self.valid = scf, key, True
        self.mf = scf.pyscf_obj()
        C = np.asarray(ci.mos)
        occ = np.asarray(ci.occ_ref, dtype=float)
        self.C, self.nmo, self.nao = C, C.shape[1], C.shape[0]
        self.Cidx = np.where(np.abs(occ - 2.0) < 1e-6)[0]
        self.Oidx = np.where(np.abs(occ - 1.0) < 1e-6)[0]
        self.Vidx = np.where(np.abs(occ) < 1e-6)[0]
        self.H = np.concatenate([self.Cidx, self.Oidx])
        self.P = np.concatenate([self.Oidx, self.Vidx])
        self.nocca, self.nvirb = len(self.H), len(self.P)
        self.naux = int(ci.naux)
        # the configuration G = C u {O1}: the first nocca - 1 hole MOs
        self.gidx = [int(c) for c in self.Cidx] + [int(self.Oidx[0])]
        self.occ_a = (occ > 0.5).astype(float)
        self.occ_b = (occ > 1.5).astype(float)
        self.is_dft = str(scf.xc).lower() != 'hf'
        self.x_ref = float(scf.hyb[2])
        self.xcgrid = None
        self.time_init = 0.0
        if self.is_dft:
            t0 = time.time()
            self.xcgrid = XCGrid(self.mf, C, self.occ_a, self.occ_b, C[:, self.H], C[:, self.P],
                                 (self.nocca, self.nvirb, self.naux), mem_budget, block_bytes, lazy=True)
            self.time_init = time.time() - t0
        self._fock_g = None
        self._fock_ens = None
        self.time_fock_g = 0.0

    @property
    def cache(self):
        return self.xcgrid is not None and self.xcgrid.cache

    def ensure_g_kernel(self, with_kxc):
        """the G-density kernel set (v, f and, when with_kxc, k) of this SCF"""
        self.xcgrid.ensure_g_kernel(self.gidx, with_kxc, g_default=self.gidx)

    def fock_g(self, ci):
        """closed-shell KS (or HF) matrix of the configuration G evaluated
        with the triplet orbitals (eq. 17 of Oh et al.), MO basis; built
        once per session from the resident DF blocks and the cached grid
        (PySCF build when the library data are not available)"""
        if self._fock_g is not None:
            return self._fock_g
        t0 = time.time()
        F = None
        ierr = mrsf_grad.grad_init(ci.eri_file())
        if ierr == 0 and (self.xcgrid is None or self.xcgrid.cache):
            F = self._fock_g_library()
        if F is None:
            from graci.methods.mrsftddft import closed_shell_fock_pyscf
            F = closed_shell_fock_pyscf(self.scf, self.mf, self.C, self.gidx)[1]
        self._fock_g = F
        self.time_fock_g = time.time() - t0
        return F

    def occ_ens(self):
        """total occupations of the ensemble density 2 P_C + P_O1 + P_O2"""
        occ = np.zeros(self.nmo)
        occ[self.Cidx] = 2.0
        occ[self.Oidx] = 1.0
        return occ

    def ensure_kernel_ens(self, with_kxc=False):
        """the closed-shell kernel set at the ensemble density rho_ens"""
        self.xcgrid.ensure_kernel_occ(self.occ_ens(), with_kxc)

    def fxc_potential_ao(self, dm):
        return self.xcgrid.fxc_potential_ao(dm)

    def fock_ens(self):
        """closed-shell KS (or HF) matrix of the ensemble density D_ens =
        2 P_C + P_O1 + P_O2 evaluated with the triplet orbitals (MO basis):
        the invariant one-body operator of EMRSF-Pi; library build from the
        resident DF blocks (J from the Coulomb vector, K from the hole-MO
        exchange getter: K[D_ens] = K[C] + K[C u O]) and the cached grid"""
        if self._fock_ens is not None:
            return self._fock_ens
        t0 = time.time()
        nmo, naux, H, P = self.nmo, self.naux, self.H, self.P
        dims = (self.nocca, self.nvirb, naux)
        if not (self.xcgrid is None or self.xcgrid.cache):
            sys.exit('\n ERROR: EMRSF-Pi requires the AO values cached on the grid '
                     '(raise mem_budget in the $mrsftddft section)')
        dq = mrsf_grad.dq(self.occ_ens(), naux)
        JHH, JHP, JPP = mrsf_grad.jblocks(dq, dims)
        J = np.zeros((nmo, nmo))
        J[np.ix_(H, H)] = JHH[:, :, 0]
        J[np.ix_(P, P)] = JPP[:, :, 0]
        J[np.ix_(H, P)] = JHP[:, :, 0]
        J[np.ix_(P, H)] = JHP[:, :, 0].T
        nC = len(self.Cidx)
        K = mrsf_grad.kg(nC, nmo) + mrsf_grad.kg(self.nocca, nmo)
        C = self.C
        F = C.T @ self.mf.get_hcore() @ C + J - 0.5*self.x_ref*K
        if self.is_dft:
            self.ensure_kernel_ens(False)
            F = F + C.T @ self.xcgrid.vxc_g_ao() @ C
        self._fock_ens = F
        self.time_fock_g += time.time() - t0
        return F

    def _fock_g_library(self):
        nmo, naux, H, P = self.nmo, self.naux, self.H, self.P
        dims = (self.nocca, self.nvirb, naux)
        occG = np.zeros(nmo); occG[self.gidx] = 2.0
        dqG = mrsf_grad.dq(occG, naux)
        JHH, JHP, JPP = mrsf_grad.jblocks(dqG, dims)
        J = np.zeros((nmo, nmo))
        J[np.ix_(H, H)] = JHH[:, :, 0]
        J[np.ix_(P, P)] = JPP[:, :, 0]
        J[np.ix_(H, P)] = JHP[:, :, 0]
        J[np.ix_(P, H)] = JHP[:, :, 0].T
        K = 2.0*mrsf_grad.kg(len(self.gidx), nmo)                 # K[D_G]: the G orbitals are doubly occupied
        C = self.C
        F = C.T @ self.mf.get_hcore() @ C + J - 0.5*self.x_ref*K
        if self.is_dft:
            self.ensure_g_kernel(False)
            F = F + C.T @ self.xcgrid.vxc_g_ao() @ C
        return F

    def free(self):
        if self.xcgrid is not None:
            self.xcgrid.free()
            self.xcgrid = None
        self._fock_g = None
        self.valid = False
