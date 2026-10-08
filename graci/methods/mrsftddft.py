"""
Module for computing MRSF-TDDFT (mixed-reference spin-flip TDDFT)
energies, densities and transition densities from an ROKS/ROHF triplet
reference
"""
import copy as copy
import sys as sys
import numpy as np
import graci.utils.timing as timing
import graci.methods.cimethod as cimethod
import graci.core.params as params
import graci.io.output as output
import graci.interfaces.mrsf.mrsf_init as mrsf_init
import graci.interfaces.mrsf.mrsf_diag as mrsf_diag
import graci.interfaces.mrsf.mrsf_density as mrsf_density
import graci.interfaces.mrsf.mrsf_overlap as mrsf_overlap

class Mrsftddft(cimethod.Cimethod):
    """Class constructor for MRSF-TDDFT objects"""
    def __init__(self, ci_obj=None):
        # parent attributes
        super().__init__()

        # user defined quantities
        # target manifold: 1 (singlets) or 3 (triplets)
        self.mult           = 1
        # extra roots tracked by the Davidson solver
        self.nextra         = 3
        # Davidson residual-norm threshold
        self.diag_tol       = 1.e-5
        # maximum number of Davidson iterations
        self.diag_iter      = 100
        # maximum subspace dimension per tracked root
        self.diag_maxvec    = 20
        # fraction of HF exchange (None: taken from the functional)
        self.hfx            = None
        # spin-pair coupling scaling factors (coco, ovov, coov),
        # None: all equal to the fraction of HF exchange
        self.spc            = None
        # amplitude threshold for printing
        self.conf_thresh    = 0.05
        # integral precision
        self.precision      = 'double'
        # memory budget (GB) for the sigma-vector work arrays
        self.mem_budget     = 1.0
        # storage of the virtual-virtual DF block: 'paired' or 'full'
        self.vv_storage     = 'paired'
        # MO energy cutoff: all MOs are kept by default
        self.mo_cutoff      = 1.e6
        # keep the integrals loaded in the library after the run
        self.keep_ints      = True
        # extended MRSF-TDDFT (EMRSF-TDDFT): add the core-to-virtual
        # configurations of the closed-shell configuration G = |C O1^2|
        self.extended       = False
        # EMRSF coupling scale c_cp (None: the fraction of HF exchange)
        self.ccp            = None

        # computed quantities
        # reference occupation vector (first nmo MOs)
        self.occ_ref        = None
        # effective fraction of HF exchange
        self.chf_eff        = None
        # effective spin-pair couplings
        self.spc_eff        = None
        # effective EMRSF coupling scale and kernel flag
        self.ccp_eff        = None
        self.use_kernel     = False
        # dimensions of the response space: the amplitude vector is the
        # (nvirb, ncol) matrix [x | y], ncol = nocca (+ nC when extended),
        # xdim = nvirb*ncol; ncv = nvirb*nC CV slots (0 if not extended)
        self.nocca          = None
        self.nvirb          = None
        self.ncol           = None
        self.ncv            = None
        self.xdim           = None
        self.naux           = None
        # weight of the CV configurations per adiabatic state (extended)
        self.gamma_cv       = None
        # compressed amplitude vectors, one (xdim, nroots) array per
        # irrep, for each representation
        self.amps           = {'adiabatic' : None}
        # Davidson convergence information per irrep
        self.diag_info      = None

        if isinstance(ci_obj, cimethod.Cimethod):
            for name,obj in ci_obj.__dict__.items():
                if hasattr(self, name):
                    if isinstance(obj, (dict, list, np.ndarray)):
                        setattr(self, name, copy.deepcopy(obj))
                    elif obj.__class__.__name__ in params.valid_objs:
                        setattr(self, name, obj.copy())
                    else:
                        setattr(self, name, obj)

# Required functions #############################################################
    def copy(self):
        """create of deepcopy of self"""
        new = Mrsftddft()

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))

        return new

    @timing.timed
    def run(self, scf, guess=None, mo_ints=None):
        """compute the MRSF-TDDFT states for all irreps"""

        # the guess object is currently ignored
        guess      = None
        scf_energy = self.set_scf(scf, ci_guess=guess, mo_ints=mo_ints)

        if scf_energy is None:
            return None

        # sanity checks
        self.check_input()

        # write the output logfile header for this run
        if self.verbose:
            output.print_mrsftddft_header(self.label)
            output.print_coords(self.scf.mol.crds, self.scf.mol.asym)

        # spin Fock matrices in the MO basis of self.mos
        fock_ao = self.scf.fock_ao
        if fock_ao is None:
            fock_ao = self.scf.build_fock()
        fock_mo = np.array([self.mos.T @ fock_ao[0] @ self.mos,
                            self.mos.T @ fock_ao[1] @ self.mos])

        # fraction of HF exchange and spin-pair couplings
        self.set_exchange()

        # initialise the library and load the integrals
        mrsf_init.set_extended(self.extended)
        mrsf_init.init(self, fock_mo)
        mrsf_init.init_ints(self, self.eri_file())
        self.nocca, self.nvirb, self.xdim, self.naux = \
            [int(x) for x in mrsf_init.dims()]
        self.ncol, self.ncv = [int(x) for x in mrsf_init.ext_dims()]

        # extended method: closed-shell KS matrix of G, kernel cache
        xcgrid = None
        if self.extended:
            fdft_mo, xcgrid = self.extended_setup()
            mrsf_init.init_ext(self, fdft_mo, self.eri_file())

        # Davidson diagonalisation, irrep by irrep
        nirr = self.n_irrep()
        maxroots = int(max(self.nstates))
        self.energies_sym = np.zeros((nirr, max(maxroots, 1)), dtype=float)
        self.amps['adiabatic'] = [np.zeros((self.xdim, 0), dtype=float)
                                  for i in range(nirr)]
        self.diag_info = [None for i in range(nirr)]

        for irrep in self.irreps_nonzero():
            nroots, ener, xvec, niter, iconv = mrsf_diag.diag(self, irrep)
            if nroots < self.nstates[irrep]:
                output.print_message('  irrep '+str(irrep)+': number of '
                                     'roots reduced to '+str(nroots))
                self.nstates[irrep] = nroots
            self.energies_sym[irrep, :nroots] = self.scf.energy + ener
            self.amps['adiabatic'][irrep] = xvec
            self.diag_info[irrep] = (niter, iconv)
            if iconv != 1:
                output.print_message('  WARNING: Davidson not converged '
                                     'for irrep '+str(irrep))

        # energies sorted by value, with the corresponding states
        self.order_energies()
        n_tot = self.n_states()

        # the grid kernel cache of the extended method is no longer needed
        if xcgrid is not None:
            xcgrid.free()

        # state density matrices: expectation values <Psi|E_pq|Psi> of the
        # (extended) configuration expansion
        dmat_sym = mrsf_density.rdm(self)

        # the library (and the loaded integrals) is kept alive so that a
        # following mrsftddft section on the same SCF object (e.g. the
        # other spin manifold) does not re-read the integrals; it is
        # released by mrsf_init.finalize() when the integrals change or
        # at the end of the run
        mrsf_init.report_timings()
        if not self.keep_ints:
            mrsf_init.finalize()

        # weight of the CV configurations per state (extended method)
        if self.extended:
            self.gamma_cv = np.zeros(n_tot, dtype=float)
            for istate in range(n_tot):
                X = self.amplitudes(istate)
                self.gamma_cv[istate] = np.sum(X[:, self.nocca:]**2)

        # store the density matrices in adiabatic energy order
        self.dmats['adiabatic'] = np.zeros((n_tot, self.nmo, self.nmo),
                                           dtype=float)
        for istate in range(n_tot):
            irr, st = self.state_sym(istate)
            self.dmats['adiabatic'][istate, :, :] = dmat_sym[irr][:, :, st]

        # print the states
        if self.verbose:
            output.print_mrsftddft_states(self)

        # build the natural orbitals in AO basis by default
        self.build_nos()

        # only print if user-requested
        if self.print_orbitals:
            self.print_nos()

        # determine promotion numbers if ref_state != -1
        if self.ref_state >= 0:
            self.print_promotion(self.ref_state)

        # print the moments
        self.print_moments()

        return True

    #
    def extended_setup(self):
        """
        data of the extended method: the closed-shell KS matrix of the
        configuration G = |C O1^2| evaluated with the triplet orbitals
        (eq. 17 of Oh et al.), in the MO basis, and the grid-kernel cache
        of the singlet CV block (None when no kernel is needed)
        """
        from pyscf import scf as pyscf_scf, dft as pyscf_dft
        import graci.interfaces.mrsf.mrsf_xc as mrsf_xc

        mf   = self.scf.pyscf_obj()
        C    = np.asarray(self.mos)
        hmap, pmap = self.orbital_classes()
        nc   = len(hmap) - 2
        gidx = list(hmap[:nc]) + [int(hmap[nc])]
        D_G  = 2.0 * C[:, gidx] @ C[:, gidx].T
        xc   = str(self.scf.xc).lower()
        pymol = mf.mol
        if xc == 'hf':
            ks = pyscf_scf.ROHF(pymol)
        else:
            ks = pyscf_dft.ROKS(pymol)
            ks.xc = mf.xc
            ks.grids = mf.grids
        if self.scf.mol.use_df:
            ks = ks.density_fit(auxbasis=self.scf.mol.ri_basis)
        # an ROHF/ROKS object with a total density gives the closed-shell
        # Fock matrix of that density (alpha and beta Fock matrices equal)
        f  = ks.get_fock(dm=D_G)
        fa = getattr(f, 'focka', None)
        fdft_ao = np.asarray(f) if fa is None else np.asarray(fa)
        fdft_mo = C.T @ fdft_ao @ C

        self.use_kernel = bool(self.mult == 1 and xc != 'hf')
        xcgrid = None
        if self.use_kernel:
            occ_g = np.zeros(self.nmo, dtype=float)
            occ_g[gidx] = 1.0
            xcgrid = mrsf_xc.XCGrid(mf, C, occ_g, occ_g, C[:, hmap], C[:, pmap],
                                    (self.nocca, self.nvirb, self.naux),
                                    self.mem_budget, block_bytes=2.0e6)
            if not xcgrid.cache:
                need = 8.0*xcgrid.ncomp*xcgrid.ngrid*(xcgrid.nao + self.nocca)/1.e9
                sys.exit('\n ERROR: the extended MRSF-TDDFT kernel requires '
                         'the AO values cached on the grid: set mem_budget '
                         'to at least '+'{:.2f}'.format(1.05*need)+' GB')
        return fdft_mo, xcgrid

    #
    def check_input(self):
        """sanity checks on the reference and the input"""

        if int(self.scf.mult) != 3:
            sys.exit('\n ERROR: MRSF-TDDFT requires a triplet (mult = 3) '
                     'SCF reference, label = '+str(self.scf.label))

        if int(self.mult) not in (1, 3):
            sys.exit('\n ERROR: MRSF-TDDFT mult must be 1 or 3')

        if self.charge != self.scf.charge:
            sys.exit('\n ERROR: MRSF-TDDFT charge must equal the SCF charge')

        self.occ_ref = np.array(self.scf.orb_occ[:self.nmo], dtype=float)
        if np.sum(np.abs(self.occ_ref - 1.) < 1.e-6) != 2:
            sys.exit('\n ERROR: the MRSF-TDDFT reference must have exactly '
                     'two singly occupied MOs within the MO space')

        nirr = self.scf.mol.n_irrep()
        if self.nstates is None:
            sys.exit('\n ERROR: nstates must be given in the mrsftddft section')
        self.nstates = np.atleast_1d(np.array(self.nstates, dtype=int))
        if len(self.nstates) != nirr:
            sys.exit('\n ERROR: nstates must have one entry per irrep ('
                     +str(nirr)+')')

        if self.precision not in ('single', 'double'):
            sys.exit('\n ERROR: precision must be single or double')
        if self.vv_storage not in ('paired', 'full'):
            sys.exit('\n ERROR: vv_storage must be paired or full')
        if self.extended and self.precision != 'double':
            sys.exit('\n ERROR: the extended MRSF-TDDFT method requires '
                     'precision = double')

        return

    #
    def set_exchange(self):
        """effective fraction of HF exchange and spin-pair couplings"""

        hyb = self.scf.hyb
        if hyb is None:
            self.scf.build_fock()
            hyb = self.scf.hyb
        omega, alpha, chf = hyb

        if abs(omega) > 1.e-12:
            sys.exit('\n ERROR: range-separated functionals are not yet '
                     'supported in MRSF-TDDFT')

        if self.hfx is not None:
            chf = float(self.hfx)
        self.chf_eff = chf
        self.ccp_eff = chf if self.ccp is None else float(self.ccp)

        if self.spc is None:
            self.spc_eff = [chf, chf, chf]
        else:
            spc = np.atleast_1d(np.array(self.spc, dtype=float))
            if spc.size == 1:
                spc = np.repeat(spc, 3)
            if spc.size != 3:
                sys.exit('\n ERROR: spc must have 1 or 3 entries')
            self.spc_eff = list(spc)

        return

    #
    def eri_file(self):
        """name of the 2-e DF integral file written by Ao2mo"""
        return '2e_eri_'+str(self.scf.label).strip()+'.h5'

    #
    def orbital_classes(self):
        """
        returns (hole_map, particle_map): the MO indices (0-based) of the
        hole set [C..., O1, O2] and the particle set [O1, O2, V...]
        """
        occ  = self.occ_ref
        cidx = np.where(np.abs(occ - 2.) < 1.e-6)[0]
        oidx = np.where(np.abs(occ - 1.) < 1.e-6)[0]
        vidx = np.where(np.abs(occ) < 1.e-6)[0]
        return np.concatenate([cidx, oidx]), np.concatenate([oidx, vidx])

    #
    def amplitudes(self, istate, rep='adiabatic'):
        """compressed amplitude array (nvirb, ncol) of adiabatic state
        istate: columns :nocca the MRSF amplitudes x(a,i), columns nocca:
        the CV amplitudes y(a,i) of the extended method"""
        irr, st = self.state_sym(istate)
        x = self.amps[rep][irr][:, st]
        ncol = self.nocca if self.ncol is None else self.ncol
        return np.reshape(x, (self.nvirb, ncol), order='F')

    #
    def dominant_amplitudes(self, istate, thresh=None, rep='adiabatic'):
        """
        list of (coefficient, label) for the amplitudes of state istate
        with |c| > thresh, sorted by decreasing |c|
        """
        if thresh is None:
            thresh = self.conf_thresh

        hmap, pmap = self.orbital_classes()
        nc = len(hmap) - 2
        X  = self.amplitudes(istate, rep)
        irrlbl = self.scf.mol.irreplbl

        def molbl(p):
            return str(p+1)+'('+str(irrlbl[self.mosym[p]])+')'

        out = []
        for i in range(X.shape[1]):
            for a in range(self.nvirb):
                c = X[a, i]
                if abs(c) < thresh:
                    continue
                if i >= self.nocca:
                    # CV configuration of the extended method
                    if a < 2:
                        continue
                    lbl = molbl(hmap[i-self.nocca])+' -> '+molbl(pmap[a])+' (CV)'
                elif a == 0 and i == nc:
                    lbl = 'O1,O2 open-shell'
                elif a == 0 and i == nc+1:
                    lbl = molbl(hmap[i])+' -> '+molbl(pmap[a])+' (O2 -> O1)'
                elif a == 1 and i == nc:
                    lbl = molbl(hmap[i])+' -> '+molbl(pmap[a])+' (O1 -> O2)'
                else:
                    lbl = molbl(hmap[i])+' -> '+molbl(pmap[a])
                out.append((c, lbl))

        out.sort(key=lambda t: -abs(t[0]))
        return out

    #
    def bitci_mrci(self):
        """MRSF objects have no bitci wave function: used by Spinorbit
        and Dyson, which are not supported for MRSF-TDDFT (Transition and
        Overlap dispatch to the MRSF library through tdm_sym/overlap_sym)"""
        sys.exit('\n ERROR: Spinorbit/Dyson sections are not yet '
                 'supported for MRSF-TDDFT objects (label = '
                 +str(self.label)+')')

    def bitci_ref(self):
        return self.bitci_mrci()

    #
    def tdm_sym(self, ket, trans_list_sym, rep='adiabatic'):
        """
        1-TDMs <bra=self|E_pq|ket> for the (bra_irrep, ket_irrep)-blocked
        list of state pairs trans_list_sym; same layout as
        graci.interfaces.bitci.mrci_1tdm.tdm
        """
        if type(ket).__name__ != 'Mrsftddft':
            sys.exit('\n ERROR: MRSF transition densities require two '
                     'Mrsftddft objects')
        if bool(getattr(self, 'extended', False)) != bool(getattr(ket, 'extended', False)):
            sys.exit('\n ERROR: MRSF transition densities require two '
                     'extended or two standard MRSF-TDDFT objects')
        if np.any(self.mos != ket.mos):
            sys.exit('\n ERROR: MRSF transition densities require the '
                     'same MOs for the bra and ket objects')

        return mrsf_density.tdm(self, ket, trans_list_sym, rep)

    #
    def overlap_sym(self, ket, trans_list_sym, method='exact', align=True, rep='adiabatic'):
        """
        overlaps <bra=self|ket> between the states of this object and those
        of another Mrsftddft object ket (in general at another geometry) for
        the (bra_irrep, ket_irrep)-blocked list of state pairs
        trans_list_sym; same layout as
        graci.interfaces.bitci.wf_overlap.overlap. method = 'exact'
        (determinant-factorised two-index determinants, no truncation) or
        'tlf0' | 'tlf1' | 'tlf2' (truncated Leibniz formula, JCTC 15, 882
        (2019)); align: align the ket MOs to the bra MOs within each orbital
        class before the evaluation. Returns (overlap list, ierr).
        """
        if type(ket).__name__ != 'Mrsftddft':
            sys.exit('\n ERROR: MRSF overlaps require two Mrsftddft objects')
        if self.extended or ket.extended:
            sys.exit('\n ERROR: overlaps are not yet available for the '
                     'extended MRSF-TDDFT method')
        for attr in ('nmo', 'nocca', 'nvirb', 'nel'):
            if getattr(self, attr) != getattr(ket, attr):
                sys.exit('\n ERROR: MRSF overlaps require bra and ket objects '
                         'with the same MO space and reference occupation '
                         'pattern (' + attr + ' differs)')
        return mrsf_overlap.overlap(self, ket, trans_list_sym, method, align, rep)
