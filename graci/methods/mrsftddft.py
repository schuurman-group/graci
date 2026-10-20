"""
Module for computing MRSF-TDDFT (mixed-reference spin-flip TDDFT)
energies, densities and transition densities from an ROKS/ROHF triplet
reference
"""
import copy as copy
import sys as sys
import os
import numpy as np
import graci.utils.timing as timing
import graci.methods.cimethod as cimethod
import graci.core.params as params
import graci.io.output as output
import graci.interfaces.mrsf.mrsf_init as mrsf_init
import graci.interfaces.mrsf.mrsf_diag as mrsf_diag
import graci.interfaces.mrsf.mrsf_density as mrsf_density
import graci.interfaces.mrsf.mrsf_overlap as mrsf_overlap
import graci.interfaces.mrsf.mrsf_stability as mrsf_stability

def closed_shell_fock_pyscf(scf, mf, C, gidx):
    """closed-shell KS (or HF) matrix of the configuration G = |C O1^2|
    (orbitals gidx doubly occupied) evaluated with the triplet orbitals C
    (eq. 17 of Oh et al.): returns (F_AO, F_MO). The same functional,
    grid and density fitting as the SCF object of the reference."""
    from pyscf import scf as pyscf_scf, dft as pyscf_dft
    D_G = 2.0 * C[:, gidx] @ C[:, gidx].T
    xc = str(scf.xc).lower()
    pymol = mf.mol
    if xc == 'hf':
        ks = pyscf_scf.ROHF(pymol)
    else:
        ks = pyscf_dft.ROKS(pymol)
        ks.xc = mf.xc
        ks.grids = mf.grids
    if scf.mol.use_df:
        ks = ks.density_fit(auxbasis=scf.mol.ri_basis)
    # an ROHF/ROKS object with a total density gives the closed-shell
    # Fock matrix of that density (alpha and beta Fock matrices equal)
    f = ks.get_fock(dm=D_G)
    fa = getattr(f, 'focka', None)
    f_ao = np.asarray(f) if fa is None else np.asarray(fa)
    return f_ao, C.T @ f_ao @ C


closed_shell_fock = closed_shell_fock_pyscf        # reference build (PySCF); the session builds it in the library


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
        # EMRSF F' correction: 'covariant' (default: full VV and CC blocks
        # of (1 - c_H)(J[D_O2] - J[D_O1]), invariant under core-core and
        # virtual-virtual rotations) or 'diagonal' (the published form,
        # energies only)
        self.fprime         = 'covariant'
        # parent configuration of the CV excitations: 'published' = the
        # closed-shell configuration G = |C O1^2| (Oh et al.; the energies
        # depend on the choice of the two SOMOs), 'pi' = EMRSF-Pi, the
        # invariant lowest pair state of the MRSF pair block (energies
        # invariant under O1-O2 rotations); 'kernel' = XC kernel at the
        # ensemble density and unscaled Coulomb term in the CV block (the
        # published structure), 'kappa' = kernel-free with the Coulomb term
        # scaled by kappa_coov
        self.ext_space      = 'published'
        self.ext_model      = 'kernel'
        # internal stability check of the ROKS triplet reference (an
        # unstable reference is re-optimised along the unstable direction)
        self.stability      = False
        # frozen core: number of doubly occupied MOs (the lowest by orbital
        # energy) excluded from the hole set of the response space; None
        # (default) = the chemical core of the molecule, 0 = none
        self.frozen_core    = None
        # exchange term of the sigma vectors: 'auto' (explicit (ij|ab)
        # tensor when it fits mem_budget, else the DF sweep), 'tensor', 'df'
        self.exchange       = 'auto'
        # Davidson correction equations: 'inner' = Jacobi-Davidson
        # corrections solved with the kernel-free operator (active for
        # singlets with the CV kernel), 'diag' = diagonal preconditioner
        self.precond        = 'inner'
        # MINRES steps per inner correction equation
        self.inner_iter     = 6

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
        # EMRSF-Pi: parent weights and pair-block eigenvalues
        self.pi_info        = None
        self.fock_mo        = None
        # compressed amplitude vectors, one (xdim, nroots) array per
        # irrep, for each representation
        self.amps           = {'adiabatic' : None}
        # Davidson convergence information per irrep
        self.diag_info      = None
        # all tracked Davidson roots per irrep (nroots + nextra Ritz
        # vectors and values): the start vectors of a warm-started run at
        # the next geometry
        self.guess_vecs     = None
        self.guess_ener     = None

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

        # guess: a previous Mrsftddft object (another geometry) whose
        # tracked Davidson vectors, aligned to the present MOs, start the
        # Davidson iterations (warm start); validated in warm_vectors
        scf_energy = self.set_scf(scf, ci_guess=None, mo_ints=mo_ints)

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

        # optional internal stability check of the reference: an unstable
        # reference is re-optimised in place (shared Scf object) and the
        # integrals are transformed again
        if self.stability:
            mrsf_stability.ensure_stable(self, mo_ints)
            fock_ao = self.scf.fock_ao
            fock_mo = np.array([self.mos.T @ fock_ao[0] @ self.mos,
                                self.mos.T @ fock_ao[1] @ self.mos])

        # initialise the library and load the integrals
        mrsf_init.set_extended(self.extended)
        mrsf_init.set_frozen_core(self.frozen_core)
        mrsf_init.init(self, fock_mo)
        mrsf_init.set_exchange(self.exchange)
        mrsf_init.init_ints(self, self.eri_file(), mo_ints=mo_ints)
        self.nocca, self.nvirb, self.xdim, self.naux = \
            [int(x) for x in mrsf_init.dims()]
        self.ncol, self.ncv = [int(x) for x in mrsf_init.ext_dims()]

        # warm start: the vectors of the guess object aligned to these MOs
        gvecs = self.warm_vectors(guess)

        # extended method: one-body operators, parent weights, kernel cache
        xcgrid = None
        if self.extended:
            self.fock_mo = fock_mo
            ext_args, xcgrid = self.extended_setup()
            mrsf_init.init_ext(self, *ext_args, self.eri_file())

        # Davidson diagonalisation, irrep by irrep
        nirr = self.n_irrep()
        maxroots = int(max(self.nstates))
        self.energies_sym = np.zeros((nirr, max(maxroots, 1)), dtype=float)
        self.amps['adiabatic'] = [np.zeros((self.xdim, 0), dtype=float)
                                  for i in range(nirr)]
        self.guess_vecs = [np.zeros((self.xdim, 0), dtype=float)
                           for i in range(nirr)]
        self.guess_ener = [np.zeros(0, dtype=float) for i in range(nirr)]
        self.diag_info = [None for i in range(nirr)]

        for irrep in self.irreps_nonzero():
            nroots, ener, xvec, niter, iconv, tracked = \
                mrsf_diag.diag(self, irrep, guess=gvecs)
            if nroots < self.nstates[irrep]:
                output.print_message('  irrep '+str(irrep)+': number of '
                                     'roots reduced to '+str(nroots))
                self.nstates[irrep] = nroots
            self.energies_sym[irrep, :nroots] = self.scf.energy + ener
            self.amps['adiabatic'][irrep] = xvec
            self.guess_ener[irrep] = self.scf.energy + tracked[0]
            self.guess_vecs[irrep] = tracked[1]
            self.diag_info[irrep] = (niter, iconv)
            if iconv != 1:
                output.print_message('  WARNING: Davidson not converged '
                                     'for irrep '+str(irrep))
        if gvecs is not None:
            output.print_message('  warm start: Davidson iterations per irrep: '
                                 + ' '.join(str(self.diag_info[i][0])
                                            for i in self.irreps_nonzero()))

        # energies sorted by value, with the corresponding states
        self.order_energies()
        n_tot = self.n_states()

        # the grid session (AO cache, kernels, F^G) stays alive for the
        # other sections and the gradients of this SCF object

        # state density matrices: expectation values <Psi|E_pq|Psi> of the
        # (extended) configuration expansion; for the extended method the
        # C->V configurations are the singles of the parent (G for
        # ext_space = published, Pi = w_G G + w_D D + w_S S for pi), whose
        # weights the density kernel receives
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

        # print the states
        if self.verbose:
            output.print_mrsftddft_states(self)

        # store the density matrices in adiabatic energy order
        self.dmats['adiabatic'] = np.zeros((n_tot, self.nmo, self.nmo),
                                           dtype=float)
        for istate in range(n_tot):
            irr, st = self.state_sym(istate)
            self.dmats['adiabatic'][istate, :, :] = dmat_sym[irr][:, :, st]

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
    def warm_vectors(self, guess):
        """
        start vectors for the Davidson iterations from a previous Mrsftddft
        object of the same manifold (in memory, or read from the checkpoint
        when guess_label names an object that is not in memory): all its
        tracked Davidson vectors (guess_vecs; the converged amplitudes of
        objects from older checkpoints), in ascending energy, aligned to the
        present MOs (orthogonal Procrustes rotations of the MO overlap
        blocks within the orbital classes, the SOMO pair only sign-aligned
        for extended objects). Returns the (xdim, nguess) array or None.
        """
        from pyscf import gto
        if guess is None and self.guess_label is not None:
            import graci.io.chkpt as chkpt
            try:
                guess = chkpt.read('Mrsftddft.' + str(self.guess_label))
            except Exception:
                guess = None
        if guess is None:
            return None
        if type(guess).__name__ != 'Mrsftddft' or guess.amps.get('adiabatic') is None:
            output.print_message('  warm start: the guess object ' + str(self.guess_label)
                                 + ' holds no MRSF-TDDFT vectors; cold start')
            return None
        for attr in ('mult', 'nmo', 'nocca', 'nvirb', 'nel', 'ncol'):
            if getattr(guess, attr, None) != getattr(self, attr):
                output.print_message('  warm start: the guess object ' + str(guess.label)
                                     + ' differs in ' + attr + '; cold start')
                return None
        # all tracked vectors of the guess object, ascending energy
        blocks, eners = [], []
        gv = getattr(guess, 'guess_vecs', None)
        for irr in range(guess.n_irrep()):
            if gv is not None and gv[irr] is not None and np.shape(gv[irr])[1] > 0:
                blocks.append(np.asarray(gv[irr]))
                eners.append(np.asarray(guess.guess_ener[irr]))
            else:
                x = np.asarray(guess.amps['adiabatic'][irr])
                if x.shape[1] == 0:
                    continue
                blocks.append(x)
                eners.append(np.asarray(guess.energies_sym[irr, :x.shape[1]]))
        if not blocks:
            return None
        xk = np.concatenate(blocks, axis=1)
        order = np.argsort(np.concatenate(eners))
        xk = xk[:, order]
        if xk.shape[0] != self.xdim:
            output.print_message('  warm start: the guess vectors have the wrong dimension; cold start')
            return None
        Sao = gto.intor_cross('int1e_ovlp', self.scf.mol.mol_obj, guess.scf.mol.mol_obj)
        xk, dmin, ierr = mrsf_overlap.align_vectors(Sao, self.mos, guess.mos, self.occ_ref,
                                                    guess.occ_ref, self.mult, self.extended, 0, xk)
        if ierr in (1, 4, 8):
            output.print_message('  warm start: MO alignment failed (code ' + str(ierr) + '); cold start')
            return None
        msg = ('  warm start from ' + str(guess.label) + ': ' + str(xk.shape[1])
               + ' vectors; smallest aligned MO overlap (core/SOMO/virtual): '
               + ' '.join('{:.4f}'.format(v) for v in dmin))
        if ierr == 5:
            msg += '  (WARNING: the SOMO pair has rotated or swapped between the geometries)'
        output.print_message(msg)
        return xk

    #
    def extended_setup(self):
        """
        data of the extended method from the per-SCF grid session (shared
        with the other sections and the gradients of the same SCF object):
        the arguments of mrsf_init.init_ext (one-body operator of the CV
        block, Fock-type coupling operators per parent, parent weights,
        shift, Coulomb scale, mode) and the grid-kernel cache of the singlet
        CV block (None when no kernel is used). Published model: the
        closed-shell KS matrix of G (eq. 17 of Oh et al.), the kernel at
        the G density; EMRSF-Pi: the invariant operators of mrsf_pi.py and
        the kernel at the ensemble density.
        """
        import graci.interfaces.mrsf.mrsf_xc as mrsf_xc
        import graci.interfaces.mrsf.mrsf_pi as mrsf_pi

        sess = mrsf_xc.get_session(self, self.mem_budget,
                                   block_bytes=getattr(self, 'kernel_block_bytes', 1.0e6))
        xc    = str(self.scf.xc).lower()
        space = str(self.ext_space).lower()
        model = str(self.ext_model).lower()
        self.use_kernel = bool(self.mult == 1 and xc != 'hf' and model == 'kernel')
        xcgrid = None
        if self.use_kernel:
            xcgrid = sess.xcgrid
            if not xcgrid.cache:
                need = 8.0*xcgrid.ncomp*xcgrid.ngrid*(xcgrid.nao + self.nocca)/1.e9
                sys.exit('\n ERROR: the extended MRSF-TDDFT kernel requires '
                         'the AO values cached on the grid: set mem_budget '
                         'to at least '+'{:.2f}'.format(1.05*need)+' GB')
        coul = 1.0 if model == 'kernel' else float(self.spc_eff[2])
        if space == 'published':
            fone = sess.fock_g(self)
            fcpl = np.stack([fone, fone, fone], axis=2)
            args = (fone, fcpl, np.array([1.0, 0.0, 0.0]), 0.0, coul, 1)
            if self.use_kernel:
                sess.ensure_g_kernel(False)
                xcgrid.install_cs_kernel(2.0)
            self.pi_info = None
        else:
            # the DF blocks of the gradient module (SOMO columns, Coulomb
            # vectors, hole exchange) serve the operator builds
            import graci.interfaces.mrsf.mrsf_grad as mrsf_grad
            ierr = mrsf_grad.grad_init(self.eri_file())
            if ierr != 0:
                sys.exit('\n ERROR: EMRSF-Pi needs the double-precision DF blocks '
                         '(gradient module initialisation failed, code '+str(ierr)+')')
            pi = mrsf_pi.pi_setup(sess, self.fock_mo[0], self.fock_mo[1], float(self.chf_eff))
            if os.environ.get('GRACI_PI_WFIX'):
                # diagnostic: the conventions of EMRSF-Pi on the pure-G parent
                wfix = np.array([1.0, 0.0, 0.0])
                pi['shift'] = float(pi['Hp'][0, 0])
                pi['fone'] = mrsf_pi.pi_fone_for_weights(sess, pi, wfix, float(self.chf_eff))
                pi['w'] = wfix
                output.print_message('  EMRSF-Pi diagnostic: parent weights forced to (1, 0, 0)')
            args = (pi['fone'], pi['fcpl'], pi['w'], pi['shift'], coul, 2)
            if self.use_kernel:
                sess.ensure_kernel_ens(False)
                xcgrid.install_cs_kernel(2.0)
            self.pi_info = {'w': pi['w'], 'pair': pi['pair']}
            output.print_message('  EMRSF-Pi: pair block eigenvalues (Ha) '
                                 + ' '.join('{:10.6f}'.format(x) for x in pi['pair'])
                                 + ', parent weights (G, D, S) '
                                 + ' '.join('{:8.5f}'.format(x) for x in pi['w']))
        return args, xcgrid

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

        ndocc = int(np.sum(np.abs(self.occ_ref - 2.) < 1.e-6))
        if self.frozen_core is None:
            # default: the chemical core (1 per Li-Ne atom, 5 per Na-Ar, ...)
            from pyscf.data import elements
            self.frozen_core = min(int(elements.chemcore(self.scf.mol.mol_obj)),
                                   max(ndocc - 1, 0))
        self.frozen_core = int(self.frozen_core)
        if self.frozen_core < 0 or self.frozen_core >= ndocc:
            sys.exit('\n ERROR: frozen_core must be >= 0 and smaller than the '
                     'number of doubly occupied MOs ('+str(ndocc)+')')

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
        self.exchange = str(self.exchange).lower()
        if self.exchange not in ('auto', 'tensor', 'df'):
            sys.exit('\n ERROR: exchange must be auto, tensor or df')
        self.precond = str(self.precond).lower()
        if self.precond not in ('diag', 'inner'):
            sys.exit('\n ERROR: precond must be diag or inner')
        self.inner_iter = int(self.inner_iter)
        if self.inner_iter < 1:
            sys.exit('\n ERROR: inner_iter must be at least 1')
        if self.extended and self.precision != 'double':
            sys.exit('\n ERROR: the extended MRSF-TDDFT method requires '
                     'precision = double')
        self.fprime = str(self.fprime).lower()
        if self.fprime not in ('diagonal', 'covariant'):
            sys.exit('\n ERROR: fprime must be diagonal or covariant')
        self.ext_space = str(self.ext_space).lower()
        if self.ext_space not in ('published', 'pi'):
            sys.exit('\n ERROR: ext_space must be published or pi')
        self.ext_model = str(self.ext_model).lower()
        if self.ext_model not in ('kernel', 'kappa'):
            sys.exit('\n ERROR: ext_model must be kernel or kappa')

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
        if (str(getattr(self, 'ext_space', 'published')).lower() !=
                str(getattr(ket, 'ext_space', 'published')).lower()):
            sys.exit('\n ERROR: MRSF transition densities require the same '
                     'extension space (ext_space) for the bra and ket objects')
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
        class before the evaluation. Standard and extended objects can be
        combined (a standard object is an extended one without CV
        amplitudes). Returns (overlap list, ierr).
        """
        for obj in (self, ket):
            if getattr(obj, 'extended', False) and str(getattr(obj, 'ext_space', 'published')).lower() == 'pi':
                sys.exit('\n ERROR: state overlaps are not yet available for ext_space = pi')

        if type(ket).__name__ != 'Mrsftddft':
            sys.exit('\n ERROR: MRSF overlaps require two Mrsftddft objects')
        for attr in ('nmo', 'nocca', 'nvirb', 'nel'):
            if getattr(self, attr) != getattr(ket, attr):
                sys.exit('\n ERROR: MRSF overlaps require bra and ket objects '
                         'with the same MO space and reference occupation '
                         'pattern (' + attr + ' differs)')
        return mrsf_overlap.overlap(self, ket, trans_list_sym, method, align, rep)
