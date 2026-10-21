"""
Analytic nuclear gradients of MRSF-TDDFT states ($mrsfgradient section)
"""
import sys
import copy
import numpy as np
import graci.core.params as params
import graci.io.output as output
import graci.utils.timing as timing
import graci.interfaces.mrsf.mrsf_gradient as mrsf_gradient


class Mrsfgradient:
    """Analytic nuclear gradients of the states of an Mrsftddft object,
    computed with the density-fitted integrals of the MRSF library and
    PySCF's derivative integrals and XC machinery"""

    def __init__(self):
        # user defined quantities
        # label of the $mrsftddft section whose states are differentiated
        self.mrsf_label    = None
        # label of a previous $mrsfgradient section (another geometry)
        # whose Z-vectors start the Z-vector iterations (warm start)
        self.guess_label   = None
        # adiabatic state indices (1-based in the input; None = all)
        self.states        = None
        # Z-vector solver: residual norm threshold (absolute 2-norm; the
        # gradient error is about 0.1-0.2 times, the coupling error about
        # 0.3-0.5 times this value), max. iterations; 'native' (default:
        # the library's shared-subspace minimal-residual solver), 'pcg'
        # (Python conjugate gradients, reference) or 'dense' (debugging,
        # small systems only); preconditioner of the native solver: 'auto'
        # (kernel-free inner steps with an XC kernel, else the diagonal),
        # 'inner' or 'diag', and the number of inner steps
        self.zvec_tol      = 1e-6
        self.zvec_iter     = 100
        self.zvec_solver   = 'native'
        self.zvec_precond  = 'auto'
        self.zvec_inner    = 2
        # include the quadrature-grid response in the reference (ROKS)
        # XC gradient (PySCF grid_response)
        self.grid_response = False
        # derivative couplings <I|d/dR J> between states of the section:
        # nested list of 1-based adiabatic state pairs, e.g.
        # [[1 2] [2 3]], one pair: [[1 2]] (None = none)
        self.nac_pairs     = None
        # memory budget (GB) for the cached AO values on the XC grid and
        # the auxiliary blocks of the derivative-integral contraction
        self.mem_budget    = 2.0
        self.verbose       = True
        self.label         = 'default'

        # computed quantities
        # the Mrsftddft object
        self.ci            = None
        # gradients (nstates, natm, 3), Hartree/Bohr, in the PySCF frame
        self.grad          = None
        # gradient of the ROKS reference (natm, 3)
        self.grad_ref      = None
        # total energies of the differentiated states
        self.energies      = None
        # symmetry labels of the states
        self.state_syms    = None
        # Z-vector iterations and residual norms
        self.zvec_niter    = None
        self.zvec_resid    = None
        # largest canonical multiplier of the frozen-core rotations per state
        self.zeta_max      = None
        # extended method: canonical multiplier of the O1-O2 rotation per
        # state and the SOMO gap (Hartree) it is divided by
        self.zeta_12       = None
        self.somo_gap      = None
        # extended method with the Pi parent: pair-block eigenvalues, parent
        # weights, pair gap, |z_w| and Q_12 per state
        self.pi_diag       = None
        # timings of the stages
        self.times         = None
        # 0-based adiabatic indices of the differentiated states
        self.state_list    = None
        # converged Z-vectors (lzdim, nstates): start vectors of a
        # warm-started gradient section at the next geometry
        self.zvec          = None
        # derivative couplings: list of dicts (pair, etf, d, dE, ...)
        self.nac           = None

    def copy(self):
        """create of deepcopy of self"""
        new = Mrsfgradient()
        var_dict = {key: value for key, value in self.__dict__.items()
                    if not key.startswith('__') and not callable(key)}
        for key, value in var_dict.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))
        return new

    @staticmethod
    def parse_pairs(value, nst_avail):
        """nac_pairs as parsed from the input ([[I J] [K L] ...], a list
        of integer arrays) -> list of 0-based (I, J) tuples. Each pair must
        hold two different existing states (1-based in the input)."""
        usage = ('$mrsfgradient: nac_pairs must be a nested list of pairs of '
                 'different existing states, e.g. nac_pairs = [[1 2] [2 3]] '
                 '(one pair: [1 2] or [[1 2]])')
        if isinstance(value, np.ndarray) and value.ndim == 1 and value.size == 2:
            value = [value]                 # a single pair written as [1 2]
        if not isinstance(value, list) or len(value) == 0:
            sys.exit(usage)
        pairs = []
        for pair in value:
            pair = np.asarray(pair).reshape(-1)
            if pair.size != 2 or not np.issubdtype(pair.dtype, np.integer) or pair[0] == pair[1]:
                sys.exit(usage)
            I, J = int(pair[0]) - 1, int(pair[1]) - 1
            if min(I, J) < 0 or max(I, J) >= nst_avail:
                sys.exit(usage)
            pairs.append((I, J))
        return pairs

    #
    @timing.timed
    def run(self, ci, guess=None):
        """compute the gradients of the requested states of the
        Mrsftddft object ci; guess: a previous Mrsfgradient object whose
        Z-vectors (aligned to the present MOs) start the Z-vector
        iterations"""

        if type(ci).__name__ != 'Mrsftddft':
            sys.exit('$mrsfgradient: mrsf_label must refer to an $mrsftddft section')
        if str(self.zvec_solver).lower() not in ('native', 'pcg', 'dense'):
            sys.exit('$mrsfgradient: zvec_solver must be native, pcg or dense')
        if str(self.zvec_precond).lower() not in ('auto', 'inner', 'diag'):
            sys.exit('$mrsfgradient: zvec_precond must be auto, inner or diag')
        self.ci = ci

        output.print_mrsfgradient_header(self.label)

        # derivative couplings (0-based pairs)
        nst_avail = ci.n_states()
        pairs = None
        if self.nac_pairs is not None:
            pairs = self.parse_pairs(self.nac_pairs, nst_avail)

        # states (0-based after parsing); with only nac_pairs given, the
        # states of the pairs (their gradients come with the couplings)
        if self.states is None and pairs is not None:
            states = sorted({s for pair in pairs for s in pair})
        elif self.states is None:
            states = list(range(nst_avail))
        else:
            states = [int(s) for s in np.atleast_1d(self.states)]
            for st in states:
                if st < 0 or st >= nst_avail:
                    sys.exit('$mrsfgradient: state %d does not exist (the %s section has %d states)'
                             % (st + 1, ci.label, nst_avail))

        # Z-vector warm start from a previous gradient section (in memory,
        # or read from the checkpoint when guess_label names one that is
        # not in memory)
        zguess = None
        if guess is None and self.guess_label is not None:
            import graci.io.chkpt as chkpt
            try:
                guess = chkpt.read('Mrsfgradient.' + str(self.guess_label))
            except Exception:
                guess = None
        if guess is not None:
            # with couplings: one start vector per state vector of the pass
            vec_states = states if pairs is None else mrsf_gradient.nac_vectors(states, pairs)[0]
            zguess = mrsf_gradient.align_zvec(guess, ci, vec_states)
            if zguess is None:
                output.print_message('  Z-vector warm start: no usable Z-vectors in '
                                     + str(getattr(guess, 'label', self.guess_label)) + '; cold start')

        if getattr(ci, 'extended', False):
            import graci.interfaces.mrsf.mrsf_egradient as mrsf_egradient
            driver = mrsf_egradient.EMRSFGradientDriver(ci, self)
        else:
            driver = mrsf_gradient.GradientDriver(ci, self)
        if pairs is None:
            res = driver.run(states, zguess=zguess)
        else:
            # gradients and couplings in one response pass
            res, self.nac = driver.gradients_and_nacs(states, pairs, zguess=zguess)
        driver.finalise()

        self.grad        = res['grad']
        self.grad_ref    = res['grad_ref']
        self.energies    = np.array([ci.energy(st) for st in states])
        self.state_syms  = [ci.scf.mol.irreplbl[ci.state_sym(st)[0]] for st in states]
        self.zvec_niter  = np.array(res['zvec_niter'])
        self.zvec_resid  = np.array(res['zvec_resid'])
        self.zeta_max    = np.array(res.get('zeta_max', []))
        self.zeta_12     = np.array(res.get('zeta_12', []))
        self.somo_gap    = res.get('somo_gap', None)
        self.pi_diag     = res.get('pi', None)
        self.times       = res['times']
        self.state_list  = np.array(states)
        self.zvec        = np.array(res['zvec'])

        if self.verbose:
            output.print_mrsfgradient_results(self, states, res)
            output.print_mrsfgradient_nacs(self)

        return True
