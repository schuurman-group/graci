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
        # adiabatic state indices (1-based in the input; None = all)
        self.states        = None
        # Z-vector solver: residual norm threshold, max. iterations,
        # 'pcg' (default) or 'dense' (debugging, small systems only)
        self.zvec_tol      = 1e-8
        self.zvec_iter     = 100
        self.zvec_solver   = 'pcg'
        # include the quadrature-grid response in the reference (ROKS)
        # XC gradient (PySCF grid_response)
        self.grid_response = False
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
        # timings of the stages
        self.times         = None
        # 0-based adiabatic indices of the differentiated states
        self.state_list    = None

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

    #
    @timing.timed
    def run(self, ci):
        """compute the gradients of the requested states of the
        Mrsftddft object ci"""

        if type(ci).__name__ != 'Mrsftddft':
            sys.exit('$mrsfgradient: mrsf_label must refer to an $mrsftddft section')
        self.ci = ci

        output.print_mrsfgradient_header(self.label)

        # states (0-based after parsing)
        nst_avail = ci.n_states()
        if self.states is None:
            states = list(range(nst_avail))
        else:
            states = [int(s) for s in np.atleast_1d(self.states)]
            for st in states:
                if st < 0 or st >= nst_avail:
                    sys.exit('$mrsfgradient: state %d does not exist (the %s section has %d states)'
                             % (st + 1, ci.label, nst_avail))

        if getattr(ci, 'extended', False):
            import graci.interfaces.mrsf.mrsf_egradient as mrsf_egradient
            driver = mrsf_egradient.EMRSFGradientDriver(ci, self)
        else:
            driver = mrsf_gradient.GradientDriver(ci, self)
        res = driver.run(states)
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
        self.times       = res['times']
        self.state_list  = np.array(states)

        if self.verbose:
            output.print_mrsfgradient_results(self, states, res)

        return True
