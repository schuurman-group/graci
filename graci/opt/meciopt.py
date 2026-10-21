"""
Minimum-energy conical intersection (MECI) optimisation, $meciopt section.

The section names a method section (method_label; at present an
$mrsftddft section) and a pair of its adiabatic states (states = [I J]).
The step-0 run of the method section, at the input geometry, provides
the first energies, gradients and coupling (the driver requests them
through prepare before the CI sections run). From there the constrained
optimiser of graci.opt.meci_step takes Cartesian steps; each new geometry
is evaluated by the method adapter of graci.opt.method_adapter on copies
of the method objects, warm-started from the previous geometry. The
step-0 objects are left unchanged, so other sections that use them see
the input geometry.
"""

import sys
import copy as copy
import numpy as np
from pyscf.lib import param

import graci.utils.timing as timing
import graci.utils.constants as constants
import graci.core.params as params
import graci.io.output as output
import graci.opt.meci_step as meci_step
import graci.opt.model_hessian as model_hessian
import graci.opt.method_adapter as method_adapter


class Meciopt:
    """Class constructor for MECI optimisation objects"""
    def __init__(self):
        # user defined quantities
        self.label          = 'meci'
        # label of the method section ($mrsftddft); None = the only one
        self.method_label   = None
        # the pair of adiabatic states [I J] (1-based, input order)
        self.states         = None
        # maximum number of optimisation steps
        self.max_iter       = 50
        # initial, smallest and largest trust radius (Bohr, norm of the
        # Cartesian step)
        self.trust_radius   = 0.3
        self.trust_min      = 0.01
        self.trust_max      = 0.5
        # convergence: energy gap (Eh), seam gradient (Eh/Bohr, largest
        # and rms Cartesian component), step (Bohr) and change of the
        # average energy (Eh), which may replace the step criteria
        self.conv_gap       = 1.e-5
        self.conv_grad_max  = 4.5e-4
        self.conv_grad_rms  = 3.e-4
        self.conv_step_max  = 1.8e-3
        self.conv_step_rms  = 1.2e-3
        self.conv_energy    = 1.e-6
        # starting Hessian: 'lindh' (model Hessian) or 'unit'
        self.hess_guess     = 'lindh'
        # Hessian update: 'bfgs' or 'none', and the number of steps used
        self.hess_update    = 'bfgs'
        self.hess_window    = 5
        # Davidson and Z-vector warm starts from the previous geometry
        self.warm_start     = True
        # full SCF and method output at every step
        self.step_output    = False
        # trajectory file (xyz, Angstrom); None = meciopt_<label>.xyz
        self.traj_file      = None
        self.verbose        = True

        # computed quantities
        # 'converged', 'max_iter', 'rejected' (too many rejected steps)
        # or 'failed' (evaluations failed)
        self.status         = None
        self.converged      = False
        self.niter          = 0
        # history of the accepted geometries (n, natm, 3) in Bohr, the
        # energies of I and J (n, 2), the gap, the average energy, the
        # multipliers (n, 2), trust radius, ratio of the merit function,
        # seam gradient and step norms, wall time per step (SCF +
        # integrals + method, and the whole evaluation) and resident memory
        # after each step (MB)
        self.geoms          = None
        self.energies       = None
        self.gaps           = None
        self.ebar           = None
        self.lam            = None
        self.trust          = None
        self.rho            = None
        self.grad_max       = None
        self.grad_rms       = None
        self.step_max       = None
        self.step_rms       = None
        self.step_wall      = None
        self.step_eval      = None
        self.step_rss       = None
        # Davidson iterations (all irreps), Z-vector iterations (both
        # states and the pair) and SCF cycles per step
        self.step_dav_iter  = None
        self.step_zvec_iter = None
        self.step_scf_cycles = None
        # final geometry in the units of the $molecule section, symbols
        self.crds_final     = None
        self.asym           = None
        self.units          = None
        # gap gradient, coupling vector and average gradient at the final
        # geometry (natm, 3), Eh/Bohr, and the intersection parameters
        self.g_final        = None
        self.h_final        = None
        self.gbar_final     = None
        self.ci_char        = None
        # wall time of the optimiser itself (s)
        self.opt_wall       = 0.
        # the method object at the final geometry
        self.final_method   = None

    def copy(self):
        """create a deepcopy of self"""
        new = Meciopt()
        for key, value in self.__dict__.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))
        return new

    def pair(self):
        """0-based (I, J), I < J"""
        st = np.atleast_1d(self.states)
        return int(min(st)) - 1, int(max(st)) - 1

    def prepare(self, method, scf, mol):
        """
        checks before the SCF and method sections run, and the requests
        of the step-0 run: the gradients of I and J and their coupling
        """
        if self.states is None or np.atleast_1d(self.states).size != 2:
            sys.exit('$meciopt: states must give the pair of states, e.g. states = [1 2]')
        st = np.atleast_1d(self.states)
        if not np.issubdtype(st.dtype, np.integer) or st[0] == st[1] or min(st) < 1:
            sys.exit('$meciopt: states must be two different 1-based states, e.g. states = [1 2]')
        I, J = self.pair()

        if method_adapter.adapter_class(method) is None:
            sys.exit('$meciopt: method_label must name an $mrsftddft section '
                     '(the only method with analytic gradients and couplings)')
        if mol.use_sym or mol.sym_grp is not None:
            sys.exit('$meciopt: MECI optimisation runs without point-group symmetry; '
                     'set use_sym = False (and no sym_grp) in the $molecule section')
        if len(mol.asym) < 3:
            sys.exit('$meciopt: at least three atoms are needed')
        if str(method.precision).lower() != 'double':
            sys.exit('$meciopt: the method section must use precision = double')
        nst = int(np.sum(method.nstates))
        if J >= nst:
            sys.exit('$meciopt: state '+str(J + 1)+' does not exist (the '
                     +str(method.label)+' section has '+str(nst)+' states)')
        if str(self.hess_guess).lower() not in ('lindh', 'unit'):
            sys.exit('$meciopt: hess_guess must be lindh or unit')
        if str(self.hess_update).lower() not in ('bfgs', 'none'):
            sys.exit('$meciopt: hess_update must be bfgs or none')

        # the step-0 run computes the pair's gradients and coupling
        method_adapter.adapter_class(method).request(method, I, J)

        # orbital noise of a loose SCF limits the gradients
        if scf.conv_tol_grad is None:
            scf.conv_tol_grad = 1.e-6
            output.print_message('$meciopt: $scf conv_tol_grad set to 1e-6 (the '
                                 'default leaves ~1e-4 Eh/Bohr of noise in the gradients)')
        elif scf.conv_tol_grad > 1.e-5:
            output.print_message('$meciopt: $scf conv_tol_grad = '+str(scf.conv_tol_grad)+
                                 ' may limit the accuracy of the gradients (1e-6 advised)')
        if nst < J + 2:
            output.print_message('$meciopt: no state above '+str(J + 1)+' is computed; '
                                 'a state crossing the pair from above would go unnoticed')
        return

    @timing.timed
    def run(self, method):
        """optimise the MECI of the pair starting from the step-0 run of
        method (an object of the method section at the input geometry)"""
        I, J = self.pair()
        mol0 = method.scf.mol
        Z = mol0.pymol().atom_charges()
        self.asym = list(mol0.asym)
        self.units = str(mol0.units)

        adapter = method_adapter.adapter_class(method)(
            method, I, J, self.label, warm_start=self.warm_start,
            step_output=self.step_output)

        if str(self.hess_guess).lower() == 'lindh':
            hess = lambda x: model_hessian.lindh(x, Z)
        else:
            hess = lambda x: 0.5*np.eye(x.size)
        conv = {'gap': self.conv_gap, 'grad_max': self.conv_grad_max,
                'grad_rms': self.conv_grad_rms, 'step_max': self.conv_step_max,
                'step_rms': self.conv_step_rms, 'energy': self.conv_energy}
        solver = meci_step.MeciSolver(hess, conv, trust=self.trust_radius,
                                      trust_min=self.trust_min,
                                      trust_max=self.trust_max,
                                      hess_update=str(self.hess_update).lower(),
                                      hess_window=self.hess_window,
                                      max_iter=self.max_iter)

        traj = self.traj_file if self.traj_file is not None \
            else 'meciopt_'+str(self.label)+'.xyz'
        last = traj[:-4]+'_last.xyz' if traj.endswith('.xyz') else traj+'_last'
        open(traj, 'w').close()

        if self.verbose:
            output.print_meciopt_header(self.label)
            output.print_meciopt_setup(self, method, I, J)
            output.print_meciopt_iter_header()

        rows = []

        def callback(it, pt, stp, flags, rho, nrej):
            row = {'it': it, 'e': pt.e, 'ebar': pt.ebar, 'de': pt.de,
                   'grad_max': float(np.max(np.abs(pt.gseam))),
                   'grad_rms': float(np.linalg.norm(pt.gseam)/np.sqrt(pt.x.size)),
                   'step_max': float(np.max(np.abs(stp.dx))),
                   'step_rms': float(np.linalg.norm(stp.dx)/np.sqrt(pt.x.size)),
                   'trust': stp.trust, 'rho': rho, 'nrej': nrej, 'flags': flags,
                   'info': pt.info, 'lam': pt.lam.copy(), 'x': pt.x.copy()}
            rows.append(row)
            if self.verbose:
                output.print_meciopt_iter(row)
            comment = ('iteration {:d}  E_{:d} = {:.10f}  E_{:d} = {:.10f}  '
                       'gap = {:.3e} Eh'.format(it, I + 1, pt.e[0], J + 1,
                                                pt.e[1], pt.de))
            X = pt.x.reshape(-1, 3)*param.BOHR
            output.write_xyz(traj, self.asym, X, comment, 'a')
            output.write_xyz(last, self.asym, X, comment, 'w')

        res = solver.run(adapter.evaluate, adapter.point0(), callback=callback,
                         on_accept=adapter.accept, on_reject=adapter.reject)

        # results
        pt = res['point']
        self.status = res['status']
        self.converged = self.status == 'converged'
        self.niter = int(res['niter'])
        self.opt_wall = float(solver.wall)
        natm = len(self.asym)
        self.geoms = np.array([r['x'].reshape(natm, 3) for r in rows])
        self.energies = np.array([r['e'] for r in rows])
        self.gaps = np.array([r['de'] for r in rows])
        self.ebar = np.array([r['ebar'] for r in rows])
        self.lam = np.array([r['lam'] for r in rows])
        self.trust = np.array([r['trust'] for r in rows])
        self.rho = np.array([np.nan if r['rho'] is None else r['rho'] for r in rows])
        self.grad_max = np.array([r['grad_max'] for r in rows])
        self.grad_rms = np.array([r['grad_rms'] for r in rows])
        self.step_max = np.array([r['step_max'] for r in rows])
        self.step_rms = np.array([r['step_rms'] for r in rows])
        self.step_wall = np.array([sum(r['info'].get(k, 0.) for k in
                                       ('scf_wall', 'ao2mo_wall', 'method_wall'))
                                   for r in rows])
        self.step_eval = np.array([r['info'].get('eval_wall', 0.) for r in rows])
        self.step_rss = np.array([np.nan if r['info'].get('rss_mb') is None
                                  else r['info']['rss_mb'] for r in rows])
        self.step_dav_iter = np.array([sum(r['info'].get('davidson_iter', []))
                                       for r in rows], dtype=int)
        self.step_zvec_iter = np.array([sum(r['info'].get('zvec_iter', []))
                                        for r in rows], dtype=int)
        self.step_scf_cycles = np.array([-1 if r['info'].get('scf_cycles') is None
                                         else r['info']['scf_cycles'] for r in rows],
                                        dtype=int)
        X = pt.x.reshape(natm, 3)
        self.crds_final = X.copy() if self.units.lower().startswith(('b', 'au')) \
            else X*param.BOHR
        self.g_final = pt.g.reshape(natm, 3)
        self.h_final = pt.h.reshape(natm, 3)
        self.gbar_final = pt.gbar.reshape(natm, 3)
        ch = meci_step.characterise(pt)
        self.ci_char = {key: float(ch[key]) for key in
                        ('pitch', 'asymmetry', 'tilt', 'rel_tilt', 'heading', 'P', 'B')}
        self.final_method = adapter.current()

        if self.verbose:
            output.print_meciopt_summary(self, I, J, constants.au2ev)

        adapter.finalise()
        return

    def result_objs(self):
        """the Molecule, Scf and method objects at the final geometry (to
        be written to the checkpoint before this object)"""
        ci = self.final_method
        if ci is None:
            return []
        return [ci.scf.mol, ci.scf, ci]
