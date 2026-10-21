"""
Interface between the MECI optimiser ($meciopt) and the method sections:
an adapter returns the energies, nuclear gradients and the interstate
coupling vector of a pair of adiabatic states I < J at a new geometry.

Each evaluation works on copies of the objects of the last accepted
geometry (Molecule, Scf, method), labelled <label>_<meciopt label>, so the
step-0 objects of the input geometry are never changed. The previous
objects are the guesses of the next evaluation: the SCF density is
projected from them, and the method's own warm starts (Davidson vectors,
Z-vectors) use them.

The evaluator output is the tuple (x, (E_I, E_J), (g_I, g_J), h, info)
of meci_step.MeciSolver, with x the geometry in Bohr (3N,), the
gradients and h in Hartree/Bohr in the same frame, and h the coupling
scaled by the gap (finite at the intersection, sign arbitrary).
"""

import sys
import time
import numpy as np
from pyscf.lib import param

import graci.core.ao2mo as ao2mo
import graci.io.output as output


def rss_mb():
    """resident memory of the process (MB; None where /proc is absent)"""
    try:
        with open('/proc/self/status') as f:
            for line in f:
                if line.startswith('VmRSS:'):
                    return float(line.split()[1])/1024.
    except OSError:
        pass
    return None


class MrsftddftAdapter:
    """
    MRSF-TDDFT (Mrsftddft objects): h = <I|dH/dR|J> + (E_J - E_I) dU
    (scaled_h_etf, translation invariant), from the gradient pass of the
    object with nac_pairs = [[I+1 J+1]] and nac_divide_gap = False
    """

    @staticmethod
    def request(method, I, J):
        """
        add the gradients of I and J (0-based) and their coupling to the
        requests of the step-0 run of method (the union with the user's
        own grad_states and nac_pairs)
        """
        if method.grad_states is not None:
            states = set(int(s) for s in np.atleast_1d(method.grad_states))
            method.grad_states = np.array(sorted(states | {I, J}), dtype=int)

        pairs = method.nac_pairs
        if pairs is None:
            pairs = []
        elif isinstance(pairs, np.ndarray) and pairs.ndim == 1:
            pairs = [pairs]
        pairs = [[int(p) for p in np.asarray(pair).reshape(-1)] for pair in pairs]
        if not any(set(p) == {I + 1, J + 1} for p in pairs):
            pairs.append([I + 1, J + 1])
        method.nac_pairs = pairs
        return

    def __init__(self, method, I, J, label, warm_start=True,
                 step_output=False):
        self.step0 = method
        self.I = I
        self.J = J
        self.warm_start = warm_start
        self.step_output = step_output
        mol = method.scf.mol
        self.labels = {'mol': str(mol.label) + '_' + str(label),
                       'scf': str(method.scf.label) + '_' + str(label),
                       'ci': str(method.label) + '_' + str(label)}

        # template of the method object at the following geometries:
        # the step-0 settings, no results, exactly the pair requested
        tmpl = method.copy()
        for attr in ('scf', 'grad', 'grad_ref', 'grad_state_list', 'zvec',
                     'zvec_niter', 'zvec_resid', 'nac', 'guess_vecs',
                     'guess_ener', 'grad_times', 'grad_pi', 'energies',
                     'energies_sym', 'diag_info', 'fock_mo'):
            setattr(tmpl, attr, None)
        tmpl.amps = {'adiabatic': None}
        for attr in ('dmats', 'natocc', 'natorb_ao'):
            setattr(tmpl, attr, {'adiabatic': None, 'diabatic': None})
        tmpl.grad_states = np.array([I, J], dtype=int)
        tmpl.nac_pairs = [[I + 1, J + 1]]
        tmpl.nac_divide_gap = False
        tmpl.print_orbitals = False
        tmpl.ref_state = -1
        tmpl.guess_label = None
        tmpl.verbose = step_output
        tmpl.label = self.labels['ci']
        self.tmpl = tmpl

        # the objects of the last accepted geometry and of the last trial
        self.base = method
        self.last = None

    def point0(self):
        """the evaluator output at the input geometry, from the step-0 run"""
        return self.extract(self.step0, {'scf_wall': 0., 'ao2mo_wall': 0.,
                                         'method_wall': 0., 'eval_wall': 0.,
                                         'rss_mb': rss_mb()})

    def evaluate(self, x):
        """SCF, integrals and MRSF-TDDFT run at the geometry x (Bohr);
        None if the SCF or the states fail"""
        t_eval = time.time()
        base = self.base
        X = np.asarray(x, dtype=float).reshape(-1, 3)
        info = {}

        # molecule at the new geometry (units of the input)
        mol = base.scf.mol.copy()
        mol.label = self.labels['mol']
        units = str(mol.units).lower()
        mol.crds = X.copy() if units.startswith(('b', 'au')) else X*param.BOHR

        # SCF from the projected density of the previous geometry
        t0 = time.time()
        scf = base.scf.copy()
        scf.label = self.labels['scf']
        scf.mol_label = mol.label
        scf.verbose = self.step_output
        scf.print_orbitals = False
        scf.restart = False
        scf.guess_label = None
        energy = scf.run(mol, base.scf)
        info['scf_wall'] = time.time() - t0
        if energy is None:
            output.print_message('$meciopt: SCF not converged at the trial geometry')
            return None

        # MO integrals, as the driver prepares them for MRSF-TDDFT
        t0 = time.time()
        eri = ao2mo.Ao2mo()
        eri.emo_cut = self.tmpl.mo_cutoff
        eri.run(scf, self.tmpl.precision, load_bitci=False,
                keep_mem=float(self.tmpl.mem_budget))
        info['ao2mo_wall'] = time.time() - t0

        # states, gradients and coupling of the pair
        t0 = time.time()
        ci = self.tmpl.copy()
        ci.scf_label = scf.label
        ok = ci.run(scf, guess=base if self.warm_start else None, mo_ints=eri)
        info['method_wall'] = time.time() - t0
        if ok is None:
            output.print_message('$meciopt: MRSF-TDDFT failed at the trial geometry')
            return None
        if ci.n_states() <= self.J:
            output.print_message('$meciopt: state '+str(self.J + 1)+' was lost '
                                 '(the number of roots was reduced)')
            return None

        self.last = ci
        info['eval_wall'] = time.time() - t_eval
        info['rss_mb'] = rss_mb()
        return self.extract(ci, info)

    def extract(self, ci, info):
        """evaluator output from a method object with the pair's gradients
        and coupling"""
        I, J = self.I, self.J
        x = ci.scf.mol.pymol().atom_coords().ravel()
        slist = [int(s) for s in np.atleast_1d(ci.grad_state_list)]
        gI = np.asarray(ci.grad[slist.index(I)]).ravel()
        gJ = np.asarray(ci.grad[slist.index(J)]).ravel()
        h = None
        for res in ci.nac:
            if set(int(p) for p in res['pair']) == {I, J}:
                h = np.asarray(res['scaled_h_etf']).ravel()
                zconv = bool(res.get('zvec_converged', True))
                zpair = int(res.get('zvec_niter') or 0)
        if h is None:
            sys.exit('$meciopt: the coupling of states '+str(I + 1)+' and '
                     +str(J + 1)+' was not computed')

        # diagnostics: convergence, iterations, neighbouring gaps
        info = dict(info)
        info['davidson_conv'] = all(d is None or d[1] == 1 for d in ci.diag_info)
        info['davidson_iter'] = [int(d[0]) for d in ci.diag_info if d is not None]
        zres = np.atleast_1d(ci.zvec_resid)
        info['zvec_conv'] = zconv and bool(np.all(zres <= ci.zvec_tol))
        info['zvec_iter'] = [int(n) for n in np.atleast_1d(ci.zvec_niter)] + [zpair]
        scf_info = getattr(ci.scf, 'scf_info', None)
        info['scf_cycles'] = None if not scf_info else scf_info.get('cycles')
        nst = ci.n_states()
        info['gap_below'] = ci.energy(I) - ci.energy(I - 1) if I > 0 else None
        info['gap_above'] = ci.energy(J + 1) - ci.energy(J) if J + 1 < nst else None
        info['grad_times'] = ci.grad_times
        return x, (ci.energy(I), ci.energy(J)), (gI, gJ), h, info

    def accept(self, pt):
        """the trial point is accepted: its objects become the base"""
        if self.last is not None:
            self.base = self.last
        self.last = None

    def reject(self, pt):
        """the trial point is rejected: its objects are dropped"""
        self.last = None

    def current(self):
        """the method object of the last accepted geometry"""
        return self.base

    def finalise(self):
        """release the library state of the last evaluation (the grid
        session and the loaded integrals)"""
        import graci.interfaces.mrsf.mrsf_init as mrsf_init
        import graci.interfaces.mrsf.mrsf_xc as mrsf_xc
        mrsf_xc.free_session()
        mrsf_init.finalize()
        return


# adapters by method class name
ADAPTERS = {'Mrsftddft': MrsftddftAdapter}


def adapter_class(method):
    """the adapter class of the method object (None if unsupported)"""
    return ADAPTERS.get(type(method).__name__)
