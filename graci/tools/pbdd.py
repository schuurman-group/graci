"""
The Pbdd object and its associated functions.

Pbdd drives a propagative block diagonalisation diabatisation: it takes a
converged DFT/MRCI(2) calculation at a reference geometry, generates the
displaced geometries around it, walks each one as a chain of diabatisation
steps, and harvests the resulting diabatic potential matrices. Given a Hessian
it goes on to fit a vibronic coupling Hamiltonian via BDDpy; given a path of
Cartesian geometries it stops at the potentials.

Note that this is distinct from graci.interfaces.overlap.bdd, which performs a
single reference/displacement block diagonalisation step. Pbdd is the object
that propagates such steps along a chain.

The full design, including the reasoning behind the keyword set, is in
BDDpy/docs/port_plan.md section 8.
"""

import sys as sys
import os as os
import copy as copy
import shutil as shutil
import numpy as np

import graci.core.params as params
import graci.core.ao2mo as ao2mo
import graci.core.molecule as molecule
import graci.io.output as output
import graci.utils.timing as timing
import graci.utils.constants as constants
import graci.interfaces.overlap.overlap as overlap


#
def load_bddpy():
    """Import BDDpy, which is only needed for the parts of a Pbdd run that
       involve normal modes. A cut job never touches it.

    Returns:
      the bddpy system, hessian, constants and symmetry-detection modules
    """

    try:
        import bddpy.system as bddpy_system
        import bddpy.hessian as bddpy_hessian
        import bddpy.constants as bddpy_constants
        import bddpy.symmetry.detect as bddpy_detect
    except ImportError:
        sys.exit('\n ERROR: a Pbdd hessian_file run requires BDDpy, which '
                 'could not be imported.\n Install it with '
                 '"pip install -e <path to BDDpy>"')

    return bddpy_system, bddpy_hessian, bddpy_constants, bddpy_detect


#
def kabsch(mobile, target):
    """Return the proper rotation R for which mobile @ R best matches
       target, in the least squares sense. Both arguments are expected to
       be centred already.

    Arguments:
      mobile: (natm, 3) coordinates to be rotated
      target: (natm, 3) coordinates to rotate onto

    Returns:
      R: (3, 3) rotation matrix, acting on row vectors from the right
    """

    u, _, vt = np.linalg.svd(np.matmul(mobile.T, target))

    # force a proper rotation: an improper one would map the molecule onto
    # its mirror image, which is not a change of frame
    d = np.sign(np.linalg.det(np.matmul(u, vt)))

    return np.matmul(u, np.matmul(np.diag([1., 1., d]), vt))


#
def rebuild_system(system_data):
    """
    Rebuild a BDDpy System from the derived state a run recorded.

    Nothing is re-derived: the modes are the symmetry-adapted ones the run
    finished with and the irreps are read rather than reassigned, so a
    collect step cannot reach a different model space than the run did.

    Arguments:
      system_data: the dict recorded by Pbdd.store_system

    Returns:
      the BDDpy System
    """

    bddpy_system, _, bddpy_constants, bddpy_detect = load_bddpy()

    system = bddpy_system.System.from_arrays(
        labels      = list(system_data['labels']),
        numbers     = np.asarray(system_data['numbers']),
        coords      = np.asarray(system_data['coords']),
        masses      = np.asarray(system_data['masses']),
        frequencies = np.asarray(system_data['frequencies']),
        modes       = np.asarray(system_data['modes']),
        source      = str(system_data.get('source', '')))

    system.program_mode_labels = list(
        system_data.get('program_mode_labels') or [])

    irreps = system_data.get('mode_irreps')
    if irreps is not None:
        system.mode_irreps = np.asarray(irreps, dtype=int)

    axes   = system_data.get('frame_axes')
    origin = system_data.get('frame_origin')

    if axes is not None:
        system.frame = bddpy_detect.frame_from_axes(
            str(system_data['point_group']),
            np.asarray(axes),
            origin = None if origin is None else np.asarray(origin))

    return system


class IrrepView:
    """
    One irrep of a symmetric CI object, presented as a single-irrep one.

    Holds references, not copies. The determinant expansions run to
    hundreds of megabytes on anything real, and overlap() needs only four
    attributes of what it is given.
    """

    def __init__(self, ci, irrep):
        self.nmo         = ci.nmo
        self.scf         = ci.scf
        self.det_strings = {'adiabatic':
                            [ci.det_strings['adiabatic'][irrep]],
                            'diabatic': None}
        self.vec_det     = {'adiabatic':
                            [ci.vec_det['adiabatic'][irrep]],
                            'diabatic': None}


class Pbdd:
    """Class constructor for the Pbdd object."""
    def __init__(self):
        # the following is determined from user input -- these are
        # keywords in the params module
        # -----------------------------------------------------------
        self.label               = 'default'
        self.verbose             = True

        # the reference calculation: a Dftmrci2 section label
        self.reference           = None

        # generate: run the reference, derive the model space, and write
        #           the displaced geometries and the reference checkpoint
        # cut:      walk one chain of geometries and record what collect
        #           needs. The geometries come from the $molecule section
        #           and the chain is seeded from origin.
        #
        # There is no mode that does both. Fitting is gkdc's job, so no
        # GRaCI run can produce an operator file -- see port_plan.md 8.11.
        self.job_type            = 'generate'
        self.hessian_file        = None

        # a checkpoint from a generate run: supplies the derived model
        # space and the wave functions a chain propagates from
        self.origin      = None

        # where generate writes the displaced geometries
        self.geom_dir            = 'pbdd'

        # which reference this run propagated from, recorded so that a
        # collect step can refuse a set of cuts that did not share one
        self.ref_source          = None

        # the reference space configuration files of the C1 reference,
        # as raw bytes keyed by file name. ref_space.propagate is handed
        # file *names* rather than unit numbers, so a cut has to have
        # these on disk -- and they are the only thing in the chain that
        # is read from the reference's bitci scratch rather than from
        # memory. Carrying them in the checkpoint is what makes
        # origin self-contained: 30 kB on a def2-TZVPD stilbene,
        # against a scratch directory that every cut would otherwise
        # have to reach.
        self.ref_confs           = {}

        # diabatisation
        self.adt_type            = 'qdpt'
        self.norm_thresh         = 0.999
        self.det_thresh          = 1e-6

        # truncation for the state symmetry mapping (section 8.5). Kept
        # separate from the diabatisation's, and far looser, because the
        # mapping only has to answer which irrep a state belongs to -- a
        # discrete label, not a number anything is computed from.
        #
        # The cost is quadratic in the number of unique beta strings:
        # liboverlap allocates betafac(nbetaB, nbetaK). On a def2-TZVPD
        # C2 system with 23 states, 0.999 asks for ~290 GB and 0.95 for
        # 30 MB. The looser value is not a compromise, it is the only one
        # that runs.
        self.map_norm_thresh     = 0.95
        self.map_det_thresh      = 1e-4

        # least fraction of a state's character that must sit in one irrep
        # for its assignment to be called. Below this the mapping is
        # ambiguous and the run stops rather than guessing.
        self.map_thresh          = 0.5
        # cut generation (hessian_file mode only)
        self.cut_scheme          = '1mode'
        self.stepsize            = 0.5
        self.npoints             = 10


        #
        # Private variables -- should not be directly referenced
        #                      outside the class
        # ------------------------------------------------------------
        # Nothing here describes a vibronic Hamiltonian. A $pbdd section
        # harvests diabatic potentials and stops; expanding them into a
        # model is gkdc's job, and every keyword controlling that fit
        # lives on gkdc's command line. The split is deliberate: it makes
        # it impossible to go from an input file to an operator file
        # without having looked at the data. See port_plan.md 8.11.
        #
        # least |S_ii|, singular value or |det S| a chain step may show
        # before collect flags it
        self.overlap_warn        = 0.7

        # allowed values for the multiple-choice keywords
        self.allowed_job_type    = ['generate', 'cut']
        self.allowed_adt_type    = ['bdd', 'qdpt']
        self.allowed_cut_scheme  = ['1mode', '2mode', '2mode_ondiag']

        # harvested data, keyed by cut name. Values are numpy arrays in
        # BDDpy's geometry-last layout: qvec is (nmodes, ngeom) and diabpot
        # is (nsta, nsta, ngeom).
        #
        # N.B. no BDDpy object may be held as an attribute of this class.
        # chkpt.write walks __dict__ and exits on anything it cannot
        # JSON-encode, which would end the run after all the expensive work
        # had already been done.
        self.qvec                = {}
        self.diabpot             = {}

        # chain health, keyed by cut name: (ngeom, 3) of
        # [min|S_ii|, smallest singular value of S, |det S|], with point 0
        # NaN since it has no predecessor. Recorded, never judged here --
        # the thresholds are applied at collect time, so revising one costs
        # a single collect run rather than re-running every cut. See
        # port_plan.md section 8.11.
        self.diagnostics         = {}

        # per-cut deviation of point 0 from the reference states, Hartree
        self.refcheck            = {}

        # the derived state: everything about the model space that comes
        # from the reference calculation and the Hessian rather than from
        # the cuts. Persisted as plain arrays so that a later collect step
        # can rebuild the System without re-deriving it -- re-deriving
        # would repeat the frame and symmetry decisions of section 8.4,
        # and the mode symmetrisation can rotate modes within a degenerate
        # block, so a second derivation is not guaranteed to agree with
        # the first. See port_plan.md section 8.11.
        self.system_data         = {}

        # state symmetries, taken from the reference calculation
        self.state_syms          = None

        # energies at the reference geometry, taken from point 0 of the
        # first chain rather than from the reference calculation itself.
        # Both sit at the same geometry, but the reference runs with
        # symmetry and the chains in C1, and their reference spaces are
        # selected differently. The expansion origin has to be consistent
        # with the data being expanded, so it comes from the same C1
        # calculations the potentials do.
        self.q0_ener             = None

        # retain each geometry's bitci scratch instead of deleting it once
        # the chain has moved past. Only for debugging a chain that has
        # gone wrong: a production run cannot afford it (see
        # release_scratch)
        self.keep_scratch        = False

        # how closely the first geometry of a cut must match the reference
        # to be recognised as it. The geometries are written to seven
        # decimals, so the copy in the file differs from the checkpoint's
        # by up to 5e-8; a real displacement is never this small.
        self.geom_match_tol      = 1e-6

        # largest RMSD, in Bohr, tolerated between the Hessian geometry and
        # the reference geometry once they have been aligned
        self.geom_tol            = 1e-3

        # largest deviation, in Hartree, tolerated between the reference
        # calculation and point 0 of a chain. The two sit at the same
        # geometry but the reference runs with symmetry and the chain in
        # C1, so their reference spaces are selected differently and the
        # energies do not agree to convergence: on CH3/cc-pVDZ the
        # difference is 5e-5. This is a check that the chain starts from
        # the same states, not a convergence test, so it is loose.
        self.ener_tol            = 1e-3

    #
    def copy(self):
        """create of deepcopy of self"""
        new = Pbdd()

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
    def run(self, ci_objs):
        """
        Run the diabatisation and, in hessian_file mode, the fit.

        Arguments:
          ci_objs: single-element list holding the reference CI object,
                   assembled by Driver.get_postscf_objs
        """

        ref_obj = ci_objs[0]

        output.print_pbdd_header(self.label, ref_obj.label)

        if self.generate_mode():
            # only a generate job maps state symmetries against the
            # reference, so only it needs the reference's determinant
            # expansions. A cut takes its wave functions from the
            # origin, or from its own point 0, and re-running
            # the reference here would be a wasted CI calculation on
            # every one of the 2 x nmodes cuts.
            self._ensure_reference_wavefunctions(ref_obj)
            self._generate(ref_obj)
            return

        # cut mode: the geometries come from the $molecule section
        self._walk_from_molecule(ref_obj)

        output.print_pbdd_summary(self.diabpot)
        output.print_pbdd_diagnostics(self.diagnostics, self.refcheck,
                                      self.overlap_warn, self.ener_tol)

        return


    #
    def generate_mode(self):
        """True if this run derives the model space and writes the
           displaced geometries"""
        return str(self.job_type).lower() == 'generate'

    #
    def _reference_frame(self, ref_obj):
        """
        The symmetry frame of the reference calculation.

        PySCF does not reorient the molecule when symmetry is on: it leaves
        the geometry in the frame it was given and records the symmetry
        frame separately. State irreps are assigned with respect to that
        frame, so the normal modes must be classified in it too, or the
        Mulliken labels of the two will not correspond.

        Arguments:
          ref_obj: the reference CI object

        Returns:
          (point_group, axes, origin): the computational point group, the
          frame's basis vectors as rows, and its origin in Bohr. axes and
          origin are None for a C1 calculation.
        """

        mol   = ref_obj.scf.mol
        pymol = mol.pymol()

        if not pymol.symmetry:
            return mol.comp_sym, None, None

        # _symm_axes rows are the frame's basis vectors in the input frame
        # and _symm_orig its origin, matching BDDpy's SymmetryFrame
        # convention, coords -> (coords - origin) @ axes.T
        axes   = getattr(pymol, '_symm_axes', None)
        origin = getattr(pymol, '_symm_orig', None)

        if axes is None or origin is None:
            # the input geometry was already in the symmetry frame, so
            # pyscf recorded no transformation
            axes   = np.eye(3)
            origin = np.zeros(3)

        axes   = np.asarray(axes, dtype=float)
        origin = np.asarray(origin, dtype=float)

        # pyscf does not guarantee a right-handed frame, and BDDpy requires
        # one. Negating an axis is safe here: every operation of the eight
        # abelian groups is diagonal in the symmetry frame, and diagonal
        # matrices commute with the sign flip, so the operations -- and
        # hence every irrep character -- are unchanged. A sign flip changes
        # the direction of an axis, not which axis is which, so it cannot
        # exchange B1 and B2.
        if np.linalg.det(axes) < 0.:
            axes = axes.copy()
            axes[0] *= -1.

        return mol.comp_sym, axes, origin

    #
    def _build_system(self, ref_obj):
        """
        Read the Hessian and return a BDDpy System whose normal modes are
        expressed in the reference calculation's frame and classified in
        its symmetry frame.

        Arguments:
          ref_obj: the reference CI object

        Returns:
          the BDDpy System object

        N.B. the System is deliberately not stored on self: chkpt.write
        walks __dict__ and exits on anything it cannot JSON-encode.
        """

        bddpy_system, bddpy_hessian, bddpy_constants, _ = load_bddpy()

        mol     = ref_obj.scf.mol
        pymol   = mol.pymol()
        # the geometry GRaCI will actually compute at, in Bohr
        ref_crds = pymol.atom_coords()
        ref_nums = np.asarray(pymol.atom_charges(), dtype=int)

        # the Hessian, in the frequency program's own frame
        data = bddpy_hessian.read(self.hessian_file)

        if len(data.numbers) != len(ref_nums):
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n the hessian holds '+str(len(data.numbers))+
                     ' atoms, the reference geometry '+str(len(ref_nums)))

        if not np.array_equal(np.asarray(data.numbers, dtype=int), ref_nums):
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n the hessian and the reference geometry list their '
                     'atoms in different orders\n hessian:   '
                     +str([int(n) for n in data.numbers])+
                     '\n reference: '+str([int(n) for n in ref_nums]))

        # align the hessian geometry onto the reference geometry, and carry
        # the normal modes with it
        hess_crds = np.asarray(data.coords, dtype=float)
        hess_cen  = hess_crds.mean(axis=0)
        ref_cen   = ref_crds.mean(axis=0)

        rot  = kabsch(hess_crds - hess_cen, ref_crds - ref_cen)
        rmsd = np.sqrt(np.mean(np.sum(
            (np.matmul(hess_crds - hess_cen, rot)
             - (ref_crds - ref_cen))**2, axis=1)))

        if rmsd > self.geom_tol:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n the hessian geometry cannot be aligned with the '
                     'reference geometry: rmsd = '+str(rmsd)+' Bohr, '
                     'tolerance = '+str(self.geom_tol)+
                     '\n the two are not the same structure. Note that a '
                     'geometry related to the reference by a reflection '
                     'rather than a rotation will also fail here.')

        # modes arrive as (ncoo, nmodes); rotate the Cartesian triples of
        # each mode by the same rotation
        natm   = len(ref_nums)
        nmodes = data.modes.shape[1]
        modes  = np.asarray(data.modes, dtype=float).reshape(natm, 3, nmodes)
        modes  = np.einsum('aim,ij->ajm', modes, rot).reshape(3*natm, nmodes)

        system = bddpy_system.System.from_arrays(
            labels      = data.labels,
            numbers     = ref_nums,
            coords      = ref_crds,
            masses      = data.masses,
            frequencies = data.frequencies,
            modes       = modes,
            source      = 'graci: '+str(self.hessian_file))

        # retained so that a modes.xyz written later can report what the
        # frequency program itself called each mode
        system.program_mode_labels = list(data.mode_labels)

        # classify the modes in the reference calculation's symmetry frame
        pt_grp, axes, origin = self._reference_frame(ref_obj)

        if axes is None:
            system.assign_symmetry(point_group='c1')
        else:
            # the symmetry analysis works in Angstrom, so use BDDpy's own
            # conversion rather than GRaCI's -- the origin has to be
            # consistent with System.coords_angstrom
            system.assign_symmetry(
                point_group = pt_grp,
                axes        = axes,
                origin      = origin * bddpy_constants.active.BOHR2ANG)

            self._check_point_group(system, ref_obj)

        self._store_system(system, ref_obj)

        return system

    #
    def _store_system(self, system, ref_obj):
        """
        Record the derived state so a later step can rebuild the System.

        The modes are taken *after* assign_symmetry, since symmetry
        adaptation may have rotated them within a degenerate block; the
        irreps are stored alongside rather than re-derived, so that
        rebuilding cannot reach a different answer than the run did.
        """

        frame = system.frame

        self.system_data = {
            'labels'      : [str(x) for x in system.labels],
            'numbers'     : np.asarray(system.numbers, dtype=int),
            'coords'      : np.asarray(system.coords, dtype=float),
            'masses'      : np.asarray(system.masses, dtype=float),
            'frequencies' : np.asarray(system.frequencies, dtype=float),
            'modes'       : np.asarray(system.modes, dtype=float),
            'mode_irreps' : (np.asarray(system.mode_irreps, dtype=int)
                             if system.mode_irreps is not None else None),
            'point_group' : str(system.point_group),
            'frame_axes'  : (np.asarray(frame.axes, dtype=float)
                             if frame is not None else None),
            'frame_origin': (np.asarray(frame.origin, dtype=float)
                             if frame is not None else None),
            'program_mode_labels' : [str(x) for x in
                                     (system.program_mode_labels or [])],
            'source'      : 'graci: '+str(self.hessian_file)}

        return

    #
    def _check_point_group(self, system, ref_obj):
        """
        Confirm that BDDpy's own view of the geometry agrees with PySCF's.

        The frame is pinned rather than detected, so this cannot catch an
        orientation disagreement -- it is there to catch a hessian computed
        for a different structure than the reference, which would otherwise
        produce normal modes for the wrong molecule.
        """

        _, _, _, bddpy_detect = load_bddpy()

        found = bddpy_detect.detect(system.coords_angstrom,
                                    system.numbers).group

        if found.lower() != system.frame.group.lower():
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n the reference calculation uses point group '
                     +str(system.frame.group)+', but the hessian geometry '
                     'is '+str(found))

        return

    #
    def _ensure_reference_wavefunctions(self, ref_obj):
        """
        Make sure the reference calculation kept its determinant
        expansions, re-running it if it did not.

        extract_wf only fires on save_wf or diabatic, so an ordinary
        reference run leaves vec_det as None and no overlap with it is
        possible. Postci objects run after every CI object, so the flag
        cannot be set retroactively -- the calculation has to be repeated.
        """

        if ref_obj.vec_det['adiabatic'] is not None:
            return

        output.print_message('Pbdd: the reference calculation did not '
                             'retain its wave functions, which the state '
                             'symmetry mapping needs. Re-running it with '
                             'save_wf = True.')

        ref_obj.save_wf = True

        # the first pass already chose an active space and wrote it onto
        # this object; re-deriving it would repeat a preliminary DFT/CIS
        # for no gain, and the C1 head then inherits from whichever pass
        # ran last rather than from the one whose energies were checked
        ref_obj.autoras = False

        eri_mo         = ao2mo.Ao2mo()
        eri_mo.emo_cut = ref_obj.mo_cutoff
        eri_mo.run(ref_obj.scf, ref_obj.precision)

        ref_obj.run(ref_obj.scf, None, mo_ints=eri_mo)

        return

    #
    def _map_state_symmetries(self, ref_obj, pt0_ci):
        """
        Assign each state of a chain the irrep of the reference states it
        overlaps with.

        Point 0 of a chain sits at the reference geometry, so the overlap
        between the two is close to a permutation matrix and the
        assignment is unambiguous wherever the states are distinguishable
        at all. Energy ordering would do the same job everywhere except
        for near-degenerate states, which is exactly where it fails.

        Done one reference irrep at a time. Nothing below overlap() knows
        about symmetry -- overlap_c takes raw determinant bit strings,
        coefficient matrices and the MO overlap -- so an irrep can be
        presented as a single-irrep object without touching its arrays.

        The earlier implementation padded every irrep into one block so a
        single call would do, which meant a dense (sum ndet, sum nstates)
        array that was half zeros by construction, and handed liboverlap
        both dimensions doubled. Since its working set goes as the square
        of the unique beta strings, that is the difference between running
        and not.

        Returns:
          (irreps, smin) -- the irrep of each chain state, and the
          smallest irrep weight any assignment relied on
        """

        npt0 = pt0_ci.vec_det['adiabatic'][0].shape[1]
        nirr = ref_obj.n_irrep()

        # MO overlaps between the two calculations, bra = the reference
        smo = pt0_ci.scf.mo_overlaps(ref_obj.scf)[:ref_obj.nmo,
                                                 :pt0_ci.nmo]

        # how much of each chain state sits in each reference irrep
        weight = np.zeros((nirr, npt0), dtype=float)

        for irrep in ref_obj.irreps_nonzero():

            view = IrrepView(ref_obj, irrep)
            nref = view.vec_det['adiabatic'][0].shape[1]

            pairs = np.array([[i, j] for i in range(nref)
                              for j in range(npt0)], dtype=int)

            Sij = overlap.overlap(view, pt0_ci, smo, pairs, 0,
                                  self.map_norm_thresh,
                                  self.map_det_thresh, False)
            Sij = np.reshape(Sij, (nref, npt0))

            weight[irrep] = np.sum(Sij**2, axis=0)

        irreps = np.argmax(weight, axis=0)
        best   = weight[irreps, np.arange(npt0)]

        # an assignment is only called when most of the state's character
        # sits in one irrep. Truncating this hard is safe precisely
        # because a threshold too loose to decide announces itself here
        # rather than quietly mislabelling a state.
        bad = np.where(best < self.map_thresh)[0]

        if bad.size:
            lines = []
            for j in bad:
                lines.append('   state %3d  weights %s'
                             % (j+1, ' '.join('%.3f' % w
                                              for w in weight[:, j])))
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n the state symmetry mapping is ambiguous for '
                     +str(bad.size)+' state(s): no irrep holds more than '
                     +str(self.map_thresh)+' of their character.\n'
                     +'\n'.join(lines)+
                     '\n Tighten map_norm_thresh (currently '
                     +str(self.map_norm_thresh)+') and try again.')

        return irreps, float(best.min())

    #
    def _check_reference_energies(self, ref_obj, ci, stem):
        """
        Confirm that point 0 of a chain reproduces the reference states.

        Point 0 sits at the reference geometry but runs in C1 with a flat
        state list, so this is the check that the flattening kept the same
        set of states -- that nothing was dropped, and that no extra state
        came in from an irrep the reference did not ask for. It is the
        precondition for mapping state symmetries onto the chain.

        Returns:
          the largest deviation, in Hartree
        """

        ref_ener = np.sort(np.asarray(ref_obj.energies).reshape(-1))
        pt0_ener = np.sort(np.asarray(ci.energies).reshape(-1))

        if ref_ener.size != pt0_ener.size:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n chain '+str(stem)+' point 0 holds '
                     +str(pt0_ener.size)+' states, the reference '
                     +str(ref_ener.size))

        dev = np.abs(ref_ener - pt0_ener).max()

        if dev > self.ener_tol:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n chain '+str(stem)+' point 0 does not reproduce '
                     'the reference states\n largest deviation = '
                     +str(dev)+' Hartree, tolerance = '
                     +str(self.ener_tol)+
                     '\n reference: '+str(ref_ener)+
                     '\n point 0  : '+str(pt0_ener))

        return dev

    #
    def _check_against_origin(self, ref_obj, stem):
        """
        Confirm a seeded cut's own q0 calculation agrees with the
        reference it propagates from.

        A cut input carries its own CI section, and the driver runs it at
        the first geometry of the cut file -- which is the reference
        geometry, since Molecule.read_xyz takes coords[0]. That
        calculation is otherwise unused: the chain propagates its
        reference space and wave functions from the origin, not
        from this section. What it is good for is a check, and it is the
        only one available -- nothing else compares a cut input to the
        generate job it belongs to, so an edited hamiltonian, functional,
        basis or state count would walk the chain with settings differing
        from the head it propagates from, and produce surfaces that look
        entirely reasonable.

        The reference space is not what is being checked. Whether the cut
        arrived at its space through autoras or through explicit
        ras1/ras3 is the user's business, and either way the space that
        propagates comes from the origin. This checks that the
        two calculations describe the same states.

        Returns:
          the largest deviation, in Hartree
        """

        cut_ener = np.sort(np.asarray(ref_obj.energies).reshape(-1))
        ref_ener = np.sort(np.asarray(self.q0_ener).reshape(-1))

        if ref_ener.size != cut_ener.size:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n this cut\'s reference calculation holds '
                     +str(cut_ener.size)+' states, the reference in\n '
                     +str(self.origin)+' holds '
                     +str(ref_ener.size)+
                     '\n The two were not set up for the same states.')

        dev = np.abs(ref_ener - cut_ener).max()

        if dev > self.ener_tol:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n this cut\'s reference calculation does not '
                     'reproduce the reference held in\n '
                     +str(self.origin)+
                     '\n largest deviation = '+str(dev)+
                     ' Hartree, tolerance = '+str(self.ener_tol)+
                     '\n reference   : '+str(ref_ener)+
                     '\n this cut    : '+str(cut_ener)+
                     '\n Both sit at the same geometry, so the CI '
                     'sections differ -- hamiltonian, functional, basis '
                     'or reference space.')

        return dev

    #
    def _release_wavefunctions(self, ci):
        """
        Drop a geometry's determinant expansions once nothing further in
        the chain needs them.

        These are the largest arrays a chain holds, and a chain of any
        length would otherwise keep every geometry's expansions alive at
        once. Only the object walking the chain knows when a geometry has
        been finished with: its successor must already have been run, and
        the chain diagnostics already computed from it.
        """

        if ci is None:
            return

        for key in ci.det_strings.keys():
            ci.det_strings[key] = None
            ci.vec_det[key]     = None

        self._release_scratch(ci)

        return

    #
    def _release_scratch(self, ci):
        """
        Delete a geometry's bitci scratch files and MO integrals.

        bitci gives every calculation its own scratch directory, named for
        the object label, and the AO -> MO integrals are written per scf
        label, so a chain of N geometries leaves N of each behind with
        nothing reclaiming them. The old multi-step workflow was bounded
        by one cut per process; a chain walked here is not, and on
        C2H4/aug-cc-pVTZ this runs to about 110 MB per geometry.

        Safe at the same moment the determinant expansions are dropped:
        the diabatisation and the overlaps both read the previous
        geometry's wave functions from memory rather than from disk, so
        nothing reopens these files once the successor has run.
        """

        if self.keep_scratch:
            return

        # take the directories from the recorded file names rather than
        # assuming where bitci puts them
        dirs = set()

        for wfn in [getattr(ci, 'mrci_wfn', None),
                    getattr(ci, 'ref_wfn', None)]:
            if wfn is None:
                continue
            for attr in ['conf_name', 'ci_name', 'avii_name']:
                names = getattr(wfn, attr, None)
                if not isinstance(names, dict):
                    continue
                for rep_names in names.values():
                    if rep_names is None:
                        continue
                    for name in rep_names:
                        if name:
                            dirs.add(os.path.dirname(str(name)))

        for path in dirs:
            if path and os.path.isdir(path):
                shutil.rmtree(path, ignore_errors=True)

        # the AO -> MO integrals are written per scf label and are equally
        # dead once the geometry is done with. Every geometry has its own
        # label, so these are never shared and never reused.
        label = str(getattr(ci, 'scf_label', '')).strip()

        if label:
            for name in ['1e_'+label+'.h5', '2e_eri_'+label+'.h5']:
                if os.path.isfile(name):
                    try:
                        os.remove(name)
                    except OSError:
                        pass

        return

    #
    def _run_point(self, ref_obj, coords, label, prev_scf=None,
                  prev_ci=None, diabatic=True, inherit_ras=False):
        """
        Run one geometry of a chain.

        Every point gets its own copies of the Molecule, Scf and CI
        objects: the chain is object-to-object in memory, so point n needs
        point n-1's CI *and* Scf objects alive and populated, and reusing
        one object would destroy the reference before it was used.

        Arguments:
          ref_obj:  the reference CI object, used as the template
          coords:   (natm, 3) geometry in Angstrom
          label:    label for this point's objects
          prev_scf: the previous point's Scf object, or None at the head
          prev_ci:  the previous point's CI object, or None at the head
          diabatic: whether to diabatise against prev_ci

        Returns:
          (mol, scf, ci)
        """

        nsta    = int(np.sum(ref_obj.nstates))
        ref_scf = ref_obj.scf
        ref_mol = ref_scf.mol

        # geometry: symmetry is switched off along a chain, since a
        # displacement generally breaks it
        mol = ref_mol.copy()
        mol.label      = label
        mol.use_sym    = False
        mol.sym_grp    = None
        mol.multi_geom = False
        mol.units      = 'angstrom'
        mol.crds       = 1. * np.asarray(coords)
        mol.run()

        # scf, started from the previous point's orbitals
        scf = ref_scf.copy()
        scf.label       = label
        scf.mol_label   = label
        scf.restart     = False
        scf.guess_label = None if prev_scf is None else prev_scf.label

        if scf.run(mol, prev_scf) is None:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n scf did not converge at point '+str(label))

        # ci: the state space is flattened, since the chain runs in C1
        ci = ref_obj.copy()
        ci.label       = label
        ci.scf_label   = label
        ci.nstates     = np.array([nsta], dtype=int)
        ci.save_wf     = True
        ci.adt_type    = self.adt_type
        ci.norm_thresh = self.norm_thresh
        ci.det_thresh  = self.det_thresh
        ci.diabatic    = diabatic and prev_ci is not None
        ci.guess_label = None if prev_ci is None else prev_ci.label

        # The active space is only ever chosen twice: at the symmetric
        # reference of a generate job, and at the head of a cut that has
        # no origin. Everywhere else a guess object is passed and
        # dftmrci/dftmrci2 build the reference space with
        # ref_space.propagate rather than ref_space.generate, so autoras
        # is not reached at all.
        #
        # The exception is the C1 reference of a generate job, which has
        # no predecessor to propagate from and must not choose for itself:
        # selection without symmetry is unreliable. On a def2-TZVPD C2
        # case the symmetric reference chose RAS1 = 37, 40-48 and RAS3 =
        # 49-65, while the C1 run at the same geometry chose RAS1 = 1-48,
        # every occupied orbital including core, with four scattered
        # virtuals for RAS3. That caller passes inherit_ras.
        #
        # ras1/ras3/nhole1/nelec3 already come across with the copy; what
        # is needed is to stop autoras overwriting them. The MO indices
        # carry over because PySCF orders MOs by energy with and without
        # symmetry.
        if inherit_ras:
            ci.autoras = False

            if self.verbose:
                output.print_pbdd_inherited_ras(label, ci.ras1, ci.ras2,
                                                ci.ras3, ci.nhole1,
                                                ci.nelec3)

        # the AO -> MO transformation lives in the driver rather than in
        # ci.run, so it has to be done here
        eri_mo         = ao2mo.Ao2mo()
        eri_mo.emo_cut = ci.mo_cutoff
        eri_mo.run(scf, ci.precision)

        ci.run(scf, prev_ci, mo_ints=eri_mo)

        return mol, scf, ci

    #
    def _chain_overlap(self, prev_ci, cur_ci):
        """
        Wave function overlaps between consecutive points of a chain.

        The ADT is built as S^-1 (S S^T)^1/2, which returns a perfectly
        well-formed matrix from a badly conditioned overlap, so a broken
        chain does not otherwise announce itself. Both diabatisation
        schemes compute an overlap internally but neither returns it, so it
        is recomputed here from the adiabatic wave functions. That is one
        extra overlap per point, small beside the CI calculation itself.

        The `bdd` scheme builds S on its way to the ADT and now hands it
        back, so for that adt_type this costs nothing. `qdpt` computes its
        overlaps inside the Fortran and does not return them, so there the
        matrix is rebuilt here.

        The health of a step is measured by the smallest singular value of
        S and by |det S|, not by the smallest diagonal element. Adiabatic
        states routinely exchange between neighbouring geometries, which
        sends min|S_ii| to zero while the state space is perfectly well
        described -- that exchange is exactly what the diabatisation
        exists to absorb. Both singular measures are invariant to such a
        permutation and fall only when a state genuinely enters or leaves
        the window.

        Returns:
          (min|S_ii|, smallest singular value of S, |det S|)
        """

        # reuse the diabatisation's own overlap where it was kept.
        # N.B. bdd returns S already transformed by the previous
        # geometry's ADT, so it is <diabatic_n-1|adiabatic_n> rather than
        # the raw adiabatic overlap. The ADT is orthogonal, so singular
        # values and |det S| -- the two measures anything is judged on --
        # are unchanged by that. min|S_ii| is basis dependent and so is not
        # strictly the same quantity as in the recomputed branch; it is
        # reported for information only and never tested against.
        if getattr(cur_ci, 'chain_smat', None) is not None:
            Sij = cur_ci.chain_smat[0]

        else:
            nstates = prev_ci.vec_det['adiabatic'][0].shape[1]
            pairs   = np.array([[i, j] for i in range(nstates)
                                for j in range(nstates)], dtype=int)

            Sij = overlap.overlap(prev_ci, cur_ci, cur_ci.smo, pairs, 0,
                                  self.norm_thresh, self.det_thresh, False)
            Sij = np.reshape(Sij, (nstates, nstates))

        return (np.min(np.abs(np.diag(Sij))),
                np.min(np.linalg.svd(Sij, compute_uv=False)),
                abs(np.linalg.det(Sij)))

    #
    @timing.timed
    def _walk_chain(self, ref_obj, coords, stem,
                   head_scf=None, head_ci=None):
        """
        Walk one chain of geometries, propagating the diabatisation.

        Point 0 is the reference geometry and is run without
        diabatisation; every later point diabatises against its
        predecessor. Each point gets its own copies of the Molecule, Scf
        and CI objects: the chain is object-to-object in memory, so point n
        needs point n-1's CI *and* Scf objects still populated, and reusing
        one object would destroy the reference before it was used.

        Arguments:
          ref_obj: the reference CI object, used as the template
          coords:  (npoints, natm, 3) geometries in Angstrom, point 0 being
                   the reference geometry
          stem:    label stem for this chain, e.g. 'q1r'

        Returns:
          (npoints, nsta, nsta) array of diabatic potential matrices,
          in Hartree
        """

        # imported here rather than at module scope: chkpt imports this
        # module in order to register the Pbdd class
        import graci.io.chkpt as chkpt

        npoints = coords.shape[0]
        nsta    = int(np.sum(ref_obj.nstates))

        ref_scf = ref_obj.scf
        ref_mol = ref_scf.mol

        diabpot = np.zeros((npoints, nsta, nsta), dtype=float)

        # point 0 has no predecessor, so its row stays NaN
        health  = np.full((npoints, 3), np.nan, dtype=float)

        # a cut is seeded from the reference calculation rather than
        # recomputing it, so that every chain starts from the same
        # wave functions and their phases cannot differ
        prev_scf = head_scf
        prev_ci  = head_ci
        seeded   = head_ci is not None

        for ipt in range(npoints):

            lbl = stem + '_' + str(ipt)

            mol, scf, ci = self._run_point(ref_obj, coords[ipt], lbl,
                                          prev_scf, prev_ci)

            if ipt == 0 and not seeded:
                dev = self._check_reference_energies(ref_obj, ci, stem)
                self.refcheck[stem] = dev
                output.print_pbdd_refcheck(stem, dev, self.ener_tol)

                if self.q0_ener is None:
                    self.q0_ener = np.asarray(ci.energies).reshape(-1).copy()

                # point 0 of every chain is the reference geometry in C1,
                # so the first one to be run supplies the mapping and no
                # extra calculation is needed
                if self.state_syms is None and ref_obj.n_irrep() > 1:
                    self.state_syms, smin = \
                        self._map_state_symmetries(ref_obj, ci)
                    output.print_pbdd_state_syms(
                        self.state_syms,
                        ref_obj.scf.mol.irreplbl, smin)

                # the ADT is the identity at the reference geometry, so the
                # diabatic potential is just the adiabatic energies
                np.fill_diagonal(diabpot[ipt],
                                 ci.energies_sym[0][:nsta])
            else:
                diabpot[ipt] = ci.diabpot[0]

                sdiag, ssvd, sdet = self._chain_overlap(prev_ci, ci)
                health[ipt] = [sdiag, ssvd, sdet]
                output.print_pbdd_step(stem, ipt, sdiag, ssvd, sdet)

            self._write_point(ci)

            # the previous point has now been used both to propagate the
            # diabatisation and to compute the diagnostics, so its
            # determinant expansions can go -- unless it is the seeded
            # head, which belongs to the generate job and is shared by
            # every other cut
            if not (seeded and ipt == 0):
                self._release_wavefunctions(prev_ci)

            prev_scf = scf
            prev_ci  = ci

        self._release_wavefunctions(prev_ci)

        self.diagnostics[stem] = health

        return diabpot

    #
    def _load_origin(self):
        """
        Read a generate run's checkpoint: the model space, and the wave
        functions a chain propagates from.

        Named explicitly rather than found at a well-known location, and
        a file that is named but cannot be used is an error rather than a
        quiet recompute. A typo or a file that had not finished staging
        would otherwise produce a job that ran fine, reported nothing, and
        used a different diabatic basis from its siblings -- and phases do
        not show up in energies, so nothing downstream could detect it.

        Returns:
          the C1 reference CI object
        """

        import graci.io.chkpt as chkpt

        path = self.origin

        if not os.path.isfile(path):
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n origin '+str(path)+' does not exist')

        groups = chkpt.contents(file_name=path) or []
        pbdd   = [g for g in groups if str(g).startswith('Pbdd.')]

        if not pbdd:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n origin '+str(path)+' holds no Pbdd '
                     'section: it is not a generate run\'s checkpoint')

        source = chkpt.read(pbdd[0], file_name=path,
                            build_subobj=False, make_mol=False)

        if not source.system_data:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n origin '+str(path)+' carries no normal '
                     'mode data')

        # the derived state is copied, not referenced, so that this cut's
        # own output records what it was computed against
        self.system_data = source.system_data
        self.q0_ener     = source.q0_ener
        self.state_syms  = source.state_syms
        self.ref_confs   = getattr(source, 'ref_confs', {})
        self.ref_source  = os.path.abspath(path)

        name = 'Dftmrci2.'+str(source.label)+'_q0'

        if name not in groups:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n origin '+str(path)+' holds no C1 '
                     'reference calculation ('+name+')')

        c1_ci = chkpt.read(name, file_name=path, build_subobj=True,
                           make_mol=True)

        self._restore_ref_confs(c1_ci)

        if c1_ci is None or c1_ci.vec_det['adiabatic'] is None:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n the C1 reference in '+str(path)+' did not '
                     'retain its wave functions, so no chain can be '
                     'propagated from it')

        return c1_ci

    #
    def _write_point(self, ci):
        """
        Write one chain point's checkpoint entry, without the bulk.

        A cut writes one entry per geometry and a suite is 2 x nmodes
        cuts, so these files have to stay small: gkdc reads exactly two
        things from a point, the diabatic potential matrix and the
        geometry, and nothing else in the object is ever read back.

        The determinant expansions, natural orbitals, density matrices,
        MO coefficients and MO overlaps are dropped for the write and put
        back afterwards -- vec_det and det_strings are still needed in
        memory, since the next point diabatises against them. On a
        def2-TZVPD stilbene those run to hundreds of MB per geometry,
        against a few kB of potential matrix.

        chkpt.write walks __dict__, so setting an attribute aside is
        enough to keep it out of the file.
        """

        import graci.io.chkpt as chkpt

        # chkpt.write recurses into sub-objects, so the point's Scf is
        # written too. It has to stay attached -- Molecule.<label> comes
        # through it and gkdc reads the geometry from there -- but its
        # arrays are as unwanted as the CI object's. orbs and rdm_ao are
        # both still needed in memory: the next point takes its MO
        # overlaps from the former and its SCF guess density from the
        # latter.
        bulky = [(ci, ['vec_det', 'det_strings', 'natorb_ao', 'natocc',
                       'dmats', 'mos', 'smo']),
                 (getattr(ci, 'scf', None), ['orbs', 'rdm_ao'])]

        stash = []
        for obj, attrs in bulky:
            if obj is None:
                continue
            for attr in attrs:
                if hasattr(obj, attr):
                    stash.append((obj, attr, getattr(obj, attr)))
                    setattr(obj, attr, None)

        try:
            chkpt.write(ci)
        finally:
            for obj, attr, val in stash:
                setattr(obj, attr, val)

        return

    #
    def _store_ref_confs(self, ci):
        """
        Copy the reference space configuration files into this object, so
        that the checkpoint carries everything a cut needs.

        ref_space.propagate is passed file names rather than unit numbers
        (ref_space.py, `confnames0 = ci_method0.ref_wfn.conf_name[rep]`),
        so the reference's configurations have to exist on disk wherever
        a cut runs. It is the only thing in a chain that comes off the
        reference's bitci scratch -- the wave functions, the MOs and the
        reference space dimensions all travel in the checkpoint already.

        Stored as raw bytes keyed by base name, so that a cut can lay
        them down in its own scratch directory under the names bitci
        expects.
        """

        self.ref_confs = {}

        names = getattr(getattr(ci, 'ref_wfn', None), 'conf_name', None)
        if not isinstance(names, dict):
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n the reference calculation recorded no reference '
                     'space configuration files,\n so no cut could '
                     'propagate from it')

        for rep, rep_names in names.items():
            for name in (rep_names or []):
                if not name:
                    continue
                if not os.path.isfile(name):
                    sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                             '\n the reference space configuration file '
                             +str(name)+'\n was not found; a cut cannot '
                             'propagate from this reference')
                with open(name, 'rb') as handle:
                    self.ref_confs[os.path.basename(name)] = \
                        np.frombuffer(handle.read(), dtype=np.uint8)

        if not self.ref_confs:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n no reference space configuration files were '
                     'stored')

        return

    #
    def _restore_ref_confs(self, ci):
        """
        Write the stored reference space configurations into this job's
        own scratch, and point the reference object at them.

        This is what lets a cut run anywhere: it needs no access to the
        directory the generate job ran in, and 2 x nmodes cuts do not
        contend on one shared scratch tree.
        """

        if not self.ref_confs:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n origin '+str(self.origin)+
                     ' carries no reference space\n configurations. It '
                     'was written by an older version of the generate '
                     'job;\n re-run it to produce a self-contained '
                     'reference.')

        root = os.path.abspath(os.path.join('bitscratch',
                                            str(ci.label)))
        os.makedirs(root, exist_ok=True)

        written = {}
        for base, data in self.ref_confs.items():
            path = os.path.join(root, str(base))
            with open(path, 'wb') as handle:
                handle.write(np.asarray(data, dtype=np.uint8).tobytes())
            written[str(base)] = path

        # rewrite the recorded names in place, preserving the per-irrep
        # ordering the reference recorded them in
        names = ci.ref_wfn.conf_name
        for rep, rep_names in names.items():
            if not rep_names:
                continue
            names[rep] = [written.get(os.path.basename(str(n)), n)
                          for n in rep_names]

        return

    #
    def _walk_from_molecule(self, ref_obj):
        """
        Walk the chain of geometries given in the $molecule section.

        A cut is identified by its geometry file and by nothing inside it:
        the label here is only for reporting, since collect recovers the
        normal coordinate from the geometry itself rather than trusting a
        name.

        With `origin` given -- the generated case -- the chain is
        seeded from that reference and the copy of it at the head of the
        geometry file is verified and skipped.

        Without it, the chain simply starts from the first geometry given,
        which is what someone wants who has a path of their own rather
        than a generated cut. There are no normal modes in that case and
        so nothing to fit; the diabatic potentials are the deliverable.
        """

        mol  = ref_obj.scf.mol
        path = mol.xyz_file

        if path is None:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n a cut job takes its geometries from the '
                     '$molecule section, which names no xyz_file')

        _, coords = molecule.read_xyz_file(path)

        if str(mol.units).lower().startswith('b'):
            coords = coords * constants.bohr2ang

        stem = os.path.basename(path)
        for suffix in ['.xyz']:
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]
        if stem.startswith('geom_'):
            stem = stem[len('geom_'):]

        c1_ci = None
        first = 0

        if self.origin is not None:

            c1_ci = self._load_origin()

            # A generated cut carries the reference as its first geometry,
            # which the origin already holds. Skipping it avoids
            # repeating the calculation and, more importantly, keeps every
            # chain propagating from the same wave functions.
            #
            # If it does not match, the geometry file and the
            # origin were not generated together, and the run
            # stops. This used to fall through to walking it as an
            # ordinary point, which silently computed something other than
            # what was asked for.
            reference = np.asarray(self.system_data['coords']) \
                * constants.bohr2ang

            if coords.shape[1:] != reference.shape:
                sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                         '\n '+str(path)+' has '+str(coords.shape[1])+
                         ' atoms, the reference '+str(reference.shape[0]))

            # the geometries are written to seven decimals, so the copy in
            # the file differs from the checkpoint's by up to 5e-8
            if np.allclose(coords[0], reference, atol=self.geom_match_tol):
                first = 1
            else:
                worst = np.abs(coords[0] - reference).max()
                sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                         '\n the first geometry of '+str(path)+' is not '
                         'the reference geometry held in\n '
                         +str(self.origin)+
                         ' (largest difference '+('%.2e' % worst)+
                         ' Angstrom).\n They were not generated together, '
                         'so this cut cannot propagate from that '
                         'reference.')

            # The geometry check above confirms the same nuclei; this one
            # confirms the same electronic structure. The driver has
            # already run this input's CI section at coords[0], so the
            # comparison is free -- and it is the only thing standing
            # between an edited cut input and a chain that propagates
            # from a head it does not match.
            dev = self._check_against_origin(ref_obj, stem)
            self.refcheck[stem] = dev
            output.print_pbdd_refcheck(stem, dev, self.ener_tol)

        walked = self._walk_chain(ref_obj, coords[first:], stem,
                                 head_scf=None if c1_ci is None
                                 else c1_ci.scf, head_ci=c1_ci)

        # stored geometry-last, matching BDDpy's DiabaticData. gkdc reads
        # the per-geometry groups rather than this, since pairing each
        # geometry with its own potential by label cannot go out of step
        # the way two parallel arrays can -- but it costs nothing to keep
        # and write_surfaces uses it.
        self.diabpot[stem] = walked.transpose(1, 2, 0)

        return

    #
    @timing.timed
    def _generate(self, ref_obj):
        """
        Derive the model space and write the displaced geometries.

        This is where everything the integration provides is worked out --
        the pinned frame, the mode symmetries, the state symmetries, the
        expansion origin. All of it comes from the reference calculation
        and the Hessian, and none of it from the cuts, which is what lets
        the cuts be farmed out without losing any of it. It is carried in
        this run's checkpoint, which becomes the cut jobs' origin.

        Two reference calculations are run: the symmetric one the user
        supplied, which fixes the frame and gives the state irreps, and a
        C1 one at the same geometry, which supplies the expansion origin,
        the wave functions every chain propagates from, and the states the
        symmetry mapping is made against.
        """

        # imported here rather than at module scope: chkpt imports this
        # module in order to register the Pbdd class
        import graci.io.chkpt as chkpt

        _, _, _, _ = load_bddpy()
        import bddpy.displace as bddpy_displace

        system = self._build_system(ref_obj)

        if self.verbose:
            output.print_pbdd_modes(system.frame.group,
                                    system.frequencies,
                                    system.mode_symmetry_labels(),
                                    system.program_mode_labels)

        # the C1 reference: the head every chain propagates from
        coords = system.coords_angstrom
        _, _, c1_ci = self._run_point(ref_obj, coords, str(self.label)+'_q0',
                                     diabatic=False, inherit_ras=True)

        dev = self._check_reference_energies(ref_obj, c1_ci, 'q0')
        self.refcheck['q0'] = dev
        output.print_pbdd_refcheck('q0', dev, self.ener_tol)

        self.q0_ener = np.asarray(c1_ci.energies).reshape(-1).copy()

        if ref_obj.n_irrep() > 1:
            self.state_syms, smin = self._map_state_symmetries(ref_obj,
                                                              c1_ci)
            output.print_pbdd_state_syms(self.state_syms,
                                         ref_obj.scf.mol.irreplbl, smin)

        # the C1 reference is written whole: a cut job reads its orbitals
        # and determinant expansions to seed its chain
        self._store_ref_confs(c1_ci)
        chkpt.write(c1_ci)

        # stepsize and npoints may be given per mode, which cannot be
        # length-checked at parse time since the mode count comes from
        # the Hessian
        try:
            cuts = bddpy_displace.generate(system,
                                           scheme  = self.cut_scheme,
                                           step    = self.stepsize,
                                           npoints = self.npoints)
        except ValueError as err:
            sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                     '\n '+str(err))

        os.makedirs(self.geom_dir, exist_ok=True)

        written = []
        for cut in cuts:
            path = os.path.join(self.geom_dir, self._cut_filename(cut))
            bddpy_displace.write_cut(cut, system, path,
                                     comment=self._cut_comments(cut))
            written.append((path, cut))

        output.print_pbdd_generated(written, self.geom_dir)

        return system

    #
    def _cut_filename(self, cut):
        """
        The xyz file name for a cut: the modes it displaces and which way.

        A cut is identified by its file, not by anything inside it, so the
        name has to carry the whole identity -- geom_q3_pos.xyz for the
        positive half-cut along mode 3, geom_q3_q5_neg.xyz for the
        negative half of the diagonal cut through modes 3 and 5. Modes are
        numbered from one, as they are everywhere the user sees them.
        """

        modes = '_'.join('q'+str(m+1) for m in cut.modes)
        way   = 'pos' if cut.direction == 'r' else 'neg'

        return 'geom_'+modes+'_'+way+'.xyz'

    #
    def _cut_comments(self, cut):
        """
        One comment line per geometry, tagging which displacement it is.

        Signed, so the file says which half-cut it belongs to as well as
        how far along: point 0 is the reference, and the negative half
        counts down. Nothing reads these -- collect recovers the normal
        coordinate from the geometry itself -- but tagging every artifact
        redundantly is cheap and makes a stray file identifiable on sight.
        """

        modes = ','.join('q'+str(m+1) for m in cut.modes)
        sign  = 1 if cut.direction == 'r' else -1

        return [modes+', displacement '+str(sign*point)
                for point in range(cut.npoints)]
