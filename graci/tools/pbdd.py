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
import copy as copy
import numpy as np

import graci.core.params as params
import graci.core.ao2mo as ao2mo
import graci.io.output as output
import graci.utils.timing as timing
import graci.utils.constants as constants
import graci.interfaces.overlap.overlap as overlap


#
def load_bddpy():
    """Import BDDpy, which is only needed for the parts of a Pbdd run that
       involve normal modes. A path_file run never touches it.

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

        # run mode: exactly one of these is set
        self.hessian_file        = None
        self.path_file           = None

        # diabatisation
        self.adt_type            = 'qdpt'
        self.norm_thresh         = 0.999
        self.det_thresh          = 1e-6
        # warn when min|S_ii| or det(S) along a chain falls below this
        self.overlap_warn        = 0.7

        # cut generation (hessian_file mode only)
        self.cut_scheme          = '1mode'
        self.stepsize            = 0.5
        self.npoints             = 10

        # the fit (hessian_file mode only)
        # [one-mode order, two-mode order], for the on- and off-diagonal
        # elements of the diabatic potential matrix respectively
        self.diag_order          = [6, 2]
        self.offdiag_order       = [6, 2]
        self.weight              = None
        self.reexpand            = None
        self.blocks              = None
        self.blockdiag_algorithm = 'svd'
        self.cartgrad            = False

        # overrides for quantities that are otherwise derived
        self.point_group         = None
        self.state_irreps        = None

        # output
        self.op_file             = None
        self.sop_file            = None
        self.opstates            = None
        self.h5_file             = None

        # text dumps of the harvested surfaces, one file per cut, for
        # plotting or for fitting functional forms of one's own
        self.print_potentials    = False
        self.print_couplings     = False

        #
        # Private variables -- should not be directly referenced
        #                      outside the class
        # ------------------------------------------------------------
        # allowed values for the multiple-choice keywords
        self.allowed_adt_type    = ['bdd', 'qdpt']
        self.allowed_cut_scheme  = ['1mode', '2mode', '2mode_ondiag']
        self.allowed_blockdiag   = ['svd', 'invsqrt']

        # the two-mode expansion order the cut schemes can constrain. A
        # diagonal 2-mode cut samples only the line Q_a = Q_b, where terms
        # of equal total degree are indistinguishable, so only the single
        # bilinear term is determined. See port_plan.md section 8.8.
        self.max_twomode_order   = 2

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
    def hessian_mode(self):
        """True if this is a Hessian-driven run, i.e. one that generates
           its own cuts and fits a vibronic Hamiltonian"""
        return self.hessian_file is not None

    #
    def reference_frame(self, ref_obj):
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
    def build_system(self, ref_obj):
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
        pt_grp, axes, origin = self.reference_frame(ref_obj)

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

            self.check_point_group(system, ref_obj)

        return system

    #
    def check_point_group(self, system, ref_obj):
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
    def read_path(self, ref_obj):
        """
        Read a path_file of Cartesian geometries.

        The reference geometry is always the head of the chain, so it is
        prepended: the path file holds the displaced points only.

        Returns:
          (npoints+1, natm, 3) geometries in Angstrom
        """

        # imported here rather than at module scope: parse imports this
        # module in order to validate a $pbdd section
        import graci.io.parse as parse

        mol = ref_obj.scf.mol.copy()
        mol.xyz_file = self.path_file

        coords = parse.parse_all_geoms(mol)

        # parse_all_geoms returns the file in its own units, and a Molecule
        # section declaring bohr means the path file is in bohr too
        if str(mol.units).lower().startswith('b'):
            coords = coords * constants.bohr2ang

        # atom_coords is always in bohr
        ref_crds = ref_obj.scf.mol.pymol().atom_coords() * constants.bohr2ang

        return np.concatenate([ref_crds.reshape(1, -1, 3), coords], axis=0)

    #
    def ensure_reference_wavefunctions(self, ref_obj):
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

        eri_mo         = ao2mo.Ao2mo()
        eri_mo.emo_cut = ref_obj.mo_cutoff
        eri_mo.run(ref_obj.scf, ref_obj.precision)

        ref_obj.run(ref_obj.scf, None, mo_ints=eri_mo)

        return

    #
    def c1_view(self, ref_obj):
        """
        A C1 view of a symmetric CI object, for overlaps against a chain.

        Nothing below overlap() knows about symmetry: overlap_c takes raw
        determinant bit strings, coefficient matrices and the MO overlap,
        and the irrep index is only a bookkeeping partition applied above
        it. So a symmetric object can be presented as a C1 one by
        concatenating its determinant lists and assembling its
        coefficients block-diagonally, padding with zeros off the blocks.
        Determinants of different irreps are disjoint by construction and
        n_int depends only on nmo, so the arrays are conformable and no
        deduplication is needed.

        Returns:
          (view, table) where view is a CI object presenting one irrep and
          table[i] is the (irrep, root) that flat state i came from
        """

        dets  = []
        vecs  = []
        table = []

        for irrep in ref_obj.irreps_nonzero():
            det = ref_obj.det_strings['adiabatic'][irrep]
            vec = ref_obj.vec_det['adiabatic'][irrep]
            dets.append(det)
            vecs.append(vec)
            table.extend([(int(irrep), root) for root in
                          range(vec.shape[1])])

        ndet    = sum(v.shape[0] for v in vecs)
        nstates = sum(v.shape[1] for v in vecs)

        block = np.zeros((ndet, nstates), dtype=float)
        idet  = 0
        ista  = 0
        for vec in vecs:
            block[idet:idet+vec.shape[0], ista:ista+vec.shape[1]] = vec
            idet += vec.shape[0]
            ista += vec.shape[1]

        view = ref_obj.copy()
        view.nstates                    = np.array([nstates], dtype=int)
        view.det_strings['adiabatic']   = [np.concatenate(dets, axis=2)]
        view.vec_det['adiabatic']       = [block]
        view.det_strings['diabatic']    = None
        view.vec_det['diabatic']        = None
        view.scf                        = ref_obj.scf

        return view, table

    #
    def map_state_symmetries(self, ref_obj, pt0_ci):
        """
        Assign each state of a chain the irrep of the reference state it
        overlaps with.

        Point 0 of a chain sits at the reference geometry, so the overlap
        between the two is close to a permutation matrix and the
        assignment is unambiguous wherever the states are distinguishable
        at all. Energy ordering would do the same job everywhere except
        for near-degenerate states, which is exactly where it fails.

        Returns:
          (irreps, smin) -- the irrep of each chain state, and the
          smallest overlap any assignment relied on
        """

        view, table = self.c1_view(ref_obj)

        nref = view.vec_det['adiabatic'][0].shape[1]
        npt0 = pt0_ci.vec_det['adiabatic'][0].shape[1]

        # MO overlaps between the two calculations, bra = the reference
        smo = pt0_ci.scf.mo_overlaps(view.scf)[:view.nmo, :pt0_ci.nmo]

        pairs = np.array([[i, j] for i in range(nref)
                          for j in range(npt0)], dtype=int)

        Sij = overlap.overlap(view, pt0_ci, smo, pairs, 0,
                              self.norm_thresh, self.det_thresh, False)
        Sij = np.reshape(Sij, (nref, npt0))

        assigned = np.argmax(np.abs(Sij), axis=0)
        irreps   = np.array([table[i][0] for i in assigned], dtype=int)
        smin     = np.min([abs(Sij[assigned[j], j]) for j in range(npt0)])

        return irreps, smin

    #
    def check_reference_energies(self, ref_obj, ci, stem):
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
    def release_wavefunctions(self, ci):
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

        return

    #
    def chain_overlap(self, prev_ci, cur_ci):
        """
        Wave function overlaps between consecutive points of a chain.

        The ADT is built as S^-1 (S S^T)^1/2, which returns a perfectly
        well-formed matrix from a badly conditioned overlap, so a broken
        chain does not otherwise announce itself. Both diabatisation
        schemes compute an overlap internally but neither returns it, so it
        is recomputed here from the adiabatic wave functions. That is one
        extra overlap per point, small beside the CI calculation itself.

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
    def walk_chain(self, ref_obj, coords, stem):
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

        prev_scf = None
        prev_ci  = None

        for ipt in range(npoints):

            lbl = stem + '_' + str(ipt)

            # geometry: symmetry is switched off along a chain, since a
            # displacement generally breaks it
            mol = ref_mol.copy()
            mol.label      = lbl
            mol.use_sym    = False
            mol.sym_grp    = None
            mol.multi_geom = False
            mol.units      = 'angstrom'
            mol.crds       = 1. * coords[ipt]
            mol.run()

            # scf, started from the previous point's orbitals
            scf = ref_scf.copy()
            scf.label       = lbl
            scf.mol_label   = lbl
            scf.restart     = False
            scf.guess_label = None if prev_scf is None else prev_scf.label

            if scf.run(mol, prev_scf) is None:
                sys.exit('\n ERROR: Pbdd, label = '+str(self.label)+
                         '\n scf did not converge at point '+str(ipt)+
                         ' of chain '+str(stem))

            # ci: the state space is flattened, since the chain runs in C1
            ci = ref_obj.copy()
            ci.label       = lbl
            ci.scf_label   = lbl
            ci.nstates     = np.array([nsta], dtype=int)
            ci.save_wf     = True
            ci.adt_type    = self.adt_type
            ci.norm_thresh = self.norm_thresh
            ci.det_thresh  = self.det_thresh
            ci.diabatic    = ipt > 0
            ci.guess_label = None if prev_ci is None else prev_ci.label

            # the AO -> MO transformation lives in the driver rather than
            # in ci.run, so it has to be done here
            eri_mo         = ao2mo.Ao2mo()
            eri_mo.emo_cut = ci.mo_cutoff
            eri_mo.run(scf, ci.precision)

            ci.run(scf, prev_ci, mo_ints=eri_mo)

            if ipt == 0:
                dev = self.check_reference_energies(ref_obj, ci, stem)
                output.print_pbdd_refcheck(stem, dev, self.ener_tol)

                if self.q0_ener is None:
                    self.q0_ener = np.asarray(ci.energies).reshape(-1).copy()

                # point 0 of every chain is the reference geometry in C1,
                # so the first one to be run supplies the mapping and no
                # extra calculation is needed
                if self.state_syms is None and ref_obj.n_irrep() > 1:
                    self.state_syms, smin = \
                        self.map_state_symmetries(ref_obj, ci)
                    output.print_pbdd_state_syms(
                        self.state_syms,
                        ref_obj.scf.mol.irreplbl, smin)

                # the ADT is the identity at the reference geometry, so the
                # diabatic potential is just the adiabatic energies
                np.fill_diagonal(diabpot[ipt],
                                 ci.energies_sym[0][:nsta])
            else:
                diabpot[ipt] = ci.diabpot[0]

                sdiag, ssvd, sdet = self.chain_overlap(prev_ci, ci)
                output.print_pbdd_step(stem, ipt, sdiag, ssvd, sdet,
                                       self.overlap_warn)

            chkpt.write(ci)

            # the previous point has now been used both to propagate the
            # diabatisation and to compute the diagnostics, so its
            # determinant expansions can go
            self.release_wavefunctions(prev_ci)

            prev_scf = scf
            prev_ci  = ci

        self.release_wavefunctions(prev_ci)

        return diabpot

    #
    def assemble_data(self, ref_obj):
        """
        Gather every chain into the single arrays BDDpy fits.

        The chains are independent, so the fit sees one flat set of points;
        which chain a point came from survives only in the origin list, for
        diagnostics.
        """

        _, _, _, _ = load_bddpy()
        import bddpy.data as bddpy_data

        stems  = sorted(self.diabpot.keys())
        qvec   = np.concatenate([self.qvec[s] for s in stems], axis=1)
        diab   = np.concatenate([self.diabpot[s] for s in stems], axis=2)
        origin = [s for s in stems
                  for _ in range(self.diabpot[s].shape[2])]

        return bddpy_data.DiabaticData(qvec=qvec, diabpot=diab,
                                       origin=origin)

    #
    def build_config(self, ref_obj, system):
        """
        Translate the Pbdd keywords into a BDDpy KdcConfig.

        The quantities kdc used to require of the user -- the reference
        energies, the point group, the state symmetries -- all come from
        the reference calculation instead.
        """

        _, _, _, _ = load_bddpy()
        import bddpy.config as bddpy_config

        state_irreps = self.state_irreps
        if state_irreps is None:
            state_irreps = self.state_syms

        return bddpy_config.KdcConfig(
            q0_energies         = self.q0_ener,
            order               = int(max(self.diag_order[0],
                                          self.offdiag_order[0])),
            diag_order          = self.diag_order,
            offdiag_order       = self.offdiag_order,
            point_group         = system.frame.group,
            state_irreps        = state_irreps,
            weight              = self.weight,
            reexpand            = self.reexpand,
            blocks              = self.blocks,
            blockdiag_algorithm = self.blockdiag_algorithm,
            cartgrad            = self.cartgrad,
            operator_states     = self.opstates,
            label_stem          = str(self.label),
            source              = 'graci: Pbdd, label = '+str(self.label))

    #
    def cut_groups(self):
        """
        Group the walked chains by the modes they displace.

        The two half-cuts through a mode are two halves of one curve, so
        they belong in one file: joining them is what makes the dumped
        surfaces plottable without further work.

        Returns:
          {(mode indices): [chain names]}, mode indices empty for a path
        """

        groups = {}

        for stem in sorted(self.diabpot.keys()):
            qvec = self.qvec.get(stem)
            if qvec is None:
                modes = ()
            else:
                modes = tuple(np.where(np.abs(qvec).max(axis=1) > 1e-12)[0])
            groups.setdefault(modes, []).append(stem)

        return groups

    #
    def write_surfaces(self, ref_obj, model=None):
        """
        Dump the harvested diabatic potentials and couplings as text.

        One file per cut, with the two half-cuts joined into a single
        curve running through the reference geometry. Energies are given
        in eV relative to the reference ground state, which is the
        convention the fitted model uses, so the ab initio and fitted
        columns are directly comparable.

        Arguments:
          ref_obj: the reference CI object, for the zero of energy
          model:   the fitted model, when there is one. Its values are
                   written alongside the ab initio ones.
        """

        _, _, bddpy_constants, _ = load_bddpy()

        eh2ev = bddpy_constants.active.EH2EV
        ezero = np.min(self.q0_ener)

        nsta  = None
        diag  = []
        offd  = []

        for modes, stems in self.cut_groups().items():

            xs = []
            ws = []
            qs = []

            for stem in stems:
                w    = self.diabpot[stem]
                qvec = self.qvec.get(stem)
                nsta = w.shape[0]

                if modes:
                    x = qvec[modes[0]]
                else:
                    # a path has no expansion coordinate, so index the
                    # geometries instead
                    x = np.arange(w.shape[2], dtype=float)

                xs.append(x)
                ws.append(np.transpose(w, (2, 0, 1)))
                if qvec is not None:
                    qs.append(qvec.T)

            x = np.concatenate(xs)
            w = np.concatenate(ws, axis=0)
            q = np.concatenate(qs, axis=0) if qs else None

            # the half-cuts share the reference geometry
            order = np.argsort(x, kind='stable')
            x, w  = x[order], w[order]
            if q is not None:
                q = q[order]
            keep  = np.concatenate([[True], np.abs(np.diff(x)) > 1e-10])
            x, w  = x[keep], w[keep]
            if q is not None:
                q = q[keep]

            w = (w - ezero * np.eye(nsta)[None, :, :]) * eh2ev

            fit = None
            if model is not None and q is not None:
                fit = model.diabatic(q)

            diag = [(i, i) for i in range(nsta)]
            offd = [(i, j) for i in range(nsta) for j in range(i+1, nsta)]

            stem = 'path' if not modes else \
                'q' + '_q'.join(str(m+1) for m in modes)

            if self.print_potentials:
                self.write_dat(str(self.label)+'_'+stem+'_pot.dat',
                               modes, x, w, fit, diag, 'V')

            if self.print_couplings and offd:
                self.write_dat(str(self.label)+'_'+stem+'_coup.dat',
                               modes, x, w, fit, offd, 'W')

        return

    #
    def write_dat(self, path, modes, x, w, fit, elements, symbol):
        """write one surface file"""

        coord = 'point' if not modes else \
            'Q' + '/Q'.join(str(m+1) for m in modes)

        head = ['%14s' % coord]
        for (i, j) in elements:
            head.append('%14s' % ('%s%d%d' % (symbol, i+1, j+1)))
            if fit is not None:
                head.append('%14s' % ('%s%d%d_fit' % (symbol, i+1, j+1)))

        with open(path, 'w') as datfile:
            datfile.write('# diabatic %s, eV, relative to the reference '
                          'ground state\n'
                          % ('couplings' if symbol == 'W' else 'potentials'))
            datfile.write('#'+''.join(head)[1:]+'\n')

            for ipt in range(len(x)):
                row = ['%14.8f' % x[ipt]]
                for (i, j) in elements:
                    row.append('%14.8f' % w[ipt, i, j])
                    if fit is not None:
                        row.append('%14.8f' % fit[ipt, i, j])
                datfile.write(''.join(row)+'\n')

        return

    #
    @timing.timed
    def fit(self, ref_obj, system):
        """
        Fit the vibronic coupling Hamiltonian and write the requested
        output files.

        Returns:
          the BDDpy VibronicModel
        """

        _, _, _, _ = load_bddpy()
        import bddpy.fitting as bddpy_fitting
        import bddpy.io.store as bddpy_store
        import bddpy.operators.mctdh as bddpy_mctdh
        import bddpy.operators.multiqd as bddpy_multiqd

        data   = self.assemble_data(ref_obj)
        config = self.build_config(ref_obj, system)
        config.validate()

        model = bddpy_fitting.fit(system, config, data)

        output.print_pbdd_fit(model, data)

        # naming a file selects its writer, so there is no separate
        # operator format keyword
        if self.op_file is not None:
            bddpy_mctdh.write(model, self.op_file)

        if self.sop_file is not None:
            bddpy_multiqd.write(model, self.sop_file, states=self.opstates)

        if self.h5_file is not None:
            bddpy_store.save(self.h5_file, model, data=data, system=system,
                             config=config)

        if self.print_potentials or self.print_couplings:
            self.write_surfaces(ref_obj, model)

        return model

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

        self.ensure_reference_wavefunctions(ref_obj)

        if self.hessian_mode():
            system = self.build_system(ref_obj)

            if self.verbose:
                output.print_pbdd_modes(system.frame.group,
                                        system.frequencies,
                                        system.mode_symmetry_labels(),
                                        system.program_mode_labels)

            bddpy_system, _, _, _ = load_bddpy()
            import bddpy.displace as bddpy_displace

            cuts = bddpy_displace.generate(system,
                                           scheme  = self.cut_scheme,
                                           step    = self.stepsize,
                                           npoints = self.npoints)

            # stored geometry-last, matching BDDpy's DiabaticData, so that
            # nothing is transposed on the way into the fit
            for cut in cuts:
                walked = self.walk_chain(ref_obj, cut.coords, cut.name)
                self.diabpot[cut.name] = walked.transpose(1, 2, 0)
                self.qvec[cut.name]    = 1. * cut.qvec

        else:
            coords = self.read_path(ref_obj)
            walked = self.walk_chain(ref_obj, coords, 'path')
            self.diabpot['path'] = walked.transpose(1, 2, 0)

        output.print_pbdd_summary(self.diabpot)

        # a path_file run has no normal modes, so there is nothing to
        # expand in and the diabatic potentials are the deliverable
        if self.hessian_mode():
            self.fit(ref_obj, system)
        elif self.print_potentials or self.print_couplings:
            self.write_surfaces(ref_obj)

        return
