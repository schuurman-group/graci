"""
Module for initialising/finalising the mrsf (MRSF-TDDFT) library
"""

import numpy as np
import graci.core.libs as libs
import graci.core.ao2mo as ao2mo


def stamp(ci_method):
    """identity of the MO integrals of a CI object (geometry + MOs)"""
    return ao2mo.mo_stamp(ci_method.scf.mol.crds, ci_method.mos)

#
def init(ci_method, fock_mo):
    """
    Initialise the mrsf library for the Mrsftddft object ci_method.
    fock_mo: (2, nmo, nmo) array of the alpha/beta spin Fock matrices
    of the reference in the MO basis of ci_method.mos
    """

    nmo   = ci_method.nmo
    nel   = ci_method.nel
    mosym = np.array(ci_method.mosym, dtype=np.int64)
    occ   = np.array(ci_method.occ_ref, dtype=float)
    moen  = np.array(ci_method.emo, dtype=float)
    fa    = np.reshape(np.asarray(fock_mo[0], dtype=float), (nmo*nmo), order='F')
    fb    = np.reshape(np.asarray(fock_mo[1], dtype=float), (nmo*nmo), order='F')
    chf   = float(ci_method.chf_eff)
    spc   = np.array(ci_method.spc_eff, dtype=float)
    escf  = float(ci_method.scf.energy)

    if ci_method.scf.mol.sym_indx <= 0:
        pgrp = 1
    else:
        pgrp = ci_method.scf.mol.sym_indx + 1

    args = (nmo, nel, int(ci_method.scf.mult), mosym, occ, moen, fa, fb,
            chf, spc, escf, pgrp, ci_method.label, ci_method.verbose)
    libs.lib_func('mrsf_initialise', args)

    return

#
def init_ints(ci_method, eri_file, mo_ints=None, mem_budget=None):
    """load the density-fitted MO integrals: from the packed tensor that
    the Ao2mo object mo_ints holds in memory when it belongs to the same
    geometry and MOs (double precision), otherwise from eri_file. The
    stamp (geometry + MOs) is part of the library's reuse key, so a file
    rewritten under the same name at another geometry is never mistaken
    for the loaded integrals. mem_budget (GB): budget for the exchange
    tensor decision and the work buffers (default: the section's)"""

    stmp = stamp(ci_method)
    budget = float(ci_method.mem_budget if mem_budget is None else mem_budget)
    eri_mem = getattr(mo_ints, 'eri_mo', None) if mo_ints is not None else None
    if eri_mem is not None and mo_ints.stamp == stmp \
            and ci_method.precision == 'double' \
            and eri_mem.shape[1] == ci_method.nmo*(ci_method.nmo + 1)//2:
        naux, npair = eri_mem.shape
        # the C-ordered (naux, npair) array is the F-contiguous (npair, naux) matrix
        args = (ci_method.precision, eri_file, ci_method.vv_storage,
                budget, stmp, int(naux), int(npair), eri_mem.T)
        libs.lib_func('mrsf_int_load_mem', args)
        if mo_ints.file_written:
            mo_ints.release()
    else:
        ao2mo.ensure_file(eri_file)
        args = ('df', ci_method.precision, eri_file, ci_method.vv_storage,
                budget, stmp)
        libs.lib_func('mrsf_int_initialise', args)

    return

#
def set_extended(flag):
    """select the extended (EMRSF) response space; must precede init()"""

    libs.lib_func('mrsf_set_extended', (bool(flag),))

    return

#
def set_frozen_core(nfc):
    """number of doubly occupied MOs (the lowest by orbital energy)
    excluded from the response space; must precede init()"""

    libs.lib_func('mrsf_set_frozen_core', (int(nfc),))

    return

#
def set_exchange(mode):
    """exchange term of the sigma vectors: 'auto' (explicit (ij|ab) tensor
    when it fits mem_budget, else the DF plane sweep), 'tensor', 'df';
    must precede init_ints()"""

    m = {'auto': 0, 'tensor': 1, 'df': 2}[str(mode).lower()]
    libs.lib_func('mrsf_set_exchange', (m,))

    return

#
def init_ext(ci_method, fone, fcpl, wpi, shift, coul, mode, eri_file):
    """extended-method data: fone (nmo, nmo) one-body operator of the CV
    block (published model: the closed-shell KS matrix of G, corrected to
    F' in the library; EMRSF-Pi: the invariant operator), fcpl (nmo, nmo,
    3) one-body operators of the Fock-type coupling terms per parent G, D,
    S, wpi (3) parent weights, shift (Pi pair energy; ignored for the
    published model), coul (scale of the singlet CV Coulomb term), mode
    (1 published, 2 Pi); coupling scale ci_method.ccp_eff and kernel flag
    ci_method.use_kernel; after init_ints() and, when the kernel is used,
    after the grid initialisation (mrsf_xc.XCGrid)"""

    nmo = int(ci_method.nmo)
    # the particle-core block is read from the file when the integrals were
    # loaded by a standard section before (the file may still be unwritten)
    ao2mo.ensure_file(eri_file)
    fone = np.asfortranarray(np.asarray(fone, dtype=float).reshape(nmo, nmo))
    fcpl = np.asfortranarray(np.asarray(fcpl, dtype=float).reshape(nmo, nmo, 3))
    wpi = np.asfortranarray(np.asarray(wpi, dtype=float).reshape(3))
    # F' correction of the published model: 1 = diagonal, 2 = covariant (Part XII)
    fprime = {'diagonal': 1, 'covariant': 2}[str(getattr(ci_method, 'fprime', 'covariant')).lower()]
    args = (nmo, fone, fcpl, wpi, float(shift), float(ci_method.ccp_eff),
            bool(ci_method.use_kernel), float(coul), int(mode), fprime, eri_file)
    libs.lib_func('mrsf_ext_initialise', args)

    return

#
def parent_tables(mult, X):
    """the library's pinned coupling tables of parent X (1 G, 2 D, 3 S) for
    the multiplicity mult: dict of F-ordered arrays Jco, Kco, Fco (2,2),
    Wco (2,3,2), Jov, Kov, Fov (2,2), Wov (2,3,2), A, B (2,2), Joo (3),
    Koo (3,4), Fblk (3) (0-based indices, see mrsf_extended.f90)"""
    shapes = [('Jco', (2, 2)), ('Kco', (2, 2)), ('Fco', (2, 2)), ('Wco', (2, 3, 2)),
              ('Jov', (2, 2)), ('Kov', (2, 2)), ('Fov', (2, 2)), ('Wov', (2, 3, 2)),
              ('A', (2, 2)), ('B', (2, 2)), ('Joo', (3,)), ('Koo', (3, 4)), ('Fblk', (3,))]
    arrs = [np.zeros(shp, dtype=np.float64, order='F') for _, shp in shapes]
    libs.lib_func('mrsf_ext_parent_tables', (int(mult), int(X)) + tuple(arrs))
    return {name: np.array(a) for (name, _), a in zip(shapes, arrs)}


def ext_dims():
    """returns (ncol, ncv): columns of the response vector and CV slots"""

    return libs.lib_func('mrsf_get_ext_dims', (0, 0))

#
def dims():
    """returns (nocca, nvirb, xdim, naux)"""

    return libs.lib_func('mrsf_get_dims', (0, 0, 0, 0))

#
def finalize():
    """finalise the mrsf library"""

    libs.lib_func('mrsf_finalise', ())

    return

#
def report_timings():
    """print the library timings of the current section and reset the
    counters (the loaded integrals are kept)"""

    libs.lib_func('mrsf_report_timings', ())

    return

#
def reset_timings():
    """reset the library timing counters without printing them"""

    libs.lib_func('mrsf_reset_timings', ())

    return
