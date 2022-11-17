"""
Module for the calculation of QDPT Type-I and Type-II diabatic states
within the DFT/MRCI(2) framework
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.output as output
import graci.io.convert as convert
import graci.core.molecule as molecule

@timing.timed
def type_ii(ci_method0, ci_method):
    """
    Type-II QDPT diabatisation within the DFT/MRCI(2) framework
    Here, ci_method0 corresponds to the previous geometry, R_n-1,
    and ci_method to the current geometry R_n
    """

    # initialise the list of ADT matrices, one per irrep
    adt_matrices = []

    # GVVPT2 regularizer index
    ireg = ci_method.allowed_regularizer.index(ci_method.regularizer)+1

    # nirr is given by the length of the nstates vector in ci obj
    nirr = ci_method.n_irrep()

    # Regularisation factor
    regfac = ci_method.regfac
    
    # the bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units, dtype=int)

    # bitci ref space eigenvector scratch file numbers
    ref_ciunits = np.array(ci_method.ref_wfn.ci_units, dtype=int)    

    # MO overlaps
    nmo0 = ci_method0.nmo
    nmo  = ci_method.nmo
    smat = np.reshape(ci_method.smo, (nmo0 * nmo), order='F')

    # frozen core orbitals
    ncore_el = np.array([molecule.atom_ncore[n] for n in
                         [molecule.atom_name.index(m)
                          for m in ci_method0.scf.mol.asym]])
    ncore_mo = int(np.sum(ncore_el/2))

    # fill in the core MO indices
    # (bitX libraries have the MOs in energy order and
    # use Fortran indexing)
    icore = np.arange(0, ncore_mo, 1) + 1
    ncore = icore.size

    # deleted core orbital flag
    delete_core = True

    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep)

        # Number of extra roots
        nextra = ci_method.nextra['pt2'][irrep]

        # A-vector scratch file number
        Aunit = ci_method.Aunits[irrep]
        
        # R0 determinant bit strings and eigenvectors
        n_int0 = ci_method0.det_strings['adiabatic'][irrep].shape[0]
        n_det0 = ci_method0.det_strings['adiabatic'][irrep].shape[2]
        n_vec0 = ci_method0.vec_det['adiabatic'][irrep].shape[1]
        dets0  = np.reshape(ci_method0.det_strings['adiabatic'][irrep],
                            (n_int0*2*n_det0), order='F')
        vec0   = np.reshape(ci_method0.vec_det['adiabatic'][irrep],
                            (n_det0*n_vec0), order='F')

        # R0 ADT matrix
        if ci_method0.adt is None:
            adt0 = np.reshape(np.eye(n_vec0), (n_vec0**2), order='F')
        else:
            adt0 = np.reshape(ci_method0.adt[irrep], (n_vec0**2), order='F')

        args = (irrep, nroots, nextra, ireg, regfac,
                n_int0, n_det0, n_vec0, dets0, vec0,
                nmo0, smat, ncore, icore, delete_core,
                ci_confunits, ref_ciunits, Aunit, adt0)

        libs.lib_func('gvvpt2_diab', args)
        
    return
