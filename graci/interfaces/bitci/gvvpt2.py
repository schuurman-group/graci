"""
Module for calculating the DFT/MRCI(2) energy and wave function corrections
via the diagonalisation of the GVVPT2 effective Hamiltonian
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
def diag_heff(ci_method):
    """
    Diagonalisation of the DFT/MRCI(2) effective Hamiltonian
    """

    # GVVPT2 regularizer index
    ireg = ci_method.allowed_regularizer.index(ci_method.regularizer)+1
    
    # Regularisation factor
    regfac = ci_method.regfac
    
    # nirr is given by the length of the nstates vector in ci obj
    nirr = ci_method.n_irrep()

    # the bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units['adiabatic'], dtype=int)

    # Bitci ref space eigenvector scratch file numbers
    ref_ciunits = np.array(ci_method.ref_wfn.ci_units['adiabatic'],
                           dtype=int)

    # Initialise the lists to hold the eigenvector and Q-space
    # scratch file numbers
    ciunit    = 0
    ciunits   = [0 for i in range(nirr)]
    qunit     = 0
    qunits    = [0 for i in range(nirr)]
    dspunit   = 0
    dspunits  = [0 for i in range(nirr)]
    aviiunit  = 0
    aviiunits = [0 for i in range(nirr)]
    
    # Loop over irreps
    for irrep in ci_method.irreps_nonzero():

        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep)

        # Number of extra roots
        nextra = ci_method.nextra['pt2'][irrep]
        
        args = (irrep, nroots, nextra, ireg, regfac,
                ci_confunits, ciunit, ref_ciunits,
                qunit, dspunit, aviiunit)

        ciunit, qunit, dspunit, aviiunit = libs.lib_func('gvvpt2', args)

        # Bitci eigenvector, Q-space, intruder state, and averaged
        # on-diagonal Hamiltonian matrix element scratch numbers
        ciunits[irrep]   = ciunit
        qunits[irrep]    = qunit
        dspunits[irrep]  = dspunit
        aviiunits[irrep] = aviiunit
    
    # Retrieve the DFT/MRCI(2) energies
    maxroots = max(ci_method.n_states_sym())
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in ci_method.irreps_nonzero():
            # Number of roots for the current irrep
            nroots = ci_method.n_states_sym(irrep)

            args = (ciunits[irrep], nroots, ener[irrep,:nroots])
            (ener[irrep,:nroots]) = \
                    libs.lib_func('retrieve_energies', args)

    # Retrieve the DFT/MRCI(2) eigenvector scratch file names
    ciname = ['' for i in range(nirr)]
    name   = ''
    for irrep in ci_method.irreps_nonzero():
        args = (ciunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        ciname[irrep] = name

    # Retrieve the spin-coupling averaged on-diagonal Hamiltonian
    # matrix element scratch file names
    aviiname = ['' for i in range(nirr)]
    for irrep in ci_method.irreps_nonzero():
        args = (aviiunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        aviiname[irrep] = name

    # Print the report of the DFT/MRCI(2) states 
    if ci_method.verbose:
        output.print_dftmrci2_states_header()
        
    ciunits = np.array(ciunits, dtype=int)
    aviiunits = np.array(aviiunits, dtype=int)
    nstates = ci_method.n_states_sym()
    args = (ci_confunits, ciunits, nstates)
    libs.lib_func('print_mrci_states', args)

    return ciunits, ciname, ener, qunits, dspunits, aviiunits, aviiname

@timing.timed
def diag_heff_follow(ci_method, ci_method0):
    """
    Diagonalisation of the DFT/MRCI(2) effective Hamiltonian
    with root following
    """

    # GVVPT2 regularizer index
    ireg = ci_method.allowed_regularizer.index(ci_method.regularizer)+1
    
    # Regularisation factor
    regfac = ci_method.regfac
    
    # nirr is given by the length of the nstates vector in ci obj
    nirr = ci_method.n_irrep()

    # the bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # Bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units['adiabatic'], dtype=int)

    # Bitci ref space eigenvector scratch file numbers
    ref_ciunits = np.array(ci_method.ref_wfn.ci_units['adiabatic'],
                           dtype=int)

    # Initialise the list to hold the eigenvector, Q-space info
    # DSP, and A-vector scratch file numbers
    ciunit    = 0
    ciunits   = [0 for i in range(nirr)]
    qunit     = 0
    qunits    = [0 for i in range(nirr)]
    dspunit   = 0
    dspunits  = [0 for i in range(nirr)]
    Aunit     = 0
    Aunits    = [0 for i in range(nirr)]
    aviiunit  = 0
    aviiunits = [0 for i in range(nirr)]
    
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
    for irrep in ci_method.irreps_nonzero():

        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep)

        # Number of extra roots
        nextra = ci_method.nextra['pt2'][irrep]

        # R0 determinant bit strings and eigenvectors
        n_int0 = ci_method0.det_strings['adiabatic'][irrep].shape[0]
        n_det0 = ci_method0.det_strings['adiabatic'][irrep].shape[2]
        n_vec0 = ci_method0.vec_det['adiabatic'][irrep].shape[1]
        dets0  = np.reshape(ci_method0.det_strings['adiabatic'][irrep],
                            (n_int0*2*n_det0), order='F')
        vec0   = np.reshape(ci_method0.vec_det['adiabatic'][irrep],
                            (n_det0*n_vec0), order='F')
        
        args = (irrep, nroots, nextra, ireg, regfac,
                n_int0, n_det0, n_vec0, dets0, vec0,
                nmo0, smat, ncore, icore, delete_core,
                ci_confunits, ciunit, ref_ciunits, qunit,
                dspunit, Aunit, aviiunit)

        ciunit, qunit, dspunit, Aunit, aviiunit = \
            libs.lib_func('gvvpt2_follow', args)
        
        
        # Bitci eigenvector, Q-space, intruder state, and A-vector
        # scratch numbers
        ciunits[irrep]  = ciunit
        qunits[irrep]   = qunit
        dspunits[irrep] = dspunit
        Aunits[irrep]   = Aunit
        aviiunits[irrep] = aviiunit
        
    # Retrieve the DFT/MRCI(2) energies
    maxroots = max(ci_method.n_states_sym())
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in ci_method.irreps_nonzero():
        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep)

        args = (ciunits[irrep], nroots, ener[irrep,:nroots])
        (ener[irrep,:nroots]) = \
            libs.lib_func('retrieve_energies', args)

    # Retrieve the DFT/MRCI(2) eigenvector scratch file names
    ciname = ['' for i in range(nirr)]
    name   = ''
    for irrep in ci_method.irreps_nonzero():
        args = (ciunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        ciname[irrep] = name

    # Retrieve the spin-coupling averaged on-diagonal Hamiltonian
    # matrix element scratch file names
    aviiname = ['' for i in range(nirr)]
    for irrep in ci_method.irreps_nonzero():
        args = (aviiunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        aviiname[irrep] = name

    # Print the report of the DFT/MRCI(2) states 
    if ci_method.verbose:
        output.print_dftmrci2_states_header()

    ciunits = np.array(ciunits, dtype=int)
    aviiunits = np.array(aviiunits, dtype=int)
    nstates = ci_method.n_states_sym()
    args = (ci_confunits, ciunits, nstates)
    libs.lib_func('print_mrci_states', args)
    
    return ciunits, ciname, ener, qunits, dspunits, Aunits, aviiunits, aviiname
