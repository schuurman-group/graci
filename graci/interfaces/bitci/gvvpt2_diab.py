"""
Module for the calculation of QDPT diabatic states within the 
DFT/MRCI(2) framework
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
def diabpot(ci_method0, ci_method):
    """
    QDPT diabatisation within the DFT/MRCI(2) framework
    Here, ci_method0 corresponds to the previous geometry, R_n-1,
    and ci_method to the current geometry R_n
    """

    # nirr is given by the length of the nstates vector in ci obj
    nirr = ci_method.n_irrep()
    
    # initialise the list of diabatic potential matrices, one per irrep
    diabpots = [None for i in range(nirr)]
    
    # initialise the list to hold the bitci diabatic configuration
    # scratch file numbers
    diab_ciunit  = 0
    diab_ciunits = [0 for i in range(nirr)]

    # initialise the list to hold the bitci diabatic state vector
    # scratch file numbers
    diab_confunit  = 0
    diab_confunits = [0 for i in range(nirr)]

    # initialise the list to hold the bitci spin-coupling averaged
    # Hii scratch file numbers
    diab_aviiunit  = 0
    diab_aviiunits = [0 for i in range(nirr)]
    
    # initialise the list to hold the numbers of diabatic confs
    diab_nconf  = 0
    diab_nconfs = [0 for i in range(nirr)]
    
    # GVVPT2 regularizer index
    ireg = ci_method.allowed_regularizer.index(ci_method.regularizer)+1

    # Regularisation factor
    regfac = ci_method.regfac
    
    # the bitci mrci wfn object
    mrci_wfn = ci_method.bitci_mrci()

    # bitci MRCI configuration scratch file numbers
    ci_confunits = np.array(mrci_wfn.conf_units['adiabatic'], dtype=int)

    # bitci ref space eigenvector scratch file numbers
    ref_ciunits = np.array(ci_method.ref_wfn.ci_units['adiabatic'],
                           dtype=int)
    
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
    if len(ci_method.icvs) > 0:
        delete_core = False
    else:
        delete_core = True

    # Loop over irreps
    for irrep in ci_method.irreps_nonzero():

        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep)
        
        # Number of extra roots
        nextra = ci_method.nextra['pt2'][irrep]

        # A-vector scratch file number
        Aunit = ci_method.Aunits[irrep]
        
        # R0 determinant bit strings and eigenvectors
        # if there are no diabatic states for the previous geometry,
        # then assume that it is the starting geometry
        if ci_method0.det_strings['diabatic'] is None:
            key = 'adiabatic'
        else:
            key = 'diabatic'
        
        n_int0 = ci_method0.det_strings[key][irrep].shape[0]
        n_det0 = ci_method0.det_strings[key][irrep].shape[2]
        n_vec0 = ci_method0.vec_det[key][irrep].shape[1]
        dets0  = np.reshape(ci_method0.det_strings[key][irrep],
                            (n_int0*2*n_det0), order='F')
        vec0   = np.reshape(ci_method0.vec_det[key][irrep],
                            (n_det0*n_vec0), order='F')

        # Output diabatic potential matrix
        diabpot = np.zeros((nroots*nroots), dtype=float)
            
        args = (irrep, nroots, nextra, ireg, regfac,
                n_int0, n_det0, n_vec0, dets0, vec0,
                nmo0, smat, ncore, icore, delete_core,
                ci_confunits, ref_ciunits, Aunit, diabpot,
                diab_ciunit, diab_confunit, diab_nconf,
                diab_aviiunit)

        diabpot, diab_ciunit, diab_confunit, diab_nconf, \
            diab_aviiunit = libs.lib_func('gvvpt2_diab', args)

        # Bitci diabatic state vector scratch file number
        diab_ciunits[irrep] = diab_ciunit

        # Bitci diabatic configuration scratch file number
        diab_confunits[irrep] = diab_confunit

        # Number of diabatic confs
        diab_nconfs[irrep] = diab_nconf

        # Bitci diabatic spin-coupling averaged Hii scratch
        # file numbers
        diab_aviiunits[irrep] = diab_aviiunit
        
        # Save the diabatic potential
        diabpots[irrep] = np.reshape(diabpot, (nroots, nroots), order='F')

    # Retrieve the DFT/MRCI(2) diabatic state vector scratch file names
    diab_cinames = ['' for i in range(nirr)]
    name = ''
    for irrep in ci_method.irreps_nonzero():
        args = (diab_ciunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        diab_cinames[irrep] = name

    # Retrieve the DFT/MRCI(2) diabatic configuration scratch file names
    diab_confnames = ['' for i in range(nirr)]
    name = ''
    for irrep in ci_method.irreps_nonzero():
        args = (diab_confunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        diab_confnames[irrep] = name

    # Retrieve the spin-coupling averaged on-diagonal Hamiltonian
    # matrix element scratch file names
    diab_aviiname = ['' for i in range(nirr)]
    for irrep in ci_method.irreps_nonzero():
        args = (diab_aviiunits[irrep], name)
        name = libs.lib_func('retrieve_filename', args)
        diab_aviiname[irrep] = name
        
    return diabpots, np.array(diab_ciunits, dtype=int), diab_cinames,\
        np.array(diab_confunits, dtype=int), diab_confnames, diab_nconfs, \
        np.array(diab_aviiunits), diab_aviiname
