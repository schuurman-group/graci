"""
Module for performing the reference space diagonalisation
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing
import graci.io.convert as convert
import graci.io.output as output
import graci.core.molecule as molecule

@timing.timed
def diag(ci_method):
    """
    Diagonalisation of the reference space Hamiltonian
    """
    
    # length of nstates vector is the number of irreps
    nirr    = ci_method.n_irrep()

    # bitci reference space wfn
    ref_wfn = ci_method.bitci_ref()

    # Print the section header
    if ci_method.verbose:
        output.print_refdiag_header()
    
    # Bitci reference configuration scratch file numbers
    confunits = np.array(ref_wfn.conf_units['adiabatic'], dtype=int)

    # Numbers of configurations
    nconf = np.array(ref_wfn.nconf['adiabatic'], dtype=int)

    # Bitci eigenvector scratch file numbers
    ciunit  = 0
    ciunits = []
    
    # Loop over irreps
    for irrep in range(nirr):
    
        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep)

        # Number of extra roots
        nextra = ci_method.nextra['max'][irrep]

        # Call to the bitci reference space diagonalisation routine
        args = (irrep, nroots+nextra, confunits, nconf, ciunit)
        (nroots, ciunit) = libs.lib_func('ref_diag_mrci',args)

        # Bitci eigenvector scratch number
        ciunits.append(ciunit)
    
        # If the number of reference space configurations for the
        # current irrep is less than the requested number of roots
        # then reset nstates accordingly
        if nroots < ci_method.n_states_sym(irrep):
            ci_method.nstates[irrep] = nroots
    
    # Retrieve the reference space energies
    maxroots = max(ci_method.n_states_sym())
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci_method.n_states_sym(irrep) > 0:
    
            # Number of roots for the current irrep
            nroots = ci_method.n_states_sym(irrep)

            # Retrieve the energies
            args = (ciunits[irrep], nroots, ener[irrep, :nroots],
                    np.array([n+1 for n in range(nroots)], dtype=int))
            
            (ener[irrep, :nroots]) = \
                    libs.lib_func('retrieve_some_energies', args)
            
    return ciunits, ener 

@timing.timed
def n_extra(ci_method):
    """
    Determination of the number of extra reference space
    eigenvectors needed
    """

    # Initialisation
    nextra = {'prune' : None,
              'guess' : None,
              'max'   : None}

    # length of nstates vector is the number of irreps
    nirr    = ci_method.n_irrep()
    
    # If pruning is off and ENPT2 guess vectors are not being used,
    # then nextra = 0
    if ci_method.prune == False and ci_method.diag_guess != 'enpt2':
        for key in nextra:
            nextra[key] = [0 for n in range(nirr)]
        return nextra
            
    # Number of extra ref space roots needed for pruning
    if ci_method.prune != 'off':
        nextra['prune'] = [ci_method.prune_extra for n in range(nirr)]
    else:
        nextra['prune'] = [0 for n in range(nirr)]
    
    # Number of extra ref space roots needed for guess vector generation
    if ci_method.diag_guess == 'enpt2':
        if len(ci_method.diag_blocksize) == 0:
            nextra['guess'] = [int(round(max(n+5, n*1.3)) - n) 
                                        for n in ci_method.nstates]
        else:
            nextra['guess'] = ci_method.diag_blocksize - ci_method.nstates
    else:
        nextra['guess'] = [0 for n in range(nirr)]

    # Number of extra ref space roots to be calculated
    nextra['max'] = []
    for n in range(nirr):
        nextra['max'].append(max([nextra[key][n] for key in nextra
                                  if key != 'max']))

    return nextra

@timing.timed
def diag_follow(ci_method, ci_method0):
    """
    Diagonalisation of the reference space Hamiltonian plus root
    following based on overlaps with the saved MRCI vectors of the
    ci_method0 object computed at a previous geometry R0
    """

     # Exit if the point groups for the two calculations are different
    if ci_method.scf.mol.comp_sym != ci_method0.scf.mol.comp_sym:
        sys.exit('\n Error in ref_space.propagate: non-equal point groups')
    
    # length of nstates vector is the number of irreps
    nirr    = ci_method.n_irrep()

    # bitci reference space wfn
    ref_wfn = ci_method.bitci_ref()

    # Print the section header
    if ci_method.verbose:
        output.print_refdiag_header()
    
    # Bitci reference configuration scratch file numbers
    confunits = np.array(ref_wfn.conf_units['adiabatic'], dtype=int)

    # Numbers of configurations
    nconf = np.array(ref_wfn.nconf['adiabatic'], dtype=int)

    # Bitci eigenvector scratch file numbers
    ciunit  = 0
    ciunits = []

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
        nextra = ci_method.nextra['max'][irrep]
        
        # R0 determinant bit strings and eigenvectors
        n_int0 = ci_method0.det_strings['adiabatic'][irrep].shape[0]
        n_det0 = ci_method0.det_strings['adiabatic'][irrep].shape[2]
        n_vec0 = ci_method0.vec_det['adiabatic'][irrep].shape[1]
        dets0  = np.reshape(ci_method0.det_strings['adiabatic'][irrep],
                            (n_int0*2*n_det0), order='F')
        vec0   = np.reshape(ci_method0.vec_det['adiabatic'][irrep],
                            (n_det0*n_vec0), order='F')
        
        # Call to the bitci root-following reference space
        # diagonalisation routine
        args = (irrep, nroots+nextra, confunits,
                n_int0, n_det0, n_vec0, dets0, vec0, nmo0, smat,
                ncore, icore, delete_core, nconf, ciunit)
        (nroots, ciunit) = libs.lib_func('ref_diag_mrci_follow',args)

        # Bitci eigenvector scratch number
        ciunits.append(ciunit)

        # If the number of reference space configurations for the
        # current irrep is less than the requested number of roots
        # then reset nstates accordingly
        if nroots < ci_method.n_states_sym(irrep):
            ci_method.nstates[irrep] = nroots
        
    # Retrieve the reference space energies
    maxroots = max(ci_method.n_states_sym())
    ener     = np.zeros((nirr, maxroots), dtype=float)
    for irrep in range(nirr):
        if ci_method.n_states_sym(irrep) > 0:
    
            # Number of roots for the current irrep
            nroots = ci_method.n_states_sym(irrep)

            # Retrieve the energies
            args = (ciunits[irrep], nroots, ener[irrep, :nroots],
                    np.array([n+1 for n in range(nroots)], dtype=int))
            
            (ener[irrep, :nroots]) = \
                    libs.lib_func('retrieve_some_energies', args)

    return ciunits, ener
