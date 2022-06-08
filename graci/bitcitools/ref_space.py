"""
Module for constructing the reference space configurations
"""

import sys
import ctypes as ctypes
import numpy as np
import graci.utils.timing as timing
import graci.core.libs as libs
import graci.io.output as output
import graci.io.convert as convert
from pyscf import gto

@timing.timed
def generate(ci_method):
    """generate the reference space object"""

    # Optional automated determination of the RAS MO
    # spaces via the analysis of the DFT/CIS eigenvectors
    if ci_method.autoras:
        autoras(ci_method)

    # number of irreps
    nirr = ci_method.n_irrep()

    # number of mos
    nmo = ci_method.scf.nmo

    # orbital occupations
    occ = ci_method.ref_occ

    # Number of orbitals in RAS1, RAS2, and RAS3
    m1 = ci_method.ras1.size
    m2 = ci_method.ras2.size
    m3 = ci_method.ras3.size

    # Number of electrons in RAS2
    n2 = 0
    for i in ci_method.ras2:
        n2 += int(occ[i-1])

    # Maximum number of holes in RAS1
    nh1 = ci_method.nhole1  
    
    # Maximum number of electrons in RAS3
    ne3 = ci_method.nelec3
    
    # Array of RAS1 orbital indices
    iras1 = np.zeros(nmo, dtype=int)
    k = -1
    for i in ci_method.ras1:
        k += 1
        iras1[k] = i
    
    # Array of RAS2 orbital indices    
    iras2 = np.zeros(nmo, dtype=int)
    k = -1
    for i in ci_method.ras2:
        k += 1
        iras2[k] = i
    
    # Array of RAS3 orbital indices
    iras3 = np.zeros(nmo, dtype=int)
    k = -1
    for i in ci_method.ras3:
        k +=1
        iras3[k] = i

    # CVS core MO flags
    cvsflag = np.zeros(nmo, dtype=int)
    for i in ci_method.icvs:
        cvsflag[i-1] = 1
    
    # Number of reference space configurations per irrep
    ref_nconf = np.zeros(nirr, dtype=int)

    # Reference space configurations scratch file number
    confunits = np.zeros(nirr, dtype=int)
    
    # Construct the reference space configurations
    args = (iras1, iras2, iras3, nh1, m1, n2, m2, ne3, m3,
            cvsflag, ref_nconf, confunits)
    (ref_nconf, confunits) = libs.lib_func('generate_ref_confs', args)

    # Retrieve the reference space configuration scratch file names
    confnames = []
    name      = ' '
    for i in range(nirr):
        args = (confunits[i], name)
        name = libs.lib_func('retrieve_filename', args)
        confnames.append(name)
    
    return ref_nconf, confunits, confnames

@timing.timed
def autoras(ci_method):
    """determination of the RAS subspaces via preliminary
    DFT/CIS calculations"""

    # Print the section header
    output.print_autoras_header()
   
    # number of irreps
    nirr = ci_method.n_irrep()

    # number of mos
    nmo = ci_method.scf.nmo

    # orbital occupations
    orb_occ = ci_method.ref_occ

    # orbital energies
    orb_ener = ci_method.scf.orb_ener

    # Number of extra roots
    n_extra = 2
    
    #
    # Diagonalisation of the DFT/CIS Hamiltonian
    #
    # Bitci eigenvector scratch file numbers
    dftcis_vec  = 0
    dftcis_unit = []

    # CVS core MO flags
    cvsflag = np.zeros(nmo, dtype=int)
    for i in ci_method.icvs:
        cvsflag[i-1] = 1

    # Aggressively loose integral screening
    loose = True
    
    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep) + n_extra

        # Call the the bitci DFT/CIS routine
        args = (irrep, nroots, cvsflag, dftcis_vec, loose)
        dftcis_vec = libs.lib_func('diag_dftcis', args)

        # Bitci eigenvector scratch number
        dftcis_unit.append(dftcis_vec)

    #
    # RAS1 and RAS3 guess
    #
    # Dominant particle/hole indices
    iph  = np.zeros(nmo, dtype=int)
    iph1 = np.zeros(nmo, dtype=int)

    # Loop over irreps
    for irrep in range(nirr):

        # Number of roots for the current irrep
        nroots = ci_method.n_states_sym(irrep) + n_extra

        # Eigenpair scratch file number
        dftcis_vec = dftcis_unit[irrep]

        # Get the particle/hole MOs corresponding to the
        # dominant DFT/CIS CSFs
        args = (irrep, nroots, cvsflag, dftcis_vec, iph1)
        (dftcis_vec, iph1) = libs.lib_func('ras_guess_dftcis', args)

        # Update the array of dominant particle/hole indices
        iph = np.maximum(iph, np.array(iph1[:], dtype=int))

    # RAS1 & RAS3 MO indices
    ras1 = []
    ras3 = []
    for n in range(nmo):
        if iph[n] == 1:
            if orb_occ[n] == 0:
                # RAS3 MO
                ras3.append(n+1)
            else:
                # RAS1 MO
                ras1.append(n+1)

    # If this is a CVS calculation, then
    # include the HOMO and HOMO-1 in the
    # RAS1 space
    if sum(cvsflag) > 0:
        # HOMO
        iend = np.nonzero(orb_occ)[0][-1]
        e = orb_ener[iend]
        j = iend
        while(True):
            j -= 1
            if (abs(orb_ener[j] - e) > 1e-4):
                istart = j+1
                break
        # HOMO-1
        istart = istart - 1
        e = orb_ener[istart]
        j = istart
        while(True):
            j -= 1
            if (abs(orb_ener[j] - e) > 1e-4):
                istart = j+1
                break
        # Update the RAS1 space
        ras1 = ras1 + [n+1 for n in range(istart, iend+1)]

    # Fill in the RAS arrays
    ci_method.ras1 = np.array(ras1, dtype=int)
    ci_method.ras3 = np.array(ras3, dtype=int)
    
    # Output the selected RAS1 & RAS3 spaces
    print('\n Selected RAS MOs:', flush=True)
    print('\n RAS1',[n for n in ci_method.ras1], flush=True)
    print(' RAS3',[n for n in ci_method.ras3], flush=True)

    # Force nhole1 = nelec3 = 2
    if ci_method.nhole1 != 2:
        print('\n Setting nhole1 = 2', flush=True) 
        ci_method.nhole1 = 2
    if ci_method.nelec3 != 2:
        print('\n Setting nelec3 = 2', flush=True) 
        ci_method.nelec3 = 2
    
    return

@timing.timed
def propagate(ci_method, ci_method0):
    """
    propagates forwards the reference space from a previous
    CI calculation
    """

    # Exit if the point groups for the two calculations are different
    if ci_method.scf.mol.comp_sym != ci_method0.scf.mol.comp_sym:
        sys.exit('\n Error in ref_space.propagate: non-equal point groups')

    # Number of irreps
    nirr = ci_method.n_irrep()
    
    # Number of reference space configurations per irrep
    ref_nconf = np.zeros(nirr, dtype=int)

    # Reference space configurations scratch file number
    confunits = np.zeros(nirr, dtype=int)

    # Names of the previous reference space configuration files
    confnames0 = ci_method0.ref_wfn.conf_name

    # MO overlaps
    mol0 = ci_method0.scf.mol.mol_obj.copy()
    mol  = ci_method.scf.mol.mol_obj.copy()
    smat = gto.intor_cross('int1e_ovlp', ci_method0.scf.mol.mol_obj,
                           ci_method.scf.mol.mol_obj)
    smat = np.matmul(np.matmul(ci_method0.scf.orbs.T, smat), ci_method.scf.orbs)
    nmo0 = ci_method0.scf.nmo
    nmo  = ci_method.scf.nmo
    smat = np.reshape(smat, (nmo0 * nmo), order='F')
    
    # Create the reference spaces for this calculation
    for irr in range(nirr):
        nconf = 0
        unit  = 0
        args  = (nmo0, nmo, smat, confnames0[irr], nconf, unit)
        (nconf, unit) = libs.lib_func('ref_space_propagate', args)
        
    # We will need to pass back the names of the ref space configuration
    # files that were created in this function...
    sys.exit()
    
    return ref_nconf, confunits
