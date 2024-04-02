"""
Orbitals module. This handles the construction and output of molecular
 orbitals from each method and class.
"""
import os as os
import numpy as np
import importlib as importlib
import graci.io.output as output
import graci.utils.timing as timing

import scipy as sp

#
@timing.timed
def build_nos(rdm, basis='AO', mos=None):
    """print the natural orbitals for each state in states list. 
    
       Arguments:
         rdms:   A reference to a function that returns the RDM of state i
         states: A list of states for which to compute the NOs
         basis:  A string argument either "AO" or "MO"
         mos:    (nao, nmo) matrix of the MOs to be used to transform
                 to the AO basis. If basis is "MO", this is ignored.

       Returns:
         natocc:  A (nmo, nst) matrix of NO occuptations
         nos:     The natural orbitals in either the AO or MO basis
    """

    # get dimensions of the NOs
    if mos is not None:
        (dim1, dim2) = mos.shape
    else:
        (dim1, dim2) = rdm.shape

    natorb = np.zeros((dim1, dim2), dtype=float)
    natocc = np.zeros((dim2), dtype=float)

    occ, nos = np.linalg.eigh(rdm)
    if basis.lower() == 'ao':
        nos = np.matmul(mos, nos)

    # we will explicitly order orbitals by occ no. here
    # (even though default eigh should be ascending), since we
    # choose more elaborate order function at later date
    ordr     = order_orbs(occ)
    occ      = occ[ordr]
    nos      = nos[:, ordr]

    return occ, nos


@timing.timed
def build_ndos(rdm, rdm_ref, basis='mo', mos=None):
    """ build natural difference orbitals"""

    # get dimensions of the NOs
    if mos is not None:
        (nbas, nmo) = mos.shape
    else:
        (nbas, nmo) = rdm.shape

    # delta is the different 1RDM between ist and
    # reference state (likely the ground state)
    delta = rdm - rdm_ref

    # form different natural orbitals and transform
    wt, orbs = np.linalg.eigh(delta)

    # transform to AO basis if requested
    if basis.lower() == 'ao':
        orbs = np.matmul(mos, orbs)

    # sort NDO wts by increasing magnitude (hole 
    # orbitals to start, then particle
    ndos = np.zeros((nbas, nmo), dtype=complex)
    wts  = np.zeros((nmo), dtype=complex)

    ordr = np.argsort(wt)
    ndos = orbs[:,ordr]
    wts  = wt[ordr]

    return wts, ndos

@timing.timed
def build_ntos(tdm, basis='ao', mos=None):
    """build the natural transition orbitals"""
    
    # get dimensions of the NOs
    if mos is not None:
        (nbas, nmo) = mos.shape
    else:
        (nbas, nmo) = tdm.shape

    # first perform SVD of 1TDMs to get hole and
    # particle orbitals and weights and convert
    # orbitals to AO basis
    part, s, hole = np.linalg.svd(tdm)
    
    if basis.lower() == 'ao':
        part = np.matmul(mos, part)
        hole = np.matmul(mos, hole.T)

    # normalize singular values
    s  = 0.5 * np.square(s)

    # sort the NTO amplitudes by decreasing magnitude
    ntos = np.zeros((2, nbas, nmo), dtype=complex)
    wts  = np.zeros((2, nmo), dtype=complex)

    ordr     = np.flip(np.argsort(s))

    # hole orbs
    ntos[0,:,:] = hole[:,ordr]
    wts[0,:]    = -s[ordr]
  
    # particle orbs
    ntos[1,:,:] = part[:,ordr]
    wts[1,:]    = s[ordr]     

    return wts, ntos 

#
def promotion_numbers(wts, ndos):
    """compute the detachment and attachment promotion numbers
       for a given set of NDOs -- assumes MO basis"""

    wt_mat   = np.diag([wt if wt < 0. else 0. for wt in wts])
    d_mat    = ndos @ wt_mat @ ndos.transpose()
    p_detach = np.trace(d_mat)

    wt_mat   = np.diag([wt if wt > 0. else 0. for wt in wts])
    d_mat    = ndos @ wt_mat @ ndos.transpose()
    p_attach = np.trace(d_mat)

    return p_detach.real, p_attach.real

#
def export_orbitals(fname, mol, orbs,
                    orb_occ=None, orb_ener=None, orb_sym=None,
                    fmt='molden', orb_dir=None, cart=None):
    """export orbitals to various file formats

       Arguments:
       fname:    the name of the file to be exported
       mol:      the corresponding GRaCI molecule object
       orbs:     an (nao, nmo) matrix correpsonding to nmo MO vectors
       orb_occ:  the (nmo) length set of occupation numbers
                 [ignored if not present or None]
       orb_ener: energies of the orbitals 
                 [ignored if not present or None]
       orb_sym:  a (nmo) length vector of symmetry indices for the NOs 
                 [ignored if not present or None]
       fmt:      the output orbital format (default: molden)
       orb_dir:  If True, outputs orbitals to an 'orbs' subdirectory.
                 The directory is created if it doesn't exist
       cart:     If True, AOs are in cartesian d/f/etc. functions


       Returns:
        None
    """

    if orb_dir is not None:
        orb_file = orb_dir+'/'+str(fname)
        # if dir_name directory doesn't exist, create it
        if not os.path.exists(orb_dir):
            os.mkdir(orb_dir)
    else:
        orb_file = str(fname)

    # if a file called fname exists, delete it
    if os.path.isfile(fname):
        os.remove(fname)

    # import the appropriate library for the file_format
    if fmt in output.orb_formats:
        orbmod = importlib.import_module('graci.io.'+fmt)
    else:
        print('orbital format type='+fmt+' not found. exiting...')
        sys.exit(1)

    orbmod.write_orbitals(orb_file, mol, orbs.real,
                           occ=orb_occ, ener=orb_ener, 
                           sym_lbl=orb_sym, cart=cart)

    return

#--------------------------------------------------------------------------
#
# Internal routines
#
#---------------------------------------------------------------------------

#
def order_orbs(vec):
    """sort the occupation vector 'occ'. This is a dedicated
       function for orbitals because I can envision more elaborate
       sorting criteria than just occupation number (i.e. core orbitals
       first, etc.)"""

    # this sorts ascending
    ordr = np.argsort(vec)

    # flip so largest occupations come first
    return np.flip(ordr)

#
def get_orb_sym(scf, orbs):
    """determine the symmetry of each orb in orbs. Right now we just use
       the symmetry of the scf/ks orbitals to determine that. We can 
       get more elaborate later if need be -- only works for orbitals
       in AO basis is primarily for printing purposes at the moment"""

    nmo      = len(orbs[0,:])
    sym_indx = np.zeros((nmo), dtype=int)

    # take the overlap of the orb with each KS orbital. We will
    # assign the symmetry of this orbital to be the same as the 
    # symmetry of the KS orbital with maximum overlap
    for orb in range(nmo):
        olap = [abs(np.dot(orbs[:,orb], scf.orbs[:,i]))
                for i in range(nmo)]
        maxo = olap.index(max(olap))
        sym_indx[orb] = scf.orb_sym[maxo]

    # make sure the symmetry orbital counts match up
    oval, ocnts = np.unique(sym_indx, return_counts=True)
    val,  cnts  = np.unique(scf.orb_sym, return_counts=True)

    #if not (ocnts==cnts).all():
    #    sys.exit('ERROR: cannot assign symmetry of natural orbs')
    return sym_indx
