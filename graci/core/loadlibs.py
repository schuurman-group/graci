"""Module for initializing compiled libraries """

import os
import sys
import numpy as np
import ctypes as ctypes
import graci.io.convert as convert

# THIS IS TEMPORARY UNTIL WE PASS THE HAMILTONIAN NAME
# DIRECTLY. THE METHOD CALLING init_bitci IS RESPONSIBLE 
# FOR "HAMILTONIAN" BEING VALID
hamiltonians   = ['canonical',
                  'grimme_standard',
                  'grimme_short',
                  'lyskov_standard',
                  'lyskov_short',
                  'heil17_standard',
                  'heil18_short']

#
def init_bitci(mol, scf, ci):
    """Initialize the bitci library"""
    
    # load the appropriate library
    bitci_path = os.environ['GRACI']+'/graci/lib/bitci/lib/libbitci.{}'.format(
        'so' if sys.platform != 'darwin' else 'dylib')
    
    if not os.path.isfile(bitci_path):
        raise FileNotFoundError('bitci library not found: '+bitci_path)
    lib_bitci = ctypes.cdll.LoadLibrary(bitci_path)

    # set all variable that have to be passed to bitci_initialise
    # (note that the pgrp and iham variables use Fortran indexing)
    imult = convert.convert_ctypes(mol.mult,               dtype='int32')
    nel   = convert.convert_ctypes(mol.nel,                dtype='int32')
    nmo   = convert.convert_ctypes(scf.nmo,                dtype='int32')
    mosym = convert.convert_ctypes(np.array(scf.orb_sym),  dtype='int64')
    moen  = convert.convert_ctypes(np.array(scf.orb_ener), dtype='double')
    isym  = mol.sym_indx + 1 if mol.sym_indx > 0 else 1
    pgrp  = convert.convert_ctypes(isym,                   dtype='int32')
    enuc  = convert.convert_ctypes(mol.enuc,               dtype='double')
    iham  = convert.convert_ctypes(hamiltonians.index(ci.hamiltonian)+1,
                                   dtype='int32')
    label = convert.convert_ctypes(ci.label, dtype='string')

    # call to bitci_initialise
    lib_bitci.bitci_initialise(ctypes.byref(imult),
                               ctypes.byref(nel),
                               ctypes.byref(nmo),
                               mosym,
                               moen,
                               ctypes.byref(pgrp),
                               ctypes.byref(enuc),
                               ctypes.byref(iham),
                               label)

    return lib_bitci

#
def init_intpyscf(mol, scf):
    """Initialize the int_pyscf library"""

    # load the appropriate library
    intpyscf_path = os.environ['GRACI']+'/graci/lib/integrals/lib/libint_pyscf.{}'. \
            format('so' if sys.platform != 'darwin' else 'dylib')

    if not os.path.isfile(intpyscf_path):
        raise FileNotFoundError('int_pyscf library not found: '+intpyscf_path)

    lib_intpyscf = ctypes.cdll.LoadLibrary(intpyscf_path)

    # set the variables that have to be passed to load_mo_integrals
    nmo     = convert.convert_ctypes(scf.nmo,    dtype='int32')
    naux    = convert.convert_ctypes(scf.naux,   dtype='int32')
    use_df  = convert.convert_ctypes(mol.use_df, dtype='logical')
    thresh  = convert.convert_ctypes(1e-14,      dtype='double')
    max_mem = convert.convert_ctypes(-1,         dtype='int32')   

    # call to load_mo_integrals
    lib_intpyscf.load_mo_integrals(ctypes.byref(nmo),
                                   ctypes.byref(naux),
                                   ctypes.byref(use_df),
                                   ctypes.byref(thresh),
                                   ctypes.byref(max_mem))
    
    return lib_intpyscf
