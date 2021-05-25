"""Module for initializing compiled libraries """

import os
import sys
import numpy as np
import copy as copy
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

# registry of bitci functions
bitci_registry = {
    'generate_ref_confs'     : ['int32','int32','int32','int32',
                                'int32','int32','int32','int32',
                                'int32','int32','int32','int32'],
    'diag_dftcis'            : ['int32','int32','int32','int32','logical'],
    'ras_guess_dftcis'       : ['int32','int32','int32','int32',
                                'int32'],
    'ref_diag_mrci'          : ['int32','int32','int32','int32',
                                'int32'],
    'retrieve_energies'      : ['int32','int32','double'],
    'retrieve_some_energies' : ['int32','int32','double','int32'],
    'generate_mrci_confs'    : ['int32','int32','int32','int32',
                                'int32', 'double','int32'],
    'mrci_prune'             : ['double','int32','int32','int32',
                                'int32', 'int32','int32'],
    'retrieve_filename'      : ['int32','string'],
    'diag_mrci'              : ['int32','int32','int32','int32',
                                'int32','double','int32','int32',
                                'logical'],
    'print_mrci_states'      : ['int32','int32','int32'],
    'refine_ref_space'       : ['int32','int32','int32','int32',
                                'double','double','int32'],
    'density_mrci'           : ['int32','int32','int32','double',
                                'int32','int32']
}

bitci_intent = {
    'generate_ref_confs'     : ['in','in','in','in','in','in','in',
                                'in','in','in','out','out'],
    'diag_dftcis'            : ['in','in','in','out','in'],
    'ras_guess_dftcis'       : ['in','in','in','out','out'],
    'ref_diag_mrci'          : ['in','out','in','in','out'],
    'retrieve_energies'      : ['in','in','out'],
    'retrieve_some_energies' : ['in','in','out','in'],
    'generate_mrci_confs'    : ['in','in','in','out','out','in','in'],
    'mrci_prune'             : ['in','in','in','in','in','in','out'],
    'retrieve_filename'      : ['in','out'],
    'diag_mrci'              : ['in','in','in','out','in','in','in',
                                'in','in'],
    'print_mrci_states'      : ['in','in','in'],
    'refine_ref_space'       : ['in','out','in','in','in','out','out'],
    'density_mrci'           : ['in','in','in','out','in','in']
}


# registry of bitsi functions
bitsi_registry = {
    'transition_density_mrci' : ['int32','int32','int32','int32','int32',
                                 'int32','double','string','string',
                                 'string','string']
}

bitsi_intent   = {
    'transition_density_mrci' : ['in','in','in','in','in','in','out',
                                 'in','in','in','in']
}

# list of existing library objects
lib_objs = {}


def lib_func(name, args):
    """call a function from a fortran library"""
    global bitci_registry, bitci_intent
    global bitsi_registry, bitsi_intent, lib_objs

    if name in bitci_registry:
        arg_list   = bitci_registry[name]
        arg_intent = bitci_intent[name]
    elif name in bitsi_registry:
        arg_list   = bitsi_registry[name]
        arg_intent = bitsi_intent[name]
        
    #arg_list   = bitci_registry[name]
    #arg_intent = bitci_intent[name]

    arg_ctype = []
    arg_ptr   = []
    for i in range(len(args)):
        
        # if argument is a string, pad to a length of 255 characters
        if isinstance(args[i], str):
            arg = args[i].ljust(255)
        else:
            arg = args[i]
        c_arg = convert.convert_ctypes(arg, dtype=arg_list[i])

        arg_ctype.append(c_arg)

        if isinstance(arg, (list, np.ndarray, str)):
            arg_ptr.append(c_arg)
        else:
            arg_ptr.append(ctypes.byref(c_arg))

    if name in bitci_registry:
        getattr(lib_objs['bitci'], name)(*arg_ptr)
    elif name in bitsi_registry:
        getattr(lib_objs['bitsi'], name)(*arg_ptr)
        
    #getattr(lib_objs['bitci'], name)(*arg_ptr)

    args_out = ()
    for i in range(len(args)):
        if arg_intent[i] == 'out':
            if isinstance(args[i], list):
                args_out += (np.ndarray((len(args[i]),), 
                          buffer=arg_ctype[i], dtype=arg_list[i]),)
            elif isinstance(args[i], np.ndarray):
                args_out += (np.ndarray((args[i].size,),
                          buffer=arg_ctype[i], dtype=arg_list[i]),)
            elif isinstance(args[i], str):
                args_out += (bytes.decode(arg_ctype[i].value),)
            else:
                args_out += (arg_ctype[i].value,)
    
    if len(args_out) == 1:
        return args_out[0]
    else:
        return args_out

#
def init_bitci(ci_method):
    """Initialize the bitci library"""
    global lib_objs
    
    # load the appropriate library
    path_str   = os.environ['GRACI']+'/graci/lib/bitci/lib/libbitci.{}'
    bitci_path = path_str.format('so' if sys.platform != 'darwin' else 'dylib')
    
    if not os.path.isfile(bitci_path):
        raise FileNotFoundError('bitci library not found: '+bitci_path)
    lib_objs['bitci'] = ctypes.cdll.LoadLibrary(bitci_path)

    # set the variables that have to be passed to intpyscf_initialise
    nmo      = convert.convert_ctypes(ci_method.scf.nmo,    
               dtype='int32')
    naux     = convert.convert_ctypes(ci_method.scf.naux,   
               dtype='int32')
    use_df   = convert.convert_ctypes(ci_method.mol.use_df, 
               dtype='logical')
    use_rrdf = convert.convert_ctypes(ci_method.mol.use_rrdf,   
               dtype='logical')
    rrdf_fac = convert.convert_ctypes(ci_method.mol.rrdf_fac,
               dtype='int32')
    thresh   = convert.convert_ctypes(1e-14,                
               dtype='double')
    max_mem  = convert.convert_ctypes(-1,                   
               dtype='int32')   

    # call to intpyscf_initialise
    lib_objs['bitci'].intpyscf_initialise(ctypes.byref(nmo),
                                  ctypes.byref(naux),
                                  ctypes.byref(use_df),
                                  ctypes.byref(use_rrdf),
                                  ctypes.byref(rrdf_fac),
                                  ctypes.byref(thresh),
                                  ctypes.byref(max_mem))
    
    # set all variable that have to be passed to bitci_initialise
    # (note that the pgrp and iham variables use Fortran indexing)
    imult = convert.convert_ctypes(ci_method.mol.mult,               
            dtype='int32')
    nel   = convert.convert_ctypes(ci_method.mol.nel,                
            dtype='int32')
    mosym = convert.convert_ctypes(np.array(ci_method.scf.orb_sym),  
            dtype='int64')
    moen  = convert.convert_ctypes(np.array(ci_method.scf.orb_ener), 
            dtype='double')

    if ci_method.mol.sym_indx <= 0:
        isym = 1
    else:
        isym = ci_method.mol.sym_indx + 1

    pgrp  = convert.convert_ctypes(isym,                   
            dtype='int32')
    enuc  = convert.convert_ctypes(ci_method.mol.enuc,               
            dtype='double')
    iham  = convert.convert_ctypes(hamiltonians.index(ci_method.hamiltonian)+1,
            dtype='int32')
    label = convert.convert_ctypes(ci_method.label, 
            dtype='string')
    
    # call to bitci_initialise
    lib_objs['bitci'].bitci_initialise(ctypes.byref(imult),
                               ctypes.byref(nel),
                               ctypes.byref(nmo),
                               mosym,
                               moen,
                               ctypes.byref(pgrp),
                               ctypes.byref(enuc),
                               ctypes.byref(iham),
                               label)
    
    return 

#
def finalise_bitci():
    """Finalises the bitci library"""
    lib_objs['bitci'].bitci_finalise()
    return

#
def init_bitsi(si_method):
    """Initialize the bitsi library"""
    global lib_objs

    # load the appropriate library
    path_str   = os.environ['GRACI']+'/graci/lib/bitci/lib/libbitsi.{}'
    bitsi_path = path_str.format('so' if sys.platform != 'darwin' else 'dylib')
    
    bra_method = si_method.final_method
    ket_method = si_method.init_method

    if not os.path.isfile(bitsi_path):
        raise FileNotFoundError('bitsi library not found: '+bitsi_path)
    lib_objs['bitsi'] = ctypes.cdll.LoadLibrary(bitsi_path)

    # if the number of mos is different between bra and ket, end
    if bra_method.scf.nmo != ket_method.scf.nmo:
        sys.exit('bra_method.scf.nmo != ket_method.scf.nmo init_bitsi')
    
    # set all variables that have to be passed to bitsi_initialise
    # (note that the pgrp uses Fortran indexing)
    imultBra = convert.convert_ctypes(bra_method.mol.mult,dtype='int32')
    imultKet = convert.convert_ctypes(ket_method.mol.mult,dtype='int32')
    nelBra   = convert.convert_ctypes(bra_method.mol.nel, dtype='int32')
    nelKet   = convert.convert_ctypes(ket_method.mol.nel, dtype='int32')
    nmo      = convert.convert_ctypes(bra_method.scf.nmo, dtype='int32')

    # not sure what this is about...
    if bra_method.mol.sym_indx <= 0:
        isym = 1
    else:
        isym = bra_method.mol.sym_indx + 1
    pgrp  = convert.convert_ctypes(isym, dtype='int32')
    
    # call to bitsi_initialise
    lib_objs['bitsi'].bitsi_initialise(ctypes.byref(imultBra),
                               ctypes.byref(imultKet),
                               ctypes.byref(nelBra),
                               ctypes.byref(nelKet),
                               ctypes.byref(nmo),
                               ctypes.byref(pgrp))
        
    return

#
def finalise_bitsi():
    """Finalises the bitsi library"""
    lib_objs['bitsi'].bitsi_finalise()
    return
