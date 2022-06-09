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
                  'heil17_short',
                  'heil18_standard',
                  'heil18_short',
                  'cvs_standard']

libraries      = ['bitci','bitsi','bitwf']

# registry of bitci functions
bitci_registry = {
    'bitci_initialise'        : ['int32','int32','int32','int64',
                                 'double','int32','double','string',
                                 'string'],
    'bitci_finalise'         : [],
    'bitci_int_initialize'   : ['string', 'string', 'string', 'string'],
    'bitci_int_finalize'     : [],
    'generate_ref_confs'     : ['int32','int32','int32','int32',
                                'int32','int32','int32','int32',
                                'int32','int32','int32','int32'],
    'ref_space_propagate'    : ['int32','int32','int32','double','string',
                                'int32','int32'],
    'diag_dftcis'            : ['int32','int32','int32','int32','logical'],
    'ras_guess_dftcis'       : ['int32','int32','int32','int32',
                                'int32'],
    'ref_diag_mrci'          : ['int32','int32','int32','int32',
                                'int32'],
    'prune_ref_space'        : ['int32','int32','int32','int32','int32'],
    'retrieve_energies'      : ['int32','int32','double'],
    'retrieve_some_energies' : ['int32','int32','double','int32'],
    'generate_mrci_confs'    : ['int32','int32','int32', 'int32',
                                'double','int32','logical'],
    'mrci_prune'             : ['double','int32','int32','int32',
                                'int32', 'int32','int32','int32'],
    'retrieve_qcorr'         : ['int32','int32','int32','int32',
                                'int32','int32','int32','double',
                                'double'],
    'retrieve_filename'      : ['int32','string'],
    'diag_mrci'              : ['int32','int32','int32','int32',
                                'int32','double','int32','int32',
                                'logical','int32','int32'],
    'print_mrci_states'      : ['int32','int32','int32'],
    'print_pmrci_states'     : ['int32','int32','int32','int32','int32',
                                'int32'],
    'refine_ref_space'       : ['int32','int32','int32','int32','int32',
                                'double','double','int32'],
    'refine_ref_space_pt2'   : ['int32','int32','int32','int32','int32',
                                'int32','int32','double','double','double',
                                'double','int32'],
    'density_mrci'           : ['int32','int32','int32','double',
                                'int32','int32'],
    'wf_mrci'                : ['int32','int32','int32','int32',
                                'int32','string','string', 'int32',
                                'int32'],
    'gvvpt2'                 : ['int32','int32','int32','double',
                                'int32','int32','int32',
                                'int32','int32'],
    'truncate_mrci_wf'       : ['int32','int32','int32','int32',
                                'double','int32']
}

bitci_intent = {
    'bitci_initialise'       : ['in','in','in','in','in','in','in',
                                'in','in','in'],
    'bitci_finalise'         : [],
    'bitci_int_initialize'   : ['in','in','in','in'],
    'bitci_int_finalize'     : [],
    'generate_ref_confs'     : ['in','in','in','in','in','in','in',
                                'in','in','in','out','out'],
    'ref_space_propagate'    : ['in','in','in','in','in','out','out'],
    'diag_dftcis'            : ['in','in','in','out','in'],
    'ras_guess_dftcis'       : ['in','in','in','out','out'],
    'ref_diag_mrci'          : ['in','out','in','in','out'],
    'prune_ref_space'        : ['in','in','in','out','in'],
    'retrieve_energies'      : ['in','in','out'],
    'retrieve_some_energies' : ['in','in','out','in'],
    'generate_mrci_confs'    : ['in','in','out','out','in','in','in'],
    'mrci_prune'             : ['in','in','in','in','in','in','out',
                                'out'],
    'retrieve_qcorr'         : ['in','in','in','in','in','in','in',
                                'out','out'],
    'retrieve_filename'      : ['in','out'],
    'diag_mrci'              : ['in','in','in','out','in','in','in',
                                'in','in','in','in'],
    'print_mrci_states'      : ['in','in','in'],
    'print_pmrci_states'     : ['in','in','in','in','in','in'],
    'refine_ref_space'       : ['in','out','in','in','in','in','out',
                                'out'],
    'refine_ref_space_pt2'   : ['in','out','in','in','in','in','in',
                                'in','in','in','out','out'],
    'density_mrci'           : ['in','in','in','out','in','in'],
    'wf_mrci'                : ['in','in','in','in','in','in','in',
                                'in','out'],
    'gvvpt2'                 : ['in','in','in','in','in','out','in',
                                'out','out'],
    'truncate_mrci_wf'       : ['in','in','in','in','in','out']
}


# registry of bitsi functions
bitsi_registry = {
    'bitsi_initialise'        : ['int32','int32','int32','int32','int32',
                                  'int32','string'],
    'bitsi_finalise'          : [],
    'transition_density_mrci' : ['int32','int32','int32','int32','int32',
                                 'int32','double','string','string',
                                 'string','string'],
    'redmat_mrci'             : ['int32','int32','int32','int32','int32',
                                 'int32','double','string','string',
                                 'string','string'],
    'cgcoeff_soc'             : ['int32','int32','double']
}

bitsi_intent   = {
    'bitsi_initialise'        : ['in','in','in','in','in','in','in'],
    'bitsi_finalise'          : [],
    'transition_density_mrci' : ['in','in','in','in','in','in','out',
                                 'in','in','in','in'],
    'redmat_mrci'             : ['in','in','in','in','in','in','out',
                                 'in','in','in','in'],
    'cgcoeff_soc'             : ['in','in','out']
}

# registry of bitwf functions
bitwf_registry = {
    'bitwf_initialise' : ['int32','int32','int32','int32','int32','int32',
                          'double','int32','string'],
    'bitwf_finalise'   : [],
    'detwf'            : ['int32','string','string','int32','string',
                          'int32'],
    'detoverlap'       : ['int32','int32','int32','int32','int32',
                          'int32','int32','double','int32','int32',
                          'logical','double']
}

bitwf_intent = {
    'bitwf_initialise' : ['in','in','in','in','in','in','in','in','in'],
    'bitwf_finalise'   : [],
    'detwf'            : ['in','in','in','in','in','out'],
    'detoverlap'       : ['in','in','in','in','in','in','in','in',
                          'in','in','in','out']
}

# list of existing library objects
lib_objs = {}

def lib_load(name):
    """load a library object and add it to lib_objs dictionary"""
    global lib_objs

    # if we haven't loaded the bitX object, do so now
    if name not in lib_objs.keys() and name in libraries:
        # load the appropriate library
        rel_path = '/graci/lib/lib/lib'+str(name)+'.{}'
        path_str = os.environ['GRACI'] + rel_path
        lib_path = path_str.format('so' if sys.platform != 'darwin' 
                                                      else 'dylib')

        if not os.path.isfile(lib_path):
            raise FileNotFoundError(str(name)+' library not found: '
                                    + lib_path)
        lib_objs[str(name)] = ctypes.cdll.LoadLibrary(lib_path)

    return

def lib_func(name, args):
    """call a function from a fortran library"""
    global bitci_registry, bitci_intent
    global bitsi_registry, bitsi_intent
    global bitwf_registry, bitwf_intent
    global lib_objs

    if name in bitci_registry:
        arg_list   = bitci_registry[name]
        arg_intent = bitci_intent[name]
    elif name in bitsi_registry:
        arg_list   = bitsi_registry[name]
        arg_intent = bitsi_intent[name]
    elif name in bitwf_registry:
        arg_list   = bitwf_registry[name]
        arg_intent = bitwf_intent[name]
    else:
        sys.exit('function: '+str(name)+' not found.') 

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
        # this is hacky -- probably a better way
        if len(args) > 0:
            getattr(lib_objs['bitci'], name)(*arg_ptr)
        else:
            getattr(lib_objs['bitci'], name)()
    elif name in bitsi_registry:
        if len(args) > 0:
            getattr(lib_objs['bitsi'], name)(*arg_ptr)
        else:
            getattr(lib_objs['bitsi'], name)()
    elif name in bitwf_registry:
        if len(args) > 0:
            getattr(lib_objs['bitwf'], name)(*arg_ptr)
        else:
            getattr(lib_objs['bitwf'], name)()
            
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
def lib_exists(name):
    """check if the library name 'name' exists in 
       lib_objs dictionary"""
    global lib_objs

    return name in lib_objs.keys()

