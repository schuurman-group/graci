"""Module for performing file operations"""

import sys
import ctypes
import numpy as np 
import graci.utils.timing as timing
import graci.molecule.molecule as molecule
import graci.methods.params as params

#
def parse_input():
    """parse the d3CI input file"""

    # Read input file. Small enough to gulp the whole thing
    with open(params.user_inp['input_file'], 'r') as infile:
        input_file = infile.readlines()

    keyword_types = {**params.mol_kword, 
                     **params.scf_kword, 
                     **params.mrci_kword,
                     **params.asci_kword,
                     **params.file_kword}

    mol    = molecule.Molecule()

    # parse geometry section
    parse_geometry(mol, input_file)

    # parse molecule section
    parse_section('molecule', params.mol_kword, 
                              params.mol_param, input_file)

    # parse SCF keywords
    parse_section('scf', params.scf_kword, 
                         params.scf_param, input_file)

    # parse MRCI keywords
    parse_section('mrci', params.mrci_kword, 
                          params.mrci_param, input_file)

    # check the input
    check_input()
    
    # update molecule object with relevant info
    mol.charge = params.mol_param['charge']
    mol.mult   = params.mol_param['mult']    

    # use the pyscf symmetry analyzer
    atms = mol.atoms
    cart = mol.geom
    mol_str = ';'.join([atms[i]+'   '+
              ' '.join([str(cart[i,j]) for j in range(3)])
                        for i in range(mol.n_atoms())])
    
    pymol = gto.M(
        dump_input = False,
        parse_arg  = False,
        verbose    = 0,
        atom       = mol_str,
        charge     = mol.charge,
        spin       = mol.spin,
        output     = params.file_param['pyscf_out'],
        basis      = params.mol_param['basis'],
        symmetry   = params.mol_param['use_sym'],
        unit       = parsms.mol_param['units'])
    pymol.build()

    mol.nel      = int(sum(pymol.nelec))
    mol.full_sym = pymol.topgroup.lower()
    mol.comp_sym = pymol.groupname.lower()
    if pymol.symmetry:
        mol.sym_indx = molecule.point_grps.index(mol.comp_sym)
    else:
        mol.sym_indx = -1
    
    return mol

#
def parse_geometry(mol, input_file):
    """extract molecular geometry"""

    nlines = len(input_file)

    # find the start of the section
    for i in range(nlines):
        if '$geometry' in input_file[i]:
            break
        elif i == nlines - 1:
            raise IOError('geometry section not found in input file')
    
    iline = i+1
    cart  = []
    atoms = []
    while iline < nlines:

        if '$end' in input_file[iline]:
            break
        line = input_file[iline].split()
        try:
            atm_indx = [molecule.atom_name[i].upper() 
              for i in range(len(molecule.atom_name))].index(line[0].upper())
        except ValueError:
            sys.exit('atom '+str(len[0])+' not found.')
        atoms.append(molecule.atom_name[atm_indx])
        try:
            cart.append([float(line[i]) for i in range(1,4)])
        except:
            sys.exit('Cannot interpret input as a geometry')

        iline += 1
        
    mol.set_geom(np.array(cart, dtype=float))
    mol.set_atoms(atoms)
    return

#
def parse_section(section, kword_types, kwords, input_file):
    """extract keywords from section"""

    nlines = len(input_file)

    # find the start of the section
    for i in range(nlines):
        if '$'+str(section) in input_file[i]:
            break
        elif i == nlines - 1:
            raise IOError(section+' section not found in input file')

    iline = i+1
    # first line is charge and multiplicity
    while iline < nlines:
        if '$end' in input_file[iline]:
            break
        (kword,value) = input_file[iline].split('=')
        kword = kword.strip().lower()

        if kword in kword_types.keys():
            val          = parse_value(value)
            correct_type = False
            if isinstance(val, np.ndarray):
                val_list = val.tolist()
                if all([isinstance(elem, kword_types[kword]) 
                                            for elem in val_list]):
                    correct_type     = True
            
            else:
                if isinstance(val, kword_types[kword]):
                    correct_type     = True                
                    
            if correct_type:
                kwords[kword] = val
            else:
                print(kword+" value is wrong type, expecting "
                            +str(kword_types[kword]),flush=True)
        iline+=1

    return

#
def parse_value(valstr):
    """Returns a value converted to the appropriate type and shape.

    By default, spaces and newlines will be treated as delimiters. The
    combination of both will be treated as a 2D array.
    """
    all_lines = valstr.split('\n')[:-1]
    split_lines = [line.split() for line in all_lines]

    if len(all_lines) > 1 and len(split_lines[0]) == 1:
        # newline and space give the same result for a vector
        split_lines = [[line[0] for line in split_lines]]

    if len(split_lines) == 1:
        if len(split_lines[0]) == 1:
            # handle single values
            return convert_value(split_lines[0][0])
        else:
            # handle vectors
            return convert_array(split_lines[0])
    else:
        # handle 2D arrays
        return convert_array(split_lines)

#
def check_input():
    """Checks on the user-supplied input"""
    
    # Make sure that nstates is an array - in the case of C1
    # symmetry this will be a single integer, but it needs
    # to be a numpy array for use later on
    if isinstance(params.user_inp['nstates'], int):
        arr = np.zeros(1, dtype=int)
        arr[0] = params.user_inp['nstates']
        params.user_inp['nstates'] = arr
    
    return
    
#
def convert_value(val):
    """Converts a string value to NoneType, bool, int, float or string."""
    if val.lower() == 'none':
        return None
    elif val.lower() == 'true':
        return True
    elif val.lower() == 'false':
        return False

    try:
        return int(val)
    except ValueError:
        pass

    try:
        return float(val)
    except ValueError:
        pass

    return val

#
def convert_array(val_list):
    """Converts a list of strings to an array of ints, floats or strings."""
    try:
        return np.array(val_list, dtype=int)
    except ValueError:
        pass

    try:
        return np.array(val_list, dtype=float)
    except ValueError:
        pass

    return np.array(val_list, dtype=str)

#
def convert_ctypes(py_val, dtype=None):
    """convert a python array into a C-data type"""

    # note: the current approach is used based on:
    # https://bugs.python.org/issue27926
    # Namely, the default Python constructor is _very_ slow.
    # if that changes, then this function can change too...

    # there are fancier ways to do this by querying both the array and
    # machine environment, but this will do for now
    if dtype == 'int32':
        type_sym = 'i';i_size = 4; ctype_sym = ctypes.c_int32
    elif dtype == 'int64':
        type_sym = 'i';i_size = 8; ctype_sym = ctypes.c_int64
    elif dtype == 'double':
        type_sym = 'd';i_size = 8; ctype_sym = ctypes.c_double
    elif dtype == 'logical':
        type_sym = 'i';i_size = 4; ctype_sym = ctypes.c_bool
    else:
        sys.exit('convert_ctypes does not recognize dtype='+str(dtype))

    if isinstance(py_val, (float, int)):
        return ctype_sym(py_val)

    if py_val.size == 1:
        c_arr = (ctype_sym * py_val.size)(py_val)
    else:
        c_arr = (ctype_sym * py_val.size)(*py_val)

    return c_arr

