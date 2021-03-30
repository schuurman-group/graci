"""Module for performing file operations"""

import sys
import numpy as np 
import graci.utils.timing as timing
import graci.methods.params as params
import graci.methods.molecule as molecule
import graci.io.output as output
from pyscf import gto

#valid input sections

#
def parse_input():
    """parse the graci input file"""

    # Read input file. Small enough to gulp the whole thing
    with open(output.file_names['input_file'], 'r') as infile:
        input_file = infile.readlines()

    run_list = []
    for valid_section in params.valid_methods:

        # parse all section keywords
        run_list.extend(parse_section(valid_section, input_file)) 

    # check the input
    check_input(run_list)

    # add the geometry to the molecule object
    for run_obj in run_list:
        if params.method_name(run_obj) == 'molecule':
            parse_geometry(run_obj, input_file)

            atms = run_obj.atoms
            cart = run_obj.geom
            mol_str = ';'.join([atms[i]+'   '+
                      ' '.join([str(cart[i,j]) for j in range(3)])
                        for i in range(run_obj.n_atoms())])

            pymol = gto.M(
                   dump_input = False,
                   parse_arg  = False,
                   verbose    = 0,
                   atom       = mol_str,
                   charge     = run_obj.charge,
                   spin       = run_obj.spin,
                   output     = output.file_names['pyscf_out'],
                   basis      = run_obj.basis,
                   symmetry   = run_obj.use_sym,
                   unit       = run_obj.units)
            pymol.build()

            run_obj.nel      = int(sum(pymol.nelec))
            run_obj.full_sym = pymol.topgroup.lower()
            run_obj.comp_sym = pymol.groupname.lower()
            if pymol.symmetry:
                run_obj.sym_indx = molecule.point_grps.index(run_obj.comp_sym)
            else:
                run_obj.sym_indx = -1

    return run_list

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
def parse_section(section, input_file):
    """extract keywords from section"""

    nlines       = len(input_file)
    section_objs = []

    # find the start of the section
    iline = 0
    while iline < nlines:
        if '$'+str(section) in input_file[iline]:
            iline += 1
            
            # section exists, create class object
            method = params.method_obj(section)

            # first line is charge and multiplicity
            while iline < nlines:
                if '$end' in input_file[iline] or iline == nlines-1:
                    # if we hit end of section, or end of file,
                    # add object to return list and continue parsing the input
                    section_objs.extend([method])
                    break

                (kword,value) = input_file[iline].split('=')
                kword = kword.strip().lower()

                if kword in params.kwords[section].keys():
                    val = parse_value(value)

                    correct_type = False
                    if isinstance(val, np.ndarray):
                        val_list = val.tolist()
                        if all([isinstance(elem, 
                                       params.kwords[section][kword])
                                       for elem in val_list]):
                            correct_type = True
                    else:
                        if isinstance(val, params.kwords[section][kword]):
                            correct_type = True

                    if correct_type:
                        setattr(method, kword, val)
                    else:
                        print(kword+" value is wrong type, expecting "
                          +str(params.kwords[section][kword]),flush=True)

                # iterate section loop
                iline += 1

        # iterate line in section search
        iline += 1

    if len(section_objs) == 0:
        print(section+' section not found in input file')

    return section_objs

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
def check_input(run_list):
    """Checks on the user-supplied input"""
    
    # Make sure that nstates is an array - in the case of C1
    # symmetry this will be a single integer, but it needs
    # to be a numpy array for use later on
    for obj in run_list:
        sec_name = params.method_name(obj)

        # Make sure that nstates is an array - in the case of C1
        # symmetry this will be a single integer, but it needs
        # to be a numpy array for use later on
        if 'nstates' in params.kwords[sec_name].keys():
            if isinstance(obj.nstates, int):
                obj.nstates = np.array([obj.nstates], dtype=int)

        # Make sure that all the RAS entries are numpy arrays.
        ras_key = ['ras1','ras2','ras3']
        for key in ras_key:
            if key in params.kwords[sec_name].keys():
                if isinstance(getattr(obj, key), int):
                    setattr(obj, key, 
                            np.array([getattr(obj, key)], dtype=int))
                elif isinstance(getattr(obj, key), list):
                    setattr(obj, key,
                            np.array(getattr(obj, key), dtype=int))

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

