"""Module for performing file operations"""

import sys
import re as re
import numpy as np 
import graci.utils.timing as timing
import graci.utils.constants as constants
import graci.core.geometry as geometry
import graci.core.params as params
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
    for obj_name in params.valid_objs:

        # parse all section keywords
        run_list.extend(parse_section(obj_name, 
                                      input_file)) 

    # check the input
    check_input(run_list)

    # create the pyscf molecule objects
    for run_obj in run_list:
        if run_obj.name() == 'molecule':
 
            # loop over geom objects and link a geometry
            # to a molecule object via the label
            for obj in run_list:

                if obj.name() == 'geometry':

                    # generate the pymol object
                    if obj.label == run_obj.label:
                        run_obj.set_pymol(obj)

    return run_list

#
def parse_section(obj_name, input_file):
    """extract keywords from section"""

    nlines       = len(input_file)
    section_objs = []

    # find the start of the section
    iline = 0
    while iline < nlines:
        if '$'+str(obj_name) in input_file[iline]:
            iline += 1
            
            # section exists, create class object
            obj = params.name2obj(obj_name)

            # if this is a geometry section, initialize cart and atoms
            # arrays. This is a bit of a hack..
            if obj_name == 'geometry':
                cart = []
                atom = []

            # first line is charge and multiplicity
            while iline < nlines:
                if '$end' in input_file[iline] or iline == nlines-1:
                    # if we hit end of section, or end of file,
                    # add object to return list and continue parsing 
                    # the input

                    # if this is geometry section and xyz file was 
                    # specified, let's load this now:
                    if obj_name == 'geometry':
                        set_geometry(obj, cart, atom)

                    section_objs.extend([obj])
                    break

                elif '=' in input_file[iline]:
                    (kword,value) = input_file[iline].split('=')
                    kword = kword.strip().lower()

                    if kword in params.kwords[obj_name].keys():
                        val      = parse_value(value)
                        expected = params.kwords[obj_name][kword]

                        if correct_type(val, expected):
                            setattr(obj, kword, val)
                        else:
                            print(kword+" is wrong type, expected: "
                              +str(expected), flush=True)

                # ...else, try to parse an a cartesian structure
                elif obj_name == 'geometry' and \
                        len(input_file[iline].strip()) > 0:

                    line = input_file[iline].split()
                    try:
                        atm_indx = geometry.atom_name.index(line[0].upper())
                        atom.append(geometry.atom_name[atm_indx])
                    except ValueError:
                        sys.exit('atom '+str(line[0])+' not found.')
                    try:
                        cart.append([float(line[i]) for i in range(1,4)])
                    except:
                        sys.exit('Cannot interpret input as a geometry')

                # iterate section loop
                iline += 1

        # iterate line in section search
        iline += 1

    return section_objs

#
def set_geometry(gm_obj, cart, atm):
    """parse xyz file, return the list of atoms and cartesian
       coordinates as lists"""

    atoms     = []
    cartesian = []

    #if gm_obj.units == 'angstrom':
    #    conv = constants.ang2bohr
    #else:
    #    conv = 1.

    # if no xyz_file is specified, use the contents of of th cart 
    # and atm arrays. If those are empty, quit since no geometry
    # is found
    if gm_obj.xyz_file is None and len(cart)==0 and len(atm)==0:
        output.print_message('No geometry file specified and no '+
             'cartesian structure found in input file. Exiting..\n')
        sys.exit()

    # if no xyz file specified, use contents of cart/atm arguments
    if gm_obj.xyz_file is None:
        atoms     = atm
        cartesian = cart

    else:
        # parse contents of xyz file
        try:
            with open(gm_obj.xyz_file, 'r') as xyzfile:
                xyz_gm = xyzfile.readlines()
        except:
            output.print_message('xyz_file: '
                  +str(gm_obj.xyz_file)+' not found.')
            sys.exit()

        # use the number of atoms rather than number of lines in file
        natm = int(xyz_gm[0].strip())

        for i in range(2, natm+2):
            line = xyz_gm[i]
            try:
                atm_indx = geometry.atom_name.index(line[0].upper())
                atoms.append(geometry.atom_name[atm_indx])
            except ValueError:
                sys.exit('atom '+str(line[0])+' not found.')
            try:
                cartesian.append([float(line[i] * conv) 
                    for i in range(1,4)])
            except:
                sys.exit('Cannot interpret input as a geometry')

    gm_obj.set_atoms(atoms)
    gm_obj.set_geom(cartesian)

    return

#
def correct_type(value, keyword_type):
    """determine if 'value' is the correct type of value specified by
       keyword_type. This function will return true if value is an array
       and all elements are of 'keyword_type'"""

    correct = False
    if isinstance(value, np.ndarray):
        val_list = value.tolist()
        if all([isinstance(elem, keyword_type) for elem in val_list]):
            correct = True
    else:
        if isinstance(value, keyword_type):
            correct = True

    return correct

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

    # we have the ability to recognize ranges, specified as either:
    # 1-4 or 1 - 4 or 1:4 or 1 : 4, etc.
    # first see if this is a range input. NOTE: we will not accept
    # ranges split over multiple lines at this time. Also, we will
    # only take the integers adjacent to the range symbol as being
    # relevant. So: 1 2 - 4 5 is intepreted as 2-4
    range_chk = []
    for i in range(len(split_lines[0][:])):
        range_chk.extend(re.split('(:)', split_lines[0][i]))
    if ':' in range_chk:
        # remove empty strings
        try:
            while True:
                range_chk.remove('')
        except ValueError:
            pass
        range_index = range_chk[:].index(':')

        if range_index == 0 or range_index == len(range_chk[:])-1:
            sys.exit(' unknown range specified: '+str(valstr))

        try:
            rstart = int(range_chk[range_index-1])
            rend   = int(range_chk[range_index+1])+1
        except:
            sys.exit('cannot convert range values: '+str(valstr))

        split_lines[0] = list(range(rstart,rend))

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
    
    for obj in run_list:
        
        # Make sure that nstates is an array - in the case of C1
        # symmetry this will be a single integer, but it needs
        # to be a numpy array for use later on
        if 'nstates' in params.kwords[obj.name()].keys():
            if isinstance(obj.nstates, int):
                obj.nstates = np.array([obj.nstates], dtype=int)

        # Make sure that all the RAS entries are numpy arrays.
        ras_key = ['ras1','ras2','ras3']
        for key in ras_key:
            if key in params.kwords[obj.name()].keys():
                if isinstance(getattr(obj, key), int):
                    setattr(obj, key, 
                            np.array([getattr(obj, key)], dtype=int))
                elif isinstance(getattr(obj, key), list):
                    setattr(obj, key,
                            np.array(getattr(obj, key), dtype=int))

        # Make sure that the icvs entry is a numpy array
        if 'icvs' in params.kwords[obj.name()].keys():
            if isinstance(obj.icvs, int):
                obj.icvs = np.array([obj.icvs], dtype=int)

        # If RR-DF is enabled, make sure that the DF flag is
        # set to True
        if 'use_rrdf' in params.kwords[obj.name()].keys():
            if (obj.use_rrdf and not obj.use_df):
                obj.use_df = True
     
        # init/final_states and i/fstate_array need to be lists, also:
        # internal state ordering is 0->n-1, vs. 1->n for input
        if obj.name() == 'transition':
            if obj.init_states is not None:
                if not isinstance(obj.init_states, (list, np.ndarray)):
                    obj.init_states = [obj.init_states]
                # shift state index to range 1 -> 0
                obj.init_states = [obj.init_states[i] - 1 
                                for i in range(len(obj.init_states))]

            if obj.final_states is not None:
                if not isinstance(obj.final_states, (list, np.ndarray)):
                    obj.final_states = [obj.final_states]
                obj.final_states = [obj.final_states[i] - 1 
                               for i in range(len(obj.final_states))]

            if obj.init_states_sym is not None:
                if not isinstance(obj.init_states_sym, (list, np.ndarray)):
                    obj.init_states_sym = [obj.init_states_sym]
                # shift irrep and state to range 1.1 -> 0.0
                obj.init_states_sym = [obj.init_states_sym[i] - 1.1
                               for i in range(len(obj.init_states_sym))]

            if obj.final_states_sym is not None:
                if not isinstance(obj.final_states_sym, (list, np.ndarray)):
                    obj.final_states_sym = [obj.final_states_sym]
                obj.final_states_sym = [obj.final_states_sym[i] - 1.1
                               for i in range(len(obj.final_states_sym))]
            
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

