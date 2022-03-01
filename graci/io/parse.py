"""Module for performing file operations"""

import sys
import re as re
import numpy as np 
import graci.utils.constants as constants
import graci.core.params as params
import graci.io.output as output

import graci.core.molecule as molecule
import graci.methods.scf as scf
import graci.methods.dftmrci as dftmrci
import graci.methods.dftmrenpt2 as dftmrenpt2
import graci.interaction.transition as transition
import graci.interaction.spinorbit as spinorbit

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

    return run_list

#
def parse_section(class_name, input_file):
    """extract keywords from section"""

    nlines       = len(input_file)
    section_objs = []
    mod_name     = str(class_name.lower())

    # find the start of the section
    iline = 0
    while iline < nlines:
        if '$'+mod_name in input_file[iline].lower():
            iline += 1
            
            # section exists, create class object
            sec_obj = getattr(globals()[mod_name], class_name)()

            # if this is a molecule section, look for an xyz geometry
            # input format
            if class_name == 'Molecule':
                cart = []
                atom = []

            # first line is charge and multiplicity
            while iline < nlines:
                if '$end' in input_file[iline] or iline == nlines-1:
                    # if we hit end of section, or end of file,
                    # add object to return list and continue parsing 
                    # the input

                    # if this is molecule section and we're exiting,
                    # set the cart/atm arrays to the cart/atom lists.
                    # Molecule will check it when run() is called. 
                    if class_name == 'Molecule':
                        sec_obj.set_geometry(atom, cart)

                    section_objs.extend([sec_obj])
                    break

                elif '=' in input_file[iline]:
                    (kword,value) = input_file[iline].split('=')
                    kword = kword.strip().lower()

                    if kword in params.kwords[class_name].keys():
                        val      = parse_value(value)
                        expected = params.kwords[class_name][kword]

                        if correct_type(val, expected):
                            setattr(sec_obj, kword, val)
                        else:
                            print(kword+" is wrong type, expected: "
                              +str(expected), flush=True)

                elif class_name == 'Molecule' and \
                        len(input_file[iline].strip()) > 0:
                    # try to interpret as a cartesian atom definition

                    line = input_file[iline].split()

                    # if line is length 4 and comprised of one string 
                    # followed by 3 numbers, we'll take it for now and
                    # check it later
                    if len(line) == 4:
                        crds = convert_array(line[1:])
                        if crds.dtype==np.int or crds.dtype==np.float:
                            cart.append(crds.tolist())
                            atom.append(line[0])
                        
                # iterate section loop
                iline += 1

        # iterate line in section search
        iline += 1

    return section_objs

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
   
    # if any section label names are repeated, append a number to 
    # ensure that each section label is unique
    section_lbls = []

    for obj in run_list:
    
        # the class name for this run object
        obj_name = type(obj).__name__

        # Make sure that nstates is an array - in the case of C1
        # symmetry this will be a single integer, but it needs
        # to be a numpy array for use later on
        if 'nstates' in params.kwords[obj_name].keys():
            if isinstance(obj.nstates, int):
                obj.nstates = np.array([obj.nstates], dtype=int)

        # Make sure that all the RAS entries are numpy arrays.
        ras_key = ['ras1','ras2','ras3']
        for key in ras_key:
            if key in params.kwords[obj_name].keys():
                if isinstance(getattr(obj, key), int):
                    setattr(obj, key, 
                            np.array([getattr(obj, key)], dtype=int))
                elif isinstance(getattr(obj, key), list):
                    setattr(obj, key,
                            np.array(getattr(obj, key), dtype=int))

        # Make sure that the icvs entry is a numpy array
        if 'icvs' in params.kwords[obj_name].keys():
            if isinstance(obj.icvs, int):
                obj.icvs = np.array([obj.icvs], dtype=int)

        # the basis set definition should actually be a dictionary
        if 'basis' in params.kwords[obj_name].keys():
            # construct basis dictionary
            if isinstance(obj.basis, (list, np.ndarray)):
                basis_str = obj.basis.copy()             
                obj.basis = {basis_str[2*i] : basis_str[2*i+1] for 
                                i in range(int(len(basis_str)/2.))}

            # if basis just a single string, apply to all atoms
            else:
                # get the list of unique atoms
                atms       = set(obj.asym)
                basis_str = obj.basis
                obj.basis = {atm : basis_str for atm in atms}

        # the basis set definition should actually be a dictionary
        if 'ri_basis' in params.kwords[obj_name].keys():
            # construct basis dictionary, but unlike obj.basis, 
            # ri_basis is allowed to be None
            if obj.ri_basis is not None:
                
                if isinstance(obj.ri_basis, (list, np.ndarray)):
                    basis_str = obj.ri_basis.copy()
                    obj.ri_basis = {basis_str[2*i] : basis_str[2*i+1] 
                             for i in range(int(len(basis_str)/2.))}

                # if basis just a single string, apply to all atoms
                else:
                    # get the list of unique atoms
                    atms       = set(obj.asym)
                    basis_str = obj.ri_basis
                    obj.ri_basis = {atm : basis_str for atm in atms}

        # If RAS spaces have not been specified, then set autoras = True
        if hasattr(obj, 'autoras'):
            obj.autoras = True
            for key in ['ras1', 'ras2', 'ras3']:
                if hasattr(obj, key):
                    if len(getattr(obj, key)) != 0:
                        obj.autoras = False
                    
        # init/final_states and i/fstate_array need to be lists, also:
        # internal state ordering is 0->n-1, vs. 1->n for input
        if type(obj).__name__ == 'Transition' or type(obj).__name__ == 'Spinorbit':
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

        # save section labels -- make sure all our unique
        section_lbls.append(obj.label)

    for obj_indx in range(len(section_lbls)):
        # save section labels -- make sure all our unique
        lbl_tst = run_list[obj_indx].label
        indices = [index for index,lbl in enumerate(section_lbls) 
                                                  if lbl == lbl_tst]
        # if multiple identical sections -- rename
        if len(indices) > 1:
            for suffix in range(len(indices)):
                run_list[indices[suffix]].label += str(suffix+1)

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

