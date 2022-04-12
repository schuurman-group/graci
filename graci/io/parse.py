"""Module for performing file operations"""

import sys
import re as re
import numpy as np 
import graci.utils.constants as constants
import graci.core.params as params
import graci.io.output as output

import graci.core.molecule as molecule
import graci.core.scf as scf
import graci.methods.dftmrci as dftmrci
import graci.methods.dftmrenpt2 as dftmrenpt2
import graci.interaction.transition as transition
import graci.interaction.spinorbit as spinorbit
import graci.interaction.overlap as overlap

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
    if isinstance(value, (np.ndarray, list)):
        if isinstance(value, np.ndarray):
            val_list = value.flatten().tolist()
        else:
            val_list = value 

        if all([isinstance(elem, keyword_type) for elem in val_list]):
            correct = True
    else:
        if isinstance(value, keyword_type):
            correct = True

    return correct

#
def parse_value(valstr):
    """Returns a value converted to the appropriate type and shape.

    By default, spaces and newlines will be treated as delimiters.
    """

    # split any braces or ':' symbols
    split_line = re.split('(:)|(\[)|(\])|\n', valstr)
    # remove instances of 'None'
    try:
        while True:
            split_line.remove(None)
    except ValueError:
        pass

    # val_list is a list of the space/newline delimited input
    val_list = []
    for line in split_line:
        val_list.extend(line.split())

    # remove empty strings
    try:
        while True:
            val_list.remove('')
    except ValueError:
        pass

    # remove commas
    try:
        while True:
           val_list.remove(',')
    except ValueError:
        pass

    # step through and replace all X:Y with range(X,Y+1)
    while ':' in val_list:
        indx = val_list.index(':')
        try:
            start = int(val_list[indx-1])
            end   = int(val_list[indx+1])+1
        except:
            sys.exit('cannot convert range values: '+str(valstr))
        num_list = list(range(start, end))
        str_list = [str(num) for num in num_list]
        new_list = val_list[:indx-1] + str_list + val_list[indx+2:]
        val_list = new_list

    # if this is a scalar: just convert the number
    if len(val_list) == 1:
        return convert_value(val_list[0])

    # if this is vector, need to determine if 1D or 2D
    # we will accept the following syntax for 1D arrays:
    #  1. kword = [X Y Z]
    #  2. kword = [X, Y, Z]
    #
    # However, we will be more strict re: 2D arrays. Need
    # to see those square braces to denote change in dim
    #  1. kword = [X Y Z] [X Y Z]
    #  2. kword = [X, Y, Z], [X, Y, Z]
    
    vec_str = []
    # scan for opening brace
    while '[' in val_list:
        start = val_list.index('[')+1
        # try to close this brace, else just take to the end
        try:
            end = val_list.index(']')
        except:
            sys.exit('missing a closing brace: '+str(valstr))
        # append the entries between the braces a new vector 
        vec_str.append(list(val_list[start:end]))
        if end == len(val_list)-1:
            break
        val_list = val_list[end+1:]

    # if just a single element of this 2D list, convert to vector
    if len(vec_str) == 1:
        return convert_array(vec_str[0])
    else:
        return convert_array(vec_str)

#
def check_input(run_list):
    """Checks on the user-supplied input"""
   
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

        # Spinorbit _has_ to enter states as a vector, just shift state
        # indices
        if type(obj).__name__ == 'Spinorbit':
            # shift statesby 1 to internal/C ordering
            obj.couple_states -= 1

        # init/final_states and i/fstate_array need to be lists, also:
        # internal state ordering is 0->n-1, vs. 1->n for input
        if (type(obj).__name__ == 'Transition' or 
            type(obj).__name__ == 'Sotransition') or
            type(obj).__name__ == 'Overlap')):
            if obj.init_states is not None:
                if not isinstance(obj.init_states, (list, np.ndarray)):
                    obj.init_states = np.array([obj.init_states])
            if obj.final_states is not None:
                if not isinstance(obj.final_states, (list, np.ndarray)):
                    obj.final_states = np.array([obj.final_states])
            # shift statesby 1 to internal/C ordering
            obj.init_states  -= 1
            obj.final_states -= 1

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

