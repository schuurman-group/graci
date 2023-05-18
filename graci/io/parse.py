"""Module for performing file operations"""

import sys
import re as re
import numpy as np 
import h5py as h5py
import graci.utils.constants as constants
import graci.core.params as params
import graci.core.hamiltonians as hamiltonians
import graci.io.output as output
import graci.utils.basis as basis
import graci.utils.rydano as rydano
import graci.core.molecule as molecule
import graci.core.scf as scf
import graci.tools.parameterize as parameterize
import graci.methods.dftmrci as dftmrci
import graci.methods.dftmrci2 as dftmrci2
import graci.interaction.transition as transition
import graci.interaction.spinorbit as spinorbit
import graci.interaction.overlap as overlap
import graci.interaction.dyson as dyson

from pyscf import gto

#valid input sections

#
def parse_input():
    """parse the graci input file"""

    # Read input file. Small enough to gulp the whole thing
    with open(output.file_names['input_file'], 'r') as infile:
        input_file = infile.readlines()

    # check for invalid sections
    check_sections(input_file)
        
    run_list = []
    for obj_name in params.valid_objs:

        # parse all section keywords
        run_list.extend(parse_section(obj_name, 
                                      input_file))

    # check for multiple-geometry molecule sections
    # and create all required replicate class objects
    # if any are found
    run_list = replicate_sections(run_list)
        
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
        if '$'+mod_name+' ' in input_file[iline].lower():

            iline += 1

            # section exists, create class object
            sec_obj = getattr(globals()[mod_name], class_name)()

            # if this is a molecule section, look for an xyz geometry
            # input format
            if class_name == 'Molecule':
                cart = []
                atom = []

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
                        expected = params.kwords[class_name][kword]
                        val      = parse_value(value, expected)

                        if correct_type(val, expected):
                            setattr(sec_obj, kword, val)
                        else:
                            sys.exit(kword+" is wrong type, expected: "
                                  +str(expected))
                    else:
                        sys.exit('Invalid keyword found in a $' \
                                 +mod_name+' section'': '+kword)
                        
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
            val_list = []
            for val in value:
                if isinstance(val, np.ndarray):
                    val_list.extend(val.flatten().tolist())
                else:
                    val_list.extend(val)

        if all([isinstance(elem, keyword_type) for elem in val_list]):
            correct = True
    else:
        if isinstance(value, keyword_type):
            correct = True

    return correct

#
def parse_value(valstr, val_type):
    """Returns a value converted to the appropriate type and shape.
    By default, spaces and newlines will be treated as delimiters.
    All keywords are converted to lower case.
    """

    # split any braces or ':' symbols
    split_line = re.split('(:)|(\[)|(\])|\n', valstr.lower())
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

    # if type is string, don't remove characters or interpret range
    # values
    if val_type is not str:

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
    if len(vec_str)==1: 
        return convert_array(vec_str[0])
    else:
        return convert_array(vec_str)
    
#
def check_sections(input_file):
    """checks for the existense of invalid sections in the input file"""

    # get the list of section headers
    header_list  = [line.strip() for line in input_file if 'section' in line]
    
    # make sure that each section corresponds to valid class object
    valid_objs = [obj.lower() for obj in params.valid_objs]
    for header in header_list:
        section_name = header[1:header.index('section')].lower().strip()
        if section_name not in valid_objs:
            sys.exit('Invalid section found: $'+section_name)

    return
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

        # Shift ref_state by -1 to get the Python indexing
        if 'ref_state' in params.kwords[obj_name].keys():
            if obj.ref_state != -1:
                obj.ref_state -= 1
                
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

        # Catch any Hamiltonian aliases
        if hasattr(obj, 'hamiltonian'):
            try:
                obj.hamiltonian = hamiltonians.aliases[obj.hamiltonian.lower()]
            except:
                pass

        # Spinorbit _has_ to enter states as a vector, just shift state
        # indices
        if type(obj).__name__ == 'Spinorbit':
            # if a single range of states is provided, rehape array to
            # to a 2D array with shape[0] = 1
            if isinstance(obj.couple_states, np.ndarray):
                obj.couple_states = [obj.couple_states]

            # couple_states is now a list of numpy arrays, shift states
            # by 1 to internal ordering
            for i in range(len(obj.couple_states)):
                obj.couple_states[i] -= 1
                 
        # init/final_states and i/fstate_array need to be lists, also:
        # internal state ordering is 0->n-1, vs. 1->n for input
        if type(obj).__name__ == 'Transition' \
           or type(obj).__name__ == 'Dyson':
            if obj.init_states is not None:
                if not isinstance(obj.init_states, (list, np.ndarray)):
                    obj.init_states = np.array([obj.init_states])
            if obj.final_states is not None:
                if not isinstance(obj.final_states, (list, np.ndarray)):
                    obj.final_states = np.array([obj.final_states])
            # shift statesby 1 to internal/C ordering
            obj.init_states  -= 1
            obj.final_states -= 1            

        if type(obj).__name__ == 'Overlap':
            if obj.bra_states is not None:
                if not isinstance(obj.bra_states, (list, np.ndarray)):
                    obj.bra_states = np.array([obj.bra_states])
            if obj.ket_states is not None:
                if not isinstance(obj.ket_states, (list, np.ndarray)):
                    obj.ket_states = np.array([obj.ket_states])
            # shift statesby 1 to internal/C ordering
            obj.bra_states -= 1
            obj.ket_states -= 1
            
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
def convert_array(arg_list):
    """Converts a list of strings to an array of ints, floats or strings."""

    # if this is a nested list, iterate over elements, else, convert to
    # nested lis
    if type(arg_list[0]) != list:
        conv_list = [arg_list]
    else:  
        conv_list = arg_list
    
    new_list = []
    for arg in conv_list:

        try:
            arr = np.array(arg).astype(int)
            new_list.append(arr)
            continue 
        except ValueError:
            pass

        try:
            arr = np.array(arg).astype(float)
            new_list.append(arr)
            continue
        except ValueError:
            pass

        arr = np.array(arg, dtype=h5py.string_dtype(encoding='utf-8'))
        new_list.append(arr)

    if type(arg_list[0]) != list:
        return new_list[0]
    else:
        return new_list
    
#
def replicate_sections(run_list):
    """
    checks for the existence of multi-geometry molecule objects
    if any are found, then replicate molecule objects are creates for
    each geometry along with any associated scf, dftmrci, si, etc.
    objects
    """

    # initialise the new list of class objects to run
    new_run_list = []

    # list of all class objects of the various different
    # types
    misc_objs   = params.valid_objs
    for g_obj in params.ci_objs + params.postci_objs + params.si_objs:
        misc_objs.remove(g_obj)
    for g_obj in ['Molecule','Scf','Parameterize']:
        misc_objs.remove(g_obj)

    mol_objs    = [obj for obj in run_list
                     if type(obj).__name__ == 'Molecule']
    scf_objs    = [obj for obj in run_list
                     if type(obj).__name__ == 'Scf']
    ci_objs     = [obj for obj in run_list
                     if type(obj).__name__ in params.ci_objs]
    postci_objs = [obj for obj in run_list
                     if type(obj).__name__ in params.postci_objs]
    si_objs     = [obj for obj in run_list
                     if type(obj).__name__ in params.si_objs]
    hparam_objs = [obj for obj in run_list
                     if type(obj).__name__ == 'Parameterize']
    graci_objs  = [obj for obj in run_list
                     if type(obj).__name__ in misc_objs]

    # check for multi-geometry xyz files
    for mol in mol_objs:

        # list of scf, ci, postci and si objects
        # associated with this molecule object
        scf_list   = [obj for obj in scf_objs
                      if obj.mol_label == mol.label]
        scf_labels = [obj.label for obj in scf_list]        
        ci_list   = [obj for obj in ci_objs
                     if obj.scf_label in scf_labels]
        ci_labels = [obj.label for obj in ci_list]
        postci_list = [obj for obj in postci_objs
                       if any(lbl in obj.couple_groups
                              for lbl in ci_labels)]
        postci_labels = [obj.label for obj in postci_list]
        all_ci_labels = list(set(ci_labels).union(set(postci_labels)))
        try:
            # overlap class objects have 'bra' and 'ket' labels
            si_list = [obj for obj in si_objs
                       if obj.bra_label in all_ci_labels
                       or obj.ket_label in all_ci_labels]
        except:
            # everything else has 'initial' and 'final' labels
            si_list = [obj for obj in si_objs
                       if obj.init_label in all_ci_labels
                       or obj.final_label in all_ci_labels]
            
            
        if mol.multi_geom:
            # Create replicate objects for all geometries
            
            # read the complete set of Cartestian coordinates
            # from the xyz file
            coords = parse_all_geoms(mol)

            # append the new run list with replicates of the
            # class objects corresponding to this molecule
            # object
            for i in range(coords.shape[0]):

                # molecule object
                new_mol       = mol.copy()
                new_mol.label = mol.label+str(i+1)
                new_mol.crds  = 1. * coords[i]
                new_run_list.append(new_mol)
                
                # scf object(s)
                for scf in scf_list:
                    new_scf           = scf.copy()
                    new_scf.label     = scf.label+str(i+1)
                    new_scf.mol_label = scf.mol_label+str(i+1)
                    if i > 0:
                        new_scf.guess_label = scf.label+str(i)
                    new_run_list.append(new_scf)
                    
                # ci object(s)
                for ci in ci_list:
                    new_ci           = ci.copy()
                    new_ci.label     = ci.label+str(i+1)
                    new_ci.scf_label = ci.scf_label+str(i+1)
                    try:
                        if (new_ci.diabatic):
                            new_ci.save_wf = True
                            if i > 0:
                                new_ci.guess_label = ci.label+str(i)
                            else:
                                new_ci.diabatic = False
                    except:
                        pass
                    new_run_list.append(new_ci)

                # post-ci objects
                for postci in postci_list:
                    new_postci           = postci.copy()
                    new_postci.label     = postci.label+str(i+1)
                    new_postci.couple_groups = [lbl+str(i+1)
                                                for lbl in postci.couple_groups]
                    new_run_list.append(new_postci)
                    
                # si object(s)
                for si in si_list:
                    new_si             = si.copy()
                    new_si.label       = ci.label+str(i+1)
                    new_si.init_label  = si.init_label+str(i+1)
                    new_si.final_label = si.final_label+str(i+1)
                    new_run_list.append(new_si)
                  
        else:
            # add the single-geometry objects to the list
            new_run_list.append(mol)
            for scf in scf_list:
                new_run_list.append(scf)
            for ci in ci_list:
                new_run_list.append(ci)
            for postci in postci_list:
                new_run_list.append(postci)
            for si in si_list:
                new_run_list.append(si)
   
    # assume graci objs are geometry indepenent
    for gobj in graci_objs:
        new_run_list.append(gobj)
 
    # might want to re-think this a bit...
    for hparam in hparam_objs:
        new_run_list.append(hparam)                

    return new_run_list

#
def parse_all_geoms(mol):
    """
    given a molecule object, reads in all geometries in mol.xyz_file file
    """

    # parse the xyz file
    with open(mol.xyz_file, 'r') as xyzfile:
        xyz = xyzfile.readlines()
        
    # remove the leading no. atom and blank lines
    xyz_clean = [string.split() for string in xyz
                if string.split() != []
                 and len(string.split()) != 1]
    n_atm  = len(mol.crds)
    n_geom = int(len(xyz_clean) / n_atm)
    
    # get the array of nuclear geometries
    coords = np.array([float(xx)
                       for x in xyz_clean
                       for xx in x[1:]]).reshape(n_geom,n_atm,3)

    return coords
