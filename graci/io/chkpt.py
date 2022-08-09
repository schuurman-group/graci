"""
Module for facilitating passing of strings between
Python and compiled dlls
"""
import os as os
import sys as sys
import base64 as base64
import h5py as h5py
import numpy as np
import struct as struct
import json as json
import graci.core.params as params
import graci.core.molecule as molecule
import graci.core.scf as scf
import graci.core.bitciwfn as bitciwfn
import graci.tools.parameterize as parameterize
import graci.methods.dftmrci as dftmrci
import graci.methods.dftmrci2 as dftmrci2
import graci.properties.moments as moments
import graci.interaction.interaction as spininfo # this makes me uncomfortable...
import graci.interaction.transition as transition
import graci.interaction.spinorbit as spinorbit
import graci.interaction.overlap as overlap
import graci.utils.rydano as rydano
import graci.io.output as output

# top level computation objects: each will have own top-level group
comp_objs = tuple([getattr(globals()[obj_class.lower()], obj_class)
                               for obj_class in params.valid_objs])

# support level computation objects: each will have own sub-group
sub_objs = tuple([getattr(globals()[obj_class.lower()], obj_class)
                               for obj_class in params.support_objs])

#
def contents(file_name=None, file_handle=None, grp_name=None):
    """return the group names present in the checkpoint file"""

    if file_handle is not None:
        chkpt = file_handle
    else:
        f_name = file_name
        if f_name is None:
            f_name = output.file_names['chkpt_file']
        try:
            chkpt = h5py.File(f_name, 'r', libver='latest')
        except:
            return None

    # if grp_name stub (i.e. class type) is not given,
    # return all top-level class objects in checkpoint file
    if grp_name is None:
        grp_lst = list(chkpt.keys())
    # else, only return objects associated with a specific
    # object
    elif grp_name in list(chkpt.keys()):
        grp_lst = list(chkpt[grp_name].keys())
    # if this is a dataset, just return the data set name
    elif grp_name in chkpt:
        grp_lst = [grp_name]
    # if neither group nor dataset, return none
    else:
        grp_lst = None

    # if we were just passed the file_name, close the chkpt
    # when done
    if file_handle is None:
        chkpt.close()

    return grp_lst

#
def write(graci_obj, file_name=None, file_handle=None, 
                                 grp_name=None, label_suffix=None):
    """write the graci_obj to the checkpoint file"""
    global sub_objs, comp_objs

    # if new_label is not none, update the object label
    l_suffix=''
    if label_suffix is not None:
        graci_obj.label += label_suffix
        l_suffix = label_suffix

    # by default, the group name will by the class type + object label
    if grp_name is None:
        grp_name = str(type(graci_obj).__name__)+'.'+str(graci_obj.label)

    # this routine is called recursively. If the file handle
    # is given, use that to access checkpoint file, else try
    # to open file given by file_name (or barring that, the chkpt_file
    # variable in the output module
    if file_handle is None:
        f_name = file_name
        if f_name is None:
            f_name = output.file_names['chkpt_file']
        chkpt = h5py.File(f_name, 'a', libver='latest')
    else:
        chkpt = file_handle

    # if grp doesn't exist, create it 
    if grp_name not in chkpt:
        obj_grp = chkpt.create_group(grp_name)
    else:
        obj_grp = chkpt[grp_name]

    # list of all objects in the class
    objs = graci_obj.__dict__

    for name, obj in objs.items():

        # if the object is a graci computation class object, save
        # the name of the object as a string -- this
        # is preferable to storing the same object
        # repeatedly
        if isinstance(obj, comp_objs):
            obj_grp.attrs[name] = dumps(link_name(obj, 
                                        suffix=l_suffix))
            continue

        # if this is a class object, create a subgroup and write it there
        if isinstance(obj, sub_objs):
            dset_name = link_name(obj, suffix=str(name))
            obj_grp.attrs[name] = dumps(dset_name)
            write(obj, file_handle=chkpt, grp_name=grp_name+'/'+dset_name)
            continue

        # if object is a numpy array, save
        # as a distinct dataset
        if isinstance(obj, np.ndarray):
            dset_name = link_name(obj, suffix=str(name))
            obj_grp.attrs[name] = dumps(dset_name)
            write_dataset(chkpt, grp_name+'/'+dset_name, obj)
            continue

        # if this is a pyscf object, simply save as 'None', this will
        # be regenerated if need be
        if str(type(obj).__name__) == 'Mole':
            obj_grp.attrs[name] = dumps(None)
            continue

        # if this is a list, replace complex objects with links
        if isinstance(obj, list):
            obj_write = [None]*len(obj) 
            for indx in range(len(obj)):
                lname = str(name)+'.'+str(indx)
                value = obj[indx]
                if isinstance(value, comp_objs):
                    dset_name = link_name(value, suffix=l_suffix)
                    obj_write[indx] = dset_name
                elif isinstance(value, sub_objs):
                    dset_name = link_name(value, suffix=lname)
                    write(obj[indx], file_handle=chkpt,
                                   grp_name=grp_name+'/'+dset_name)
                    obj_write[indx] = dset_name
                elif isinstance(value, np.ndarray):
                    dset_name = link_name(value, suffix=lname)
                    write_dataset(chkpt, grp_name+'/'+dset_name, value)
                    obj_write[indx] = dset_name
                else:
                    obj_write[indx] = obj[indx]

            # the list has been rendered safe for writing
            try:
                attr_str = dumps(obj_write)
                obj_grp.attrs[name] = attr_str
            except:
                sys.exit('Could not write list: '+str(name))
            continue
        
        # if this is a dictionary, check if any
        # values are numpy arrays or class objects:
        # things that we don't want to store as attributes
        if isinstance(obj, dict):
            obj_write = obj.copy()
            for key,value in obj_write.items():
                lname = str(name)+'.'+str(key)
                if isinstance(value, comp_objs):
                    dset_name = link_name(value, suffix=l_suffix)
                    obj_write[key] = dset_name
                elif isinstance(value, sub_objs):
                    dset_name = link_name(value, suffix=lname)
                    write(value, file_handle=chkpt, 
                                 grp_name=grp_name+'/'+dset_name)
                    obj_write[key] = dset_name
                elif isinstance(value, np.ndarray):
                    dset_name = link_name(value, suffix=lname)
                    write_dataset(chkpt, grp_name+'/'+dset_name, value)
                    obj_write[key] = dset_name

            # the dictionary has been rendered safe for writing
            try:
                obj_grp.attrs[name] = dumps(obj_write)
            except:
                sys.exit('Could not write dictionary: '+str(name))
            continue

        # everything else should be handled by json
        try:
            attr_str = dumps(obj)
            obj_grp.attrs[name] = attr_str
        except TypeError:
            sys.exit('Writing attribute = '+str(name) + 
                     ' to chkpt file failed')

    # close the chkpt file if the file_name was passed
    if file_handle is None:
        chkpt.close()

    return

#
def read(grp_name, build_subobj=True, make_mol=True, 
         file_name=None, file_handle=None):
    """return a GRaCI object, with name 'obj_name', after parsing 
       checkpoint file. If build_subobj is True, also build GRaCI 
       objects contained in 'obj_name'. If the file_handle is passed
       directly, use that instead of opening file. If make_mol is 
       'True', call run() routine in GRaCI mol object to regenerate
       the pyscf mol object (which is not written to chkpt file)"""
    global sub_objs

    # if file doesn't exist, just return None
    if file_handle is None:
        f_name = file_name
        if f_name is None:
            f_name = output.file_names['chkpt_file']
        try:
            chkpt = h5py.File(f_name, 'r', libver='latest')
        except:
            return None
    else:
        chkpt = file_handle

    # if requested item doesn't exist, just return None
    if grp_name not in chkpt:
        return None
    else:
        class_grp = chkpt[grp_name]

    # group names must have the form ClassName.label
    obj_class = os.path.basename(grp_name).split('.')[0]

    # create the object -- assumes module has already been 
    # imported above
    Chkpt_obj = getattr(globals()[obj_class.lower()], obj_class)()

    # first set the class variables that are saved as attributes
    for name, attr in class_grp.attrs.items():

        if not hasattr(Chkpt_obj, name):
            # going to fail silently on this one for now, since
            # users can set attributes that are not class variables
            # using glabel
            continue
            #raise AttributeError(
            #'Class '+str(obj_class)+' has no attribute '+str(name))

        # run the parser through JSON
        try:
            setattr(Chkpt_obj, name, loads(attr))
        except:
            sys.exit('Cannot set class object, name='+str(name))

        # Scroll through lists looking for links 
        if isinstance(getattr(Chkpt_obj, name), list):
            list_obj = getattr(Chkpt_obj, name)
            for indx in range(len(list_obj)):
                # this is a link to a top level object
                value = list_obj[indx]
                gobj, sub_name = data_name(grp_name, value)
                if sub_name is not None:
                    if gobj:
                        obj = read(sub_name,
                               build_subobj = build_subobj,
                               make_mol = make_mol,
                               file_handle = chkpt)
                        list_obj[indx] = obj
                    else:
                        list_obj[indx] = np.array(chkpt[sub_name])
            continue

        # if this is a dictionary, run through the keys and load
        # any non-json serializable objects
        if isinstance(getattr(Chkpt_obj, name), dict):
            dict_obj = getattr(Chkpt_obj, name)
            for key, value in dict_obj.items():
                gobj, sub_name = data_name(grp_name, value)
                if sub_name is not None:
                    if gobj:
                        obj = read(sub_name,
                               build_subobj = build_subobj,
                               make_mol = make_mol,
                               file_handle = chkpt)
                        dict_obj[key] = obj
                    else:
                        dict_obj[key] = np.array(chkpt[sub_name])
            continue

        # check the for links to complex objects
        name_str = str(getattr(Chkpt_obj, name))

        # attribute holds a link to a complex object
        gobj, sub_name = data_name(grp_name, name_str)

        if sub_name is not None:
            if gobj:
                obj = read(sub_name,
                           build_subobj = build_subobj,
                           make_mol = make_mol,
                           file_handle = chkpt)
            else:
                obj = np.array(chkpt[sub_name])
            setattr(Chkpt_obj, name, obj)
            continue

    # don't close the file if it was called by handle -- if it's called
    # recursively, may trip up higher levels
    if file_handle is None:
        chkpt.close()

    # this is a bit of hack -- generate Molecule class if need be
    if make_mol and obj_class == 'Molecule':
        Chkpt_obj.run()

    return Chkpt_obj

#----------------------------------------------------------------------------

#
def write_dataset(chkpt_handle, dset_name, dset):
    """
    Write a dataset to the grp_handle

    Arguments:
    grp_handle:   handle to checkpoint file group
    dset_name:    name of the dataset (a string)
    dset:         the dataset, typically, a numpy array
    """

    # if dataset exists, delete it
    if dset_name in chkpt_handle:
        del chkpt_handle[dset_name]

    try:
        chkpt_handle.create_dataset(dset_name, data=dset, dtype=dset.dtype,
                                    compression="gzip", compression_opts=9)
    # if non-native HDf5 type, DIE
    except TypeError:
        sys.exit('Could not create dataset: '+str(dset_name))

    return

#
def read_dataset(chkpt_handle, dset_name):
    """
    Read a dataset from chkpt_handle

    Arguments;
    chkpt_handle:   file handle for the chkpt file
    dset_name:      name of the dataset (str)

    Returns:
    dset:           a numpy array with same shape as dset_name
    """

    if dset_name not in chkpt_handle:
        sys.exit('Could not locate '+str(dset_name)+' in chkpt file')

    try:
        dset = np.ndarray(chkpt_handle[dset_name])
    except TypeError:
        sys.exit('Could not read '+str(dset_name)+' as numpy array')

    return dset

# 
def write_attribute(chkpt_handle, obj_name, name, value):
    """
    Write an attribute to a dataset or group

    Arguments:
    chkpt_handle:   file handle for the chkpt file
    obj_name:       group or dataset name
    name:           name of attribute
    value:           attribute to write. attribute will be run through
                    overloaded json method dumps()
    """

    if obj_name not in chkpt_handle:
        sys.exit('Could not locate '+str(obj_name)+' in chkpt file')

    try:
        chkpt_handle[obj_name].attrs[name] = dumps(value)
    except:
        sys.exit('Could not write attribute: '+str(name))

    return

#
def read_attribute(chkpt_handle, obj_name, name):
    """
    Read an attribute from a dataset or group

    Arguments:
    chkpt_handle:   file handle for the chkpt file
    obj_name:       group or dataset name
    name:           name of attribute

    Returns
    value:           attribute to write. attribute will be run through
                     overloaded json method dumps()
    """

    if obj_name not in chkpt_handle:
        sys.exit('Could not locate '+str(obj_name)+' in chkpt file')

    try:
        value = loads(chkpt_handle[obj_name].attrs[name])
    except:
        value = None

    return value

#
def data_name(grp, link_name):
    """
    Generate the group name given the link string
    """

    lstr      = str(link_name)
    gobj      = False
    data_name = None
    if '.OBJ' in lstr:
        gobj = True
        if '.OBJC' in lstr:
            data_name = grp + '/'+lstr
        else:
            data_name = lstr.replace('.OBJ','')
    elif 'NUMPY' in lstr:
        data_name = grp + '/' + lstr

    return gobj, data_name

#
def link_name(obj, suffix=''):
    """
    Generate a string file name that points to an object in checkpoint file
    """
    global comp_objs, sub_objs

    if isinstance(obj, comp_objs):
        lstr = '.' + str(obj.label) + str(suffix) + '.OBJ'
        return str(type(obj).__name__) + lstr 

    if isinstance(obj, sub_objs):
        return str(type(obj).__name__)+'.'+str(suffix)+'.OBJC'

    if isinstance(obj, np.ndarray):
        return 'NUMPY.'+str(suffix)

    return None

#
class GraciEncoder(json.JSONEncoder):
    def default(self, obj):
        """
        if input object is a ndarray it will be converted into a dict holding dtype,
        shape and the data base64 encoded
        """
        if isinstance(obj, bytes):
            return obj.decode()
        # json doesn't like numpy data types
        if isinstance(obj, np.int64):
            return int(obj)
        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, obj)

# Overload dump/load to default use this behavior.
def dumps(*args, **kwargs):
    kwargs.setdefault('cls', GraciEncoder)
    return json.dumps(*args, **kwargs)

def loads(*args, **kwargs):
    return json.loads(*args, **kwargs)

#-----------------------------------------------------------
# targeted for removal
def read_wfn(obj_name, wfn_indx, cutoff, file_name=None):
    """reads the wfn object from the obj_name group. Returns
       a numpy array with the determinant coefficients and a 
       numpy array with the orbital occupations"""

    # using 64-bit integers to store det occupations
    STR_LEN = 64

    f_name = file_name
    if f_name is None:
        f_name = output.file_names['chkpt_file']
    try:
        chkpt = h5py.File(f_name, 'r', libver='latest')
    except:
        return None

    grp_name = '/'+str(obj_name)

    cf_name  = 'wfn_cf'+str(wfn_indx)
    det_name = 'wfn_dets'

    grp_dsets = list(chkpt[grp_name].keys())

    if cf_name not in grp_dsets or det_name not in grp_dsets:
        print(cf_name + ' and/or ' + det_name +
                            ' not found in group '+grp_name)
        return None

    cf_name  = grp_name + '/' + cf_name
    det_name = grp_name + '/' + det_name
    (ndet_cf,   dum)  = chkpt[cf_name].shape
    (ndet_wfn, n_tot) = chkpt[det_name].shape

    if ndet_wfn != ndet_cf:
        print("WARNING: determinant and coefficient lists are not"+
                " the same length: "+str(ndet_wfn)+"!="+str(ndet_cf))

    # some preliminaries so we can allocated the output coefficient
    # and det list arrays
    ndet  = min(ndet_wfn, ndet_cf)
    n_int = int(n_tot / 2)
    nmo   = int(chkpt[det_name].attrs['nmo'])

    cfs_raw  = np.array(chkpt[cf_name][:, 0], dtype=float)
    det_ints = chkpt[det_name][...]

    # next, we sort the coefficients, and only convert those dets
    # corresponding to contributions > the cutoff.
    # we flip array to get values in descending order
    ordr = np.flip(np.argsort(np.absolute(cfs_raw)))

    # now we step through the cf array 
    cfs  = np.zeros((ndet,), dtype=float)
    dets = np.zeros((nmo, ndet), dtype=np.short)

    for idet in range(ndet):

        indx = ordr[idet]

        if abs(cfs_raw[indx]) < cutoff:
            break

        cfs[idet] = cfs_raw[indx]

        alpha = []
        beta  = []

        # use numpy intrinsic binary_repr -- if width specified,
        # returns two's complement value
        for i_int in range(n_int):
            alpha  += list(np.binary_repr(det_ints[indx, i_int],
                                           width=STR_LEN))[::-1]
            beta   += list(np.binary_repr(det_ints[indx, n_int+i_int],
                                           width=STR_LEN))[::-1]

        # only take first nmo entries -- the rest are padding
        dabs = [int(alpha[i]) + int(beta[i]) for i in range(nmo)]
        rdet = np.array(dabs, dtype=np.short)
        for i in range(nmo):
            if int(alpha[i]) == 0:
                rdet[i] = -dabs[i]

        if idet == 0:
            ecnt = sum(dabs)
        elif sum(dabs) != ecnt:
            print('Wrong number of electrons in determinant. '+
                  '\nidet=' + str(idet) +
                  '\nints=' + str(det_ints[indx, :]) +
                  '\n alpha = '+str(alpha) +
                  '\n beta = '+str(beta) +
                  '\nrdet='+str(rdet))
            sys.exit(1)

        dets[:, idet] = rdet

    return cfs[:idet], dets[:,:idet]
