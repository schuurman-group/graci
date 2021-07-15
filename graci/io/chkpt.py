"""
Module for facilitating passing of strings between
Python and compiled dlls
"""
import sys as sys
import h5py as h5py
import numpy as np
import struct as struct

import graci.core.params as params
import graci.core.molecule as molecule
import graci.core.parameterize as parameterize
import graci.methods.scf as scf
import graci.methods.dftmrci as dftmrci
import graci.methods.dftmrenpt2 as dftmrenpt2
import graci.methods.dftcis as dftcis
import graci.properties.transition as transition
import graci.properties.spinorbit as spinorbit
import graci.utils.timing as timing

#
def contents(file_name, grp_name=None):
    """return the group names present in the checkpoint file"""

    try:
        chkpt = h5py.File(file_name, 'r', libver='latest')
    except:
        return None

    if grp_name is None:
        grp_lst = list(chkpt.keys())
    elif grp_name in list(chkpt.keys()):
        grp_lst = list(chkpt[grp_name].keys())
    else:
        grp_lst = None

    chkpt.close()

    return grp_lst

#
def write(file_name, graci_obj):
    """write the graci_obj to the checkpoint file"""

    # get the data sets and attributes
    dsets, attrs = pack(graci_obj)

    if len(dsets)+len(attrs) == 0:
        return

    # the group name will by the class type + object label
    grp_name = str(type(graci_obj).__name__)+'.'+str(graci_obj.label)

    chkpt = h5py.File(file_name, 'a', libver='latest')

    # create the group if doesn't already exist
    if grp_name not in list(chkpt.keys()):
        obj_grp = chkpt.create_group(grp_name)
    else:
        obj_grp = chkpt[grp_name]

    # first write the datasets
    for name, dset in dsets.items():

        # if the dataset already exists, resize if necessary
        if name in obj_grp.keys():
            # get shape and resize
            if obj_grp[name].shape != dset.shape:
                try:
                    obj_grp[name].resize(dset.shape)
                    continue
                except:
                    del obj_grp[name]

        # if we're here, need to create a new dataset
        try:
            obj_grp.create_dataset(name, data=dset, dtype=dset.dtype, 
                                compression="gzip", compression_opts=9)
            # if non-native HDf5 type, just convert to a string
        except TypeError:
            dset_str = str(dset)
            obj_grp.create_dataset(name, data=dset_str, 
                                compression="gzip", compression_opts=9)

    # next write the attributes
    for name, attr in attrs.items():
        try:
            obj_grp.attrs[name] = attr
        # if attribute is non-native HDF5 type, convert to string
        except TypeError:
            attr_str = str(attr)
            obj_grp.attrs[name] = attr_str
            
    chkpt.close()
    return

#
def read(file_name, obj_name, build_subobj=True, make_mol=True,
         file_handle=None):
    """return a GRaCI object, with name 'obj_name', after parsing 
       checkpoint file. If build_subobj is True, also build GRaCI 
       objects contained in 'obj_name'. If the file_handle is passed
       directly, use that instead of opening file. If make_mol is 
       'True', call run() routine in GRaCI mol object to regenerate
       the pyscf mol object (which is not written to chkpt file"""

    # if file doesn't exist, just return None
    if file_handle is not None:
        chkpt = file_handle
    else:
        try:
            chkpt = h5py.File(file_name, 'r', libver='latest')
        except:
            return None

    # if requested item doesn't exist, just return None
    if obj_name not in list(chkpt.keys()):
        return None
    else:
        class_grp = chkpt[obj_name]

    # group names must have the form ClassName.label
    obj_class = obj_name.split('.')[0]

    # create the object -- assumes module has already been 
    # imported above
    Chkpt_obj = getattr(globals()[obj_class.lower()], obj_class)()

    # populate with the save datasets (we know these 
    # are N-rank numpy arrays)
    for name, dset in class_grp.items():

        # only load if class has attribute. 
        if hasattr(Chkpt_obj, name):
            nparray = np.empty(dset.shape, dset.dtype)
            nparray = dset[...]
            setattr(Chkpt_obj, name, nparray)

    # populate with the saved attributes (these can
    # be more complicated)
    for name, attr in class_grp.attrs.items():
   
        if not hasattr(Chkpt_obj, name):
            raise AttributeError(
            'Class '+str(obj_class)+' has no attribute '+str(name))

        # if this is a graci object, read this object from 
        # this file
        if 'GRACI_OBJ.' in str(attr) and build_subobj:
            subgrp_name = str(attr).replace('GRACI_OBJ.', '')
            subobj = read(file_handle, subgrp_name, 
                          build_subobj = build_subobj, 
                          make_mol = make_mol, 
                          file_handle = chkpt)
            setattr(Chkpt_obj, name, subobj)
            continue

        # if this is a dictonary, need to convert
        # from string
        if isinstance(getattr(Chkpt_obj, name), dict):
            try:
                setattr(Chkpt_obj, name, ast.literal.eval(attr))
            except:
                print("Chkpt error reading dictionary: "+str(name))
            continue

        # if this is a string 'None', convert to NoneType. This 
        # feels hacky, but, None is not a native type, so...
        if str(attr) == 'None':
            setattr(Chkpt_obj, name, None)
            continue

        setattr(Chkpt_obj, name, attr) 

    # don't close the file if it was called by handle -- if it's called
    # recursively, may trip up higher levels
    if file_handle is None:
        chkpt.close()

    if make_mol and obj_class == 'Molecule':
        Chkpt_obj.run()

    return Chkpt_obj

#
def read_wfn(file_name, obj_name, wfn_indx, cutoff):
    """reads the wfn object from the obj_name group. Returns
       a numpy array with the determinant coefficients and a 
       numpy array with the orbital occupations"""

    # using 64-bit integers to store det occupations
    STR_LEN = 64

    try:
        chkpt = h5py.File(file_name, 'r', libver='latest')
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

#--------------------------------------------------------

#
def pack(graci_obj):
    """packs the contents of the current object into two
       dictionaries: datasets and attributes"""

    graci_objs = tuple([getattr(globals()[obj_class.lower()], obj_class) 
                                    for obj_class in params.valid_objs])

    # list of all objects in the class
    objs = graci_obj.__dict__
    obj_attr = {}
    obj_dset = {}

    for name, obj in objs.items():
        # if object is a numpy array, save
        # as a distinct dataset
        if isinstance(obj, np.ndarray):
            obj_dset[name] = obj
            continue

        # if the object is a graci class object, save
        # the name of the object as a string -- this
        # is preferable to storing the same object
        # repeatedly
        if isinstance(obj, graci_objs):
            obj_attr[name] = 'GRACI_OBJ.'+ str(type(obj).__name__) + \
                             '.' + str(obj.label)
            continue

        obj_attr[name] = obj

    return obj_dset, obj_attr


