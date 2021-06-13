"""
Module for facilitating passing of strings between
Python and compiled dlls
"""
import sys as sys
import h5py as h5py
import numpy as np
import graci.core.params as params
import graci.core.molecule as molecule
import graci.core.geometry as geometry

orb_types  = ['scf', 'no', 'ndo', 'nto']
calc_types = params.method_objs


#
def write_chkpt(file_name, graci_obj):
    """write the graci_obj to the checkpoint file"""

    # get the data sets and attributes
    dsets, attrs = pack(graci_obj)

    if len(dsets)+len(attrs) == 0:
        return

    # the group name will by the class type + object label
    grp_name = str(type(graci_obj).__name__)+'.'+str(graci_obj.label)

    chkpt = h5py.File(file_name, 'a', libver='latest')
    
    # create the group if doesn't already exist
    if grp_name not in chkpt.keys():
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
            obj_grp.create_dataset(name, data=dset, dtype=dset.dtype)
            # if non-native HDf5 type, just convert to a string
        except TypeError:
            dset_str = str(dset)
            obj_grp.create_dataset(name, data=dset_str)

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
def pack(graci_obj):
    """packs the contents of the current object into two 
       dictionaries: datasets and attributes"""
 
    # list of all objects in the class
    objs = graci_obj.__dict__
    obj_attr = {}
    obj_dset = {}

    for name, obj in objs.items():
        # if object is a numpy array of rank >=2, save
        # as a distinct dataset
        if isinstance(obj, np.ndarray):
            if len(obj.shape) >= 2:
                obj_dset[name] = obj
                continue

        obj_attr[name] = obj

    return obj_dset, obj_attr

#
def chkpt_contents(file_name):
    """return the group names present in the checkpoint file"""

    return

#---------------------------------
#
def read_chkpt(file_name, obj_name):
    """return a molecule object after parsing checkpoint file. If label 
       is 'None', then return a molecule object only if there is just 
       one in file, else fail due to ambiguous request"""

    chkpt = h5py.File(file_name, 'a', libver='latest')

