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
import graci.core.parameterize as parameterize
import graci.methods.scf as scf
import graci.methods.dftmrci as dftmrci
import graci.methods.dftcis as dftcis
import graci.properties.moments as moments
import graci.properties.transition as transition
import graci.properties.spinorbit as spinorbit

#
def contents(file_name):
    """return the group names present in the checkpoint file"""

    try:
        chkpt = h5py.File(file_name, 'r', libver='latest')
    except:
        return None

    grp_lst = list(chkpt.keys())

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
            obj_grp.create_dataset(name, data=dset, dtype=dset.dtype, 
                                                   compression="gzip")
            # if non-native HDf5 type, just convert to a string
        except TypeError:
            dset_str = str(dset)
            obj_grp.create_dataset(name, data=dset_str, 
                                                   compression="gzip")

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

        if not hasattr(Chkpt_obj, name):
            raise AttributeError(
            'Class '+str(obj_class)+' has no attribute '+str(name))
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




