"""
module that determines what types of calculations run
and in what order
"""
import graci.methods.params as params

def run(calc_array):
    """Determine how to run the calculation given the array
       of method objects in the array argument"""


    # this is currently a very simple implementation of this method.
    # eventually there will have to be a lot of logic involved in 
    # determining how to run the various sections.

    # right now, it just looks for a mol and scf objects, and then
    # runs each (non-scf) method using the mol and scf objects 
    # in the list
    mol_arr    = []
    scf_arr    = []
    method_arr = []
    for obj in calc_array:
        if params.method_name(obj) == 'molecule':
            mol_arr.extend([obj])
        elif params.method_name(obj) == 'scf':
            scf_arr.extend([obj])
        else:
            method_arr.extend([obj])

    for mol in mol_arr:
        for scf in scf_arr:
            for method in method_arr:
                method.run(mol, scf)

    return


#######################################################

