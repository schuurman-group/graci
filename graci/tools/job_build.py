"""
module that determines what types of calculations run
and in what order
"""
import graci.methods.params as params

def run(calc_array):
    """Determine how to run the calculation given the array
       of method objects in the array argument"""

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

