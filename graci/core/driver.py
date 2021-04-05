"""
module that determines what types of calculations run
and in what order
"""
import graci.core.params as params

class Driver:
    def __init__(self):
        """constructor for Driver object. Not sure what parameters
           will be needed to be passed to the driver at the moment"""
        self.label = ''

    def name(self):
        """ return the name of the class object as a string"""
        return 'driver'

    def run(self, calc_array):
        """Determine how to run the calculation given the array
           of method objects in the array argument"""


        # this is currently a very simple implementation of this method.
        # eventually there will have to be a lot of logic involved in 
        # determining how to run the various sections.

        # right now, it just looks for a mol and scf objects, and then
        # runs each (non-scf) method using the mol and scf objects 
        # in the list
        gm_arr     = []
        mol_arr    = []
        scf_arr    = []
        method_arr = []
        for obj in calc_array:
            # identify the geometries in the run_list
            if obj.name() == 'geometry':
                gm_arr.extend([obj])
            elif obj.name() == 'molecule':
                mol_arr.extend([obj])
            elif obj.name() == 'scf':
                scf_arr.extend([obj])
            elif obj.name() in params.valid_objs:
                method_arr.extend([obj])

        # initialize each of the molecule objects with
        # the appropriate geometry object if they don't
        # already exit
        for mol_obj in mol_arr:
            if mol_obj.pymol_exists() is False:
                for gm_obj in gm_arr:
                    if mol_obj.label == gm_obj.label:
                        mol_obj.set_pymol(gm_obj)
                        break
            if mol_obj.pymol_exists() is False:
                output.print_message('molecule section, label='+
                        str(mol_obj.label)+
                        ' has no geometry set. Removing from run list')
                mol_arr.remove(mol_obj)

        for mol in mol_arr:
            for scf in scf_arr:
                for method in method_arr:
                    method.run(mol, scf)

        return


#######################################################

