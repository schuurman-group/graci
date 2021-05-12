"""
module that determines what types of calculations run
and in what order
"""
import graci.core.params as params

class Driver:
    def __init__(self):
        """constructor for Driver object. Not sure what parameters
           will be needed to be passed to the driver at the moment"""
        self.label = 'default'

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
        gm_grp     = []
        mol_grp    = []
        scf_grp    = []
        method_grp = []
        si_grp     = []
        for obj in calc_array:
            # identify the geometries in the run_list
            if obj.name() == 'geometry':
                gm_grp.extend([obj])
            elif obj.name() == 'molecule':
                mol_grp.extend([obj])
            elif obj.name() == 'scf':
                scf_grp.extend([obj])
            elif obj.name() in params.method_objs:
                method_grp.extend([obj])
            elif obj.name() in params.si_objs:
                si_grp.extend([obj])

        # initialize each of the molecule objects with
        # the appropriate geometry object if they don't
        # already exit
        for mol_obj in mol_grp:
            if mol_obj.pymol_exists() is False:
                for gm_obj in gm_grp:
                    if mol_obj.label == gm_obj.label:
                        mol_obj.set_pymol(gm_obj)
                        break
            if mol_obj.pymol_exists() is False:
                output.print_message('molecule section, label='+
                        str(mol_obj.label)+
                        ' has no geometry set. Removing from run list')
                mol_grp.remove(mol_obj)

        # run each of the methods in turn
        for mol in mol_grp:
            for scf in scf_grp:
                for method in method_grp:
                    method.set_mol(mol)
                    method.set_scf(scf)
                    method.run()

        # run the state interaction computations
        for si in si_grp:

            final_method = None
            init_method  = None

            if len(method_grp) == 1:
                final_method = method_grp[0]
                init_method  = method_grp[0]

            else:
                # if there's only one method, set the default values
                # of si.final_method and si.init_method to be the same
                for method in method_grp:
                    if method.label == si.final_method:
                        final_method = method
                    if method.label == si.init_method:
                        init_method = method

            if final_method is not None and init_method is not None:
                si.set_init_method(init_method)
                si.set_final_method(final_method)
                si.run()


#######################################################

