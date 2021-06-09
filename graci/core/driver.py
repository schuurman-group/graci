"""
module that determines what types of calculations run
and in what order
"""
import graci.core.params as params
import graci.core.libs as libs
import graci.io.output as output

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

        # the first step is to match geometries to molecule sections
        gm_objs      = []
        mol_objs     = []
        scf_objs     = []
        postscf_objs = []
        si_objs      = []
        for obj in calc_array:
            # identify the geometries in the run_list
            if obj.name() == 'geometry':
                gm_objs.append(obj)
            elif obj.name() == 'molecule':
                mol_objs.append(obj)
            elif obj.name() == 'scf':
                scf_objs.append(obj)
            elif obj.name() in params.method_objs:
                postscf_objs.append(obj)
            elif obj.name() in params.si_objs:
                si_objs.append(obj)

        # Molecle sections 
        # ----------------------------------------------------

        # match geometry objects to molecule objects
        for mol_obj in mol_objs:

            for gm_obj in gm_objs:
                # if labels match -- set the geometry
                # to the molecule object
                if mol_obj.label == gm_obj.label:
                    mol_obj.set_pymol(gm_obj)
                    break

            # if we didn't match labels, but there is a single
            # geometry object, 
            if mol_obj.pymol_exists() is False and len(gm_objs)==1:
                mol_obj.set_pymol(gm_obj)

            # if we have a label problem, then we should exit
            # with an error
            if mol_obj.pymol_exists() is False:
                output.print_message('molecule section, label='+
                        str(mol_obj.label)+
                        ' has no geometry set. Please check input')
                sys.exit(1)

        # print output file header
        output.print_header(calc_array)

        # SCF Sections 
        # -----------------------------------------------------

        # match geometry objects to molecule objects
        for scf_obj in scf_objs:

            for mol_obj in mol_objs:
                # if labels match -- set the geometry
                # to the molecule object
                if scf_obj.label == mol_obj.label:
                    scf_obj.set_mol(mol_obj)
                    break

            # if we didn't match labels, but there is a single
            # molecule object, 
            if scf_obj.mol_exists() is False and len(mol_objs)==1:
                scf_obj.set_mol(mol_obj)

            # if we have a label problem, then we should exit
            # with an error
            if scf_obj.mol_exists() is False:
                output.print_message('scf section, label='+
                        str(scf_obj.label)+
                        ' has no molecule section. Please check input')
                sys.exit(1)
        
            scf_obj.run()

            # initialize the MO integrals following the SCF run, but
            # finalise previous integrals if they exist
            if libs.lib_exists('bitci'):
                libs.finalise_intpyscf()
            libs.init_intpyscf(scf_obj.mol, scf_obj)

            # Post-SCF Sections 
            #-----------------------------------------------------
            # run all the post-scf routines that map to the current
            # scf object.

            for postscf in postscf_objs:

                # if labels match -- set the geometry
                # to the molecule object
                if postscf.label == scf_obj.label:
                    postscf.set_scf(scf_obj)
                    break

                # if we didn't match labels, but there is a single
                # molecule object, 
                if postscf.scf_exists() is False and len(scf_objs)==1:
                    postscf.set_scf(scf_obj)

                # if we have a label problem, then we should exit
                # with an error
                if postscf.scf_exists() is False:
                    output.print_message(postscf.name() +
                            ' section, label=' + str(postscf.label) +
                            ' has no scf section. Please check input')
                    sys.exit(1)

                postscf.run()

        # State interaction Sections 
        # ----------------------------------------------------

        for si_obj in si_objs:

            # first check if init/final_method is set. If not
            # and there is a single postscf method, we can
            # determine a sensible default
            if (si_obj.init_method is None and si_obj.final_method 
                    is None and len(postscf_objs) == 1):
                si_obj.set_bra(postscf_objs[0])
                si_obj.set_ket(postscf_objs[0])

            # else, we use the label names to determine bra
            # and ket states. If things don't match up, fail
            # with an error message
            else:
                for postscf in postscf_objs:
                    if si_obj.init_method == postscf.label:
                        si_obj.set_ket(postscf)
                    if si_obj.final_method == postscf.label:
                        si_obj.set_bra(postscf)

            if (si_obj.bra_exists() is False or 
                si_obj.ket_exists() is False):
                output.print_message(si_obj.name()+' section, label='+
                        str(si_obj.label)+
                        ' has no bra/ket defined. Please check input')
                sys.exit(1)

            si_obj.run()


        return


