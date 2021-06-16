"""
module that determines what types of calculations run
and in what order
"""
import sys as sys
import graci.core.params as params
import graci.core.libs as libs
import graci.io.output as output
import graci.io.chkpt as chkpt
import graci.core.molecule as molecule

class Driver:
    def __init__(self):
        """constructor for Driver object. Not sure what parameters
           will be needed to be passed to the driver at the moment"""
        self.label = 'default'

    def run(self, calc_array):
        """Determine how to run the calculation given the array
           of postscf objects in the array argument"""

        # the first step is to match geometries to molecule sections
        mol_objs     = []
        scf_objs     = []
        postscf_objs = []
        si_objs      = []
        for obj in calc_array:
            # identify the geometries in the run_list
            if type(obj).__name__ == 'Molecule':
                mol_objs.append(obj)
            elif type(obj).__name__ == 'Scf':
                scf_objs.append(obj)
            elif type(obj).__name__ in params.method_objs:
                postscf_objs.append(obj)
            elif type(obj).__name__ in params.si_objs:
                si_objs.append(obj)

        # Molecule sections 
        # ----------------------------------------------------

        # generate the pyscf GTO Mole objects
        for mol_obj in mol_objs:
            mol_obj.run()
            chkpt.write(output.file_names['chkpt_file'], mol_obj)

        # print output file header
        output.print_header(calc_array)

        # SCF Sections 
        # -----------------------------------------------------

        # match scf objects to molecule objects
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
            chkpt.write(output.file_names['chkpt_file'], scf_obj)

            #TESTING
            #scf_test = chkpt.read(output.file_names['chkpt_file'], scf_obj, build_objs=True)

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
                    output.print_message(type(postscf).__name__ +
                            ' section, label=' + str(postscf.label) +
                            ' has no scf section. Please check input')
                    sys.exit(1)

                postscf.run()
                chkpt.write(output.file_names['chkpt_file'], postscf)

        # State interaction Sections 
        # ----------------------------------------------------

        for si_obj in si_objs:

            # first check if init/final_method is set. If not
            # and there is a single postscf method, we can
            # determine a sensible default
            if (si_obj.init_label is None and si_obj.final_label
                    is None and len(postscf_objs) == 1):
                si_obj.set_bra(postscf_objs[0])
                si_obj.set_ket(postscf_objs[0])

            # else, we use the label names to determine bra
            # and ket states. If things don't match up, fail
            # with an error message
            else:
                for postscf in postscf_objs:
                    if si_obj.init_label == postscf.label:
                        si_obj.set_ket(postscf)
                    if si_obj.final_label == postscf.label:
                        si_obj.set_bra(postscf)

            if (si_obj.bra_exists() is False or 
                si_obj.ket_exists() is False):
                output.print_message(type(si_obj).__name__+' section, '+
                        'label='+str(si_obj.label)+
                        ' has no bra/ket defined. Please check input')
                sys.exit(1)

            si_obj.run()
            chkpt.write(output.file_names['chkpt_file'], si_obj)

        return


