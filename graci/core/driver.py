"""
module that determines what types of calculations run
and in what order
"""
import sys as sys
import os as os
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
        overlap_objs = []
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
            elif type(obj).__name__ in params.overlap_objs:
                overlap_objs.append(obj)
                
        # Load required libraries
        #-----------------------------------------------------
        # for now, assume postscf will require bitci
        if len(postscf_objs) > 0:
            libs.lib_load('bitci')

        if len(si_objs) > 0:
            libs.lib_load('bitsi')

        if len(overlap_objs) > 0:
            libs.lib_load('bitwf')
            
        # Molecule sections 
        # ----------------------------------------------------

        # generate the pyscf GTO Mole objects
        for mol_obj in mol_objs:
            mol_obj.run()
            chkpt.write(mol_obj)

        # print output file header
        output.print_header(calc_array)

        # SCF Sections 
        # -----------------------------------------------------

        # match scf objects to molecule objects
        for scf_obj in scf_objs:

            # if restarting, load KS orbitals from chkpt, but
            # still do AO -> MO transformation (via call to
            # scf_obj.load())
            if scf_obj.restart:
                # this should be changed: we should only set
                # scf_obj to the object read from the chkpt file
                # after we've confirmed they're the same..
                scf_obj = chkpt.read('Scf.' + scf_obj.label, 
                                     build_subobj = True,
                                     make_mol = True)

                if scf_obj is None:
                    sys.exit('Cannot restart Scf, section = Scf.' + 
                              str(scf_obj.label) + 
                             ' not found in chkpt file = ' + 
                              str(output.file_names['chkpt_file']))

                # evidence that this is imperfect:
                scf_obj.restart = True
                # call load to ensure AO -> MO transformation is run
                # and that orbitals are printed, etc.
                scf_obj.load()

            # else assign molecule object and call run() routine
            else:
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
                chkpt.write(scf_obj)

            # initialize the MO integrals following the SCF run, but
            # finalise previous integrals if they exist. Not sure if this
            # should ultimately go here (probably not), but fine for now
            if len(postscf_objs) > 0:
                libs.lib_func('bitci_int_finalize', [])
                type_str = 'exact'
                if scf_obj.mol.use_df:
                    type_str = 'df'
                libs.lib_func('bitci_int_initialize', 
                              ['pyscf', type_str, scf_obj.moint_1e, 
                                scf_obj.moint_2e_eri])

            # Post-SCF Sections 
            #-----------------------------------------------------
            # run all the post-scf routines that map to the current
            # scf object.

            for postscf in postscf_objs:

                # if labels match -- set the geometry
                # to the molecule object
                if postscf.label == scf_obj.label or len(scf_objs)==1:
                    postscf.set_scf(scf_obj)

                    postscf.run()
                    chkpt.write(postscf)

        # All SCF + POST-SCF objects are run first before state interaction
        # objects are evaluated

        # State interaction Sections 
        # ----------------------------------------------------

        for si_obj in si_objs:

            obj_list = self.extract_si_obj_list(si_obj, postscf_objs)
            if None in obj_list:
                output.print_message(type(si_obj).__name__+' section, '+
                        'label='+str(si_obj.label)+
                        ' is missing an interaction object. '+
                        ' Please check input')
                sys.exit(1)

            si_obj.run(obj_list)
            chkpt.write(si_obj)

        return

    # 
    def extract_si_obj_list(si_obj, postscf_objs):
        """extract the list of required objects for an state 
           iteraction object based on user input"""

        # list of objects to be passed to interaction class
        if type(si_obj).__name__ == 'Transition':
            lbls = [si_obj.final_label, si_obj.init_label]
        elif type(si_obj).__name__ == 'Spinorbit':
            lbls = si_obj.couple_groups
        elif type(si_obj).__name__ == 'Overlap':
            lbls = [si_obj.bra_label, si_obj.ket_label]

        obj_list = [None]*len(lbls)

        # objects stored as bra/ket
        for postscf in postscf_objs:
            if postscf.label in lbls:
                obj_list[lbl.index(postscf.label)] = postscf

        # if user labels are  not set (i.e. None) and there's 
        # only one postscf object, set label to that object
        indices = [i for i, j in enumerate(lbls) if j == None]
        if len(postscf_objs) == 1:
            for indx in indices:
                lbls[indx] = postscf_objs[0]
                
        return obj_list
