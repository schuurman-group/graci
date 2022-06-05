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

class Driver:
    def __init__(self):
        """constructor for Driver object. Not sure what parameters
           will be needed to be passed to the driver at the moment"""
        self.label = 'default'

    def run(self, calc_array):
        """Determine how to run the calculation given the array
           of postscf objects in the array argument"""

        # the first step is to match geometries to molecule sections
        mol_objs    = []
        scf_objs    = []
        ci_objs     = []
        postci_objs = []
        si_objs     = []
        for obj in calc_array:
            # identify the geometries in the run_list
            if type(obj).__name__ == 'Molecule':
                mol_objs.append(obj)
            elif type(obj).__name__ == 'Scf':
                scf_objs.append(obj)
            elif type(obj).__name__ in params.ci_objs:
                ci_objs.append(obj)
            elif type(obj).__name__ in params.postci_objs:
                postci_objs.append(obj)
            elif type(obj).__name__ in params.si_objs:
                si_objs.append(obj)
                
        # Load required libraries
        #-----------------------------------------------------
        # for now, assume postscf will require bitci
        if len(ci_objs) > 0:
            libs.lib_load('bitci')

        if len(si_objs) or len(postci_objs) > 0:
            libs.lib_load('bitsi')
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

                mol_lbls = [obj.label for obj in mol_objs]
                if scf_obj.label in mol_lbls:
                    mol_obj = mol_objs[mol_lbls.index(scf_obj.label)]
                elif len(mol_objs) == 1:
                    mol_obj = mol_objs[0]
                else:
                    mol_obj = None

                # if we have a label problem, then we should exit
                # with an error
                if mol_obj is None:
                    output.print_message('scf section, label='+
                            str(scf_obj.label)+
                            ' has no molecule section. Please check input')
                    sys.exit(1)

                scf_obj.run(mol_obj)
                chkpt.write(scf_obj)

            # initialize the MO integrals following the SCF run, but
            # finalise previous integrals if they exist. Not sure if this
            # should ultimately go here (probably not), but fine for now
            if len(ci_objs) > 0:
                libs.lib_func('bitci_int_finalize', [])
                if scf_obj.mol.use_df:
                    type_str = 'df'
                else:
                    type_str = 'exact'
                libs.lib_func('bitci_int_initialize', 
                              ['pyscf', type_str, scf_obj.moint_1e, 
                                scf_obj.moint_2e_eri])

            # CI Sections 
            #-----------------------------------------------------
            # run all the post-scf routines that map to the current
            # scf object.

            for ci_obj in ci_objs:

                # if labels match -- set the geometry
                # to the molecule object
                if ci_obj.label == scf_obj.label or len(scf_objs)==1:
                    ci_obj.run(scf_obj)
                    chkpt.write(ci_obj)

        # All SCF + CI objects are created and run() called before 
        # PostCI and subsequently SI objects are run()

        # PostCI Sections 
        # -- these can take ci_objects as arguments
        # ----------------------------------------------------
        for postci_obj in postci_objs:
            arg_list = self.get_postscf_objs(postci_obj, ci_objs)

            postci_obj.run(arg_list)
            chkpt.write(postci_obj)

        # State Interaction sections
        # -- these can take ci_objects or postci_objects as arguments
        #------------------------------------------------------------
        for si_obj in si_objs:
            arg_list = self.get_postscf_objs(si_obj, 
                                             ci_objs + postci_objs)            
            si_obj.run(arg_list)
            chkpt.write(si_obj)

        return

    #
    def get_postscf_objs(self, run_obj, avail_objs):
        """scan the run_obj to determine which postscf objects
           it requires in order to execute, and return them
           as a list
        """

        if type(run_obj).__name__ in params.postci_objs:
            lbls = list(run_obj.couple_groups)
        elif type(run_obj).__name__ in params.si_objs:
            lbls = [run_obj.final_label, run_obj.init_label]
        else:
            print('Cannot construct argument list for run() '+
                  ' method for object: '+str(run_obj.label))
            sys.exit(1)

        arg_list = [None]*len(lbls)

        # objects stored as bra/ket
        for chk_obj in avail_objs:
            indxs = [i for i, x in enumerate(lbls) if x==chk_obj.label]
            for indx in indxs:
                arg_list[indx] = chk_obj

        # if user labels are  not set (i.e. None) and there's 
        # only one postscf object, set label to that object
        indxs = [i for i, j in enumerate(lbls) if j == None]
        if len(avail_objs) == 1:
            for indx in indxs:
                arg_list[indx] = avail_objs[0]

        # we are agressive with exit calls. If there is a problem
        # with the requested coupling of sections, exit here
        if None in arg_list:
            output.print_message(type(run_obj).__name__+
                    ' section, label='+str(run_obj.label)+
                    ' is missing an interaction object. '+
                    ' Please check input')
            sys.exit(1)

        return arg_list

