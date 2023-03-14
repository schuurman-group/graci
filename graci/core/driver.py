"""
module that determines what types of calculations run
and in what order
"""
import sys as sys
import os as os
import graci.utils.basis as basis
import graci.utils.rydano as rydano
import graci.core.params as params
import graci.core.libs as libs
import graci.core.ao2mo as ao2mo
import graci.io.output as output
import graci.io.chkpt as chkpt

import numpy as np

class Driver:
    def __init__(self):
        """constructor for Driver object. Not sure what parameters
           will be needed to be passed to the driver at the moment"""
        self.label = 'default'

    def run(self, calc_array, save_to_chkpt=True):
        """Determine how to run the calculation given the array
           of postscf objects in the array argument"""

        # the first step is to match geometries to molecule sections
        mol_objs    = []
        scf_objs    = []
        ci_objs     = []
        postci_objs = []
        si_objs     = []
        param_objs  = []

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
            elif type(obj).__name__ == 'Parameterize':
                param_objs.append(obj)

        # Sanity check that sections of the same type have 
        # distinct label identifiers
        #-----------------------------------------------------
        self.check_labels(mol_objs)
        self.check_labels(scf_objs)
        self.check_labels(ci_objs)
        self.check_labels(postci_objs)
        self.check_labels(si_objs)

        # Load required libraries
        #-----------------------------------------------------
        # for now, assume postscf will require the bitci and
        # overlap libraries
        if len(ci_objs) or len(param_objs) > 0:
            libs.lib_load('bitci')
            libs.lib_load('overlap')
            
        if len(si_objs) or len(postci_objs) or len(param_objs) > 0:
            libs.lib_load('bitsi')
            libs.lib_load('bitwf')

        # Generate PySCF objects 
        # ----------------------------------------------------
        # generate the pyscf GTO Mole objects
        for mol_obj in mol_objs:
            mol_obj.run()
            if save_to_chkpt:
                chkpt.write(mol_obj)

        # print output file header
        output.print_header(calc_array)

        # Molecule Sections
        # ----------------------------------------------------
        #
        for mol_obj in mol_objs:
            # if we need to perform run-time basis set modifications
            if mol_obj.add_rydberg is not None:
                mol_obj = self.modify_basis(calc_array, mol_obj)

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
                scf_load = chkpt.read('Scf.' + scf_obj.label, 
                                      build_subobj = True,
                                      make_mol = True)
                
                if scf_load is None:
                    sys.exit('Cannot restart Scf, section = Scf.' + 
                              str(scf_obj.label) + 
                             ' not found in chkpt file = ' + 
                              str(output.file_names['chkpt_file']))

                # evidence that this is imperfect:
                scf_load.restart = True
                # call load to ensure AO -> MO transformation is run
                # and that orbitals are printed, etc.
                scf_load.load()

                # overwrite the original scf object in scf_objs
                # with the 'filled in' object
                # (required in case this scf object is going to be
                # used as the guess for another)
                scf_objs[scf_objs.index(scf_obj)] = scf_load.copy()
                
            # else assign molecule object and call run() routine
            else:

                # grab the corresponding molecule section
                mol_obj = self.match_sections(scf_obj.mol_label,'label',
                                              mol_objs, match_all=False) 

                # if we can't match a mol object, exit
                if mol_obj is None:
                    ostr = '\nCannot find molecule object for Scf ' + \
                    str(scf_obj.label) + ': Exiting...'
                    output.print_message(ostr)
                    sys.exit(1)

                # guess SCF object
                scf_guess = self.match_sections(scf_obj.guess_label, 
                                                'label', scf_objs, 
                                                 match_all=False)

                # run the SCF calculation
                scf_obj.run(mol_obj, scf_guess)
                
                # write scf object to checkpoint file
                if save_to_chkpt:
                    chkpt.write(scf_obj)

            # CI Sections 
            #-----------------------------------------------------
            # run all the post-scf routines that map to the current
            # scf object.
            # if there is a single scf object, ignore labels
            ci_calcs = self.match_sections(scf_obj.label, 
                                          'scf_label', ci_objs, 
                                           match_all=True)
            eri_mo   = ao2mo.Ao2mo()
            for ci_calc in ci_calcs:

                if ci_calc is None:
                    continue

                # perform AO -> MO integral transformation
                if eri_mo.emo_cut is None or \
                        eri_mo.emo_cut < ci_calc.mo_cutoff or \
                          eri_mo.precision != ci_calc.precision:
                    eri_mo.emo_cut = ci_calc.mo_cutoff
                    eri_mo.run(scf_obj, ci_calc.precision)

                # update ci object with results of ao2mo
                ci_calc.update_eri(eri_mo)

                # guess CI object
                ci_guess = self.match_sections(ci_calc.guess_label, 
                                               'label', ci_objs, 
                                                match_all=False)
               
                ci_calc.run(scf_obj, ci_guess)
                chkpt.write(ci_calc)

        # All SCF + CI objects are created and run() called before 
        # PostCI and subsequently SI objects are run()

        # PostCI Sections 
        # -- these can take ci_objects as arguments
        # ----------------------------------------------------
        for postci_obj in postci_objs:
            obj_list = self.get_postscf_objs(postci_obj, ci_objs)

            postci_obj.run(obj_list)
            if save_to_chkpt:
                chkpt.write(postci_obj)

        # State Interaction sections
        # -- these can take ci_objects or postci_objects as arguments
        #    ...but expect exactly two arguments: a bra object and
        #    a ket object
        #------------------------------------------------------------
        for si_obj in si_objs:
            [bra, ket] = self.get_postscf_objs(si_obj, 
                                               ci_objs + postci_objs)            
            si_obj.run(bra, ket)
            if save_to_chkpt:
                chkpt.write(si_obj)

        # Hamiltonian Parameterization
        # -- this is a special case for now
        #-------------------------------------------------------------
        for param_obj in param_objs:
            param_obj.run()

        return
 
    #
    def match_sections(self, label, sec_lbl, sec_lst, match_all=False):
        """Match the sections in sec_lst that have the variable 
           sec_lbl equal to label. If lst is True, return all
           cases that match, else, just return the first match

        Arguments:
           label:   the string label to match
           sec_lbl: the variable name that is being comared to label
           sec_lst: list of candidate objects
           match_all: If False, return first true match, if true, 
                    return all matches

        Returns:
           match:   list (if lst==True) of matches, else, just
                    the matching object
        """

        sec_lbls = [getattr(sec_lst[i], sec_lbl) 
                            for i in range(len(sec_lst))]
        secs     = []

        for i in range(len(sec_lst)):
            if label == sec_lbls[i]:
                secs.append(sec_lst[i])

        # If multiple sections match and match_all is false, this
        # likely constitutes a labelling problem, and we're going
        # to exit the program, rather than just give a warning
        if not match_all and len(secs) > 1:
            output.print_message('Multiple objects matched to label=' +
                    str(label)+': '+str([sec.label for sec in secs]))
            sys.exit(1)

        # if we still don't have any matches, return None 
        if len(secs) == 0:
            secs = [None]

        if match_all:
            return secs
        else:
            return secs[0]

    # 
    def check_labels(self, obj_lst):
        """Check that each label in object list is unique

        Arguments:
          obj_list: list of objects, each possessing 'label' variable

        Returns:
          None
        """

        if len(obj_lst) == 0:
            return 

        lbls = [obj.label for obj in obj_lst]
        uniq = list(set(lbls))

        if len(uniq) != len(lbls):
            c_name = type(obj_lst[0]).__name__
            ostr = '\nSection type: '+str(c_name) + '\n' + \
                   'Section labels: '+str(lbls) + '\n' + \
                   '--> Sections of the same type require unique labels'
            output.print_message(ostr)
            sys.exit(1)

        return

    #
    def modify_basis(self, calc_array, mol_obj):
        """
        modify the atomic basis set at run-time if need be
        """
 
        # right now this is just to add Rydberg functions
        if mol_obj.add_rydberg is not None:

            # if this is a valid contraction, construct a default
            # rydano object
            ncon = basis.str_to_contract(mol_obj.add_rydberg)
            if ncon is not None:
                ryd_basis = rydano.Rydano()
                ryd_basis.label = 'auto-generated with defaults'
                ryd_basis.contract = mol_obj.add_rydberg
                ryd_basis.run(mol_obj)

            # else, load the rydano object defined in input file
            else:
                ryd_basis = None
                for obj in calc_array:
                    if (type(obj).__name__  == 'Rydano' and
                              obj.label == mol_obj.add_rydberg):
                        ryd_basis = obj
                        break
                if ryd_basis is None:
                    print('Molecule.rydano = ' +
                           str(mol_obj.add_rydberg) + ' not found.')
                    sys.exit(1)
                ryd_basis.run(mol_obj)

        return mol_obj


    #
    def get_postscf_objs(self, run_obj, avail_objs):
        """scan the run_obj to determine which postscf objects
           it requires in order to execute, and return them
           as a list
        """

        if type(run_obj).__name__ in params.postci_objs:
            lbls = list(run_obj.couple_groups)
        elif type(run_obj).__name__ in params.si_objs:
            if hasattr(run_obj, 'final_label'):
                lbls = [run_obj.final_label, run_obj.init_label]
            elif hasattr(run_obj, 'bra_label'):
                lbls = [run_obj.bra_label, run_obj.ket_label]
            else:
                print('Cannot find init/ket, final/bra states in '+
                      'state interaction object: '+str(run_obj.label))
                sys.exit(1)
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

