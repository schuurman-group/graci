"""
The Parameterize object and its associated functions.
"""
import sys as sys
import numpy as np
import mpi4py as mpi4py
import os as os

import graci.core.libs as libs
import graci.core.driver as driver
import graci.io.chkpt as chkpt
import graci.io.parse as parse

class Parameterize:
    """Class constructor for the Parameterize object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.algorithm   = None
        self.label       = 'parameterize'

        self.hamiltonian    = ''
        self.graci_ref_data = ''
        self.target_file    = ''
        self.pthresh        = 1.e-5

        # -------------------------------
        self.target_data    = {}
        self.ndata          = 0
        self.iter           = 0

    def run(self):
        """re-parameterize the Hamiltonian"""
 
        # parse the target data file
        target_data = parse_target_file()

        # parse the graci reference file
        scf_data, ci_data = parse_graci_file(target_data)

        # pull the reference CI states
        ref_wf = extract_ci_states(ref_chkpt,ci_data[molecule],data)

        # print the parameterize header
        # show the mappings between the target data and the data in the
        # reference checkpoint file
        output.print_param_header(target_data, ci_data)

        # set up the directory structure
        for target in target_data.keys():
            try:
                os.mkdir(target)
            except FileExistsError:
                os.rmdir(target)
                os.mkdir(target)

        # run the first pass
        eval_energies(target_data, scf_data, ci_data, run_scf=True)

        # set the initial values of the parameters using the
        # named Hamiltonian
        args = (self.hamiltonian)
        p0 = libs.lib_func('get_params', args)

        lsq_out = scipy.optimize.leastsq(err_func, p0, 
                             args = (target_data, scf_data, ci_data),
                             ftol = self.pthresh, full_output=True)
        (fitp, cov_p, info, mesg, ierr) = lsq_out  
    

        output.print_param_results(fitp, info)

        return

    #
    def err_func(params, target_data, scf_data, ci_data):
        """
        Evaluate the error function. In this case, simple RMSD

        Arguments:

        targets:      target values
        values:       current values

        Returns:
        norm of the different between target and current values
        """

        print(' ITERATION '+self.iter+' params='+str(params))
        self.iter += 1

        # update the Hamiltonian
        libs.lib_func('update_params', params)

        iter_ener = eval_energies(target_data, scf_data, ci_data)
        dif_vec   = np.zeros(self.ndata, dtype=float)

        # this approach assumes dict is ordered! Only true from Python
        # 3.6 onwards...
        ide = 0
        for molecule, states in target_data.items():
            for trans, ener in target_data[molecule].items():
                init, final  = trans.split().strip()
                de_iter      = iter_ener[molecule][final] - \
                               iter_ener[molecule][init]
                dif_vec[ide] = de_iter - ener

        return dif_vec


    #
    def eval_energies(target_data, scf_data, ci_data, run_scf=False):
        """
        evaluate all the energies in the graci data set

        """

        # open the reference data file and get the contents
        if os.path.isFile(self.graci_ref_data):
            ref_chkpt = h5py.File(self.graci_ref_data, 'r',
                                                   libver='latest')
        else:
            print('File: '+str(self.graci_ref_data)+' not found. '+
                  ' Exiting...')
            sys.exit(1)

        energies = {}

        # run through all the reference data objects
        for molecule, states in target_data.items():

            # pull the reference SCF and CI objects
            scf_name = scf_data[molecule]
            scf_obj  = chpt.read(scf_name, file_handle=ref_chkpt)

            ci_names = ci_data[molecule]
            ci_objs  = [chkpt.read(ci, file_handle=ref_chkpt) 
                                                for ci in ci_names]

            # pull the reference CI states
            ref_det, ref_cf = extract_ci_states(ref_chkpt, 
                                                ci_data[molecule], 
                                                states)

            # change to the appropriate sub-directory
            os.chdir(molecule)

            # run the scf (and generate integral files if need be)
            if run_scf:
                scf_obj.load()

            # run all the CI objects 
            for ci_obj in ci_objs:
                ci_obj.run(scf_obj)

            # use overlap with ref states to identify relevant states
            ener_match = identify_states(ref_det, ref_cf, 
                                         ci_objs, states)

            energies[molecule] = {}
            for state, ener in ener_match.items():
                energies[molecule][state] = ener

        return energies
 
    #
    def identify_states(ref_det, ref_cf, ci_objs, states):
        """
        compute overlaps to identify states
        """

        eners = {}

        for ci_obj in ci_objs:
            ci_dets = ci_obj.det_strings
            ci_cf   = ci_obj.vec_det

            for lbl in ref_det.keys():
                Sij = compute_overlaps(ref_det[lbl], ref_cf[lbl], 
                                            ci_dets, ci_cf, 
                                            ci_obj.scf.orbs)
                if max(Sij) >= 0.9:
                    st = maxloc(Sij)
                    eners[lbl] = ci_obj.energies[st]

        if list(eners.keys()).sort() != list(ref_det.keys()).sort():
            sys.exit('Could not identify states: ' + str(ci_objs) + 
                     ' lbl:'+str(states))

        return eners

    # 
    def extract_ci_states(chkpt, obj_names, states):
        """
        Extract det list expansions for the reference states

        Arguments:
        chkpt:        file handle for checkpoint file
        ci_objs:      CI objects to check for needed states
        states:       labels for the states to extract

        Returns:
        ref_wf:       dictionary of bitstring wfs
        """

        ref_dets = {}
        ref_cf   = {}

        # run through the ci objects looking for states with the
        # correct label
        for obj_name in obj_names:

            # variables names are stored as attributes
            var_names = chkpt[obj_name].attrs
            ener_dset = var_names['energies']
            vec_dset  = var_names['det_strings']
            cf_dset   = var_names['vec_det']

            attrs = chkpt[ci_obj+'/'+ener_dset].attrs

            if 'labels' not in attrs.keys():
                continue

            for lbl,st in attrs['labels'].items():
                if lbl in states.keys(): 
                    st = states[lbl]
                    ref_vec[lbl] = chkpt[obj_name+'/'+vec_dset][st]
                    ref_cf[lbl]  = chkpt[obj_name+'/'+cf_dset][:,st]

        if list(states.keys()).sort() != list(ref_vec.keys()).sort():
            sys.exit('Could not find all requested states: ' + 
                      str(obj_names)+' lbls='+str(states))

        return ref_vec, ref_cf

    # 
    def compute_overlaps(ref_det, ref_cf, test_dets, test_cf, orbs):
        """
        compute overlaps
        """
       
        Sij = np.zeros(len(test_dets.keys()), dtype=float)

        # number of integers to represent bit string
        nmo      = orbs.shape[1]
        nint     = ref_det.shape[1]
        ndetb    = ref_det.shape[0]
        ndetk    = test_dets.shape[1]
        nrtb     = 1
        nrtk     = test_dets.shape[0]
        detb     = [ref_det]
        detk     = test_dets
        vecb     = [ref_cf]
        veck     = test_cf
        smo      = np.identity(nmo)
        thrsh    = 0.9
        ncore    = 0
        icore    = []
        lfrzcore = True
        npairs   = nrtk
        ipairs   = [[1,i] for i in range(nrtk)]
        args   = (nmo, nint, ndetb, ndetk, nrtb, nrtk, detb, detk, 
                  vecb, veck, smo, thrsh, ncore, icore, lfrzcore, 
                  npairs, Sij, ipairs)

        libs.lib_func('overlap_wf', args)


        return Sij

    #
    def parse_graci_file(target_vals):
        """
        Parse the GRaCI reference data and confirm what data is present
        that can be compared to the target values

        Arguments: 
        target_vals:    dictionary containing the reference data

        Returns:
        graci_vals:     a dictionary containing the corresponding GRaCI
                        objects
        """

       # open the reference data file and get the contents
        if os.path.isFile(self.graci_ref_data):
            ref_chkpt = h5py.File(self.graci_ref_data, 'r',
                                                   libver='latest')
        else:
            print('File: '+str(self.graci_ref_data)+' not found. '+
                  ' Exiting...')
            sys.exit(1)

        # get the top-level contents of the checkpoint file
        ref_contents = chkpt.contents(file_name = self.graci_ref_data)

        ci_objs  = {}
        scf_objs = {}
        for molecule in target_vals.keys():
            ci_objs[molecule]  = []
            scf_objs[molecule] = None

            # if the molecule string is in the reference 
            # object name, add it ref_objs dict
            for name in ref_contents:
                ci_type = any([ci in name for ci in params.ci_objs])
                if molecule in name:
                    if ci_type:
                        ci_objs[molecule].append(name)
                    elif 'Scf' in name:
                        scf_objs[molecule] = name

        return scf_objs, ci_objs

    #
    def parse_target_file():
        """
        Parse the file containing the data we're going to parameterize
        wrt to

        Arguments: None
        Returns:   A dictionary containing the target data
        """

        with open(self.target_file, 'r') as t_file:
            t_file_lines = t_file.readlines()

        self.ndata  = 0
        target_dict = {}
        for line in tfile_lines:
            molecule              = line[0]
            target_dict[molecule] = {} 
            # format is: state1 state2 energy
            for i in range(1,len(line),3):
                try:
                    states = line[i].strip()+' '+line[i+1].strip()
                    val    = parse.convert_value(line[i+2])
                    target_dict[molecule][states] = val
                    self.ndata += 1
                except:
                    print('Error parsing value as str/float: '+str(val))
                    print('line = '+str(line))
                    sys.exit(1)
        
        return target_dict

