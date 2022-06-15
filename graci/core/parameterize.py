"""
The Parameterize object and its associated functions.
"""
import sys as sys
import numpy as np
import h5py as h5py
import mpi4py as mpi4py
import os as os
import graci.core.libs as libs
import graci.core.driver as driver
import graci.core.params as params
import graci.io.chkpt as chkpt
import graci.io.parse as parse
import graci.overlaptools.overlap as overlap

class Parameterize:
    """Class constructor for the Parameterize object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.algorithm   = None
        self.label       = 'parameterize'

        self.hamiltonian    = ''
        self.graci_ref_file = ''
        self.target_file    = ''
        self.pthresh        = 1.e-5

        # -------------------------------
        self.target_data    = {}
        self.ref_states     = {}
        self.ndata          = 0
        self.iter           = 0

    def run(self):
        """re-parameterize the Hamiltonian"""
 
        # parse the target data file
        target_data = self.parse_target_file()

        # parse the graci reference file
        scf_data, ci_data, self.ref_states = \
                               self.parse_graci_file(target_data)

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
        self.eval_energies(target_data, scf_data, ci_data, run_scf=True)

        # set the initial values of the parameters using the
        # named Hamiltonian
        args = (self.hamiltonian)
        p0 = libs.lib_func('get_params', args)

        lsq_out = scipy.optimize.leastsq(self.err_func, p0, 
                             args = (target_data, scf_data, ci_data),
                             ftol = self.pthresh, full_output=True)
        (fitp, cov_p, info, mesg, ierr) = lsq_out  
    

        output.print_param_results(fitp, info)

        return

    #
    def err_func(self, hparams, target_data, scf_data, ci_data):
        """
        Evaluate the error function. In this case, simple RMSD

        Arguments:

        targets:      target values
        values:       current values

        Returns:
        norm of the different between target and current values
        """

        print(' ITERATION '+self.iter+' params='+str(hparams))
        self.iter += 1

        # update the Hamiltonian
        libs.lib_func('update_params', hparams)

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
    def eval_energies(self, target_data, scf_names, ci_names, 
                                                    run_scf=False):
        """
        evaluate all the energies in the graci data set

        """

        # open the reference data file and get the contents
        if os.path.isFile(self.graci_ref_file):
            ref_chkpt = h5py.File(self.graci_ref_file, 'r',
                                                   libver='latest')
        else:
            print('File: '+str(self.graci_ref_file)+' not found. '+
                  ' Exiting...')
            sys.exit(1)

        energies = {}

        # run through all the reference data objects
        for molecule, states in target_data.items():

            # pull the reference SCF and CI objects
            scf_name = scf_names[molecule]
            scf_obj  = chpt.read(scf_name, file_handle=ref_chkpt)

            ci_name = ci_names[molecule]
            ci_objs  = [chkpt.read(ci, file_handle=ref_chkpt) 
                                                for ci in ci_name]

            # pull the reference CI states
            #ref_det, ref_cf = extract_ci_states(ref_chkpt, 
            #                                    ci_data[molecule], 
            #                                    states)

            # copy the scf and CI objects and re-run them
            scf_run = scf_obj.copy()
            ci_runs = [ci_obj.copy() for ci_obj in ci_objs]
         
            # change to the appropriate sub-directory
            os.chdir(molecule)

            # run the scf (and generate integral files if need be)
            if run_scf:
                scf_run.load()

            # run all the CI objects 
            for ci_run in ci_runs:
                ci_run.run(scf_run)

            # use overlap with ref states to identify relevant states
            ener_match = identify_states(molecule, ci_objs, ci_runs) 

            energies[molecule] = {}
            for state, ener in ener_match.items():
                energies[molecule][state] = ener

        return energies
 
    #
    def identify_states(self, molecule, ref_ci, new_ci):
        """
        compute overlaps to identify states
        """

        eners  = {}
        bra_st = []
        smo    = np.identity(ref_ci.scf.nmo, dtype=float)
        # iterate over the reference states
        for ref_obj in ref_ci:
            bra_st  = []
            bra_lbl = ref_obj.label
            for lbl, indx in self.ref_states[moleclue][bra_lbl].items():
                ipairs = np.asarray([[indx, st] for st in 
                                  range(new_ci.n_states())], dtype=int)
                Sij = overlap.overlap(ref_ci, new_ci, smo, pairs, 
                                                               1, 0.95)
                if np.amax(Sij) >= 0.9:
                    st = np.argmax(Sij)
                    eners[lbl] = new_ci.energies[st]

        if list(eners.keys()).sort() != list(ref_det.keys()).sort():
            sys.exit('Could not identify states: ' + str(ci_objs) + 
                     ' lbl:'+str(states))

        return eners

    #
    def parse_graci_file(self, target_vals):
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
        if os.path.isfile(self.graci_ref_file):
            ref_chkpt = h5py.File(self.graci_ref_file, 'r',
                                                   libver='latest')
        else:
            print('File: '+str(self.graci_ref_file)+' not found. '+
                  ' Exiting...')
            sys.exit(1)

        # get the top-level contents of the checkpoint file
        ref_contents = chkpt.contents(file_name = self.graci_ref_file)

        ci_objs  = {}
        scf_objs = {}
        st_indxs = {}
        for molecule in target_vals.keys():
            ci_objs[molecule]  = []
            scf_objs[molecule] = None
            st_indxs[molecule] = {}
            st_strs            = []
            for key in target_vals[molecule].keys():
                st_strs.extend(key.strip().split())
            st_lst = list(set(st_strs))

            # if the molecule string is in the reference 
            # object name, add it ref_objs dict
            for name in ref_contents:
                ci_type = any([ci in name for ci in params.ci_objs])
                if molecule in name:
                    if ci_type:
                        ci_objs[molecule].append(name)
                        # parse energy array to find target lbls
                        ci_lbl, st_lbls = self.read_state_labels(
                                                 ref_chkpt, name)
                        st_indxs[molecule][ci_lbl] = st_lbls
                    elif 'Scf' in name:
                        scf_objs[molecule] = name

        ref_ckhpt.close()

        return scf_objs, ci_objs, st_indxs

    #
    def read_state_labels(self, ref_chkpt, ci_name):
        """
        Read state labels for checkpoint file
        """

        # variables names are stored as attributes
        var_names = ref_chkpt[ci_name].attrs
        ener_dset = var_names['energies']
        attrs = ref_chkpt[ci_name+'/'+ener_dset].attrs

        ci_lbl = ref_chkpt[ci_name].attrs['label']

        if 'labels' not in attrs.keys():
            lbl_dic = dict()
        else:
            lbl_dic = dict(attrs['labels'])

        return ci_lbl, lbl_dic

    #
    def parse_target_file(self):
        """
        Parse the file containing the data we're going to parameterize
        wrt to

        Arguments: None
        Returns:   A dictionary containing the target data
        """

        with open(self.target_file, 'r') as t_file:
            t_file_lines = t_file.readlines()

        print('tfile='+str(t_file_lines))

        self.ndata  = 0
        target_dict = {}
        for line in t_file_lines:
            str_arr               = line.strip().split()
            molecule              = str_arr[0]
            target_dict[molecule] = {} 
            # format is: state1 state2 energy
            for i in range(1,len(str_arr),3):
                try:
                    states = str_arr[i].strip()+' '+str_arr[i+1].strip()
                    val    = parse.convert_value(str_arr[i+2])
                    target_dict[molecule][states] = val
                    self.ndata += 1
                except:
                    print('Error parsing value as str/float: '+str(val))
                    print('line = '+str(str_arr))
                    sys.exit(1)
        
        print('target_dict='+str(target_dict))
        return target_dict

