"""
The Parameterize object and its associated functions.
"""
import sys as sys
import scipy.optimize as sp_opt
import numpy as np
import h5py as h5py
import mpi4py as mpi4py
import os as os
import graci.core.libs as libs
import graci.core.driver as driver
import graci.core.params as params
import graci.io.chkpt as chkpt
import graci.io.parse as parse
import graci.io.output as output
import graci.utils.constants as constants
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
        self.verbose        = 1

        # -------------------------------
        self.target_data    = {}
        self.ref_states     = {}
        self.ndata          = 0
        self.iter           = 0
        self.current_h      = 0

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
        p0 = np.asarray([0.500779, 0.356986, 0.573523, 1.9266], 
                                            dtype=float)

        output.print_param_header(target_data, ci_data, 
                                            self.ref_states, p0)

        # set up the directory structure
        for target in target_data.keys():
            try:
                os.mkdir(target)
            except FileExistsError:
                os.rmdir(target)
                os.mkdir(target)

        # run the first pass
        ener_init = self.eval_energies(p0, target_data, scf_data, 
                                           ci_data, run_scf=True)

        self.iter      = 1
        self.current_h = p0
        res = sp_opt.minimize(self.err_func, p0, 
                              args = (target_data, scf_data, ci_data),
                              method = 'Nelder-Mead',
                              tol = self.pthresh, 
                              callback = self.status_func)

        ener_final = self.eval_energies(res.x, target_data, scf_data, 
                                               ci_data, run_scf=True)

        output.print_param_results(res, target_data, 
                                   ener_init, ener_final )

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

        iter_ener = self.eval_energies(hparams, target_data, 
                                         scf_data, ci_data)
        dif_vec   = np.zeros(self.ndata, dtype=float)

        # this approach assumes dict is ordered! Only true from Python
        # 3.6 onwards...
        ide = 0
        for molecule, states in target_data.items():
            for trans, ener in target_data[molecule].items():
                init, final  = trans.strip().split()
                de_iter      = iter_ener[molecule][final] - \
                               iter_ener[molecule][init]
                dif_vec[ide] = de_iter*constants.au2ev - ener
                ide         += 1

        return np.linalg.norm(dif_vec)

    #
    def status_func(self, xk):
        """
        Give status of optimization
        """

        output.print_param_iter(self.iter, list(xk), dif)

        self.iter     += 1
        self.current_h = xk

        return

    #
    def eval_energies(self, hparams, target_data, scf_names, ci_names, 
                                                       run_scf=False):
        """
        evaluate all the energies in the graci data set

        """

        # open the reference data file and get the contents
        if os.path.isfile(self.graci_ref_file):
            ref_chkpt = h5py.File(self.graci_ref_file, 'r',
                                                   libver='latest')
        else:
            print('File: '+str(self.graci_ref_file)+' not found. '+
                  ' Exiting...')
            sys.exit(1)

        topdir = os.getcwd()
        energies = {}

        # run through all the reference data objects
        for molecule, states in target_data.items():

            # pull the reference SCF and CI objects
            scf_name = scf_names[molecule]
            scf_obj  = chkpt.read(scf_name, file_handle=ref_chkpt)

            ci_name = ci_names[molecule]
            ci_objs  = [chkpt.read(ci, file_handle=ref_chkpt) 
                                                for ci in ci_name]

            # copy the scf and CI objects and re-run them
            scf_run = scf_obj.copy()
            ci_runs = [ci_obj.copy() for ci_obj in ci_objs]
         
            # change to the appropriate sub-directory
            os.chdir(topdir+'/'+str(molecule))

            # run the scf (and generate integral files if need be)
            if run_scf:
                # run silent
                scf_run.verbose = 0
                scf_run.load()

            # initialize the integrals
            if scf_run.mol.use_df:
                type_str = 'df'
            else:
                type_str = 'exact'
            libs.lib_func('bitci_int_initialize', ['pyscf', type_str, 
                            scf_run.moint_1e, scf_run.moint_2e_eri])

            # run all the CI objects 
            for ci_run in ci_runs:
                # run silent
                ci_run.verbose = 0
                ci_run.update_hparam(np.asarray(hparams))
                ci_run.run(scf_run, None)

            # deallocate the integral arrays
            libs.lib_func('bitci_int_finalize', [])

            # use overlap with ref states to identify relevant states
            ener_match = self.identify_states(molecule, ci_objs, ci_runs) 

            energies[molecule] = {}
            for state, ener in ener_match.items():
                energies[molecule][state] = ener

        os.chdir(topdir)
        ref_chkpt.close()

        return energies
 
    #
    def identify_states(self, molecule, ref_ci, new_ci):
        """
        compute overlaps to identify states
        """

        eners  = {}
        bra_st = []
        smo    = np.identity(ref_ci[0].scf.nmo, dtype=float)
        # iterate over the reference states
        for iobj in range(len(ref_ci)):
            bra_lbl = ref_ci[iobj].label
            st_dict = self.ref_states[molecule][bra_lbl]
            bra_st.extend(list(st_dict.keys()))

            for lbl, indx in st_dict.items():
                ipairs = np.asarray([[indx, st] for st in 
                            range(new_ci[iobj].n_states())], dtype=int)
                Sij = overlap.overlap(ref_ci[iobj], new_ci[iobj],
                            smo, ipairs, 0, 0.95)
                if np.amax(np.absolute(Sij)) >= 0.9:
                    st = np.argmax(np.absolute(Sij))
                    eners[lbl] = new_ci[iobj].energies[st]

        
        if list(eners.keys()).sort() != bra_st.sort():
            sys.exit('Molecule: ' + str(molecule) + 
                     ' -- Could not identify states: ' + str(bra_st))

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
                        lbl, states = self.read_state_labels(ref_chkpt,
                                                             name)
                        st_indxs[molecule][lbl] = states
                    elif 'Scf' in name:
                        scf_objs[molecule] = name

        ref_chkpt.close()

        return scf_objs, ci_objs, st_indxs

    #
    def read_state_labels(self, ref_chkpt, ci_name):
        """
        Read state labels for checkpoint file
        """

        ci_lbl = chkpt.read_attribute(ref_chkpt, ci_name, 'label')
        dset_path = chkpt.read_attribute(ref_chkpt, ci_name, 'energies')

        dset_name = ci_name+'/'+dset_path
        lbl_dic = chkpt.read_attribute(ref_chkpt, dset_name, 'label')

        # convert the values in the lbl_dic to integers
        for lbl,indx in lbl_dic.items():
            # convert from range 1,n to 0,n-1
            lbl_dic[lbl] = parse.convert_value(indx) - 1

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
        
        return target_dict

