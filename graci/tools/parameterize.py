"""
The Parameterize object and its associated functions.
"""
import sys as sys
import scipy.optimize as sp_opt
import numpy as np
import h5py as h5py
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor
import os as os
import shutil 
import graci.core.libs as libs
import graci.core.params as params
import graci.io.chkpt as chkpt
import graci.io.parse as parse
import graci.io.output as output
import graci.utils.constants as constants
import graci.interfaces.overlaptools.overlap as overlap

class Parameterize:
    """Class constructor for the Parameterize object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.algorithm      = 'Nelder-Mead'
        self.label          = 'parameterize'

        self.hamiltonian    = 'heil17_short'
        self.graci_ref_file = ''
        self.target_file    = ''
        self.pthresh        = 1.e-4
        self.verbose        = True
        self.max_iter       = 1000

        # -------------------------------
        self.ndata          = 0
        self.iiter          = 0
        self.current_h      = 0
        self.error          = 0

    #
    def run(self):
        """re-parameterize the Hamiltonian"""

        # sanity check the input
        self.sanity_check()

        # parse the target data file
        target_data = self.parse_target_file()

        # parse the graci reference file
        scf_data, ci_data, ref_states = \
                               self.parse_graci_file(target_data)

        # print the parameterize header
        # show the mappings between the target data and the data in the
        # reference checkpoint file
        np   = libs.lib_func('retrieve_nhpar', (self.hamiltonian, 0))

        p0   = np.zeros(np, dtype=float)
        args = (self.hamiltonian, np, p0)
        p0   = libs.lib_func('retrieve_hpar', args)

        output.print_param_header(target_data, ci_data, ref_states, p0)

        # set up the directory structure
        for target in target_data.keys():
            try:
                os.mkdir(target)
            except FileExistsError:
                os.rmdir(target)
                os.mkdir(target)
         
        # for time being, we're going to copy a separate reference file 
        # for each process. This is not optimal and could probably be
        # improved at a later date.
        if params.nproc > 1:
            for i in range(params.nproc):
                rfile = self.graci_ref_file
                shutil.copyfile(rfile, rfile+'.'+str(i))

        # the first pass sets up the reference objects
        ener_init = self.evaluate_energies(p0, target_data, ref_states, 
                                    scf_data, ci_data, runscf=True)

        self.iter      = 1
        self.current_h = p0
        res = sp_opt.minimize(self.err_func, p0, 
                              args = (target_data, ref_states, 
                                                  scf_data, ci_data),
                              method = self.algorithm,
                              tol = self.pthresh,
                              callback = self.status_func)

        ener_final = self.evaluate_energies(res.x, target_data, 
                                      ref_states, scf_data, ci_data)
        output.print_param_results(res, target_data, ener_init, 
                                                       ener_final)

        return

    #
    def sanity_check(self):
        """
        sanity check the input
        """

        algos = ['Nelder-Mead']
        if self.algorithm not in algos:
            sys.exit(str(self.algorithm)+' not in '+str(algos))


        if not os.path.isfile(self.graci_ref_file):
            print('graci_ref_file: ' + str(self.graci_ref_file) + 
                                             ' not found. Exiting')
            sys.exit(1)

        if not os.path.isfile(self.target_file):
            print('target_file: ' + str(self.target_file) +
                                             ' not found. Exiting')
            sys.exit(1)

        return

    #
    def err_func(self, hparams, target_data, ref_states, scf_data, 
                                                         ci_data):
        """
        Evaluate the error function. In this case, simple RMSD

        Arguments:

        targets:      target values
        values:       current values

        Returns:
        norm of the different between target and current values
        """

        iter_ener = self.evaluate_energies(hparams, 
                                         target_data, ref_states,
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
 
        self.error = np.linalg.norm(dif_vec)

        return self.error

    #
    def status_func(self, xk):
        """
        Give status of optimization
        """

        dif = np.linalg.norm(xk - self.current_h)
        output.print_param_iter(self.iter, list(xk), self.error)

        self.iiter     += 1
        self.current_h = xk

        if self.iiter >= self.max_iter:
            output.print_message('Max. number of iterations reached.')
            sys.exit(1)

        return

    #
    def evaluate_energies(self, hparams, target, ref_states, scf_names,
                                              ci_names, runscf = False):
        """
        evaluate all the energies in the graci data set

        """

        topdir = os.getcwd()
        energies = {}

        if params.nproc > 1:
            with MPIPoolExecutor(max_workers=params.nproc,
                             root=0) as executor:
                if executor is not None:

                    args = ((hparams, molecule, topdir, 
                            self.graci_ref_file, ref_states[molecule], 
                            scf_names[molecule], ci_names[molecule], 
                            True, runscf) for molecule in target.keys())

                    for results in executor.starmap(self.eval_energy, args):
                        energies.update(results)

        else:
            for molecule, states in target.items():
                mol_results = self.eval_energy(hparams, molecule, 
                                  topdir, self.graci_ref_file, 
                                  ref_states[molecule], 
                                  scf_names[molecule], 
                                  ci_names[molecule], 
                                  params.nproc, False, runscf)

                energies.update(mol_results)

        return energies

    #
    def eval_energy(self, hparams, molecule, topdir, ref_file, ref_state, 
                            scf_name, ci_name, parallel, runscf):
        """
        eval_energy
        """

        # if parallel, this is run on a worker and libraries need to be 
        if parallel:
            libs.lib_load('bitci')
            libs.lib_load('overlap')
            rank =  MPI.Comm.Get_parent().Get_rank()
            local_file = ref_file + '.' + str(rank)
        else:
            local_file = ref_file

        if os.path.isfile(local_file):
            ref_chkpt  = h5py.File(local_file, 'r', libver='latest')
        else:
            print('File: '+str(local_file)+' not found in ' + 
                                            str(topdir) + ' Exiting...')
            sys.exit(1)

        # pull the reference SCF and CI objects
        scf_obj  = chkpt.read(scf_name, file_handle=ref_chkpt)
        ci_refs  = [chkpt.read(ci, file_handle=ref_chkpt)
                                                for ci in ci_name]
        ref_chkpt.close()

        # change to the appropriate sub-directory
        os.chdir(topdir+'/'+str(molecule))

        # run the scf (and generate integral files if need be)
        if runscf:
            # run silent
            scf_obj.verbose = False
            scf_obj.load()

        # initialize the integrals
        if scf_obj.mol.use_df:
            type_str = 'df'
        else:
            type_str = 'exact'
        libs.lib_func('bitci_int_initialize', ['pyscf', type_str,
                        scf_obj.moint_1e, scf_obj.moint_2e_eri])

        ci_runs = [ci_ref.copy() for ci_ref in ci_refs]

        # run all the CI objects
        for ci_run in ci_runs:
            # run silent
            ci_run.verbose = False
            ci_run.update_hparam(np.asarray(hparams))
            ci_run.run(scf_obj, None)

        # deallocate the integral arrays
        libs.lib_func('bitci_int_finalize', [])

        # use overlap with ref states to identify relevant states
        ener_match = self.identify_states(molecule, ref_state, ci_refs, 
                                                                ci_runs)

        energies = {molecule : {}}
        for state, ener in ener_match.items():
            energies[molecule][state] = ener
       
        os.chdir(topdir)
         
        return energies


    #
    def identify_states(self, molecule, ref_state, ref_ci, new_ci):
        """
        compute overlaps to identify states
        """

        eners  = {}
        smo    = np.identity(ref_ci[0].scf.nmo, dtype=float)
        # iterate over the reference states
        for iobj in range(len(ref_ci)):
            bra_lbl = ref_ci[iobj].label
            st_dict = ref_state[bra_lbl]
 
            bra_st  = list(st_dict.values())
            ket_st  = list(range(new_ci[iobj].n_states()))
            Smat = overlap.overlap_st(ref_ci[iobj], new_ci[iobj], 
                                    bra_st, ket_st, smo, 0.95, False)

            for lbl, ref_st in st_dict.items():
                Sij = np.absolute(Smat[bra_st.index(ref_st),:])
                if np.amax(Sij) <= 0.9:
                    output.print_message('MAX overlap for molecule=' +
                            str(molecule) + ' state=' + str(lbl) +
                            ' is ' + str(np.amax(Sij)))
                st = ket_st[np.argmax(Sij)]
                eners[lbl] = new_ci[iobj].energies[st]
                
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

