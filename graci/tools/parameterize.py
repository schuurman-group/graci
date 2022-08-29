"""
The Parameterize object and its associated functions.
"""
import sys as sys
import math
import scipy.optimize as sp_opt
import numpy as np
import h5py as h5py
from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor
import os as os
import shutil as shutil
import graci.core.libs as libs
import graci.core.params as params
import graci.io.chkpt as chkpt
import graci.io.parse as parse
import graci.io.output as output
import graci.utils.constants as constants
import graci.utils.timing as timing
import graci.interfaces.overlap.overlap as overlap

class Parameterize:
    """Class constructor for the Parameterize object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.algorithm      = 'Nelder-Mead'
        self.label          = 'parameterize'

        self.hamiltonian    = 'heil17_standard'
        self.init_params    = None
        self.graci_ref_file = ''
        self.target_file    = ''
        self.pthresh        = 1.e-4
        self.verbose        = False 
        self.max_iter       = 1000
        self.bounds         = None

        # -------------------------------
        self.ndata          = 0
        self.iiter          = 0
        self.current_h      = 0
        self.error          = 0
        self.de_thr         = 0.5
        self.log_file        = None

    #
    def run(self):
        """re-parameterize the Hamiltonian"""

        # sanity check the input
        hbounds = self.sanity_check()

        # parse the target data file
        target_data = self.parse_target_file()

        # parse the graci reference file
        scf_objs, ci_objs, ref_states = \
                               self.parse_graci_file(target_data)

        # set initial parameter set
        pref, p0 = self.set_init_params()
 
        # print header info to output file defined in module output
        output.print_param_header(target_data, ci_objs, ref_states, p0)

        # set up directory structure: each molecule gets a subdirectory
        scf_dirs = self.create_dirs(scf_objs)

        # the first pass sets up the orbitals and integrals -- no need
        # to recompute these every time
        ener_init = self.evaluate_energies(pref, ref_states, scf_dirs, 
                                         scf_objs, ci_objs, runscf=True)

        # save the default logfile name for writing updates of the
        # reparam procedure
        self.logfile = output.file_names['out_file']

        # optimize parameters using scipy routines
        self.iiter      = 1
        self.current_h = p0
        res = sp_opt.minimize(self.err_func, p0, 
                              args = (target_data, ref_states, 
                                      scf_dirs, scf_objs, ci_objs),
                              bounds = hbounds,
                              method = 'Nelder-Mead',
                              tol = self.pthresh,
                              callback = self.status_func)

        # one final eval_energy call with the converged params
        ener_final = self.evaluate_energies(res.x, ref_states, scf_dirs,
                                                      scf_objs, ci_objs)

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
            msg = str(self.algorithm)+' not in '+str(algos)
            self.hard_exit(msg)

        if not os.path.isfile(self.graci_ref_file):
            msg = 'graci_ref_file: ' + str(self.graci_ref_file) + \
                                             ' not found. Exiting'
            self.hard_exit(msg)

        if not os.path.isfile(self.target_file):
            msg = 'target_file: ' + str(self.target_file) + \
                                             ' not found. Exiting'
            self.hard_exit(msg)

        if self.bounds is not None:
            npar = 0
            args = (self.hamiltonian, npar)
            npar = libs.lib_func('retrieve_nhpar', args)

            if self.bounds.shape[0] != npar:
                msg = 'Length of bounds array != number of H parameters'
                self.hard_exit(msg)

            if self.bounds.shape[1] != 2:
                msg = 'bounds.shape[1] != 2: specify upr and lwr bounds'
                self.hard_exit(msg)

            hmin = np.min([self.bounds[i,0] for i in range(npar)])
            if hmin <= 0.:
                msg = 'lower bound for H[param] < 0: '+str(hmin)
                output.print_message(msg)

            bounds = sp_opt.Bounds(lb=self.bounds[:,0],
                                   ub=self.bounds[:,1])
        else:
            bounds = None

        return bounds

    #
    def set_init_params(self):
        """
        set the initial parameter values, either using default 
        Hamiltonian values, or, user supplied values
        """

        npar = 0
        args = (self.hamiltonian, npar)
        npar = libs.lib_func('retrieve_nhpar', args)

        # set initial parameter set
        pref = np.zeros(npar, dtype=float)
        args = (self.hamiltonian, npar, pref)
        pref = libs.lib_func('retrieve_hpar', args)

        if self.init_params is None:
            p0 = pref
        else:
            p0 = self.init_params
            if (type(p0).__name__ != 'ndarray' or p0.shape[0] != npar):
                msg = 'Cannot use '+str(p0)+' as initial param set'
                self.hard_exit(msg)

        return pref, p0

    #
    def err_func(self, hparams, target_data, ref_states, scf_dirs, 
                                                 scf_objs, ci_objs):
        """
        Evaluate the error function. In this case, simple RMSD

        Arguments:

        targets:      target values
        values:       current values

        Returns:
        norm of the different between target and current values
        """

        iter_ener = self.evaluate_energies(hparams, ref_states,
                                           scf_dirs, scf_objs, ci_objs)
        dif_vec   = np.zeros(self.ndata, dtype=float)

        # this approach assumes dict is ordered! Only true from Python
        # 3.6 onwards...
        ide        = 0
        for molecule in target_data.keys():
            for trans, ener in target_data[molecule].items():
                init, final  = trans.strip().split()
                de_iter      = iter_ener[molecule][final] - \
                               iter_ener[molecule][init]
                dif_vec[ide] = de_iter*constants.au2ev - ener
                ide         += 1
 
        self.error    = np.linalg.norm(dif_vec)

        return self.error

    #
    def status_func(self, xk):
        """
        Give status of optimization
        """

        dif = np.linalg.norm(xk - self.current_h)
        output.print_param_iter(self.iiter, list(xk), self.error)

        self.iiter     += 1
        self.current_h = xk

        if self.iiter >= self.max_iter:
            msg = 'Max. number of iterations completed.'
            self.hard_exit(msg) 

        return

    @timing.timed
    def evaluate_energies(self, hparams, ref_states, scf_dirs,
                                   scf_names, ci_names, runscf=False):
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
                         scf_dirs[molecule], scf_names[molecule], 
                         ci_names[molecule], True, runscf)
                                   for molecule in ref_states.keys())

                    for res in executor.starmap(self.eval_energy, args):
                        energies.update(res)

        else:

            for molecule in ref_states.keys():
                mol_results = self.eval_energy(hparams, molecule,
                                  topdir, self.graci_ref_file, 
                                  ref_states[molecule],
                                  scf_dirs[molecule], 
                                  scf_names[molecule], 
                                  ci_names[molecule], 
                                  False, runscf)

                energies.update(mol_results)

        return energies

    def eval_energy(self, hparams, molecule, topdir, ref_file, ref_state, 
                           scf_dirs, scf_name, ci_name, parallel, runscf):
        """
        eval_energy
        """

        # if parallel, this is run on a worker and libraries need to be 
        if parallel:
            libs.lib_load('bitci')
            libs.lib_load('overlap')

        ref_chkpt  = h5py.File(ref_file, 'r', libver='latest')

        # pull the reference SCF and CI objects
        scf_objs = {}
        for scf_link in scf_name:
            if scf_link not in scf_objs.keys():
                scf_objs[scf_link] = chkpt.read(scf_link, 
                                                file_handle=ref_chkpt)
 
        # these are reference ci objects: they don't change and
        # are used for state comparison
        ci_refs   = [chkpt.read(ci, file_handle=ref_chkpt)
                                                for ci in ci_name]

        # these are the ci objects that get run with the current
        # values of the Hamiltonian parameters
        ci_objs = [ci_ref.copy() for ci_ref in ci_refs]
        ref_chkpt.close()

        # set the output file stream if we're going to do any printing
        output.file_names['out_file'] = topdir+'/'+str(molecule)+'.log'

        # change to the appropriate sub-directory
        os.chdir(topdir+'/'+str(molecule))

        curr_dir  = ''
        energies  = {molecule : {}}

        for i in range(len(ci_objs)):
        
            if curr_dir != scf_dirs[i]:

                # change to the appropriate sub-directory
                os.chdir(topdir+'/'+str(molecule)+'/'+scf_dirs[i])

                if runscf:
                    scf_objs[scf_name[i]].verbose = self.verbose 
                    scf_objs[scf_name[i]].load()

                if scf_objs[scf_name[i]].mol.use_df:
                    type_str = 'df'
                else:
                    type_str = 'exact'

                if curr_dir != '':
                    libs.lib_func('bitci_int_finalize', [])

                libs.lib_func('bitci_int_initialize', ['pyscf', 
                             type_str, scf_objs[scf_name[i]].moint_1e, 
                             scf_objs[scf_name[i]].moint_2e_eri])

                curr_dir = scf_dirs[i]

            ci_objs[i].verbose = self.verbose 
            roots_found = False
            i_add       = 0
            n_add       = 5
            add_states  = np.ceil(2*ci_objs[i].nstates / n_add).astype(int)  
 
            while not roots_found and i_add < n_add:
                # only update the Hamiltonian parameters if the CI object
                # is employing the Hamiltonian we're optimizing
                if ci_objs[i].hamiltonian == self.hamiltonian:
                    ci_objs[i].update_hparam(np.asarray(hparams))
                ci_objs[i].run(scf_objs[scf_name[i]], None)

                # use overlap with ref states to identify relevant states
                roots_found, eners = self.identify_states(molecule, 
                                   ref_state[i], ci_refs[i], ci_objs[i])

                if not roots_found:
                    ci_objs[i].nstates += add_states
                    i_add              += 1
            
            # if all roots found, update the energy directory 
            if roots_found:        
                energies[molecule].update(eners) 

            # else, hard exit
            else:
                msg = str(molecule)+' failed to match states. Current ' + \
                      ' hparam values = ' + str(hparams)
                msg += '\n ' + str(ci_objs[i].label) + ' ' + \
                        str(ci_objs[i].energies)
                self.hard_exit(msg)

        libs.lib_func('bitci_int_finalize', [])

        # always end in initial directory
        os.chdir(topdir)
       
        # reset the output to the default logfile
        output.file_names['out_file'] = self.log_file
 
        return energies


    @timing.timed
    def identify_states(self, molecule, ref_st, ref_ci, new_ci):
        """
        compute overlaps to identify states
        """
        #s_thrsh = 1./np.sqrt(2.)
        s_thrsh = 0.5

        eners  = {}
        smo    = np.identity(ref_ci.scf.nmo, dtype=float)

        # iterate over the reference states
        bra_st  = list(ref_st.values())
        ket_st  = list(range(new_ci.n_states()))
        Smat = overlap.overlap_st(ref_ci, new_ci, bra_st, ket_st, smo,
                                                         0.975, False)

        roots_found = True
        for lbl, bst in ref_st.items():
            Sij   = np.absolute(Smat[bra_st.index(bst),:])
            s_max = np.amax(Sij)

            if s_max <= s_thrsh:
                roots_found = False

            kst = ket_st[np.argmax(Sij)]
            eners[lbl] = new_ci.energies[kst]

        return roots_found, eners

    #
    def create_dirs(self, scf_objs):
        """
        Create a directory structure for running computations
        """

        scf_dirs = {}

        # set up the directory structure
        for molecule, scf_names in scf_objs.items():

            if os.path.isdir(molecule):
                shutil.rmtree(molecule)
            os.mkdir(molecule)

            idir               = -1
            name_set           = []
            scf_dirs[molecule] = [] 

            for iscf in range(len(scf_names)):

                if scf_names[iscf] in name_set:
                    indx = scf_names.index(scf_names[iscf])
                    scf_dirs[molecule].append(scf_dirs[molecule][indx])

                else:
                    idir += 1
                    dir_name = 'scf' + str(idir)
                    os.mkdir(molecule + '/' + dir_name)
                    name_set.append(scf_names[iscf])
                    scf_dirs[molecule].append(dir_name)          
                    
        return scf_dirs

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
        ref_file = h5py.File(self.graci_ref_file, 'r', libver='latest')

        # get the top-level contents of the checkpoint file
        ref_contents = chkpt.contents(file_handle = ref_file)

        ci_objs    = {}
        scf_objs   = {}
        ref_states = {}
        for molecule in target_vals.keys():

            # get the list of state labels we're looking for
            target_lbls = []
            for lbls,ener in target_vals[molecule].items():
                target_lbls.extend(lbls.strip().split())            
            target_lbls = list(set(target_lbls))
            target_lbls.sort()

            ci_objs[molecule]    = []
            scf_objs[molecule]   = []
            ref_states[molecule] = []
            states_found         = []

            # if the molecule string is in the reference 
            # object name, add it ref_objs dict
            for name in ref_contents:
                is_ci  = any([ci in name for ci in params.ci_objs])
                st_map = chkpt.read_attribute(ref_file, name, 
                                                       'state_map')
                
                # name of section is 'original.name.molecule'
                m_name = name.strip().split('.')[-1]
                if molecule == m_name and is_ci and st_map is not None:
                    ci_objs[molecule].append(name)
                    ref_states[molecule].append(st_map)
                    states_found.extend(list(st_map.keys()))
                    states_found.sort()

                    # get the name of corresponding scf object
                    scf_link = chkpt.read_attribute(ref_file, name, 
                                                             'scf')
                    is_scf, scf_name = chkpt.data_name('', scf_link)
                    if is_scf:
                        scf_objs[molecule].append(scf_name)

            if target_lbls != states_found:
                msg = 'ERROR - molecule: ' + str(molecule) + ' -- ' +\
                       str(target_lbls) + ' != ' + str(states_found)
                self.hard_exit(msg)
                    
        ref_file.close()

        return scf_objs, ci_objs, ref_states

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
                    msg = 'Error parsing value as str/float: ' + \
                           str(val) + '\nline = '+str(str_arr)
                    self.hard_exit(msg)
        
        return target_dict

    #
    def hard_exit(self, message):
        """
        exit the program with a sys.exit call
        """
        print("\n"+str(message)+"\n", flush=True)
        sys.exit(1)

