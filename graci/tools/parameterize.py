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

        self.hamiltonian    = 'heil17_short'
        self.init_params    = None
        self.graci_ref_file = ''
        self.target_file    = ''
        self.pthresh        = 1.e-4
        self.verbose        = False 
        self.max_iter       = 1000

        # -------------------------------
        self.ndata          = 0
        self.iiter          = 0
        self.current_h      = 0
        self.error          = 0
        self.de_thr         = 0.5

    #
    def run(self):
        """re-parameterize the Hamiltonian"""

        # sanity check the input
        self.sanity_check()

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
        self.create_dirs(scf_objs)

        # the first pass sets up the orbitals and integrals -- no need
        # to recompute these every time
        ener_init = self.evaluate_energies(pref, ref_states, scf_objs,
                                                 ci_objs, runscf=True)

        # optimize parameters using scipy routines
        self.iiter      = 1
        self.current_h = p0
        res = sp_opt.minimize(self.err_func, p0, 
                              args = (target_data, ref_states, 
                                      scf_objs, ci_objs),
                              method = 'Nelder-Mead',
                              tol = self.pthresh,
                              callback = self.status_func)

        # one final eval_energy call with the converged params
        ener_final = self.evaluate_energies(res.x, ref_states,
                                                      scf_objs, ci_objs)

        output.print_param_results(res, target_data, ener_init,
                                                             ener_final)

        return

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

        return

    #
    def err_func(self, hparams, target_data, ref_states, scf_objs, 
                                                          ci_objs):
        """
        Evaluate the error function. In this case, simple RMSD

        Arguments:

        targets:      target values
        values:       current values

        Returns:
        norm of the different between target and current values
        """

        iter_ener = self.evaluate_energies(hparams, ref_states,
                                           scf_objs, ci_objs)
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
    def evaluate_energies(self, hparams, ref_states, scf_names,
                                              ci_names, runscf=False):
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
                        scf_names[molecule], ci_names[molecule], True,
                        runscf) for molecule in ref_states.keys())

                    for res in executor.starmap(self.eval_energy, args):
                        energies.update(res)

        else:

            for molecule in ref_states.keys():
                mol_results = self.eval_energy(hparams, molecule,
                                  topdir, self.graci_ref_file, 
                                  ref_states[molecule], 
                                  scf_names[molecule], 
                                  ci_names[molecule], 
                                  False, runscf)

                energies.update(mol_results)

        return energies

    def eval_energy(self, hparams, molecule, topdir, ref_file, ref_state, 
                            scf_name, ci_name, parallel, runscf):
        """
        eval_energy
        """

        # if parallel, this is run on a worker and libraries need to be 
        if parallel:
            libs.lib_load('bitci')
            libs.lib_load('overlap')

        ref_chkpt  = h5py.File(ref_file, 'r', libver='latest')

        # pull the reference SCF and CI objects
        scf_objs = []
        for iscf in range(len(scf_name)):
            if iscf == 0 or scf_name[iscf] != scf_name[iscf-1]:
                scf_objs.append(chkpt.read(scf_name[iscf], 
                                           file_handle=ref_chkpt))
 
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

        curr_scf = ''
        iscf     = -1
        for i in range(len(ci_objs)):
        
            if curr_scf != scf_name[i]:
                iscf += 1
                # change to the appropriate sub-directory
                os.chdir(topdir+'/'+str(molecule)+'/scf'+str(iscf))
 
                if runscf:
                    scf_objs[iscf].verbose = self.verbose 
                    scf_objs[iscf].load()

                if scf_objs[iscf].mol.use_df:
                    type_str = 'df'
                else:
                    type_str = 'exact'

                if i != 0:
                    libs.lib_func('bitci_int_finalize', [])
                libs.lib_func('bitci_int_initialize', ['pyscf', 
                             type_str, scf_objs[iscf].moint_1e, 
                             scf_objs[iscf].moint_2e_eri])

                curr_scf = scf_name[i]

            ci_objs[i].verbose = self.verbose 
            ci_objs[i].update_hparam(np.asarray(hparams))
            ci_objs[i].run(scf_objs[iscf], None)
            
        libs.lib_func('bitci_int_finalize', [])

        # use overlap with ref states to identify relevant states
        ener_match = self.identify_states(molecule, ref_state, ci_refs, 
                                                               ci_objs)
 
        # if we're unable to match the states output some info that
        # may help diagnose the problem
        if ener_match is None:
            msg = str(molecule)+' failed to match states. Current ' + \
                  ' hparam values = ' + str(hparams)
            for i in range(len(ci_objs)):
                msg += '\n ' + str(ci_objs[i].label) + ' ' + \
                       str(ci_objs[i].energies)
            self.hard_exit(msg)

        energies = {molecule : {}}
        for state, ener in ener_match.items():
            energies[molecule][state] = ener

        # always end in initial directory
        os.chdir(topdir)
        
        return energies


    @timing.timed
    def identify_states(self, molecule, ref_state, ref_ci, new_ci):
        """
        compute overlaps to identify states
        """
        #s_thrsh = 1./np.sqrt(2.)
        s_thrsh = 0.5

        eners  = {}
        smo    = np.identity(ref_ci[0].scf.nmo, dtype=float)
        # iterate over the reference states
        for iobj in range(len(ref_ci)):
            ref_st  = ref_state[iobj]

            bra_st  = list(ref_st.values())
            ket_st  = list(range(new_ci[iobj].n_states()))
            Smat = overlap.overlap_st(ref_ci[iobj], new_ci[iobj], 
                                    bra_st, ket_st, smo, 0.975, False)

            for lbl, bst in ref_st.items():
                Sij   = np.absolute(Smat[bra_st.index(bst),:])
                s_max = np.amax(Sij)

                if s_max <= s_thrsh:
                    print('ERROR: overlap for ' + str(molecule) +
                          ' state=' + str(lbl) +' is ' + str(Sij))
                    return None

                kst = ket_st[np.argmax(Sij)]
                eners[lbl] = new_ci[iobj].energies[kst]

        return eners

    #
    def create_dirs(self, scf_objs):
        """
        Create a directory structure for running computations
        """
        # set up the directory structure
        for molecule in scf_objs.keys():
            nscf = len(list(set(scf_objs[molecule])))

            if os.path.isdir(molecule):
                shutil.rmtree(molecule)

            os.mkdir(molecule)
            for iscf in range(nscf):
                os.mkdir(molecule+'/scf'+str(iscf))

        return

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

