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
        scf_objs, ci_objs, ref_states = \
                               self.parse_graci_file(target_data)

        # print the parameterize header
        # show the mappings between the target data and the data in the
        # reference checkpoint file
        nparam   = 0
        args     = (self.hamiltonian, nparam)
        nparam   = libs.lib_func('retrieve_nhpar', args)

        p0   = np.zeros(nparam, dtype=float)
        args = (self.hamiltonian, nparam, p0)
        p0   = libs.lib_func('retrieve_hpar', args)

        output.print_param_header(target_data, ci_objs, ref_states, p0)

        # set up directory structure: each molecule gets a subdirectory
        self.create_dirs(scf_objs)

        # the first pass sets up the reference objects
        ener_init = self.evaluate_energies(p0, ref_states, scf_objs,
                                                 ci_objs, runscf=True)

        self.iiter      = 1
        self.current_h = p0
        res = sp_opt.minimize(self.err_func, p0, 
                              args = (target_data, ref_states, 
                                      scf_objs, ci_objs),
                              method = 'Nelder-Mead',
                              tol = self.pthresh,
                              callback = self.status_func)

        ener_final = self.evaluate_energies(res.x, ref_states,
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
        ide = 0
        for molecule in target_data.keys():
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
        output.print_param_iter(self.iiter, list(xk), self.error)

        self.iiter     += 1
        self.current_h = xk

        if self.iiter >= self.max_iter:
            output.print_message('Max. number of iterations reached.')
            sys.exit(1)

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
                            scf_names[molecule], ci_names[molecule], 
                            True, runscf) for molecule in ref_states.keys())

                    for results in executor.starmap(self.eval_energy, 
                                                                  args):
                        energies.update(results)

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

        # change to the appropriate sub-directory
        os.chdir(topdir+'/'+str(molecule))

        curr_scf = ''
        iscf     = -1
        for i in range(len(ci_objs)):
        
            if curr_scf != scf_name[i]:
                iscf += 1
                os.chdir(topdir+'/'+str(molecule)+'/scf'+str(iscf))
 
                if runscf:
                    scf_objs[iscf].verbose = False
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

            ci_objs[i].verbose = False
            ci_objs[i].update_hparam(np.asarray(hparams))
            ci_objs[i].run(scf_objs[iscf], None)

        libs.lib_func('bitci_int_finalize', [])

        # use overlap with ref states to identify relevant states
        ener_match = self.identify_states(molecule, ref_state, ci_refs, 
                                                               ci_objs)

        energies = {molecule : {}}
        for state, ener in ener_match.items():
            energies[molecule][state] = ener
      
        os.chdir(topdir)
         
        return energies


    @timing.timed
    def identify_states(self, molecule, ref_state, ref_ci, new_ci):
        """
        compute overlaps to identify states
        """

        eners  = {}
        smo    = np.identity(ref_ci[0].scf.nmo, dtype=float)
        # iterate over the reference states
        for iobj in range(len(ref_ci)):
            ref_st  = ref_state[iobj]

            bra_st  = list(ref_st.values())
            ket_st  = list(range(new_ci[iobj].n_states()))
            Smat = overlap.overlap_st(ref_ci[iobj], new_ci[iobj], 
                                    bra_st, ket_st, smo, 0.95, False)

            for lbl, bst in ref_st.items():
                Sij = np.absolute(Smat[bra_st.index(bst),:])
                if np.amax(Sij) <= 0.9:
                    output.print_message('MAX overlap for molecule=' +
                            str(molecule) + ' state=' + str(lbl) +
                            ' is ' + str(np.amax(Sij)))
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
                os.shutil.rmtree(molecule)

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

            # if the molecule string is in the reference 
            # object name, add it ref_objs dict
            for name in ref_contents:
                is_ci  = any([ci in name for ci in params.ci_objs])
                st_map = chkpt.read_attribute(ref_file, name, 
                                                       'state_map')
                 
                if molecule in name and is_ci and st_map is not None:
                    ci_objs[molecule].append(name)
                    ref_states[molecule].append(st_map)

                    # get the name of corresponding scf object
                    scf_link = chkpt.read_attribute(ref_file, name, 
                                                             'scf')
                    is_scf, scf_name = chkpt.data_name('', scf_link)
                    if is_scf:
                        scf_objs[molecule].append(scf_name)

            # check that we found all the reference states
            states_found = []
            for st_dic in ref_states[molecule]:
                states_found.extend(list(st_dic.keys()))
            states_found.sort()

            if target_lbls != states_found:
                print('molecule: ' + str(molecule) + ' -- ' + 
                       str(target_lbls) + ' != ' + str(states_found), 
                       flush=True)
                sys.exit(1)
                    
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
                    print('Error parsing value as str/float: '+str(val))
                    print('line = '+str(str_arr))
                    sys.exit(1)
        
        return target_dict
