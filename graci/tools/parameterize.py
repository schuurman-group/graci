"""
The Parameterize object and its associated functions.
"""
import sys as sys
import math
import scipy.optimize as sp_opt
import numpy as np
import h5py as h5py
import os as os
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

import shutil as shutil
import graci.core.libs as libs
import graci.core.params as params
import graci.core.ao2mo as ao2mo
import graci.io.chkpt as chkpt
import graci.io.parse as parse
import graci.io.output as output
import graci.utils.constants as constants
import graci.utils.timing as timing
import graci.interfaces.overlap.overlap as overlap
import graci.methods.dftmrci2 as dftmrci2

class Parameterize:
    """Class constructor for the Parameterize object."""
    def __init__(self):
        # the following is determined from user input 
        # (or subject to user input) -- these are keywords
        # in params module
        self.label           = 'default'
        self.job_type        = 'opt'
        self.wfn_lib         = ''
        self.energy_lib      = ''
        self.verbose         = False
        self.opt_algorithm   = 'nelder-mead'
        self.opt_target      = 'rmsd'
        self.conv            = 0.01 
        self.max_iter        = 1000
        self.method          = 'dftmrci'
        self.xc              = 'qtp17'

        #  hamiltonian to use for valence state calculation
        self.hamiltonians    = []
        #  parameters to use for valence state hamiltonian
        self.params          = []
        #  do we optimize the parameters for of the hamiltonian
        self.opt             = []
        # which parameters to freeze during optimization
        self.freeze          = [[]]
        # bounds for ham parameters
        self.bounds          = [] 

        # ----------------------------------------------------------
        # this is the internal dictionary parameterize uses to keep
        # track of optimization options
        self.opt_options   = {}
        # bounds put into format the optimizer will accept
        self.opt_bnds      = None
        # initial parameters -- saved for finall comparison
        self.p_0           = None
        # list of parameters to optimize
        self.p_n           = None

        # scan run options
        self.ngrid           = None
        self.scan_var        = []

        # number of worker proceses for parallel runs
        self.max_workers     = 1

        self.n_opt         = 0
        self.n_ref         = 0
        self.iiter         = 0
        self.error         = 0
        self.de_thr        = 0.5
        self.log_file      = None
        self.valid_algos   = ['nelder-mead','differentialevolution']
        self.valid_opt_targ = ['mae','rmsd']

    #
    def run(self):
        """re-parameterize the Hamiltonian"""

        print('cpu count='+str(mp.cpu_count()))

        # parse the target data file
        exc_ref, states  = self.parse_ref_file()

        # parse the graci reference file
        scf_objs, ci_objs = self.parse_wfn_file(states)

        # sanity check the input
        self.sanity_check(ci_objs)

        # set initial parameter set
        self.n_opt, self.opt_options, self.opt_bnds = self.set_init_params()
        self.p_0 = self.extract_opt_param()

        # print header info to output file defined in module output
        output.print_param_header(self.opt_options, exc_ref, ci_objs)

        # allocate number of workers/threads
        self.allocate_workers(ci_objs)

        # set up directory structure: each molecule gets a subdirectory
        # and each scf calculation gets a sub-sub-directory
        scf_dirs = self.create_dirs(scf_objs)

        # the first pass sets up the orbitals and integrals -- no need
        # to recompute these every time
        exc_init = self.evaluate_energies(scf_dirs, scf_objs, ci_objs,
                                                        gen_orbs=True)

        # save the default logfile name for writing updates of the
        # reparam procedure
        self.logfile = output.file_names['out_file']

        if self.job_type == 'opt':

            # optimize parameters using scipy routines
            self.iiter = 1
            self.p_n   = self.extract_opt_param()

            if self.opt_algorithm != 'differentialevolution':
                res = sp_opt.minimize(self.err_func, self.p_n, 
                         args = (exc_ref, scf_dirs, scf_objs, ci_objs),
                         bounds   = self.opt_bnds,
                         method   = self.opt_algorithm,
                         tol      = self.conv,
                         callback = self.status_func)
            else:
                res = sp_opt.differential_evolution(
                         self.err_func, self.opt_bnds,
                         args = (exc_ref, scf_dirs, scf_objs, ci_objs),
                         callback = self.status_func,
                         polish   = False,
                         tol      = self.conv,
                         x0       = self.p_n)

            self.update_opt_params(res.x)
            # one final eval_energy call with the converged params
            exc_final = self.evaluate_energies(scf_dirs, scf_objs,
                                                           ci_objs)

            output.print_param_results(self.opt_options, res, exc_ref, 
                                                    exc_init, exc_final)
        elif self.job_type == 'scan':
            args = (exc_ref, scf_dirs, scf_objs, ci_objs)
            self.scan(self.p_0, self.scan_var, self.ngrid, args)

        elif self.job_type == 'analysis':
            output.print_param_analysis(self.p_0, exc_ref, exc_init)

        return

    #
    def set_init_params(self):
        """
        set the initial parameter values, either using default 
        Hamiltonian values, or, user supplied values
        """
        n_opt = 0
        opt_options = {}

        # Load the parameter dictionaries 
        #-----------------------------------------------------------
        for index in range(len(self.hamiltonians)):
            ham = self.hamiltonians[index]

            opt_options[ham] = {}
            opt_options[ham]['optimize'] = self.opt[index]
            opt_options[ham]['params']   = self.params[index]
            opt_options[ham]['freeze']   = self.freeze[index]
            opt_options[ham]['bounds']   = self.bounds[index]
            if opt_options[ham]['optimize']:
                n_opt += len(opt_options[ham]['params']) - \
                         len(opt_options[ham]['freeze'])

        # Set the optimization bounds
        # ---------------------------------------------------------
        opt_bnds = np.zeros( (n_opt, 2), dtype=float)
        n_set = 0
        for index in range(len(self.hamiltonians)):
            ham = self.hamiltonians[index]

            if not opt_options[ham]['optimize']:
                continue

            bnd = [val for i, val in
                       enumerate(opt_options[ham]['bounds'])
                       if i not in opt_options[ham]['freeze']]
            opt_bnds[n_set:n_set + len(bnd), :] = np.array(bnd, 
                                                    dtype=float)
            n_set += len(bnd)

        # Make sure CI method is lower case
        # ----------------------------------------------------------
        self.method = self.method.lower()

        return n_opt, opt_options, opt_bnds 

    #
    def scan(self, p_init, scan_var, ngrid, args):
        """
        scan the parameter values 
        """

        bounds = []
        labels = []
        p_scan = []
        n_scan = 0

        # start by freezing all coordiantes
        for ham in self.hamiltonians:
            self.opt_options[ham]['freeze'] = \
                    list(range(len(self.opt_options[ham]['params'])))

        # loop over strings in scan_var
        rm_ind = []
        for p_str in self.scan_var:
            p_index = int(p_str[-1])
            rm_ind.append(p_index)
            ham     = p_str[:-1]
      
            if ham not in self.opt_options.keys():
                msg = 'Hamiltonian: '+str(ham)+' not recognized.'
                self.hard_exit(msg)

            #self.opt_options[ham]['freeze'].pop(p_index)
            labels.append(ham)
            bounds.append(self.opt_options[ham]['bounds'][p_index])
            p_scan.append(self.opt_options[ham]['params'][p_index])
            n_scan += 1

        # unfreeze the scan indices
        rm_ind.sort(reverse=True)
        for i in range(len(rm_ind)):
            self.opt_options[ham]['freeze'].pop(rm_ind[i])

        delta = [(bounds[i][1] - bounds[i][0]) / (ngrid[i]-1)
                  if ngrid[i] > 1 else 0. for i in range(n_scan)]

        output.print_param_scan_head(labels, p_scan, bounds, ngrid)
       
        hscan    = np.zeros(n_scan, dtype=float)
        step     = [0] * n_scan
        step[-1] = -1
        done = False
        while not done:
            
            param = n_scan - 1
            while step[param] == (ngrid[param]-1):
                step[param] = 0
                param      -= 1
                if param < 0:
                    done = True
                    break

            if done:
                break

            step[param] += 1
            for i in range(n_scan):
                hscan[i] = bounds[i][0] + step[i]*delta[i]

            err = self.err_func(hscan, *args) 
            output.print_param_scan_iter(hscan, step, err)

        return

    #
    def err_func(self, p_opt, exc_ref, scf_dirs, scf_objs, ci_objs):
        """
        Evaluate the error function. In this case, simple RMSD

        Arguments:

        targets:      target values
        values:       current values

        Returns:
        norm of the different between target and current values
        """

        # set the hamiltonian parameters: both frozen and optimized 
        #p_full  = self.to_full_param_set(p_opt)
        self.update_opt_params(p_opt)
        ener_i  = self.evaluate_energies(scf_dirs, scf_objs, ci_objs)
        dif_vec = np.zeros(self.n_ref, dtype=float)

        # this approach assumes dict is ordered! Only true from Python
        # 3.6 onwards...
        nde        = 0
        for molecule in exc_ref.keys():
            for trans, ener in exc_ref[molecule].items():
                init, final  = trans.strip().split()
                de_iter      = ener_i[molecule][final] - \
                               ener_i[molecule][init]
                dif_vec[nde] = de_iter*constants.au2ev - ener
                nde         += 1

        if self.opt_target == 'rmsd':
            self.error = np.linalg.norm(dif_vec) / np.sqrt(nde)
        elif self.opt_target == 'mae':
            self.error = 0.
            if nde > 0:
                self.error = np.sum(np.absolute(dif_vec)) / nde

        return self.error
 
    #
    @timing.timed
    def evaluate_energies(self, scf_dirs, scf_names, ci_objs,
                                                      gen_orbs=False):
        """
        evaluate all the energies in the graci data set

        """
        # init_ci and final_ci should have same key list
        mol_names = ci_objs.keys()
        topdir    = os.getcwd()
        energies  = {}

        if self.max_workers > 1:

            with ProcessPoolExecutor(max_workers=self.max_workers,
                      mp_context=mp.get_context("spawn")) as executor:

                if executor is not None:

                    args = ((molecule, topdir, self.wfn_lib,
                        scf_dirs[molecule], scf_names[molecule],         
                        ci_objs[molecule], True, gen_orbs)
                        for molecule in mol_names)

                    for res in executor.map(self.eval_energy, *zip(*args)):
                        energies.update(res)

        else:

            for molecule in mol_names:
                mol_results = self.eval_energy(
                             molecule, topdir, self.wfn_lib, 
                             scf_dirs[molecule], scf_names[molecule], 
                             ci_objs[molecule], False, gen_orbs)

                energies.update(mol_results)

        return energies

    # 
    def eval_energy(self, molecule, topdir, wfn_file, scf_dirs,
                         scf_name, ci_names, parallel, gen_orbs):
        """
        eval_energy
        """

        #print(molecule+' energy on rank='+str(MPI.COMM_WORLD.Get_rank()),flush=True)
        # if parallel, this is run on a worker and libraries need to be 
        if parallel:
            libs.lib_load('bitci')
            libs.lib_load('overlap')

        wfn_chkpt  = h5py.File(wfn_file, 'r', libver='latest')
        output.file_names['out_file'] = topdir+'/'+str(molecule)+'.log'
        mol_dir = topdir+'/'+str(molecule) 

        os.chdir(mol_dir)

        # set PYSCF_TMPDIR to current directory
        os.environ['PYSCF_TMPDIR'] = mol_dir

        energies  = {molecule : {}}
        mo_ints   = ao2mo.Ao2mo()
        
        # identify the CI objects
        scf_objs  = {ci_name:None for ci_name in ci_names.keys()}
 
        # loop over unique ci state objects
        for ci_name in ci_names.keys():     

            # ci_ref contains the wfns corresponding to the target
            # states
            ci_ref = chkpt.read(ci_name, file_handle=wfn_chkpt)
            
            # check the type of CI object. If same as method, no 
            # need to create new one
            ci_type = str(ci_ref.__class__.__name__).lower()

            if ci_type == self.method:
                ci_opt = ci_ref.copy()
            else:
                ci_class = self.method.capitalize()
                ci_opt = getattr(globals()[ci_class.lower()], ci_class)(ci_ref)

            # determine if we need to update, or set, hamiltonian parameters
            ham_name = ci_opt.hamiltonian
            if ham_name in self.opt_options.keys():
                params = np.array(self.opt_options[ham_name]['params'], 
                                                           dtype=float)
                ci_opt.update_hparam(params)

            # set verbosity to match requested output level
            ci_opt.verbose = self.verbose

            # move into the directory with the correspoinding scf
            # object (and MO integrals)
            os.chdir(mol_dir+'/'+scf_dirs[ci_name])

            # if scf_objs is None that means we need to generate it...
            if scf_objs[ci_name] is None:

                #...either by re-running it b/c it's the first time
                # function is called, or, b/c we're optimizing the
                # functional
                if gen_orbs:
                    scf_obj = chkpt.read(scf_name[ci_name],
                                         file_handle=wfn_chkpt)
                    scf_obj.verbose = self.verbose
                    scf_obj.load()
                    scf_obj.xc = self.xc
                    scf_obj.run(scf_obj.mol, None)
                    scf_objs[ci_name] = scf_obj

                    mo_ints.emo_cut      = ci_opt.mo_cutoff
                    mo_ints.precision_2e = ci_opt.precision
                    mo_ints.run(scf_obj)
    
                    # if only running scf once, save the orbs to file
                    fname = 'TMP_'+scf_name[ci_name]+'.chkpt.h5'
                    chkpt.write(scf_obj, file_name=fname, 
                                        grp_name = scf_name[ci_name])

                #...or by loading the scf object from a temporary
                # chkpt file
                else:
                    fname = 'TMP_'+scf_name[ci_name]+'.chkpt.h5'
                    scf_obj = chkpt.read(scf_name[ci_name],
                                                   file_name=fname)
                    scf_obj.verbose = self.verbose
                    scf_obj.load()
                    scf_objs[ci_name] = scf_obj

                    # set the mo_cutoff to ensure orb count is correct
                    mo_ints.emo_cut      = ci_opt.mo_cutoff
                    mo_ints.precision_2e = ci_opt.precision
                    mo_ints.load_bitci(scf_obj)

            ci_opt.update_eri(mo_ints)
            all_found  = False
            i_add      = 0
            n_add      = 5 
            add_states = np.asarray([1]*len(ci_opt.nstates), dtype=int)

            # iterate a couple times, in case we need to add roots to
            # find the states of interest
            while not all_found and i_add < n_add:

                ci_opt.run(scf_objs[ci_name], None, mo_ints = mo_ints)

                # use overlap with ref states to identify t states
                roots_found, eners = self.identify_states(molecule, 
                                     ci_ref, ci_names[ci_name], ci_opt)

                all_found = all(roots_found.values())
                if not all_found:
                    ci_opt.nstates += add_states
                    i_add          += 1
         
            energies[molecule].update(eners)

            # give a heads up that we failed to match roots:
            if not all_found:
                fail = [st for st,fnd in roots_found.items() if not fnd]
                msg  = str(molecule) + ' failed to match states: '
                msg += str(fail)  + ' h_param = '+str(ci_opt.hparam)
                output.print_message(msg)

        libs.lib_func('bitci_int_finalize', [])
        wfn_chkpt.close()

        # always end in initial directory
        os.chdir(topdir)
       
        # reset the output to the default logfile
        output.file_names['out_file'] = self.log_file
 
        return energies


    @timing.timed
    def identify_states(self, molecule, ci_ref, ref_states, ci_new):
        """
        compute overlaps to identify states
        """
        #s_thrsh = 1./np.sqrt(2.)
        s_thrsh = 0.6

        eners  = {}

        # in functional optimization runs, this will be something
        # other than identity matrix
        smo = ci_new.scf.mo_overlaps(ci_ref.scf)[:ci_ref.nmo,:ci_new.nmo]

        # iterate over the reference states
        bra_st  = list(ref_states.values())
        ket_st  = list(range(ci_new.n_states()))
        Smat = overlap.overlap_st(ci_ref, ci_new, bra_st, ket_st, smo,
                                  0.975, 1e-6, self.verbose) 

        root_found   = {st:True for st in ref_states.keys()}
        states_found = []
        for lbl, bst in ref_states.items():
            Sij     = np.absolute(Smat[bra_st.index(bst),:])
            s_srt   = np.flip(np.argsort(Sij))
            ind     = 0
            max_ind = s_srt[ind]
            while (ind < len(s_srt)-1 and Sij[max_ind] >= s_thrsh and 
                                            max_ind in states_found):
                ind    += 1
                max_ind = s_srt[ind]

            # if we can't find the state, set the energy to
            # 0.95*E[0]. This is undeniably arbitrary/hacky,
            # but empirically, this works reasaonbly well as a 
            # penalty function value for the optimizers 
            if Sij[max_ind] <= s_thrsh:
                root_found[lbl] = False
                eners[lbl]      = 0.95*ci_new.energies[0]
                # global optimizers will go to weird places,
                # gradient descent, not so much. If grad.
                # descent, or equiv., print a warning if we
                # can't match the state
                if self.opt_algorithm != 'differentialevolution':
                    msg = molecule+' '+str(lbl)+': Sij='
                    msg += ''.join(['{:8.4f}']*len(Sij)).format(*Sij)
                    print(msg) 
            else:
                kst        = ket_st[max_ind]
                eners[lbl] = ci_new.energies[kst]
                states_found.append(kst)
 
        return root_found, eners

    #
    def create_dirs(self, scf_objs):
        """
        Create a directory structure for running computations
        """

        scf_dirs = {}

        # set up the directory structure
        for molecule, ci2scf in scf_objs.items():

            if os.path.isdir(molecule):
                shutil.rmtree(molecule)
            os.mkdir(molecule)

            idir               = -1
            name_set           = {}
            scf_dirs[molecule] = {} 

            for ci_name, scf_name in ci2scf.items():

                if scf_name not in name_set.keys():
                    idir += 1
                    dir_name = 'scf' + str(idir)
                    os.mkdir(molecule + '/' + dir_name)
                    name_set[scf_name] = dir_name
                    scf_dirs[molecule][ci_name]  = dir_name

                else :
                    scf_dirs[molecule][ci_name] = name_set[scf_name]

        return scf_dirs

    #
    def parse_wfn_file(self, states):
        """
        Parse the GRaCI reference data and confirm what data is present
        that can be compared to the target values

        Arguments: 
        target_vals:    dictionary containing the reference data

        Returns:
        graci_vals:     a dictionary containing the corresponding GRaCI
                        objects
        """

        if not os.path.isfile(self.wfn_lib):
            msg = 'wfn_lib: ' + str(self.wfn_lib) + ' not found.'
            self.hard_exit(msg)

        # open the reference data file and get the contents
        wfn_file = h5py.File(self.wfn_lib, 'r', libver='latest')

        # get the top-level contents of the checkpoint file
        wfn_contents = chkpt.contents(file_handle = wfn_file)

        scf_objs  = {}
        ci_objs   = {}

        for molecule in states.keys():

            st_lst = set(states[molecule])

            scf_objs[molecule]   = {} 
            ci_objs[molecule]    = {} 
            states_found         = []

            # if the molecule string is in the reference 
            # object name, add it ref_objs dict
            for ci_name in wfn_contents:
                is_ci  = any([ci in ci_name for ci in params.ci_objs])
                ci_lbls = chkpt.read_attribute(wfn_file, ci_name, 
                                                         'state_map')

                # name of section is 'original.name.molecule'
                m_name = ci_name.strip().split('.')[-1]
                if (molecule == m_name and is_ci 
                                       and isinstance(ci_lbls, dict)):

                    # check if any of the required state labels are
                    # in this CI object
                    if any([lbl in ci_lbls.keys() 
                                     for lbl in states[molecule]]):

                        ci_map = {st: ci_lbls[st] 
                                  for st in states[molecule] 
                                  if st in ci_lbls.keys()}
                        ci_objs[molecule][ci_name] = ci_map
                        states_found.extend(list(ci_map.keys()))

                    # if this object contains ci states, identify the
                    # relevant scf object
                    if ci_name in ci_objs[molecule].keys(): 
 
                        # get the name of corresponding scf object
                        scf_link = chkpt.read_attribute(wfn_file, 
                                                         ci_name, 'scf')
                        is_scf, scf_name = chkpt.data_name('', scf_link)
                        if is_scf:
                            scf_objs[molecule][ci_name] = scf_name

            if len(st_lst - set(states_found)) != 0:
                msg = 'ERROR - molecule: ' + str(molecule) + ' -- ' +\
                       str(st_lst) + ' != ' + str(states_found)
                self.hard_exit(msg)
                    
        wfn_file.close()

        return scf_objs, ci_objs

    #
    def parse_ref_file(self):
        """
        Parse the file containing the data we're going to parameterize
        wrt to

        Arguments: None
        Returns:   A dictionary containing the target data
        """

        if not os.path.isfile(self.energy_lib):
            msg = 'energy_lib: ' + str(self.energy_lib) + ' not found.'
            self.hard_exit(msg)

        with open(self.energy_lib, 'r') as t_file:
            t_file_lines = t_file.readlines()

        self.n_ref = 0
        trans_ener = {}
        states     = {}

        for line in t_file_lines:
            str_arr                = line.strip().split()
            if len(str_arr) == 0:
                continue
            molecule               = str_arr[0]
            trans_ener[molecule]   = {} 
            states[molecule]       = []
            # format is: init_state final_state energy
            for i in range(1,len(str_arr),3):
                try:
                    si     = str_arr[i].strip()
                    sf     = str_arr[i+1].strip()

                    if si not in states[molecule]:
                        states[molecule].append(si)
                    if sf not in states[molecule]:
                        states[molecule].append(sf)

                    st_pair     = si+' '+sf
                    val         = parse.convert_value(str_arr[i+2])
                    trans_ener[molecule][st_pair] = val
                    self.n_ref += 1

                except:
                    msg = 'Error parsing value as str/float: ' + \
                           str(val) + '\nline = '+str(str_arr)
                    self.hard_exit(msg)

        return trans_ener, states

    #
    def sanity_check(self, ci_objs):
        """
        sanity check the input
        """

        # general options that apply to all job_types:
        if self.job_type not in ['opt', 'scan', 'analysis']:
            msg = 'job_type: '+str(self.jobtype)+' not recognized.'
            self.hard_exit(msg)

        if self.xc is None:
            msg = 'No functional specified.'
            self.hard_exit(msg)

        # this options only matter if we're doing a parameter
        # optimization run
        if self.job_type == 'opt':

            # ensure lengths of arrays are appropriate
            nham   = len(self.hamiltonians)
            nparam = len(self.params)
            nopt   = len(self.opt)
            nbnds  = len(self.bounds) 

            if len(set([nham, nparam, nopt, nbnds])) != 1:
                msg = 'Length of list arguments must be the same: '+ \
                      'len(hamiltonians)={:>1d}, '+                  \
                      'len(params)={:>1d}, '+                        \
                      'len(opt)={:>1d}, '+                           \
                      'len(bounds)={:>1d} must be the same: '.format(
                          nham, nparam, nopt, nbnds)
                self.hard_exit(msg)

            self.opt_algorithm = self.opt_algorithm.lower()
            if self.opt_algorithm not in self.valid_algos:
                msg = str(self.opt_algorithm) + ' not in ' + \
                      str(self.valid_algos)
                self.hard_exit(msg)

            self.opt_target = self.opt_target.lower()
            if self.opt_target.lower() not in self.valid_opt_targ:
                msg = str(self.opt_target) + ' not in ' + \
                      str(self.valid_opt_targ)
                self.hard_exit(msg)

        if self.job_type == 'scan' and self.ngrid is None:
            msg = 'job_type=scan, ngrid='+str(self.ngrid)+', error'
            self.hard_exit(msg)

        return

    #
    def update_opt_params(self, p_opt):
        """
        to come
        """

        n = 0

        for index in range(len(self.hamiltonians)):
            ham = self.hamiltonians[index]

            if not self.opt_options[ham]['optimize']:
                continue
          
            n_p = len(self.opt_options[ham]['params']) - \
                    len(self.opt_options[ham]['freeze'])
            params = p_opt[n:n+n_p].tolist()
            n += n_p

            m = 0
            for i in range(len(self.opt_options[ham]['params'])):
                if i not in self.opt_options[ham]['freeze']:
                    self.opt_options[ham]['params'][i] = params[m]
                    m += 1

        return 

    #
    def extract_opt_param(self):
        """
        to come
        """
       
        p_opt = np.zeros(self.n_opt, dtype=float)
        n     = 0

        for index in range(len(self.hamiltonians)):
            ham = self.hamiltonians[index]

            # skip if no parameters to optimize
            if not self.opt_options[ham]['optimize']:
                continue

            params = np.array([val for i, val in 
                           enumerate(self.opt_options[ham]['params']) 
                           if i not in self.opt_options[ham]['freeze']], 
                           dtype=float)
            p_opt[n:n+params.shape[0]] = params
            n += params.shape[0]

        return p_opt

    #
    def status_func(self, xk, convergence=None):
        """
        Give status of optimization
        """

        dif     = np.linalg.norm(xk - self.p_n)

        # some optimization methods want an optional
        # 'convergence' argument. If set, pass that
        # to print status function, else, last value
        # of the error function
        if convergence is not None:
            status_err = convergence
        else:
            status_err = self.error

        output.print_param_iter(self.iiter, xk, status_err)

        self.iiter += 1
        self.p_n    = xk.copy()

        if self.iiter >= self.max_iter:
            msg = 'Max. number of iterations completed.'
            self.hard_exit(msg)

        return

    #
    def allocate_workers(self, objs):
        """
        allocate the number of worker threads, and number of
        threads per worker. The latter will simply be set
        via the environment variable, OMP_NUM_THREAD
        """
        
        if params.nproc == 1:
            self.max_workers = 1
        else:
            self.max_workers = params.nproc
      
        return

    #
    def hard_exit(self, message):
        """
        exit the program with a sys.exit call
        """
        print("\n"+str(message)+"\n", flush=True)
        sys.exit(1)

