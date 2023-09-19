"""
The Parameterize object and its associated functions.
"""
import sys as sys
import math
import scipy.optimize as sp_opt
import numpy as np
import h5py as h5py
import os as os
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

        #  hamiltonian to use for valence state calculation
        self.valence         = None
        # optimize hamiltonian
        self.val_opt         = False
        # cimethod for valence states
        self.val_method      = None
        # keyword arguments for valence ci method
        self.val_args        = []
        #  parameters to use for valence state hamiltonian
        self.val_params      = None
        # which parameters to freeze during optimization
        self.val_freeze      = []
        # bounds for ham parameters
        self.val_bounds      = None

        #  hamiltonian to use the final state calculation 
        self.cvs             = None
        # optimize hamiltonian
        self.cvs_opt         = False
        # cimethod for initil states
        self.cvs_method      = None
        # keyword arguments for initial ci method
        self.cvs_args        = []
        #  parameters to use for initial state hamiltonian
        self.cvs_params      = None 
        # which parameters to freeze during optimization
        self.cvs_freeze      = []
        # bounds for ham parameters
        self.cvs_bounds      = None

        # which xc function - pyscf resolvable name
        self.xc              = None
        # whether to optimize xc functoin
        self.xc_opt          = False
        # variables that hold parameters 
        self.xc_params       = None
        # how the coefficients are defined in terms of parameters
        self.xc_coef         = None
        # the X and C component functionals
        self.xc_func         = None
        # which parameters to freeze during optimization
        self.xc_freeze       = []
        # bounds for xc parameters
        self.xc_bounds       = None

        # scan run options
        self.ngrid           = None
        self.scan_var        = []

        # -------------------------------
        self.opt           = {'valence': False, 'cvs':False, 'xc':False}
        self.labels        = {'valence': None, 'cvs': None, 'xc': None}
        self.bounds        = {'valence': None, 'cvs': None, 'xc': None}
        self.freeze        = {'valence': [], 'cvs': [], 'xc': []}
        self.params        = {'valence': None, 'cvs': None, 'xc': None}
        self.kwords        = {'valence': {}, 'cvs': {}}
        self.method        = {'valence': None, 'cvs': None}
 
        self.n_opt         = 0
        self.n_ref         = 0
        self.iiter         = 0
        self.p_n           = None
        self.error         = 0
        self.de_thr        = 0.5
        self.log_file      = None
        self.valid_algos   = ['nelder-mead','differentialevolution']
        self.valid_opt_targ = ['mae','rmsd']

    #
    def run(self):
        """re-parameterize the Hamiltonian"""

        # parse the target data file
        exc_ref, states  = self.parse_ref_file()

        # parse the graci reference file
        scf_objs, ci_objs = self.parse_wfn_file(states)

        # sanity check the input
        self.sanity_check(ci_objs)

        # set initial parameter set
        p_init = self.set_init_params()

        # print header info to output file defined in module output
        output.print_param_header(p_init, exc_ref, ci_objs)

        # set up directory structure: each molecule gets a subdirectory
        # and each scf calculation gets a sub-sub-directory
        scf_dirs = self.create_dirs(scf_objs)

        # the first pass sets up the orbitals and integrals -- no need
        # to recompute these every time
        exc_init = self.evaluate_energies(p_init, scf_dirs, scf_objs,
                                                  ci_objs, gen_orbs=True)
 
        # save the default logfile name for writing updates of the
        # reparam procedure
        self.logfile = output.file_names['out_file']

        if self.job_type == 'opt':

            # optimize parameters using scipy routines
            self.iiter       = 1
            self.p_n, p_bnds = self.extract_opt_param(p_init)

            if self.opt_algorithm != 'differentialevolution':
                res = sp_opt.minimize(self.err_func, self.p_n, 
                         args = (exc_ref, scf_dirs, scf_objs, ci_objs),
                         bounds   = p_bnds,
                         method   = self.opt_algorithm,
                         tol      = self.conv,
                         callback = self.status_func)
            else:
                res = sp_opt.differential_evolution(
                         self.err_func, p_bnds,
                         args = (exc_ref, scf_dirs, scf_objs, ci_objs),
                         callback = self.status_func,
                         polish   = False,
                         tol      = self.conv,
                         x0       = self.p_n)

            p_final = self.to_full_param_set(res.x)
            # one final eval_energy call with the converged params
            exc_final = self.evaluate_energies(p_final, scf_dirs,
                                                     scf_objs, ci_objs)

            output.print_param_results(p_final, res, exc_ref, 
                                                    exc_init, exc_final)
        elif self.job_type == 'scan':
            args = (exc_ref, scf_dirs, scf_objs, ci_objs)
            self.scan(p_init, self.scan_var, self.ngrid, args)

        elif self.job_type == 'analysis':
            output.print_param_analysis(p_init, exc_ref, exc_init)

        return

    #
    def set_init_params(self):
        """
        set the initial parameter values, either using default 
        Hamiltonian values, or, user supplied values
        """

        # Load the parameter dictionaries 
        #-----------------------------------------------------------

        # going to set initial and final ci object keywords here as
        # well. This should be moved eventually.
        hams = [self.valence, self.cvs]
        ci   = [self.val_method, self.cvs_method]
        args = [self.val_args, self.cvs_args]
        keys = ['valence', 'cvs']
        for i in range(len(keys)):
            # check that method is sensible
            if ci[i] is not None:
                ci_class = ci[i].lower().capitalize()
                if ci_class not in params.ci_objs:
                    msg = ci_class+' not a valid CI method.'
                    self.hard_exit(msg)
                self.method[keys[i]] = ci[i]

            # set the method keywords
            arg = args[i]
            if arg is not None and isinstance(arg, (list, np.ndarray)):
                for k in range(0, len(arg),2):
                    [key, val_str] = arg[k:k+2]
                    if '[' in val_str and ']' in val_str:
                        self.kwords[keys[i]][key] = \
                                           parse.convert_array(val_str)
                    else:
                        self.kwords[keys[i]][key] = \
                                           parse.convert_value(val_str)

            # set the hamiltonian name
            self.kwords[keys[i]]['hamiltonian'] = hams[i]

        lbls = [self.valence,         self.cvs,        self.xc]
        op   = [self.val_opt,     self.cvs_opt,    self.xc_opt]
        bnds = [self.val_bounds,  self.cvs_bounds, self.xc_bounds]
        frz  = [self.val_freeze,  self.cvs_freeze, self.xc_freeze]
        param = [self.val_params, self.cvs_params, self.xc_params]

        keys = ['valence', 'cvs', 'xc']
        for i in range(len(keys)):
            self.labels[keys[i]] = lbls[i]
            self.opt[keys[i]]    = op[i]
            self.bounds[keys[i]] = bnds[i]
            self.freeze[keys[i]] = frz[i]
            self.params[keys[i]] = param[i]

        # Identify those parameters that are to be optimized
        #----------------------------------------------------------
        self.n_opt = 0
        pfull      = {}
        pbounds    = {}

        for i in range(len(keys)):

            if self.labels[keys[i]] is not None:

                # if this is for XC functional, generate functional string
                if keys[i] == 'xc':
                    pfull['xc'] = self.xc_func_str(self.params['xc'])
                    self.params[keys[i]] = pfull['xc']

                else:

                    # pull the initial parameters from input, 
                    # if they're set
                    if self.params[keys[i]] is not None:
                        pfull[keys[i]] = self.params[keys[i]]
                    
                    # if initial params not set by user, pull them from
                    # bitci
                    else:
                        nh   = 0
                        args = (self.labels[keys[i]], nh)
                        nh   = libs.lib_func('retrieve_nhpar', args)
                        hp   = np.zeros(nh, dtype=float)
                        args = (self.labels[keys[i]], nh, hp)
                        hp   = libs.lib_func('retrieve_hpar', args)
                        pfull[keys[i]]       = hp
                        self.params[keys[i]] = hp

                if self.opt[keys[i]]:
                    self.n_opt += len(self.params[keys[i]]) - \
                                  len(self.freeze[keys[i]])

        # Set up optimization bounds
        #-------------------------------------------------------------
        for i in range(len(keys)):

            if self.labels[keys[i]] is not None:
            
                # if bounds not set, set bounds list to proper length
                if self.bounds[keys[i]] is None:
                    self.bounds[keys[i]] = [None]*len(pfull[keys[i]])

                # check that bounds are correct length
                elif (len(self.bounds[keys[i]]) != 
                                         len(self.params[keys[i]])):
                    msg = 'type='+str(keys[i])+' bounds shape mismatch:'
                    msg += str(len(self.bounds[keys[i]].shape)) + ' != '     
                    msg += str(len(pfull[keys[i]].shape))
                    self.hard_exit(msg)
 
        return pfull

    #
    def scan(self, p_vals, scan_var, ngrid, args):
        """
        scan the parameter values 
        """

        keys   = ['xc', 'valence', 'cvs']
        bounds = []
        label  = []
        p0     = []
        n_scan  = 0

        # start by freezing all coordiantes
        for key in keys:
            if self.params[key] is not None and len(self.params[key])>0:
                self.freeze[key] = list(range(len(self.params[key])))

        # loop over strings in scan_var
        for p_str in self.scan_var:
            p_val = int(p_str[-1])
 
            i = 0
            while i<3 and keys[i] not in p_str:
                i += 1

            if i==3:
                continue

            self.opt[keys[i]] = True
            self.freeze[keys[i]].remove(p_val)
            label.append(keys[i]) 
            bounds.append(self.bounds[keys[i]][p_val])
            p0.append(self.params[keys[i]][p_val])
            n_scan += 1

        delta = [(bounds[i][1] - bounds[i][0]) / (ngrid[i]-1)
                  if ngrid[i] > 1 else 0. for i in range(n_scan)]

        output.print_param_scan_head(label, p0, bounds, ngrid)
        
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
            hscan = [bounds[i][0] + step[i]*delta[i] 
                                          for i in range(n_scan)]

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
        p_full  = self.to_full_param_set(p_opt)
        ener_i  = self.evaluate_energies(p_full, scf_dirs,
                                                 scf_objs, ci_objs)
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
            self.error = np.linalg.norm(dif_vec)
        elif self.opt_target == 'mae':
            self.error = 0.
            if nde > 0:
                self.error = np.sum(np.absolute(dif_vec)) / nde

        return self.error
 
    #
    @timing.timed
    def evaluate_energies(self, p_full, scf_dirs, scf_names, ci_objs,
                                                      gen_orbs=False):
        """
        evaluate all the energies in the graci data set

        """
        # init_ci and final_ci should have same key list
        mol_names = ci_objs.keys()
        topdir    = os.getcwd()
        energies  = {}

        if params.nproc > 1:
            from mpi4py.futures import MPIPoolExecutor

            with MPIPoolExecutor(max_workers=params.nproc,
                                                 root=0) as executor:
                if executor is not None:

                    args = ((p_full, molecule, topdir, self.wfn_lib,
                             scf_dirs[molecule], scf_names[molecule], 
                             ci_objs[molecule], True, gen_orbs)
                             for molecule in mol_names)

                    for res in executor.starmap(self.eval_energy, args):
                        energies.update(res)
        else:

            for molecule in mol_names:
                mol_results = self.eval_energy(
                             p_full, molecule, topdir, self.wfn_lib, 
                             scf_dirs[molecule], scf_names[molecule], 
                             ci_objs[molecule], False, gen_orbs)

                energies.update(mol_results)

        return energies

    # 
    def eval_energy(self, p_full, molecule, topdir, wfn_file, scf_dirs,
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
            
            # determine if this is a CVS Hamiltonian or
            # valence Hamiltonian
            if len(ci_ref.icvs) > 0:
                ci_typ = 'cvs'
            else:
                ci_typ = 'valence'           

            # create the CI object 
            if self.method[ci_typ] is None:
                ci_opt = ci_ref
            else:
                ci_class = self.method[ci_typ].lower().capitalize()
                ci_opt = getattr(globals()[ci_class.lower()], ci_class)(ci_ref)

            # and update the keywords as requested
            for kword,val in self.kwords[ci_typ].items():
                if hasattr(ci_opt, kword):
                    setattr(ci_opt, kword, val)
            ci_opt.update_hparam(np.asarray(p_full[ci_typ]))
            #print('computing states: '+str(ci_names[ci_name])+' using hamiltonian: '+str(self.kwords[ci_typ]['hamiltonian']),flush=True)

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
                if gen_orbs or self.opt['xc']:
                    scf_obj = chkpt.read(scf_name[ci_name],
                                         file_handle=wfn_chkpt)
                    scf_obj.verbose = self.verbose
                    scf_obj.load()
                    scf_obj.xc = p_full['xc']
                    scf_obj.run(scf_obj.mol, None)
                    scf_objs[ci_name] = scf_obj

                    mo_ints.emo_cut = ci_opt.mo_cutoff
                    mo_ints.run(scf_obj)
    
                    # if only running scf once, save the orbs to file
                    if not self.opt['xc']:
                        fname = 'TMP_'+scf_name[ci_name]+'.chkpt.h5'
                        chkpt.write(scf_obj,
                         file_name=fname, grp_name = scf_name[ci_name])

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
                    mo_ints.emo_cut = ci_opt.mo_cutoff
                    mo_ints.load_bitci(scf_obj)

            ci_opt.update_eri(mo_ints)
            all_found  = False
            i_add      = 0
            n_add      = 5 
            add_states = np.asarray([1]*len(ci_opt.nstates), dtype=int)

            # iterate a couple times, in case we need to add roots to
            # find the states of interest
            while not all_found and i_add < n_add:

                ci_opt.run(scf_objs[ci_name], None)

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
                                                   0.975, self.verbose) 

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

        scf_objs = {}
        ci_objs  = {}

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

        if self.xc is None and self.xc_params is None:
            msg = 'No functional specified.'
            self.hard_exit(msg)

        if self.val_opt and self.valence is None:
            msg =  'Hamiltonian optimization requested, '
            msg += 'but no value of ham specified'
            self.hard_exit(msg)

        if self.cvs_opt and self.cvs is None:
            msg =  'CVS Hamiltonian optimization requested, '
            msg += 'but no value of cvs specified'
            self.hard_exit(msg)

        # this options only matter if we're doing a parameter
        # optimization run
        if self.job_type == 'opt':

            if not any([self.val_opt, self.cvs_opt, self.xc_opt]):
                msg = 'job_type=op, but no parameters to be optimized'
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

        # check that we've specified all the hamiltonians we need to
        # determined the necessary states
        wfn_file = h5py.File(self.wfn_lib, 'r', libver='latest')
        obj_types = {'valence': False, 'cvs': False}
        for molecule,mdict in ci_objs.items():
            for ci_name, states in mdict.items():
                ci_obj = chkpt.read(ci_name, file_handle=wfn_file)
                if len(ci_obj.icvs) > 0:
                    obj_types['cvs'] = True
                else:
                    obj_types['valence'] = True
            
            # if we've already found both types of objects, exit,
            # no need to continue
            if all(obj_types.values()):
                break

        wfn_file.close()

        if obj_types['valence'] and self.valence is None:
            msg =  'Valence CI objects to be computed, but valence '
            msg += 'Hamiltonian not specified.'
            self.hard_exit(msg)

        if obj_types['cvs'] and self.cvs is None:
            msg =  'CVS CI objects to be computed, but CVS '
            msg += 'Hamiltonian not specified.'
            self.hard_exit(msg)
          
        return

    #
    def xc_func_str(self, xc_param):
        """
        build an xc functional string given the current values 
        in xc_param
        """

        # if functional name specified, simply return it
        if self.labels['xc'] is not None:
            return self.labels['xc']      

        #...else construct functional string
        # first indx for exchange part, second indx for correlation part
        n_xc   = [len(self.xc_coef[i]) for i in range(2)]
        xc_str = ''
        cfmt   = '{:7.5f}'

        for ixc in range(2):
            for ifunc in range(n_xc[ixc]):
                cstr = self.xc_coef[ixc][ifunc]

                if 'xc_' in cstr:
                    fstr = self.xc_func[ixc][ifunc]

                    ind  = cstr.find('xc_')+3
                    arg  = int(cstr[ind])
                    cf   = eval( cstr.replace('xc_'+str(arg),
                                               str(xc_param[arg])))
                else:
                    cf  = eval(cstr)

                xc_str += cfmt.format(cf) + '*' + fstr
                if ifunc != n_xc[ixc]-1:
                    xc_str += ' + '

            if ixc == 0:
                xc_str += ','      

        return xc_str

    #
    def to_full_param_set(self, p_opt):
        """
        to come
        """
        p_full = {}
        n      = 0

        keys = ['xc', 'valence', 'cvs']
        for i in range(len(keys)):

            if self.labels[keys[i]] is None:       
                p_full[keys[i]] = None
                continue

            if (self.params[keys[i]] is None or 
                     isinstance(self.params[keys[i]], str)):
                p_vals = self.params[keys[i]]
            else:
                p_vals = self.params[keys[i]].copy()

            if self.opt[keys[i]]:
                nall  = len(p_vals)
                ni    = nall - len(self.freeze[keys[i]])
                opt_p = np.setdiff1d(np.array(range(nall)), 
                                     self.freeze[keys[i]])
                p_vals[opt_p] = p_opt[n : n + ni]
                n    += ni

            if keys[i] == 'xc':
                p_vals = self.xc_func_str(p_vals)

            p_full[keys[i]] = p_vals

        return p_full

    #
    def extract_opt_param(self, p_full):
        """
        to come
        """
       
        p_opt = np.zeros(self.n_opt, dtype=float)
        b_opt = np.zeros((self.n_opt,2), dtype=float)
        n     = 0

        keys    = ['xc', 'valence', 'cvs']
        for i in range(len(keys)):

            if self.opt[keys[i]]:
                nall = len(self.params[keys[i]])
                nopt = nall - len(self.freeze[keys[i]])
                pval = [self.params[keys[i]][j] for j in range(nall) 
                                 if j not in self.freeze[keys[i]]]
                bval = [self.bounds[keys[i]][j] for j in range(nall)
                                 if j not in self.freeze[keys[i]]]
              
                p_opt[n : n+nopt] = pval
                b_opt[n : n+nopt] = bval
                n += nopt 

        return p_opt, b_opt

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
    def hard_exit(self, message):
        """
        exit the program with a sys.exit call
        """
        print("\n"+str(message)+"\n", flush=True)
        sys.exit(1)

