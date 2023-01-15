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
        #  hamiltonian to use for initial state calculation
        self.init_state_h    = None
        # cimethod for initial states
        self.init_ci_method  = None
        # keyword arguments for initial ci method
        self.init_ci_args    = None
        #  parameters to use for initial state hamiltonian
        self.h_init_params   = None
        #  hamiltonian to use the final state calculation 
        self.final_state_h   = None
        # ci method type for final states
        self.final_ci_method = None
        # keyword arguments for final ci method
        self.final_ci_args   = None
        #  parameters to use for the final state hamiltonian
        self.h_final_params  = None
        # couple h_init and h_final, i.e. same hamiltonian
        self.single_hamiltonian = False
        # XC functional to use: assumes in libxc format to
        # facilitate using internal functions
        self.xc              = None
        # numerical parameters for user defined xc functional
        self.xc_params       = None
        # how the coefficients are defined in terms of parameters
        self.xc_coef         = None
        # the X and C component functionals
        self.xc_func         = None

        # reparam optimization options
        self.opt_algorithm   = 'Nelder-Mead'
        self.opt_h_init      = False
        self.opt_h_final     = False
        self.opt_xc          = False
        self.conv            = 1.e-3
        self.max_iter        = 1000
        self.bounds_h_init   = None
        self.bounds_h_final  = None
        self.bounds_xc       = None
        self.freeze_h_init   = []
        self.freeze_h_final  = []
        self.freeze_xc       = []

        # scan run options
        self.ngrid           = None
        self.scan_var        = []

        # -------------------------------
        self.init_kwords    = {}
        self.final_kwords   = {}
        self.n_opt          = 0
        self.p_ref          = None
        self.n_ref          = 0
        self.iiter          = 0
        self.p_n            = None
        self.error          = 0
        self.de_thr         = 0.5
        self.log_file       = None
        self.valid_algos    = ['Nelder-Mead','DifferentialEvolution']

    #
    def run(self):
        """re-parameterize the Hamiltonian"""

        # sanity check the input
        self.sanity_check()

        # parse the target data file
        exc_ref, init_state, final_state = self.parse_ref_file()

        # parse the graci reference file
        scf_objs, init_ci, final_ci = self.parse_wfn_file(init_state, 
                                                          final_state)

        # set initial parameter set
        self.p_ref, p_bounds  = self.set_init_params()

        # print header info to output file defined in module output
        output.print_param_header(self.p_ref,exc_ref, init_ci, final_ci)

        # set up directory structure: each molecule gets a subdirectory
        # and each scf calculation gets a sub-sub-directory
        scf_dirs = self.create_dirs(scf_objs)

        # the first pass sets up the orbitals and integrals -- no need
        # to recompute these every time
        exc_init = self.evaluate_energies(self.p_ref, scf_dirs, scf_objs,
                                        init_ci, final_ci, gen_orbs=True)

        # save the default logfile name for writing updates of the
        # reparam procedure
        self.logfile = output.file_names['out_file']

        if self.job_type == 'opt':

            # optimize parameters using scipy routines
            self.iiter       = 1
            self.p_n, p_bnds = self.extract_opt_param(self.p_ref, 
                                                              p_bounds)

            if self.opt_algorithm != 'DifferentialEvolution':
                res = sp_opt.minimize(self.err_func, self.p_n, 
                              args = (exc_ref, scf_dirs, scf_objs,
                                                init_ci, final_ci),
                              bounds   = p_bnds,
                              method   = self.opt_algorithm,
                              tol      = self.conv,
                              callback = self.status_func)
            else:
                res = sp_opt.differential_evolution(
                                   self.err_func, p_bnds,
                                   args = (exc_ref, scf_dirs, scf_objs,
                                                    init_ci, final_ci),
                                   callback = self.status_func,
                                   polish   = False,
                                   tol      = self.conv,
                                   x0       = self.p_n)

            p_final = self.to_full_param_set(res.x)
            # one final eval_energy call with the converged params
            exc_final = self.evaluate_energies(p_final, scf_dirs,
                                            scf_objs, init_ci, final_ci)

            output.print_param_results(p_final, res, exc_ref, 
                                                    exc_init, exc_final)
        elif self.job_type == 'scan':
            args = (exc_ref, scf_dirs, scf_objs, init_ci, final_ci)
            self.scan(self.p_ref,p_bounds,self.scan_var,self.ngrid,args)

        elif self.job_type == 'analysis':
            output.print_param_analysis(self.p_ref, exc_ref, exc_init)

        return

    #
    def set_init_params(self):
        """
        set the initial parameter values, either using default 
        Hamiltonian values, or, user supplied values
        """

        self.n_opt = 0
        pfull      = {}
        pbounds    = {}

        #--------------------------------------------------------------
        # XC functional

        # if we're optimizing the functional parameters, add
        # them to the parameter list. If not, set this value
        # to the xc string passed to pyscf
        if self.opt_xc:
            self.n_opt += len(self.xc_params) - len(self.freeze_xc)
        pfull['xc'] = self.xc_func_str(self.xc_params)

        #--------------------------------------------------------------
        # Initial state Hamiltonian

        # hamiltonian for initial states
        if self.h_init_params is not None:
            pfull['h_init'] = self.h_init_params
        else:
            nh   = 0
            args = (self.init_state_h, nh)
            nh   = libs.lib_func('retrieve_nhpar', args)
            hp   = np.zeros(nh, dtype=float)
            args = (self.init_state_h, nh, hp)
            hp   = libs.lib_func('retrieve_hpar', args)
            pfull['h_init'] = hp

        if self.opt_h_init:
            self.n_opt += len(pfull['h_init']) - len(self.freeze_h_init)

        #---------------------------------------------------------------
        # Final state Hamiltonian

        # hamiltonian for final states
        if self.h_final_params is not None:
            pfull['h_final'] = self.h_final_params
        else:
            nh   = 0
            args = (self.final_state_h, nh)
            nh   = libs.lib_func('retrieve_nhpar', args)
            hp   = np.zeros(nh, dtype=float)
            args = (self.final_state_h, nh, hp)
            hp   = libs.lib_func('retrieve_hpar', args)
            pfull['h_final'] = hp

        if self.opt_h_final and self.final_state_h != self.init_state_h:
            self.n_opt += len(pfull['h_final']) - len(self.freeze_h_final)

        #---------------------------------------------------------------
        # Set up optimization bounds
        
        # bounds: only applies to 'opt' and 'scan' runs
        if self.xc_params is None:
           pbounds['xc'] = None

        elif self.bounds_xc is None:
           pbounds['xc'] = [None] * len(self.xc_params)

        elif len(self.bounds_xc) == len(self.xc_params):
           pbounds['xc'] = self.bounds_xc
     
        else:
            msg  = 'bounds_xc/xc_params shape mismatch: '
            msg += str(len(self.bounds_xc)) + ' != ' 
            msg += str(len(self.xc_params))
            self.hard_exit(msg)

        # bounds work a little differently for hamiltonian params
        if self.bounds_h_init is None:
           pbounds['h_init'] = [None] * len(pfull['h_init'])

        elif len(self.bounds_h_init) == len(pfull['h_init']):
           pbounds['h_init'] = self.bounds_h_init

        else:
            msg  = 'bounds_h_init/h_init_params shape mismatch: '
            msg += str(len(self.bounds_h_init.shape)) + ' != '     
            msg += str(len(pfull['h_init'].shape))
            self.hard_exit(msg)

        # bounds work a little differently for hamiltonian params
        if self.bounds_h_final is None:
           pbounds['h_final'] = [None] * len(pfull['h_final'])

        elif len(self.bounds_h_final) == len(pfull['h_final']):
           pbounds['h_final'] = self.bounds_h_final

        else:
            msg  = 'bounds_h_final/h_final_params shape mismatch: '
            msg += str(len(self.bounds_h_final)) + ' != ' 
            msg += str(len(pfull['h_final']))
            self.hard_exit(msg)

        return pfull, pbounds

    #
    def err_func(self, p_opt, exc_ref, scf_dirs, scf_objs, 
                                               init_ci, final_ci):
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
        ener_i  = self.evaluate_energies(p_full, scf_dirs, scf_objs, 
                                                    init_ci, final_ci)
        dif_vec = np.zeros(self.n_ref, dtype=float)

        # this approach assumes dict is ordered! Only true from Python
        # 3.6 onwards...
        ide        = 0
        for molecule in exc_ref.keys():
            for trans, ener in exc_ref[molecule].items():
                init, final  = trans.strip().split()
                de_iter      = ener_i[molecule][final] - \
                               ener_i[molecule][init]
                dif_vec[ide] = de_iter*constants.au2ev - ener
                ide         += 1

        self.error    = np.linalg.norm(dif_vec)

        return self.error

    @timing.timed
    def scan(self, p_ref, p_bounds, scan_var, ngrid, args):
        """
        scan the parameter values 
        """

        # number of coordinates to scan over
        n_scan = len(scan_var)

        # freeze all other parameters other than the ones we
        # are going to scan
        if self.xc_params is not None:
            self.freeze_xc  = list(range(len(self.xc_params)))
        self.freeze_h_init  = list(range(len(p_ref['h_init'])))
        self.freeze_h_final = list(range(len(p_ref['h_final'])))        
        bounds              = []
        p0                  = []
        label               = []

        # identify the scan coordinates
        for p_str in self.scan_var:
            p_val = int(p_str[-1])
            if 'xc' in p_str:
                self.opt_xc = True
                self.freeze_xc.remove(p_val)
                p0.append(self.xc_params[p_val])
                bounds.append(self.bounds_xc[p_val])
                lbl = 'xc'
                label.append('xc')
            if 'h_init' in p_str:
                self.opt_h_init = True
                self.freeze_h_init.remove(p_val)
                if self.single_hamiltonian:
                    self.freeze_h_final.remove(p_val)
                    lbl = 'h_init/final'
                else:
                    lbl = 'h_init'
                p0.append(self.h_init_params[p_val])
                bounds.append(self.bounds_h_init[p_val])
            if 'h_final' in p_str:
                self.opt_h_final = True
                self.freeze_h_final.remove(p_val)
                if self.single_hamiltonian:
                    self.freeze_h_init.remove(p_val)
                    lbl = 'h_init/final'
                else:
                    lbl = 'h_final'
                p0.append(self.h_final_params[p_val])
                bounds.append(self.bounds_h_final[p_val])
            label.append(lbl)

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

    @timing.timed
    def evaluate_energies(self, p_full, scf_dirs, scf_names, 
                                   init_ci, final_ci, gen_orbs=False):
        """
        evaluate all the energies in the graci data set

        """

        topdir = os.getcwd()
        energies = {}

        if params.nproc > 1:

            with MPIPoolExecutor(max_workers=params.nproc,
                                                 root=0) as executor:
                if executor is not None:

                    args = ((p_full, molecule, topdir, self.wfn_lib,
                             scf_dirs[molecule], scf_names[molecule], 
                             init_ci[molecule], final_ci[molecule], 
                             True, gen_orbs)
                             for molecule in init_ci.keys())

                    for res in executor.starmap(self.eval_energy, args):
                        energies.update(res)
        else:

            for molecule in init_ci.keys():
                mol_results = self.eval_energy(p_full, molecule,
                                  topdir, self.wfn_lib, 
                                  scf_dirs[molecule], 
                                  scf_names[molecule], 
                                  init_ci[molecule],
                                  final_ci[molecule], 
                                  False, gen_orbs)

                energies.update(mol_results)

        return energies

    # 
    def eval_energy(self, p_full, molecule, topdir, wfn_file, scf_dirs,
                         scf_name, init_ci, final_ci, parallel, gen_orbs):
        """
        eval_energy
        """

        # if parallel, this is run on a worker and libraries need to be 
        if parallel:
            libs.lib_load('bitci')
            libs.lib_load('overlap')

        wfn_chkpt  = h5py.File(wfn_file, 'r', libver='latest')
        output.file_names['out_file'] = topdir+'/'+str(molecule)+'.log'
        mol_dir = topdir+'/'+str(molecule) 

        os.chdir(mol_dir)
        energies  = {molecule : {}}
        mo_ints   = ao2mo.Ao2mo()
        
        # identify the unique CI objects from the initial and 
        # final state dictionaries   
        ci_unique = list( set( list(init_ci.keys()) + 
                               list(final_ci.keys())))
        scf_objs  = {ci_unique[i]:None for i in range(len(ci_unique))}
 
        # loop over unique ci state objects
        for ci_name in ci_unique:     

            # ci_ref contains the wfns corresponding to the target
            # states
            ci_ref = chkpt.read(ci_name, file_handle=wfn_chkpt)
           
            ref_states = {}
            # if this hamiltonian is for the initial states, apply
            # arguments that apply to the initial states
            if ci_name in init_ci.keys():

                if self.init_ci_method is None:
                    ci_opt = ci_ref.copy()
                else:
                    ci_class = self.init_ci_method.lower().capitalize()
                    ci_opt = getattr(globals()[ci_class.lower()], ci_class)(ci_ref)               

                for kword,val in self.init_kwords.items():
                    if hasattr(ci_opt, kword):
                        setattr(ci_opt, kword, val)
                ci_opt.update_hparam(np.asarray(p_full['h_init']))
                ref_states.update(init_ci[ci_name])

            elif ci_name in final_ci.keys():

                if self.final_ci_method is None:
                    ci_opt = ci_ref.copy()
                else:
                    ci_class = self.final_ci_method.lower().capitalize()
                    ci_opt = getattr(globals()[ci_class.lower()], ci_class)(ci_ref)

                for kword,val in self.final_kwords.items():
                    if hasattr(ci_opt, kword):
                        setattr(ci_opt, kword, val)
                ci_opt.update_hparam(np.asarray(p_full['h_final']))

            # this feels a little hacky: if same hamiltonian for 
            # init and final states -- make sure we add final 
            # states to the list we search for
            if ci_name in final_ci.keys():
                ref_states.update(final_ci[ci_name])

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
                if gen_orbs or self.opt_xc:
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
                    if not self.opt_xc:
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
                #print(molecule+' identify start',flush=True)
                roots_found, eners = self.identify_states(molecule, 
                                             ci_ref, ref_states, ci_opt)

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

            # if all roots found, update the energy directory 
            #if roots_found:       
            #    energies[molecule].update(eners) 

            # if we fail to match states, set energies to zero, thereby
            # resulting in very large errors
            #else:
            #    msg = str(molecule)+' failed to match states: ' + \
            #          str(ref_states) + ' hparam=' + str(ci_opt.hparam)
            #         
            #    self.hard_exit(msg)

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
        s_thrsh = 1./np.sqrt(2.)
        #s_thrsh = 0.6

        eners  = {}

        # in functional optimization runs, this will be something
        # other than identity matrix
        smo = ci_new.scf.mo_overlaps(ci_ref.scf)[:ci_ref.nmo,:ci_new.nmo]
        #if ref_ci.nmo != new_ci.nmo:
        #    print(molecule+' ref.nmo='+str(ref_ci.nmo)+','+str(new_ci.nmo),flush=True)

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
    def parse_wfn_file(self, init_states, final_states):
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
        wfn_file = h5py.File(self.wfn_lib, 'r', libver='latest')

        # get the top-level contents of the checkpoint file
        wfn_contents = chkpt.contents(file_handle = wfn_file)

        scf_objs = {}
        init_ci  = {}
        final_ci = {}

        for molecule in init_states.keys():

            st_lst = set(init_states[molecule] + 
                         final_states[molecule])

            scf_objs[molecule]   = {} 
            init_ci[molecule]    = {} 
            final_ci[molecule]   = {}
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
                                     for lbl in init_states[molecule]]):

                        ci_map = {st: ci_lbls[st] 
                                  for st in init_states[molecule] 
                                  if st in ci_lbls.keys()}
                        init_ci[molecule][ci_name] = ci_map
                        states_found.extend(list(ci_map.keys()))

                    # check if any of the required state labels are
                    # in this CI object
                    if any([lbl in ci_lbls.keys() 
                                    for lbl in final_states[molecule]]):

                        ci_map = {st: ci_lbls[st] 
                                  for st in final_states[molecule] 
                                  if st in ci_lbls.keys()}
                        final_ci[molecule][ci_name] = ci_map
                        states_found.extend(list(ci_map.keys()))
        
                    # if this object contains ci states, identify the
                    # relevant CI object
                    if (ci_name in init_ci[molecule].keys() or 
                           ci_name in final_ci[molecule].keys()):
 
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

        return scf_objs, init_ci, final_ci

    #
    def parse_ref_file(self):
        """
        Parse the file containing the data we're going to parameterize
        wrt to

        Arguments: None
        Returns:   A dictionary containing the target data
        """

        with open(self.energy_lib, 'r') as t_file:
            t_file_lines = t_file.readlines()

        self.n_ref   = 0
        trans_ener   = {}
        init_states  = {}
        final_states = {}
        for line in t_file_lines:
            str_arr                = line.strip().split()
            molecule               = str_arr[0]
            trans_ener[molecule]   = {} 
            init_states[molecule]  = []
            final_states[molecule] = []
            # format is: init_state final_state energy
            for i in range(1,len(str_arr),3):
                try:
                    si     = str_arr[i].strip()
                    sf     = str_arr[i+1].strip()

                    if si not in init_states[molecule]:
                        init_states[molecule].append(si)
                    if sf not in final_states[molecule]:
                        final_states[molecule].append(sf)

                    st_pair     = si+' '+sf
                    val         = parse.convert_value(str_arr[i+2])
                    trans_ener[molecule][st_pair] = val
                    self.n_ref += 1

                except:
                    msg = 'Error parsing value as str/float: ' + \
                           str(val) + '\nline = '+str(str_arr)
                    self.hard_exit(msg)

            # check that the intersection of initial and final
            # states is zero
            set_i = set(init_states[molecule])
            set_f = set(final_states[molecule])
            if len(set_i.intersection(set_f)) > 0:
                msg = 'warning - '+molecule+' init/final state overlap'
                msg += '\ninit states='+str(init_states[molecule]) +'\n'
                msg += '\nfinal states='+str(final_states[molecule])
                print(msg)
                    
        return trans_ener, init_states, final_states

    #
    def sanity_check(self):
        """
        sanity check the input
        """

        # general options that apply to all job_types:

        if self.job_type not in ['opt', 'scan', 'analysis']:
            msg = 'job_type: '+str(self.jobtype)+' not recognized.'
            self.hard_exit(msg)

        if not os.path.isfile(self.wfn_lib):
            msg = 'wfn_lib: ' + str(self.wfn_lib) + ' not found.'
            self.hard_exit(msg)

        if not os.path.isfile(self.energy_lib):
            msg = 'energy_lib: ' + str(self.energy_lib) + ' not found.'
            self.hard_exit(msg)

        if self.init_state_h is None:
            msg = 'initial state hamiltonian not specified'
            self.hard_exit(msg)

        if self.final_state_h is None:
            msg = 'final state hamiltonian not specified'
            self.hard_exit(msg)

        if self.xc is None and self.xc_params is None:
            msg = 'No functional specified.'
            self.hard_exit(msg)

        if self.init_ci_method is not None:
            ci_class = self.init_ci_method.lower().capitalize()
            if ci_class not in params.ci_objs:
                msg = ci_class+' not a valid CI method.'
                self.hard_exit(msg)

        if self.final_ci_method is not None:
            ci_class = self.final_ci_method.lower().capitalize()
            if ci_class not in params.ci_objs:
                msg = ci_class+' not a valid CI method.'
                self.hard_exit(msg)

        # if single hamiltonian is requested, we're going to reset
        # some values to ensure internal consistency
        if self.single_hamiltonian:
            if self.init_state_h != self.final_state_h:
                msg  = 'single hamiltonian requested, '
                msg += 'but init_state_h != final_state_h'
                self.hard_exit(msg)

            # optimize both or neither
            (opti,optf) = (self.opt_h_init, self.opt_h_final)
            self.opt_h_init = self.opt_h_final = opti or optf
            
            # set the parameters the value of the init_param set
            if self.h_init_params is not None:
                self.h_final_params = self.h_init_params
            elif self.h_final_param is not None:
                self.h_init_params = self.h_final_params
            else:
                msg = 'either h_init/final_params need to be specified'
                self.hard_exit(msg)

            # frozen parameters union of initial and final sets
            freeze = set(self.freeze_h_init).union(
                                     set(self.freeze_h_final))
            self.freeze_h_init = self.freeze_h_final = freeze
 
            # unify the bounds for initial and final hamiltonian
            bounds = None
            if self.bounds_h_init is not None:
                bounds = self.bounds_h_init
            if self.bounds_h_final is not None:
                if bounds is None:
                    bounds = self.bounds_h_final
                else:
                    bounds = [[ max(bounds[i][0],
                                    self.bounds_h_final[i][0]),
                                min(bounds[i][1],
                                    self.bounds_h_final[i][1])] 
                                for i in range(len(bounds))]
            self.bounds_h_init = self.bounds_h_final = bounds

        # if user wants to optimize functional, confirm that
        # input makes sense
        if (self.xc_params is not None and
              self.xc_coef is not None and
                 self.xc_func is not None):
            shp_cf = [len(self.xc_coef[i])
                      for i in range(len(self.xc_coef))]
            shp_fn = [len(self.xc_func[i])
                      for i in range(len(self.xc_func))]

        nhpar = 0
        args = (self.init_state_h, nhpar)
        nhpar = libs.lib_func('retrieve_nhpar', args)
        if self.h_init_params is not None and \
          len(self.h_init_params) != nhpar:
            msg = '# h_init_params != nhpar: ' + \
                  str(len(self.h_init_params)) + '!= ' + str(nhpar)
            self.hard_exit(msg)

        nhpar = 0
        args = (self.final_state_h, nhpar)
        nhpar = libs.lib_func('retrieve_nhpar', args)
        if self.h_final_params is not None and \
          len(self.h_final_params) != nhpar:
            msg = '# h_final_params != nhpar: ' + \
                  str(len(self.h_final_params)) + '!= ' + str(nhpar)
            self.hard_exit(msg)

        # this options only matter if we're doing a parameter
        # optimization run
        if self.job_type == 'opt':

            if not any([self.opt_h_init, self.opt_h_final,self.opt_xc]):
                msg = 'job_type=op, but no parameters to be optimized'
                self.hard_exit(msg)

            if self.opt_algorithm not in self.valid_algos:
                msg = str(self.opt_algorithm) + ' not in ' + \
                      str(self.valid_algos)
                self.hard_exit(msg)

        if self.job_type == 'scan' and self.ngrid is None:
            msg = 'job_type=scan, ngrid='+str(self.ngrid)+', error'
            self.hard_exit(msg)

        # going to set initial and final ci object keywords here as
        # well. This should be moved eventually.
        if (self.init_ci_args is not None and 
             isinstance(self.init_ci_args, (list, np.ndarray))):
            for k in range(0,len(self.init_ci_args),2):
                [key, val_str] = self.init_ci_args[k:k+2]
                if '[' in val_str and ']' in val_str:
                    self.init_kwords[key] = parse.convert_array(val_str)
                else:
                    self.init_kwords[key] = parse.convert_value(val_str) 
        self.init_kwords['hamiltonian'] = self.init_state_h
 
        # going to set initial and final ci object keywords here as
        # well. This should be moved eventually.
        if (self.final_ci_args is not None and
             isinstance(self.final_ci_args, (list, np.ndarray))):
            for k in range(0,len(self.final_ci_args),2):
                [key, val_str] = self.final_ci_args[k:k+2]
                if '[' in val_str and ']' in val_str:
                    self.final_kwords[key] = parse.convert_array(val_str)
                else:
                    self.final_kwords[key] = parse.convert_value(val_str)
        self.final_kwords['hamiltonian'] = self.final_state_h

        # if there is a single hamiltonian for both initial and 
        # final states, the args dicts are made identifical and 
        # are the union of the two arg lists
        if self.single_hamiltonian:
            self.init_kwords.update(self.final_kwords)
            self.final_kwords = self.init_kwords

        return


    #
    def xc_func_str(self, xc_param):
        """
        build an xc functional string given the current values 
        in xc_param
        """

        if self.xc is not None:
            return self.xc        

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

        # set the hamiltonian parameters: both frozen and optimized 
        p_vals = self.xc_params
        if self.opt_xc:
            nall  = len(p_vals)
            n     = nall - len(self.freeze_xc)
            opt_p = np.setdiff1d(np.array(range(nall)), self.freeze_xc)
            p_vals[opt_p] = p_opt[:n]
        p_full['xc']  = self.xc_func_str(p_vals)       

        p_vals = self.p_ref['h_init'].copy()
        if self.opt_h_init:
            nall  = len(p_vals)
            ni    = nall - len(self.freeze_h_init)
            opt_p = np.setdiff1d(np.array(range(nall)), 
                                 self.freeze_h_init)
            p_vals[opt_p] = p_opt[n : n + ni]
            n    += ni
        p_full['h_init']  = p_vals 

        # if initial final hamiltonian are the same, set
        # h_final to h_initial
        if self.single_hamiltonian:
            p_full['h_final'] = p_full['h_init']

        else:
            p_vals = self.p_ref['h_final'].copy()
            if self.opt_h_final:
                nall  = len(p_vals)
                nf    = nall - len(self.freeze_h_final)
                opt_p = np.setdiff1d(np.array(range(nall)), 
                                     self.freeze_h_final)
                p_vals[opt_p] = p_opt[n : n + nf]

            p_full['h_final'] = p_vals


        return p_full

    #
    def extract_opt_param(self, p_full, bnd_full):
        """
        to come
        """
       
        p_opt   = np.zeros(self.n_opt, dtype=float)
        bnd_opt = np.zeros((self.n_opt,2), dtype=float)
        n       = 0

        if self.opt_xc:
            nx            = len(self.xc_params)
            n            += nx - len(self.freeze_xc)
            p_opt[:n]     = [self.xc_params[i] for i in range(nx)
                                             if i not in self.freeze_xc]
            bnd_opt[:n,:] = [bnd_full['xc'][i] for i in range(nx)             
                                             if i not in self.freeze_xc]

        if self.opt_h_init:
            nh    = len(p_full['h_init'])
            h_opt = [p_full['h_init'][i] for i in range(nh) 
                                         if i not in self.freeze_h_init]
            b_opt = [bnd_full['h_init'][i] for i in range(nh) 
                                         if i not in self.freeze_h_init]
            nh_opt    = len(h_opt)
            p_opt[n:n + nh_opt]     = h_opt
            bnd_opt[n:n + nh_opt,:] = b_opt
            n        += nh_opt 

        if self.opt_h_final and not self.single_hamiltonian:
            nh    = len(p_full['h_final'])
            h_opt = [p_full['h_final'][i] for i in range(nh) 
                                        if i not in self.freeze_h_final]
            b_opt = [bnd_full['h_final'][i] for i in range(nh)
                                        if i not in self.freeze_h_final]
            nh_opt        = len(h_opt)
            p_opt[n:]     = h_opt
            bnd_opt[n:,:] = b_opt
            n    += nh_opt

        return p_opt, bnd_opt

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

