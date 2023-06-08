"""Module for storing global parameters"""

# molecule section input keywords and data types
molecule_kword =     {'xyz_file'    : str,
                      'units'       : str,
                      'basis'       : str,
                      'ri_basis'    : str,
                      'ao_cart'     : bool,
                      'use_sym'     : bool,
                      'sym_grp'     : str,
                      'use_df'      : bool,
                      'add_rydberg' : str,
                      'label'       : str}

rydano_kword =       {'xc'         : str,
                      'charge'     : int,
                      'mult'       : int,
                      'origin'     : float,
                      'contract'   : str,
                      'verbose'    : bool,
                      'nprimitive' : int,
                      'print_ano'  : bool,
                      'label'      : str}

#----------------------------------------------

# DFT section input keywords and data types
scf_kword      = {'xc'             : str,
                  'restart'        : bool,
                  'print_orbitals' : bool,
                  'verbose'        : bool,
                  'charge'         : int,
                  'mult'           : int,
                  'max_iter'       : int,
                  'init_guess'     : str,
                  'diag_method'    : str,
                  'lvl_shift'      : float,
                  'grid_level'     : int,
                  'conv_tol'       : float,
                  'cosmo'          : bool,
                  'solvent'        : str,
                  'mol_label'      : str,
                  'guess_label'    : str,
                  'label'          : str}

# MRCI section input keywords and data types
dftmrci_kword  = {'mult'           : int,
                  'charge'         : int,
                  'nstates'        : int,
                  'hamiltonian'    : str,
                  'hparam'         : float,
                  'ras1'           : int,
                  'ras2'           : int,
                  'ras3'           : int,
                  'nhole1'         : int,
                  'nelec3'         : int,
                  'icvs'           : int,
                  'autoras'        : bool,
                  'refiter'        : int,
                  'ref_prune'      : bool,
                  'guess_label'    : str,
                  'precision'      : str,
                  'prune'          : bool,
                  'prune_thresh'   : float,
                  'prune_qcorr'    : bool,
                  'prune_extra'    : int,
                  'regfac'         : float,
                  'diag_guess'     : str,
                  'diag_method'    : str,
                  'diag_tol'       : float,
                  'diag_iter'      : int,
                  'diag_blocksize' : int,
                  'diag_deflate'   : bool,
                  'print_orbitals' : bool,
                  'save_wf'        : bool,
                  'ddci'           : bool,
                  'ref_state'      : int,
                  'verbose'        : bool,
                  'mo_cutoff'      : float,
                  'scf_label'      : str,
                  'label'          : str}

# DFT/MRCI(2) section input keywords and data types
dftmrci2_kword  = {'mult'           : int,
                   'charge'         : int,
                   'truncate'       : bool,
                   'regularizer'    : str,
                   'regfac'         : float,
                   'nstates'        : int,
                   'hamiltonian'    : str,
                   'hparam'         : float,
                   'ras1'           : int,
                   'ras2'           : int,
                   'ras3'           : int,
                   'nhole1'         : int,
                   'nelec3'         : int,
                   'icvs'           : int,
                   'autoras'        : bool,
                   'refiter'        : int,
                   'ref_prune'      : bool,
                   'guess_label'    : str,
                   'precision'      : str,
                   'diabatic'       : bool,
                   'adt_type'       : str,
                   'print_orbitals' : bool,
                   'save_wf'        : bool,
                   'ddci'           : bool,
                   'ref_state'      : int,
                   'nbuffer'        : int,
                   'verbose'        : bool,
                   'refsel'         : str,
                   'mo_cutoff'      : float,
                   'scf_label'      : str,
                   'label'          : str}

#---------------------------------------------------

# spinorbit section input keywords and data typess
spinorbit_kword = {'couple_groups'     : str,
                   'couple_states'     : int,
                   'print_soc_thresh'  : float,
                   'mf2e'              : str,
                   'representation'    : str,
                   'verbose'           : bool,
                   'label'             : str}

# transition moments input keywords and data types
transition_kword = {'init_states'      : int, 
                    'init_label'       : str,
                    'final_states'     : int,
                    'final_label'      : str,
                    'print_orbitals'   : bool,
                    'print_orb_thresh' : float,
                    'representation'   : str,
                    'verbose'          : bool,
                    'label'            : str}

overlap_kword   = {'bra_states'        : int,
                   'bra_label'         : str,
                   'ket_states'        : int,
                   'ket_label'         : str,
                   'norm_thresh'       : float,
                   'verbose'           : bool,
                   'label'             : str}

dyson_kword     = {'init_states'       : int,
                   'init_label'        : str,
                   'final_states'      : int,
                   'final_label'       : str,
                   'norm_thresh'       : float,
                   'print_orbitals'    : bool,
                   'representation'    : str,
                   'verbose'           : bool,
                   'label'             : str}

#----------------------------------------------------

# used to parameterize new Hamiltonians
parameterize_kword = {'label'          : str,     
                      'job_type'       : str,
                      'wfn_lib'        : str, 
                      'energy_lib'     : str, 
                      'verbose'        : bool, 
                      'valence'            : str,
                      'val_opt'        : bool,
                      'val_method'     : str,
                      'val_args'       : str,
                      'val_params'     : float,
                      'val_bounds'     : float,
                      'val_freeze'     : int,
                      'cvs'            : str,
                      'cvs_opt'        : bool,
                      'cvs_method'     : str,
                      'cvs_args'       : str,
                      'cvs_params'     : float,
                      'cvs_bounds'     : float,
                      'cvs_freeze'     : int,
                      'xc'             : str, 
                      'xc_opt'         : bool,
                      'xc_params'      : float, 
                      'xc_coef'        : str, 
                      'xc_func'        : str, 
                      'xc_bounds'      : float,
                      'xc_freeze'      : int,
                      'opt_algorithm'  : str,
                      'conv'           : float, 
                      'max_iter'       : int, 
                      'ngrid'          : int, 
                      'scan_var'       : str} 

######################################################################
## only the following are directly accessed from outside the module ##
######################################################################

# these are the valid computation classes. This is somewhat
# inartful.
base_objs    = ['Molecule', 'Rydano', 'Scf', 'Parameterize']
ci_objs      = ['Dftmrci', 'Dftmrci2']
postci_objs  = ['Spinorbit']
si_objs      = ['Transition', 'Overlap', 'Dyson']
support_objs = ['Bitciwfn','Moments', 'Cigroup']
valid_objs   = base_objs + ci_objs + postci_objs + si_objs

# maximum number of processors: we set this as a distinct variable, since
# we can't depend on mpirun/mpiexec -n X, since spawn-based parallelism
# expects -n 1 in these cases
#nproc        = 1

##############################################
kwords = {'Molecule'     : molecule_kword,
          'Rydano'       : rydano_kword,
          'Parameterize' : parameterize_kword,
          'Scf'          : scf_kword,
          'Dftmrci'      : dftmrci_kword,
          'Dftmrci2'     : dftmrci2_kword,
          'Spinorbit'    : spinorbit_kword,
          'Transition'   : transition_kword,
          'Overlap'      : overlap_kword,
          'Dyson'        : dyson_kword}


