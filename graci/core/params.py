"""Module for storing global parameters"""

import graci.core.molecule as molecule
import graci.core.parameterize as parameterize
import graci.methods.scf as scf
import graci.methods.dftmrci as dftmrci
import graci.methods.dftmrenpt2 as dftmrenpt2
import graci.properties.moments as moments
import graci.interaction.transition as transition
import graci.interaction.spinorbit as spinorbit

# molecule section input keywords and data types
molecule_kword =     {'xyz_file' : str,
                      'units'    : str,
                      'basis'    : str,
                      'ri_basis' : str,
                      'ao_cart'  : bool,
                      'use_sym'  : bool,
                      'use_df'   : bool,
                      'mult'     : int,
                      'charge'   : int,
                      'label'    : str
                     }

#----------------------------------------------

# DFT section input keywords and data types
scf_kword      = {'xc'             : str,
                  'restart'        : bool,
                  'print_orbitals' : bool,
                  'verbose'        : int,
                  'charge'         : int,
                  'mult'           : int,
                  'label'          : str}

# MRCI section input keywords and data types
dftmrci_kword  = {'mult'           : int,
                  'charge'         : int,
                  'nstates'        : int,
                  'hamiltonian'    : str,
                  'ras1'           : int,
                  'ras2'           : int,
                  'ras3'           : int,
                  'nhole1'         : int,
                  'nelec3'         : int,
                  'icvs'           : int,
                  'autoras'        : bool,
                  'refiter'        : int,
                  'ref_prune'      : bool,
                  'prune'          : bool,
                  'prune_thresh'   : float,
                  'prune_qcorr'    : bool,
                  'prune_extra'    : int,
                  'diag_guess'     : str,
                  'diag_method'    : str,
                  'diag_tol'       : float,
                  'diag_iter'      : int,
                  'diag_blocksize' : int,
                  'diag_deflate'   : bool,
                  'print_orbitals' : bool,
                  'save_wf'        : bool,
                  'label'          : str
                 }

# DFT/MR-ENPT2 section input keywords and data types
dftmrenpt2_kword  = {'mult'           : int,
                     'charge'         : int,
                     'nstates'        : int,
                     'hamiltonian'    : str,
                     'ras1'           : int,
                     'ras2'           : int,
                     'ras3'           : int,
                     'nhole1'         : int,
                     'nelec3'         : int,
                     'icvs'           : int,
                     'autoras'        : bool,
                     'refiter'        : int,
                     'ref_prune'      : bool,
                     'print_orbitals' : bool,
                     'save_wf'        : bool,
                     'label'          : str
                     }

#---------------------------------------------------

# transition moments input keywords and data types
transition_kword = {'init_states'      : int, 
                    'final_states'     : int,
                    'init_states_sym'  : float,
                    'final_states_sym' : float,
                    'all_final_states' : bool,
                    'init_label'       : str,
                    'final_label'      : str,
                    'print_orbitals'   : bool,
                    'label'            : str
                   }

# spinorbit section input keywords and data typess
spinorbit_kword = {'init_states'      : int, 
                   'final_states'     : int,
                   'init_states_sym'  : float,
                   'final_states_sym' : float,
                   'all_final_states' : bool,
                   'init_label'       : str,
                   'final_label'      : str,
                   'print_thresh'     : float,
                   'mf2e'             : str,
                   'label'            : str
                   }

# overlap section input keywords and data types
overlap_kword = {'bra_states' : int,
                 'ket_states' : int,
                 'bra_label'  : str,
                 'ket_label'  : str,
                 'calc'       : str
                 }

#----------------------------------------------------

# used to parameterize new Hamiltonians
parameterize_kword = {'algorithm' : str,
                      'label'     : str}

######################################################################
## only the following are directly accessed from outside the module ##
######################################################################

# these are the valid computation classes. This is somewhat
# inartful.
valid_objs   = ['Molecule', 'Parameterize','Scf', 'Dftmrci',
                'Dftmrenpt2', 
                'Transition', 'Spinorbit',
                'Overlap']
method_objs  = ['Scf', 'Dftmrci', 'Dftmrenpt2']
si_objs      = ['Transition', 'Spinorbit']
overlap_objs = ['Overlap']

##############################################
kwords = {'Molecule'     : molecule_kword,
          'Parameterize' : parameterize_kword,
          'Scf'          : scf_kword,
          'Dftmrci'      : dftmrci_kword,
          'Dftmrenpt2'   : dftmrenpt2_kword,
          'Transition'   : transition_kword,
          'Spinorbit'    : spinorbit_kword,
          'Overlap'      : overlap_kword
         }

