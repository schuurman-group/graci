"""Module for storing global parameters"""

import graci.core.molecule as molecule
import graci.core.geometry as geometry
import graci.core.parameterize as parameterize
import graci.methods.scf as scf
import graci.methods.dftmrci as dftmrci
import graci.methods.dftcis as dftcis
import graci.properties.moments as moments
import graci.properties.transition as transition
import graci.properties.spinorbit as spinorbit

# geometry section, input keywords and data types
geometry_kword =     {'xyz_file' : str,
                      'units'    : str,
                      'label'    : str}

# molecule section input keywords and data types
molecule_kword =     {'charge'   : int,
                      'mult'     : int,
                      'basis'    : str,
                      'ri_basis' : str,
                      'use_sym'  : bool,
                      'use_df'   : bool,
                      'use_rrdf' : bool,
                      'rrdf_fac' : int,
                      'label'    : str
                     }

#----------------------------------------------

# DFT section input keywords and data types
scf_kword      = {'xc'          : str,
                  'label'       : str}

# MRCI section input keywords and data types
dftmrci_kword  = {'nstates'        : int,
                  'hamiltonian'    : str,
                  'de_select'      : float,
                  'ras1'           : int,
                  'ras2'           : int,
                  'ras3'           : int,
                  'nhole1'         : int,
                  'nelec3'         : int,
                  'icvs'           : int,
                  'autoras'        : bool,
                  'ciorder'        : int,
                  'refiter'        : int,
                  'ref_prune'      : bool,
                  'prune'          : str,
                  'prune_extra'    : int,
                  'diag_guess'     : str,
                  'diag_method'    : str,
                  'diag_tol'       : float,
                  'diag_iter'      : int,
                  'diag_blocksize' : int,
                  'diag_deflate'   : bool,
                  'label'          : str
                 }

# DFT/CIS section input keywords and data types
dftcis_kword   = {'nstates'        : int,
                  'hamiltonian'    : str,
                  'de_select'      : float,
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
                    'final_label'     : str
                   }

# molecule section input keywords and data typess
spinorbit_kword = {'bra_states' : int,
                   'ket_states' : int,
                   'bra_label'  : str,
                   'ket_label'  : str
                   }

#----------------------------------------------------

# used to parameterize new Hamiltonians
parameterize_kword = {'algorithm' : str,
                      'label'     : str}


######################################################################
## only the following are directly accessed from outside the module ##
######################################################################

def name2obj(name):
    if name == 'Geometry':
        return geometry.Geometry()
    elif name == 'Molecule':
        return molecule.Molecule()
    elif name == 'Bitciwfn':
        return bitciwfn.Bitciwfn()
    elif name == 'Parameterize':
        return parameterize.Parameterize()
    elif name == 'Scf':
        return scf.Scf()
    elif name == 'Dftcis':
        return dftcis.Dftcis()
    elif name == 'Dftmrci':
        return dftmrci.Dftmrci()
    elif name == 'Moments':
        return moments.Moments()
    elif name == 'Transition':
        return transition.Transition()
    elif name == 'Spinorbit':
        return spinorbit.Spinorbit()
    else:
        print('obj: '+str(name)+' not recognized')
        return None

# these are the valid computation classes. This is somewhat
# inartful.
valid_objs = ['Geometry', 'Molecule', 'Bitciwfn', 'Parameterize', 
              'Scf', 'Dftmrci', 'Dftmrcis', 
              'Transition', 'Spinorbit']
method_objs = ['Scf', 'Dftmrci', 'Dftmrcis']
si_objs     = ['Transition', 'Spinorbit']

##############################################
kwords = {'Geometry'     : geometry_kword,
          'Molecule'     : molecule_kword,
          'Parameterize' : parameterize_kword,
          'Scf'          : scf_kword,
          'Dftcis'       : dftcis_kword,
          'Dftmrci'      : dftmrci_kword,
          'Transition'   : transition_kword,
          'Spinorbit'    : spinorbit_kword
         }

