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
                      'rrdf'     : bool,
                      'label'    : str
                     }

# used to parameterize new Hamiltonians
parameterize_kword = {'algorithm' : str,
                      'label'     : str}

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
                  'prune'          : str,
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

# moments section, input keywords and data types
moments_kword   = {'states'        : int
                  }

# transition moments input keywords and data types
transition_kword = {'init'         : int,
                    'final'        : int
                   }

# molecule section input keywords and data typess
spinorbit_kword = {'bra_states' : int,
                   'ket_states' : int
                  }


######################################################################
## only the following are directly accessed from outside the module ##
######################################################################

def name2obj(name):
    if name == 'geometry':
        return geometry.Geometry()
    elif name == 'molecule':
        return molecule.Molecule()
    elif name == 'parameterize':
        return parameterize.Parameterize()
    elif name == 'scf':
        return scf.Scf()
    elif name == 'dftcis':
        return dftcis.Dftcis()
    elif name == 'dftmrci':
        return dftmrci.Dftmrci()
    elif name == 'moments':
        return moments.Moments()
    elif name == 'transition':
        return transition.Transition()
    elif name == 'spinorbit':
        return spinorbit.Spinorbit()
    else:
        print('obj: '+str(name)+' not recognized')
        return None

# these are the valid computation classes. This is somewhat
# inartful.
valid_objs = ['geometry', 'molecule', 'parameterize', 
              'scf', 'dftmrci', 'dftmrcis', 
              'moments', 'transition', 'spinorbit']


##############################################
kwords = {'geometry'     : geometry_kword,
          'molecule'     : molecule_kword,
          'parameterize' : parameterize_kword,
          'scf'          : scf_kword,
          'dftcis'       : dftcis_kword,
          'dftmrci'      : dftmrci_kword,
          'moments'      : moments_kword,
          'transition'   : transition_kword,
          'spinorbit'    : spinorbit_kword
         }

