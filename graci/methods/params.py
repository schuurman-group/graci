"""Module for storing global parameters"""

import graci.methods.molecule as molecule
import graci.methods.scf as scf
import graci.methods.dftcis as dftcis
import graci.methods.dftmrci as dftmrci
import graci.methods.cvsdftmrci as cvsdftmrci

# molecule section input keywords and data types
molecule_kword = {'units'   : str,
                  'charge'   : int,
                  'mult'     : int,
                  'basis'    : str,
                  'ri_basis' : str,
                  'use_sym'  : bool,
                  'use_df'   : bool,
                  'name'     : str
                  }

# DFT section input keywords and data types
scf_kword      = {'xc'          : str,
                  'name'        : str}

# MRCI section input keywords and data types
dftmrci_kword  = {'nstates'        : int,
                  'hamiltonian'    : str,
                  'de_select'      : float,
                  'ras1'           : int,
                  'ras2'           : int,
                  'ras3'           : int,
                  'nhole1'         : int,
                  'nelec3'         : int,
                  'autoras'        : bool,
                  'ciorder'        : int,
                  'refiter'        : int,
                  'asci'           : str,
                  'diag_method'    : str,
                  'diag_tol'       : float,
                  'diag_iter'      : int,
                  'diag_blocksize' : int,
                  'diag_deflate'   : bool,
                  'name'           : str
                 }

cvsdftmrci_kword =  {'nstates'        : int,
                  'hamiltonian'    : str,
                  'de_select'      : float,
                  'ras1'           : int,
                  'ras2'           : int,
                  'ras3'           : int,
                  'nhole1'         : int,
                  'nelec3'         : int,
                  'autoras'        : bool,
                  'ciorder'        : int,
                  'refiter'        : int,
                  'asci'           : str,
                  'diag_method'    : str,
                  'diag_tol'       : float,
                  'diag_iter'      : int,
                  'diag_blocksize' : int,
                  'diag_deflate'   : bool,
                  'name'           : str
                 }

dftcis_kword   = {'nstates'        : int,
                  'hamiltonian'    : str,
                  'de_select'      : float,
                 }
######################################################################
## only the following are directly accessed from outside the module ##
######################################################################

def method_obj(name):
    if name == 'molecule':
        return molecule.Molecule()
    elif name == 'scf':
        return scf.Scf()
    elif name == 'dftcis':
        return dftcis.Dftcis()
    elif name == 'dftmrci':
        return dftmrci.Dftmrci()
    elif name == 'cvsdftmrci':
        return cvsdftmrci.Cvsdftmrci()
    else:
        print('method_obj: '+str(name)+' not recognized')
        return None

def method_name(obj):
    if isinstance(obj, molecule.Molecule):
        return 'molecule'
    elif isinstance(obj, scf.Scf):
        return 'scf'
    elif isinstance(obj, dftcis.Dftcis):
        return 'dftcis'
    elif isinstance(obj, dftmrci.Dftmrci):
        return 'dftmrci'
    elif isinstance(obj, cvsdftmrci.Cvsdftmrci):
        return 'cvsdrtmrci'
    else:
        print(str(obj)+' not recognized as a method object')
        return None

# these are the valid computation classes. This is somewhat
# inartful.
valid_methods = ['molecule', 'scf', 'dftmrci', 'cvsdftmrci', 'dftmrcis']

##############################################
kwords = {'molecule'   : molecule_kword,
          'scf'        : scf_kword,
          'dftcis'     : dftcis_kword,
          'dftmrci'    : dftmrci_kword,
          'cvsdftmrci' : cvsdftmrci_kword,
         }

# MRCI and DFT/MRCI Hamiltonian labels
hamiltonians   = ['canonical',
                  'grimme_standard',
                  'grimme_short',
                  'lyskov_standard',
                  'lyskov_short',
                  'heil17_standard',
                  'heil18_short']

