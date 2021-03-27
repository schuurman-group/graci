"""Module for storing global parameters"""

# molecules section input keywords and data types
mol_kword =  {'units'       : str,
              'charge'      : int,
              'mult'        : int,
              'basis'       : str,
              'ri_basis'    : str,
              'use_sym'     : bool,
              'use_df'      : bool
             }

mol_param = {'units'          : 'bohr',
              'charge'         : 0,
              'mult'           : 1,
              'basis'          : 'sto-3g',
              'ri_basis'       : '',
              'use_sym'        : True,
             }

# DFT section input keywords and data types
scf_kword =  {'xc'          : str}

scf_param =  {'xc'          : 'hf'}

# MRCI section input keywords and data types
mrci_kword = {'nstates'        : int,
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
              'interface_only' : bool,
              'diag_algorithm' : str,
              'diag_tol'       : float,
              'diag_iter'      : int,
              'diag_blocksize' : int,
              'diag_deflate'   : bool
             }

mrci_param = {'nstates'        : [],
              'hamiltonian'    : 'canonical',
              'de_select'      : 0.8,
              'ras1'           : [],
              'ras2'           : [],
              'ras3'           : [],
              'nhole1'         : 0,
              'nelec3'         : 0,
              'autoras'        : False,
              'ciorder'        : 2,
              'refiter'        : 5,
              'asci'           : 'off',
              'interface_only' : False,
              'diag_algorithm' : 'gendav',
              'diag_tol'       : 0.0001,
              'diag_iter'      : 50,
              'diag_blocksize' : [],
              'diag_deflate'   : False
             }


# ASCI-style configuration selection thresholds
asci_kword = {'off'   : bool,
              'tight'  : float,
              'normal' : float,
              'loose'  : float}

asci_param = {'off'   : False, 
              'tight'  : 1e-5,
              'normal' : 1e-4,
              'loose'  : 1e-3}

# keywords related to file names 
file_kword  = {'input_file'    : str,
               'out_file'       : str,
               'pyscf_out'      : str,
               '1ei'            : str,
               '2ei'            : str
              }

file_param = {'input_file'     : '',
               'out_file'       : '',
               'pyscf_out'      : '',
               '1ei'            : '',
               '2ei'            : ''
              }

# MRCI and DFT/MRCI Hamiltonian labels
hamiltonians = ['canonical',
                'grimme_standard',
                'grimme_short',
                'lyskov_standard',
                'lyskov_short',
                'heil17_standard',
                'heil18_short']
