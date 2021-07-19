"""
Module for computing HF or KS-DFT orbitals and energies
Initial pass, at least, will employ PySCF
"""

import os
import sys
import h5py
import math
import numpy as np
import importlib
import graci.io.output as output
import graci.utils.timing as timing
from pyscf.lib import logger
from pyscf import gto, scf, dft, symm, ao2mo, df
from pyscf.tools import molden


class Scf:
    """Class constructor for SCF object"""
    def __init__(self):
        # user defined input paramaters
        self.xc             = 'hf'
        self.print_orbitals = False
        self.label          = 'Scf'
 
        # computed quantities
        self.mol          = None 
        self.energy       = None
        self.orbs         = None
        self.orb_occ      = None
        self.orb_ener     = None
        self.orb_irrep    = []
        self.orb_sym      = []
        self.nmo          = 0
        self.nao          = 0
        self.naux         = 0
        self.rdm_ao       = None
        self.moint_1e     = None
        self.moint_2e_eri = None
 

# Required functions #######################################################

    def set_mol(self, mol):
        """set the mol object for the scf object"""
        self.mol = mol

    def mol_exists(self):
        """return True is self.mol is not None"""
        try:
            return type(self.mol).__name__ is 'Molecule'
        except:
            return False

    @timing.timed
    def run(self):
        """compute the DFT energy and KS orbitals"""
    
        # construct the molecule object
        pymol = self.mol.pymol()

        output.print_scf_header()

        # if var.d3_inp['xc']='hf', use canonical hf orbitals
        if self.xc == 'hf': 
            if self.mol.use_df:
                if self.mol.spin == 0.:
                    mf = scf.RHF(pymol).density_fit()
                else:
                    mf = scf.ROHF(pymol).density_fit()
            else: 
                if self.mol.spin == 0.:
                    mf = scf.RHF(pymol)
                else:
                    mf = scf.ROHF(pymol)
                
        # this is a DFT computation
        else:
            if self.mol.use_df:
                if self.mol.spin == 0.:
                    mf = dft.RKS(pymol).density_fit()
                else:
                    mf = dft.ROKS(pymol).density_fit()
            else:
                if self.mol.spin == 0.:
                    mf = dft.RKS(pymol)
                else:
                    mf = dft.ROKS(pymol)
            # set the XC functional to BHLYP
            mf.xc = self.xc
        
        # save integrals -- tie them to the mol object for
        # this file
        self.moint_1e     = '1e_'+str(self.mol.label).strip()+'.h5'
        self.moint_2e_eri = '2e_eri_'+str(self.mol.label).strip()+'.h5'

        if self.mol.use_df:
            mf.with_df.auxbasis       = self.mol.ri_basis
            # this will be generated during the calculation,
            # saved upon completion
            mf.with_df._cderi_to_save = self.moint_2e_eri
        else:
            if self.xc != 'hf' and hasattr(mf, 'with_df'):
                mf.with_df = False    

        # run the dft computation
        self.energy = mf.kernel()

        # if not converged, kill things
        if not mf.converged:
            sys.exit('Reference SCF computation did not converge.')
        
        # orb_sym are the symmetry labels 
        if pymol.symmetry: 
            orb_sym = symm.label_orb_symm(pymol, pymol.irrep_name, 
                                   pymol.symm_orb, mf.mo_coeff)
            orb_id  = symm.label_orb_symm(pymol, pymol.irrep_id,         
                                   pymol.symm_orb, mf.mo_coeff)
        else:
            orb_sym = ['a'] * len(mf.mo_coeff)
            orb_id  = [0] * len(mf.mo_coeff)    

        # Do the AO -> MO transformation
        if self.mol.use_df:
            ij_trans = np.concatenate(([mf.mo_coeff], [mf.mo_coeff]))
            df.outcore.general(pymol, ij_trans, 
                    self.moint_2e_eri, 
                    self.mol.ri_basis, dataname='eri_mo')
        else:
            eri_ao = pymol.intor('int2e_sph', aosym='s8')
            eri_mo = ao2mo.incore.full(eri_ao, mf.mo_coeff)
            with h5py.File(self.moint_2e_eri, 'w') as f:
                f['eri_mo'] = eri_mo
    
        # Construct the core Hamiltonian
        one_nuc_ao = pymol.intor('int1e_nuc')
        one_kin_ao = pymol.intor('int1e_kin')
        hcore_ao   = one_nuc_ao + one_kin_ao

        # transform the core Hamiltonian to MO basis explicitly 
        # and write to file
        ref_orbs   = mf.mo_coeff
        h1_mo      = np.matmul(np.matmul(ref_orbs.T, hcore_ao), ref_orbs)
        with h5py.File(self.moint_1e, 'w') as f:
            f['hcore_mo'] = h1_mo

        # extract orbitals, occupations and energies
        self.orbs      = ref_orbs
        # orb_occ are the MO occupation numbers
        self.orb_occ   = mf.mo_occ
        # orb_ener is the array of MO energies
        self.orb_ener  = mf.mo_energy
        # orb_irrep is the irrep label -- this dirty: easier for
        # checkpointing if this is NOT a numpy array of strings
        self.orb_irrep = [orb_sym[i] for i in range(len(orb_sym))]
        # orb_sym is the irrep index
        self.orb_sym   = orb_id
        # number of MOs
        self.nmo       = len(mf.mo_occ)

        if (self.mol.use_df):
            self.naux      = int(mf.with_df.auxmol.nao_nr())

        # construct density matrix
        self.rdm_ao = np.zeros((self.nmo, self.nmo), dtype=float)
        for i in range(self.nmo):
            self.rdm_ao += self.orb_occ[i]*np.outer(self.orbs[:,i], \
                                                    self.orbs[:,i])

        # print the summary of the output to file
        output.print_scf_summary(self)

        # write the Molden file if requested
        if self.print_orbitals:
            self.export_orbitals()

        return

    #
    def n_state(self):
        """return number of states optmizied, i.e. 1"""
        return 1

    #
    def n_state_sym(self, irrep):
        """number of states computed"""

        if irrep == self.state_sym:
            return 1
        else:
            return 0

    #
    def energy(self):
        """return the energy of state 'state'"""

        return self.energy

    #
    def rdm1(self, basis='mo'):
        """return the density matrix for state i"""
        
        if basis == 'mo':
            return np.diag(self.orb_occ) 
        else:
            return self.rdm_ao

    #
    def natural_orbs(self, basis='ao'):
        """return hf orbitals"""

        if basis == 'ao':
            return self.orbs
        else:
            return np.identity(self.nmo)

    # 
    def export_orbitals(self, file_format='molden', 
                                   orb_dir=True, cart=True):
        """export orbitals to molden format"""

        # append the file_format label to file name
        fname = 'mos_' + str(file_format).lower() + '_' + str(self.label)

        if orb_dir:
            fname = 'orbs/'+str(fname)
            if not os.path.exists('orbs'):
                os.mkdir('orbs')

        # if a file called fname exists, delete it
        if os.path.isfile(fname):
            os.remove(fname)

        # import the appropriate library for the file_format
        if file_format in output.orb_formats:
            orbtype = importlib.import_module('graci.io.'+file_format)
        else:
            print('orbital format type=' + file_format +
                                        ' not found. exiting...')
            sys.exit(1)

        orbtype.write_orbitals(fname, 
                               self.mol, 
                               self.orbs, 
                               sym_lbl = self.orb_irrep, 
                               ener = self.orb_ener, 
                               cart = cart)

        return

    #
    def slater_dets(self, state):
        """return the slater determinant list for state 'state'"""

        return

    #
    def csfs(self, state):
        """csf list for state 'state'"""

        return

