"""
Module for computing HF or KS-DFT orbitals and energies
Initial pass, at least, will employ PySCF
"""

import sys
import h5py
import math
import numpy as np
import graci.io.output as output
import graci.utils.timing as timing
from pyscf.lib import logger
from pyscf import gto, scf, dft, symm, ao2mo, df
from pyscf.tools import molden


class Scf:
    """Class constructor for SCF object"""
    def __init__(self):
        # user defined input paramaters
        self.xc        = 'hf'
        self.label     = 'default'

        # computed quantities
        self.energy    = None
        self.orbs      = None
        self.orb_occ   = None
        self.orb_ener  = None
        self.orb_irrep = []
        self.orb_sym   = []
        self.nmo       = 0
        self.naux      = 0

    def name(self):
        """ return the name of the class object as a string"""
        return 'scf'

    def run(self, mol):
        """compute the DFT energy and KS orbitals"""
    
        # construct the molecule object
        timing.start('scf.run')
        pymol = mol.pymol()

        output.print_scf_header()
   
        # if var.d3_inp['xc']='hf', use canonical hf orbitals
        if self.xc == 'hf': 
            if mol.use_df:
                if mol.spin == 0.:
                    mf = scf.RHF(pymol).density_fit()
                else:
                    mf = scf.ROHF(pymol).density_fit()
            else: 
                if mol.spin == 0.:
                    mf = scf.RHF(pymol)
                else:
                    mf = scf.ROHF(pymol)
                
        # this is a DFT computation
        else:
            if mol.use_df:
                if mol.spin == 0.:
                    mf = dft.RKS(pymol).density_fit()
                else:
                    mf = dft.ROKS(pymol).density_fit()
            else:
                if mol.spin == 0.:
                    mf = dft.RKS(pymol)
                else:
                    mf = dft.ROKS(pymol)
            # set the XC functional to BHLYP
            mf.xc = self.xc
        
        # save integrals
        if mol.use_df:
            mf.with_df.auxbasis       = mol.ri_basis
            # this will be generated during the calculation, 
            # saved upon completion
            mf.with_df._cderi_to_save = output.file_names['2ei']
        else:
            if self.xc != 'hf':
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
        if mol.use_df:
            ij_trans = np.concatenate(([mf.mo_coeff], [mf.mo_coeff]))
            df.outcore.general(pymol, ij_trans, output.file_names['2ei'], 
                     mol.ri_basis, dataname='eri_mo')        
        else:
            eri_ao = pymol.intor('int2e_sph', aosym='s8')
            eri_mo = ao2mo.incore.full(eri_ao, mf.mo_coeff)        
            with h5py.File(output.file_names['2ei'], 'w') as f:
                f['eri_mo'] = eri_mo
    
        # Construct the core Hamiltonian
        one_nuc_ao = pymol.intor('int1e_nuc')
        one_kin_ao = pymol.intor('int1e_kin')
        hcore_ao   = one_nuc_ao + one_kin_ao

        # transform the core Hamiltonian to MO basis explicitly 
        # and write to file
        ref_orbs   = mf.mo_coeff
        h1_mo      = np.matmul(np.matmul(ref_orbs.T, hcore_ao), ref_orbs)
        with h5py.File(output.file_names['1ei'], 'w') as f:
            f['hcore_mo'] = h1_mo

        # extract orbitals, occupations and energies
        self.orbs      = ref_orbs
        self.orb_occ   = mf.mo_occ
        self.orb_ener  = mf.mo_energy
        self.orb_irrep = orb_sym
        self.orb_sym   = orb_id
        self.nmo       = len(mf.mo_occ)
        if (mol.use_df):
            self.naux      = int(mf.with_df.auxmol.nao_nr())

        # print the summary of the output to file
        output.print_scf_summary(mol, self)

        # write the Molden file
        with open('mos.molden', 'w') as f1:
            molden.header(pymol, f1)
            molden.orbital_coeff(pymol, f1, ref_orbs, ene=mf.mo_energy,
                                 occ=mf.mo_occ)

        timing.stop('scf.run')

        return

    def density(self, mol):
        """return the density matrix for state i"""


        return

    def energy(self, state=0):
        """return the energy of state 'state'"""

        if state != 0:
            output.print_message('Calling scf.energy() for state != 0')

        return self.energy


    def slater_dets(self, state):
        """return the slater determinant list for state 'state'"""

        return

    def csfs(self, state):
        """csf list for state 'state'"""

        return

