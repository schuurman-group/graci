"""
The Ao2mo object and its associated functions.
"""
import sys as sys
import os as os
import numpy as np
import h5py as h5py
import scipy.io as sp_io
import graci.core.libs as libs
import graci.utils.timing as timing
from pyscf import gto, ao2mo, df

#

def moints_exist(scf):
    """
    Return true if the 1-e and 2-e integrals exist for correspnding
    scf object
    """
    mo_1e_eri = '1e_'+str(scf.label).strip()+'.h5'
    mo_2e_eri = '2e_eri_'+str(scf.label).strip()+'.h5'

    return os.path.isfile(mo_2e_eri) and os.path.isfile(mo_1e_eri)

class Ao2mo:
    """Class constructor for ao2mo object"""

    def __init__(self):
        self.precision_2e = None
        self.moint_2e_eri = None
        self.moint_1e     = None
        self.nmo          = None
        self.emo          = None
        self.mosym        = None
        self.emo_cut      = None
        self.orbs         = None
        self.label        = 'default'
        self.allowed_precision = ['single', 'double']

    @timing.timed
    def run(self, scf, int_precision='single'):
        """perform AO to MO integral transformation using the current
           orbitals"""

        self.load_scf(scf)

        # by default, reload bitci using newly generate MO integral
        # files
        if int_precision.strip() in self.allowed_precision:
            self.precision_2e = int_precision.strip()
        else:
            print('Integral precision not recognized: '+
                   str(int_precision.strip())+
                  ', proceeding with single precision')
            self.precision_2e = 'single'

        # Do the AO -> MO transformation
        if scf.mol.use_df:
            ij_trans = np.concatenate(([self.orbs], 
                                       [self.orbs]))
            df.outcore.general(scf.mol.pymol(), 
                                ij_trans,
                                'tmp_eri',
                                auxbasis = scf.mol.ri_basis,
                                dataname='eri_mo')

            eri    = h5py.File('tmp_eri', 'r')
            eri_mo = np.array(eri.get('eri_mo'))
            os.remove('tmp_eri')

            #df.outcore.general(scf.mol.pymol(), ij_trans, 
            #                   self.moint_2e_eri,
            #                   auxbasis = scf.mol.ri_basis, 
            #                   dataname='eri_mo')

        else:
            eri_ao = scf.mol.pymol().intor('int2e_sph', aosym='s8')
            eri_mo = ao2mo.incore.full(eri_ao, self.orbs)
            del(eri_ao)
            #eri_mo = ao2mo.incore.full(eri_ao, self.orbs)
            #with h5py.File(self.moint_2e_eri, 'w') as f:
            #    f['eri_mo'] = eri_mo

        self.write_integrals(eri_mo, self.precision_2e, 
                                                 self.moint_2e_eri)
        del(eri_mo)

        # Construct the core Hamiltonian
        one_nuc_ao = scf.mol.pymol().intor('int1e_nuc')
        one_kin_ao = scf.mol.pymol().intor('int1e_kin')
        hcore_ao   = one_nuc_ao + one_kin_ao

        # transform the core Hamiltonian to MO basis explicitly 
        # and write to file
        h1_mo = np.matmul(np.matmul(self.orbs.T, hcore_ao), self.orbs)
        self.write_integrals(h1_mo, 'double', self.moint_1e)
        #with h5py.File(self.moint_1e, 'w') as f:
        #    f['hcore_mo'] = h1_mo

        self.load_bitci(scf)

        return

    #
    def load_bitci(self, scf):
        """
        Reload bitci with the current integral files
        """

        self.load_scf(scf)

        libs.lib_func('bitci_int_finalize', [])
        if scf.mol.use_df:
            type_str = 'df'
        else:
            type_str = 'exact'

        libs.lib_func('bitci_int_initialize',
                ['pyscf', type_str, self.precision_2e, 
                           self.moint_1e, self.moint_2e_eri])

        return

    def load_scf(self, scf):
        """
        load an scf object and set internal variables based on
        scf variables
        """

        # if mo energy cutoff not set, include all orbitals
        if self.emo_cut is None:
            self.emo_cut = scf.orb_ener[-1]

        self.nmo = sum(map(lambda x : x <= self.emo_cut, scf.orb_ener))
        self.orbs    = scf.orbs[:,:self.nmo]
        self.emo     = scf.orb_ener[:self.nmo]
        self.mosym   = scf.orb_sym[:self.nmo]

        # set default file names
        self.moint_2e_eri = '2e_eri_'+str(scf.label).strip()+'.h5'
        self.moint_1e     = '1e_'+str(scf.label).strip()+'.h5'

        return

    def write_integrals(self, int_tensor, precision, file_name):
        """
        write integrals to a Fortran binary file with name
        file_name 
        """
        f = sp_io.FortranFile(file_name, 'w')

        for i in range(len(int_tensor.shape)):
            f.write_record(int_tensor.shape[i])

        out_tensor = np.transpose(int_tensor)
        if precision == 'single':
            f.write_record(np.float32(out_tensor))
        else:
            f.write_record(out_tensor) 

        f.close()

        return

