"""
The Ao2mo object and its associated functions.
"""
import sys as sys
import numpy as np
import h5py as h5py
import graci.core.libs as libs
import graci.utils.timing as timing
from pyscf import gto, ao2mo, df

#
class Ao2mo:
    """Class constructor for ao2mo object"""

    def __init__(self, precision='double'):
        self.precision    = precision
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
    def run(self, scf, int_precision=None):
        """perform AO to MO integral transformation using the current
           orbitals"""

        self.load_scf(scf)

        # Do the AO -> MO transformation
        if scf.mol.use_df:
            ij_trans = np.concatenate(([self.orbs], 
                                       [self.orbs]))
            df.outcore.general(scf.mol.pymol(), ij_trans, 
                               self.moint_2e_eri,
                               auxbasis = scf.mol.ri_basis, 
                               dataname='eri_mo')
        else:
            eri_ao = scf.mol.pymol().intor('int2e_sph', aosym='s8')
            eri_mo = ao2mo.incore.full(eri_ao, self.orbs)
            with h5py.File(self.moint_2e_eri, 'w') as f:
                f['eri_mo'] = eri_mo
            del(eri_mo)

        # Construct the core Hamiltonian
        one_nuc_ao = scf.mol.pymol().intor('int1e_nuc')
        one_kin_ao = scf.mol.pymol().intor('int1e_kin')
        hcore_ao   = one_nuc_ao + one_kin_ao

        # transform the core Hamiltonian to MO basis explicitly 
        # and write to file
        h1_mo = np.matmul(np.matmul(self.orbs.T, hcore_ao), self.orbs)
        with h5py.File(self.moint_1e, 'w') as f:
            f['hcore_mo'] = h1_mo

        # by default, reload bitci using newly generate MO integral
        # files
        if (int_precision is not None and 
                         int_precision in self.allowed_precision):
            self.precision = int_precision
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
                ['pyscf', type_str, self.precision, 
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
