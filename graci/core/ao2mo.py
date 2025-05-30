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
 
        # we will keep records to <= 2GB in order to maintain compatability
        # with FortranFile, as it does not support subrecords for standard
        # compilers
        nrecs, cprec = self.n_records(int_tensor.shape, precision)

        # write dimensions of the tensor
        for i in range(len(int_tensor.shape)):
            f.write_record(int_tensor.shape[i])

        # write information regarding record sizes
        f.write_record(nrecs)
        f.write_record(cprec)

        # always write entire columns of data (rend constant)
        rend = int_tensor.shape[0]

        # for single precision we want to ensure we do not inadvertently
        # copy data. So: we create a float32 view, then step through 
        # every other column index since the data now occupies half the 
        # space. This would be ::4 if we implement half precision. NOTE: 
        # the transpose view must be create AFTER float32 view, else
        # scipy/numpy squawks about Fortran ordering
        if precision == 'single':
            out_tensor        = int_tensor.view(np.float32)
            out_tensor[:,::2] = int_tensor
            print_tensor      = out_tensor[:,::2].T
            for j in range(nrecs):
                cend = min((j+1)*cprec, int_tensor.shape[1])
                f.write_record(print_tensor[j*cprec:cend, 0:rend])

        # double precision straightforward: note that taking transpose 
        # just creates a new view and no data is copied.
        else:
            out_tensor = int_tensor.T
            for j in range(nrecs):
                cend = min((j+1)*cprec, int_tensor.shape[1])
                f.write_record(out_tensor[j*cprec:cend, 0:rend])

        f.close()

        return

    def n_records(self, tensor_dims, precision):
        """
        determine the number of columns to write at once to ensure the
        record length remains < 2GB
        """

        # we should revisit this in the future: the maximum record size
        # is currently, arbitrarily, set to 268435456*0.9 = 241591910
        #  double precision numbers
        #nfp = 268435456*0.9
        nfp = 241591910

        if precision == 'single':
            nfp *= 2
           
        nelem = np.prod(tensor_dims)
        nrec  = int(np.ceil(nelem / nfp))
        cpr   = int(np.ceil(tensor_dims[1] / nrec))

        return nrec, cpr








