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
from pyscf import lib as pyscf_lib

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
        self.precision_2e = 'single' 
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
            # A fixed temp file name, an h5py handle that is never
            # closed, and os.remove() called while that handle is still
            # open. On a second SCF object in the same process PySCF
            # writes a new 'tmp_eri' while a live HDF5 handle still
            # refers to the previous, unlinked inode, and the read back
            # returns a mixture: on a def2-TZVPD stilbene the first
            # 14194 of 34453 columns were correct and the rest garbage.
            # Unique name per call, and close before unlinking.
            tmp_eri = 'tmp_eri_%s_%d' % (str(scf.label), os.getpid())

            # [AO2MO-IN] is the INPUT to the transformation sound? If the
            # orbitals and the molecule check out here but eri_mo comes
            # back corrupt, the fault is inside df.outcore.general or the
            # h5py read; if orbs is already wrong, it is upstream in
            # load_scf / the Scf object.
            _pm = scf.mol.pymol()
            print(' [AO2MO-IN] scf=%-12s orbs%s sum|orbs|=%.17e'
                  ' emo_cut=%s nmo=%s | pymol id=%s natm=%d nbas=%d'
                  ' nelec=%s charge=%s'
                  % (str(scf.label), str(self.orbs.shape),
                     float(np.abs(self.orbs).sum()), str(self.emo_cut),
                     str(self.nmo), hex(id(_pm)), _pm.natm, _pm.nbas,
                     str(_pm.nelec), str(_pm.charge)), flush=True)

            # The DF transformation is a RACE once a CI calculation has
            # run in this process: two back-to-back calls with identical
            # inputs give different wrong answers (6.91e4 vs 8.09e4
            # against a correct 3.10e4), while the same call before any
            # CI is bit-reproducible, and the same call in a standalone
            # script is reproducible 3/3 even with libbitci loaded.
            # bitci is built against libiomp5 and pyscf's libao2mo/
            # libnp_helper link BOTH libgomp and libiomp5, so once
            # bitci's OpenMP pool is live the transform's threaded
            # regions are no longer safe.
            #
            # Serialise the transformation, as mkl_compat.f90 already
            # does for the overlap code (commit 80f0239). Costs
            # transform wall time; buys a correct answer.
            _nthr = pyscf_lib.num_threads()
            pyscf_lib.num_threads(1)
            try:
                df.outcore.general(scf.mol.pymol(), 
                                    ij_trans,
                                    tmp_eri,
                                    auxbasis = scf.mol.ri_basis,
                                    dataname='eri_mo')
            finally:
                pyscf_lib.num_threads(_nthr)

            # [AO2MO-RD] the inputs are sound and eri_mo comes back
            # corrupt, so the fault is in df.outcore.general or in this
            # read. On stilbene the dataset is 349 MB and is pulled in a
            # single np.array() call; read it BOTH ways and compare. If
            # the blocked sum is right and the bulk sum is wrong, the
            # single large read is at fault, not the transformation.
            with h5py.File(tmp_eri, 'r') as eri:
                dset = eri['eri_mo']
                nrow = dset.shape[0]
                blk  = max(1, nrow // 16)
                s_blk = 0.0
                m_blk = 0.0
                for i0 in range(0, nrow, blk):
                    chunk = dset[i0:min(i0+blk, nrow)]
                    s_blk += float(np.abs(chunk).sum())
                    m_blk  = max(m_blk, float(np.abs(chunk).max()))
                eri_mo = np.array(dset)

            s_bulk = float(np.abs(eri_mo).sum())
            print(' [AO2MO-RD] scf=%-12s dset%s %.1f MB | blocked sum=%.17e'
                  ' max=%.6e | bulk sum=%.17e max=%.6e | %s'
                  % (str(scf.label), str(dset.shape),
                     eri_mo.nbytes/1024**2, s_blk, m_blk,
                     s_bulk, float(np.abs(eri_mo).max()),
                     'AGREE' if abs(s_blk-s_bulk) <= 1e-6*abs(s_blk)
                     else '*** BULK READ DIFFERS ***'), flush=True)

            # [AO2MO-2X] TEMPORARY: repeat the transformation with
            # byte-identical inputs and compare. Same wrong answer twice
            # => corrupted process state (deterministic). Different wrong
            # answers => a race. df.outcore.general is deterministic in
            # isolation (dftest.py, 3/3 SAME on hartree at 8 threads),
            # so whatever breaks it is something GRaCI's process carries.
            tmp2 = tmp_eri + '_2x'
            df.outcore.general(scf.mol.pymol(), ij_trans, tmp2,
                               auxbasis=scf.mol.ri_basis,
                               dataname='eri_mo')
            with h5py.File(tmp2, 'r') as _e2:
                _a2 = np.array(_e2['eri_mo'])
            os.remove(tmp2)
            _s2 = float(np.abs(_a2).sum())
            print(' [AO2MO-2X] scf=%-12s pass1=%.17e pass2=%.17e  %s'
                  % (str(scf.label), s_bulk, _s2,
                     'IDENTICAL -> deterministic (process state)'
                     if abs(_s2-s_bulk) <= 1e-12*abs(s_bulk)
                     else '*** DIFFERENT -> race ***'), flush=True)
            del _a2

            os.remove(tmp_eri)

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

        # [AO2MO] TEMPORARY DIAGNOSTIC -- remove with [INTCHK]/[READCHK].
        # bitci reads sum|bra_ket| = 3.10e4 for the first CI calculation
        # and 8.83e4 for the second, from files with identical dimensions.
        # This prints the tensor as PySCF produced it, BEFORE it is
        # written, so the corruption can be placed on one side or the
        # other of write_integrals without keeping the files.
        print(' [AO2MO] scf=%-12s shape=%-16s dtype=%-9s'
              ' sum|eri_mo|=%.17e max=%.6e nonfinite=%d'
              % (str(scf.label), str(eri_mo.shape), str(eri_mo.dtype),
                 float(np.abs(eri_mo).sum()), float(np.abs(eri_mo).max()),
                 int((~np.isfinite(eri_mo)).sum())), flush=True)

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








