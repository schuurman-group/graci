"""
The Ao2mo object and its associated functions.
"""
import sys as sys
import os as os
import hashlib
import numpy as np
import h5py as h5py
import scipy.io as sp_io
import graci.core.libs as libs
import graci.utils.timing as timing
from pyscf import gto, ao2mo, df, lib

# Ao2mo objects holding the packed MO tensor in memory, keyed on the
# 2-e integral file name (the file is written lazily by ensure_file)
_live = {}


def mo_stamp(crds, orbs):
    """hash of a geometry and an MO coefficient array: the identity of a
    set of MO integrals (used as the reuse key of the MRSF library and of
    the per-SCF grid session)"""
    h = hashlib.blake2b(digest_size=16)
    h.update(np.ascontiguousarray(crds, dtype=float).tobytes())
    h.update(np.ascontiguousarray(orbs, dtype=float).tobytes())
    return h.hexdigest()


def ensure_file(fname):
    """write the 2-e integral file fname now if an Ao2mo object holds its
    tensor in memory without having written the file yet"""
    obj = _live.get(fname)
    if obj is not None and not obj.file_written:
        obj.write_file()
    return os.path.isfile(fname)


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
        # packed MO DF tensor (naux, nmo(nmo+1)/2) kept in memory for the
        # MRSF library (double precision only), its stamp (geometry + MOs)
        # and the state of the file/bitci consumers
        self.eri_mo       = None
        self.stamp        = None
        self.file_written = False
        self.bitci_loaded = False

    @timing.timed
    def run(self, scf, int_precision='single', load_bitci=True, write_file=True,
            keep_mem=0.0):
        """perform AO to MO integral transformation using the current
           orbitals. load_bitci: reload libbitci with the new integral
           files (not needed by MRSF-TDDFT consumers); write_file: write
           the 2-e Fortran file now (otherwise lazily, see ensure_file);
           keep_mem: keep the packed double-precision MO tensor in memory
           for the MRSF library when it is smaller than keep_mem GB"""

        import time
        self.load_scf(scf)
        self.times = {}
        t0 = time.time()

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
            cderi = scf.df_tensor()
            nao = self.orbs.shape[0]
            if cderi is not None and cderi.ndim == 2 and cderi.shape[1] == nao*(nao+1)//2:
                # the AO DF tensor of the SCF is in memory: transform it
                # (same lower-Cholesky convention as df.outcore.general)
                eri_mo = self.transform_df(cderi)
                self.times['source'] = 'SCF DF tensor'
            else:
                self.times['source'] = 'new DF build'
                ij_trans = np.concatenate(([self.orbs], 
                                           [self.orbs]))
                df.outcore.general(scf.mol.pymol(), 
                                    ij_trans,
                                    'tmp_eri',
                                    auxbasis = scf.mol.ri_basis,
                                    dataname='eri_mo')

                eri    = h5py.File('tmp_eri', 'r')
                eri_mo = np.array(eri.get('eri_mo'))
                eri.close()
                os.remove('tmp_eri')

        else:
            eri_ao = scf.mol.pymol().intor('int2e_sph', aosym='s8')
            eri_mo = ao2mo.incore.full(eri_ao, self.orbs)
            del(eri_ao)

        self.times['transform'] = time.time() - t0
        t0 = time.time()

        # the tensor stays in memory for the MRSF library when requested
        # (double precision only: the single-precision writer converts the
        # array in place); the file is written now or on demand
        self.stamp = mo_stamp(scf.mol.crds, self.orbs)
        self.eri_mo = None
        self.file_written = False
        if self.precision_2e == 'double' and eri_mo.nbytes <= keep_mem*1.0e9:
            self.eri_mo = np.ascontiguousarray(eri_mo, dtype=np.float64)
        _live[self.moint_2e_eri] = self
        if write_file or self.eri_mo is None:
            self.write_integrals(eri_mo, self.precision_2e, 
                                                     self.moint_2e_eri)
            self.file_written = True
        del(eri_mo)
        self.times['write'] = time.time() - t0
        t0 = time.time()

        # Construct the core Hamiltonian
        one_nuc_ao = scf.mol.pymol().intor('int1e_nuc')
        one_kin_ao = scf.mol.pymol().intor('int1e_kin')
        hcore_ao   = one_nuc_ao + one_kin_ao

        # transform the core Hamiltonian to MO basis explicitly 
        # and write to file
        h1_mo = np.matmul(np.matmul(self.orbs.T, hcore_ao), self.orbs)
        self.write_integrals(h1_mo, 'double', self.moint_1e)

        self.times['1e'] = time.time() - t0
        t0 = time.time()
        self.bitci_loaded = False
        if load_bitci:
            self.load_bitci(scf)
        self.times['bitci'] = time.time() - t0
        if scf.verbose:
            import graci.io.output as output
            output.print_message('MO integrals: ' + self.times['source']
                                 + ', transformation {:.2f} s, file {:.2f} s, 1e {:.2f} s, bitci {:.2f} s'.format(
                                     self.times['transform'], self.times['write'], self.times['1e'], self.times['bitci'])
                                 + (' (tensor kept in memory)' if self.eri_mo is not None else ''))

        return

    #
    def transform_df(self, cderi):
        """packed MO tensor (naux, nmo(nmo+1)/2) from the packed AO DF
        tensor cderi (naux, nao(nao+1)/2) of the SCF, in blocks of
        auxiliary functions (PySCF's C transformation, s2 -> s2)"""
        from pyscf.ao2mo import _ao2mo
        naux = cderi.shape[0]
        nmo = self.orbs.shape[1]
        out = np.empty((naux, nmo*(nmo+1)//2), dtype=np.float64)
        blk = max(1, int(2.0e8 // (8*max(cderi.shape[1], out.shape[1]))))
        for p0, p1 in lib.prange(0, naux, blk):
            _ao2mo.nr_e2(np.ascontiguousarray(cderi[p0:p1]), self.orbs, (0, nmo, 0, nmo),
                         aosym='s2', mosym='s2', out=out[p0:p1])
        return out

    #
    def write_file(self):
        """write the 2-e integral file from the tensor held in memory"""
        if self.file_written or self.eri_mo is None:
            return
        self.write_integrals(self.eri_mo, self.precision_2e, self.moint_2e_eri)
        self.file_written = True

    #
    def release(self):
        """drop the tensor held in memory (the file is written first if
        it has not been)"""
        self.write_file()
        self.eri_mo = None
        if _live.get(self.moint_2e_eri) is self:
            del _live[self.moint_2e_eri]

    #
    def load_bitci(self, scf):
        """
        Reload bitci with the current integral files
        """

        self.load_scf(scf)
        ensure_file(self.moint_2e_eri)
        self.bitci_loaded = True

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








