#!/usr/bin/env python
# Copyright 2014-2022 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

'''
Generalized Kohn-Sham
'''


import numpy
from functools import reduce
import scipy.linalg
from pyscf import lib
from pyscf.lib import logger
from pyscf.scf import ghf
#from pyscf.dft import rks
from graci.pdft import rks
from graci.pdft import project
from graci.pdft.rks import prune_small_rho_grids_
from pyscf.dft.numint2c import NumInt2C


def get_veff(ks, mol=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
    '''Coulomb + XC functional

    .. note::
        This function will change the ks object.

    Args:
        ks : an instance of :class:`RKS`
            XC functional are controlled by ks.xc attribute.  Attribute
            ks.grids might be initialized.
        dm : ndarray or list of ndarrays
            A density matrix or a list of density matrices

    Kwargs:
        dm_last : ndarray or a list of ndarrays or 0
            The density matrix baseline.  If not 0, this function computes the
            increment of HF potential w.r.t. the reference HF potential matrix.
        vhf_last : ndarray or a list of ndarrays or 0
            The reference Vxc potential matrix.
        hermi : int
            Whether J, K matrix is hermitian

            | 0 : no hermitian or symmetric
            | 1 : hermitian
            | 2 : anti-hermitian

    Returns:
        matrix Veff = J + Vxc.  Veff can be a list matrices, if the input
        dm is a list of density matrices.
    '''
    if mol is None: mol = ks.mol
    if dm is None: dm = ks.make_rdm1().real

    #ks.initialize_grids(mol, dm)
    ext_basis = ks.ext_basis

    t0 = (logger.process_clock(), logger.perf_counter())

    ground_state = isinstance(dm, numpy.ndarray) and dm.ndim == 2

    # Build the AO projection operators and projected 1pdm
    if ks._SQQS is None:
        ks._build_proj() ## sets ks.SQQS object
    SQ = ks._SQQS[0]; QS = ks._SQQS[1]

    pdm = reduce(numpy.dot, (QS,dm,SQ)).real

    if ks.grids.coords is None:
        ks.grids.build(with_non0tab=True)
        if ks.small_rho_cutoff > 1e-20 and ground_state:
            # Filter grids the first time setup grids
            ks.grids = prune_small_rho_grids_(ks, mol, dm, ks.grids)
    #    #t0 = logger.timer(ks, 'setting up grids', *t0)

    if ks.nlc != '':
        if ks.nlcgrids.coords is None:
            ks.nlcgrids.build(with_non0tab=True)
            if ks.small_rho_cutoff > 1e-20 and ground_state:
                # Filter grids the first time setup grids
                ks.nlcgrids = prune_small_rho_grids_(ks, mol, dm, ks.nlcgrids)
            #t0 = logger.timer(ks, 'setting up nlc grids', *t0)

    ni = ks._numint
    # Enable Range-Separated Hybrids
    omega, alpha, hyb = ni.rsh_and_hybrid_coeff(ks.xc, spin=mol.spin)

    if hermi == 2:  # because rho = 0
        n, exc, vxc = 0, 0, 0
    else:
        ## Regular XC contributions.
        ##        i.e., vxc excludes fraction of EEX.
        max_memory = ks.max_memory - lib.current_memory()[0]
        n, exc, vxc = ni.get_vxc(mol, ks.grids, ks.xc, dm,
                                 hermi=hermi, max_memory=max_memory)
        #print("VXC:",vxc.dtype)
        if vxc.dtype is not numpy.complex128:
            vxc = vxc.astype(numpy.complex128)
        logger.debug(ks, 'nelec by numeric integration = %s', n)

        if(abs(ks.phyb)>1e-10):
          pxc = ks.xc

          #### Identify XC functional components
          if ks.xcstr is not None:
              px, pc = ks.xcstr
          else:
              px = pxc

          ## Projected VXC in Projected RDM1.
          np, excp, vxcp0 = ni.get_vxc(mol, ks.grids, px, pdm,
                                       hermi=hermi, max_memory=max_memory)
          vxcp = numpy.einsum('ik,kj->ij',SQ,numpy.einsum('ik,kj->ij',vxcp0,QS))

          #print("VXCP:",vxcp.dtype)
          ## Version 2
          #vxc -= ks.phyb*vxcp.real/(1 - hyb)
          #exc -= ks.phyb*excp.real/(1 - hyb)
          vxc -= ks.phyb*vxcp/(1 - hyb)
          exc -= ks.phyb*excp/(1 - hyb)

        #if ks.do_nlc():
        if ks.nlc != '':
            if ni.libxc.is_nlc(ks.xc):
                xc = ks.xc
            else:
                assert ni.libxc.is_nlc(ks.nlc)
                xc = ks.nlc
            n, enlc, vnlc = ni.nr_nlc_vxc(mol, ks.nlcgrids, xc, dm,
                                          hermi=hermi, max_memory=max_memory)
            exc += enlc
            vxc += vnlc
            logger.debug(ks, 'nelec with nlc grids = %s', n)
        t0 = logger.timer(ks, 'vxc', *t0)

    # Add EEX from the projected density matrix
    if(abs(ks.phyb)>1e-10):
      vxxp0 = ks.get_k(mol,pdm,hermi)
      vxxp  = numpy.einsum('ik,kj->ij', SQ, numpy.einsum('ik,kj->ij', vxxp0, QS))
      exxp  = numpy.einsum('ij,ji', pdm, vxxp0).real * .5
      #print("VXXP:",vxxp.dtype)
      ## Version 2
      #vxc -= (ks.phyb) * vxxp.real
      #exc -= (ks.phyb) * exxp.real
      vxc -= (ks.phyb) * vxxp
      exc -= (ks.phyb) * exxp

    if not ni.libxc.is_hybrid_xc(ks.xc):
        vk = None
        if (ks._eri is None and ks.direct_scf and
            getattr(vhf_last, 'vj', None) is not None):
            ddm = numpy.asarray(dm) - numpy.asarray(dm_last)
            vj = ks.get_j(mol, ddm, hermi)
            vj += vhf_last.vj
        else:
            vj = ks.get_j(mol, dm, hermi)
        vxc += vj
    else:
        #omega, alpha, hyb = ni.rsh_and_hybrid_coeff(ks.xc, spin=mol.spin)
        if (ks._eri is None and ks.direct_scf and
            getattr(vhf_last, 'vk', None) is not None):
            ddm = numpy.asarray(dm) - numpy.asarray(dm_last)
            vj, vk = ks.get_jk(mol, ddm, hermi)
            vk *= hyb
            if omega != 0:
                vklr = ks.get_k(mol, ddm, hermi, omega=omega)
                vklr *= (alpha - hyb)
                vk += vklr
            vj += vhf_last.vj
            vk += vhf_last.vk
        else:
            vj, vk = ks.get_jk(mol, dm, hermi)
            vk *= hyb
            if omega != 0:
                vklr = ks.get_k(mol, dm, hermi, omega=omega)
                vklr *= (alpha - hyb)
                vk += vklr

        ## If Hybrid
        if(abs(ks.phyb)>1e-10):
            #vkp0 = ks.get_k(mol, pdm, hermi)
            #vkp0 *= hyb
            if abs(omega) > 1e-10:
                vklrp0 = ks.get_k(mol, pdm, hermi, omega=omega)
                vklrp0 *= (alpha - hyb)
                #vkp0 += vklrp0
                vkp0 = vklrp0
                vkp = numpy.einsum('ik,kj->ij', SQ, numpy.einsum('ik,kj->ij', vkp0, QS))
                vk -= ks.phyb * vkp

        #vxc += (vj - vk).real
        vxc += vj - vk

        if ground_state:
            exc -= numpy.einsum('ij,ji', dm, vk).real * .5

    if ground_state:
        ecoul = numpy.einsum('ij,ji', dm, vj).real * .5
    else:
        ecoul = None

    vxc = lib.tag_array(vxc, ecoul=ecoul, exc=exc, vj=vj, vk=vk)
    return vxc

energy_elec = rks.energy_elec


class GKS(rks.KohnShamPDFT, ghf.GHF):
    '''Generalized Kohn-Sham'''

    get_veff = get_veff
    energy_elec = energy_elec

    def __init__(self, mol, xc='LDA,VWN', phyb=0, paos=None, ext_basis = '3-21G', use_ext_basis=True):
        ghf.GHF.__init__(self, mol)
        rks.KohnShamPDFT.__init__(self, xc, phyb, paos, ext_basis, use_ext_basis)
        self._numint = NumInt2C()
        self.__dict__['with_spin'] = True

    def dump_flags(self, verbose=None):
        ghf.GHF.dump_flags(self, verbose)
        rks.KohnShamPDFT.dump_flags(self, verbose)
        logger.info(self, 'collinear = %s', self._numint.collinear)
        if self._numint.collinear[0] == 'm':
            logger.info(self, 'mcfun spin_samples = %s', self._numint.spin_samples)
            logger.info(self, 'mcfun collinear_thrd = %s', self._numint.collinear_thrd)
            logger.info(self, 'mcfun collinear_samples = %s', self._numint.collinear_samples)
        return self

    def _build_proj(self):
        '''
        Function to build projector.
        '''
        if self.use_ext_basis:
            SQQS = project.build_spin_proj_in_ext_basis(self, ext_basis=self.ext_basis)
        else:
            #SQQS = project.build_proj_in_basis(self)
            raise NotImplementedError()
        self._SQQS = SQQS
        return

    @property
    def collinear(self):
        return self._numint.collinear
    @collinear.setter
    def collinear(self, val):
        self._numint.collinear = val

    @property
    def spin_samples(self):
        return self._numint.spin_samples
    @spin_samples.setter
    def spin_samples(self, val):
        self._numint.spin_samples = val

    def nuc_grad_method(self):
        raise NotImplementedError

    def to_hf(self):
        '''Convert to GHF object.'''
        return self._transfer_attrs_(self.mol.GHF())

    to_gpu = lib.to_gpu


if __name__ == '__main__':
    from pyscf import gto
    mol = gto.Mole()
    mol.verbose = 3
    mol.atom = 'H 0 0 0; H 0 0 1; O .5 .6 .2'
    mol.basis = 'ccpvdz'
    mol.build()

    mf = GKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()

    dm = mf.init_guess_by_1e(mol)
    dm = dm + 0j
    nao = mol.nao_nr()
    numpy.random.seed(12)
    dm[:nao,nao:] = numpy.random.random((nao,nao)) * .1j
    dm[nao:,:nao] = dm[:nao,nao:].T.conj()
    mf.kernel(dm)
    mf.canonicalize(mf.mo_coeff, mf.mo_occ)
    mf.analyze()
    print(mf.spin_square())
    print(mf.e_tot - -76.2760115704274)
