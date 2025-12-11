#!/usr/bin/env python
# Copyright 2014-2020 The PySCF Developers. All Rights Reserved.
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
Non-Relativistic Restricted Projected Kohn-Sham
'''

#import time
import numpy
import warnings
from functools import reduce
from scipy import linalg
from pyscf import lib
from pyscf.lib import logger
from pyscf import gto, scf
from pyscf.scf import hf
from pyscf.scf import _vhf
from pyscf.scf import jk
from pyscf.dft import gen_grid
from pyscf.dft import numint
from pyscf.dft import libxc as lxc
#from pdft import numint
from graci.pdft import project, pscf
from pyscf import __config__
## NOTE: currently, paos not initialized (saved as attribute).
## paos (i.e., caos) called within function to build operator; implicit.
## if you set paos manually, it doesn't do anything.

def get_veff(ks, mol=None, dm=None, dm_last=0, vhf_last=0, hermi=1):
    '''Coulomb + XC functional

    .. note::
        This function will modify the input ks object.

    Args:
        ks : an instance of :class:`RKS`
            XC functional are controlled by ks.xc attribute.  Attribute
            ks.grids might be initialized.
            Projected hybrid is controlled by fraction ks.phyb  (analogous to hyb) 
            which controls use of projection operators ks.SQ  and ks.QS 
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
    if dm is None: dm = ks.make_rdm1()
    ext_basis = ks.ext_basis

    #t0 = (logger.process_clock(), logger.perf_counter())

    ground_state = (isinstance(dm, numpy.ndarray) and dm.ndim == 2)

    # Build the AO projection operators and projected 1pdm
    # 1. Evaluate an XC contribution from DM
    # 2. Subtract off the corresponding contribution from pdm
    # 3. Add EEX from PDM.
    assert (ks.SQQS is not None) ## maybe add this to ks.build() ??
    SQ = ks.SQQS[0]; QS = ks.SQQS[1]
    D = len(ks.phyb) ## number of projector operators.

    #pdm = numpy.einsum('ik,kj->ij',ks.QS,numpy.einsum('ik,kj->ij',dm,ks.SQ))
    pdm = []
    for i in range(0, D):
      pdm_i = reduce(numpy.dot, (QS[i],dm,SQ[i]))
      pdm.append(pdm_i)

    if ks.grids.coords is None:
        ks.grids.build(with_non0tab=True)
        if ks.small_rho_cutoff > 1e-20 and ground_state:
            # Filter grids the first time setup grids
            ks.grids = prune_small_rho_grids_(ks, mol, dm, ks.grids)
        #t0 = logger.timer(ks, 'setting up grids', *t0)
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
        max_memory = ks.max_memory - lib.current_memory()[0]
        ## Regular XC contributions.
        ##        i.e., vxc excludes fraction of EEX.
        n, exc, vxc = ni.nr_rks(mol, ks.grids, ks.xc, dm, max_memory=max_memory)
        if (D!=1) or (abs(ks.phyb[0])>1e-10): #if more than one projection defined OR phyb>0.
          pxc = ks.xc

          #### Identify XC functional components
          if ks.xcstr is not None:
            px, pc = ks.xcstr
          else:
            px = pxc

          ## Projected VXC in Projected RDM1.
          for i in range(0, D):
            np, excp, vxcp0 = ni.nr_rks(mol, ks.grids, px, pdm[i], max_memory=max_memory)
            vxcp = numpy.einsum('ik,kj->ij',SQ[i],numpy.einsum('ik,kj->ij',vxcp0,QS[i]))

            ## Version 1
            #vxc -= ks.phyb*vxcp
            #exc -= ks.phyb*excp
            ## Version 2
            vxc -= ks.phyb[i]*vxcp/(1 - hyb)
            exc -= ks.phyb[i]*excp/(1 - hyb)

        if ks.nlc != '':
          assert('VV10' in ks.nlc.upper())
          _, enlc, vnlc = ni.nr_rks(mol, ks.nlcgrids, ks.xc+'__'+ks.nlc, dm,
                                      max_memory=max_memory)
          exc += enlc
          vxc += vnlc
        logger.debug(ks, 'nelec by numeric integration = %s', n)
        #t0 = logger.timer(ks, 'vxc', *t0)

    # Add EEX from the projected density matrix
    if (D!=1) or (abs(ks.phyb[0])>1e-10):
      for i in range(0, D):
        vxxp0 = ks.get_k(mol,pdm[i],hermi)
        vxxp  = numpy.einsum('ik,kj->ij', SQ[i], numpy.einsum('ik,kj->ij', vxxp0, QS[i]))
        exxp  = numpy.einsum('ij,ji', pdm[i], vxxp0).real * .5

        ## Version 1
        #vxc -= (ks.phyb) * (1 - hyb) * vxxp * .5
        #exc -= (ks.phyb) * (1 - hyb) * exxp * .5

        ## Version 2
        vxc -= (ks.phyb[i]) * vxxp * .5
        exc -= (ks.phyb[i]) * exxp * .5

    if abs(hyb) < 1e-10 and abs(alpha) < 1e-10:
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
      if (ks._eri is None and ks.direct_scf and
        getattr(vhf_last, 'vk', None) is not None):
        ddm = numpy.asarray(dm) - numpy.asarray(dm_last)
        vj, vk = ks.get_jk(mol, ddm, hermi)
        vk *= hyb
        if abs(omega) > 1e-10:  # For range separated Coulomb operator
          vklr = ks.get_k(mol, ddm, hermi, omega=omega)
          vklr *= (alpha - hyb)
          vk += vklr
        vj += vhf_last.vj
        vk += vhf_last.vk
      else:
        vj, vk = ks.get_jk(mol, dm, hermi)
        vk *= hyb
        #print("retrieved (vj - hyb * vk)")
        if abs(omega) > 1e-10:
          vklr = ks.get_k(mol, dm, hermi, omega=omega)
          vklr *= (alpha - hyb)
          vk += vklr

      ## If Hybrid (Range-Sep.)
      if (D!=0) or (abs(ks.phyb[0])>1e-10):
        #vkp0 = ks.get_k(mol, pdm, hermi)
        #vkp0 *= hyb
        if abs(omega) > 1e-10:
          for i in range(0, D):
            vklrp0 = ks.get_k(mol, pdm[i], hermi, omega=omega)
            vklrp0 *= (alpha - hyb)
            #vkp0 += vklrp0
            vkp0 = vklrp0
            vkp = numpy.einsum('ik,kj->ij', SQ[i], numpy.einsum('ik,kj->ij', vkp0, QS[i]))
            vk -= ks.phyb[i] * vkp

        vxc += vj - (vk * .5)

        if ground_state:
            exc -= numpy.einsum('ij,ji', dm, vk).real * .5 * .5

    if ground_state:
        ecoul = numpy.einsum('ij,ji', dm, vj).real * .5
    else:
        ecoul = None

    vxc = lib.tag_array(vxc, ecoul=ecoul, exc=exc, vj=vj, vk=vk)
    return vxc

def get_vsap(ks, mol=None):
    '''Superposition of atomic potentials

    S. Lehtola, Assessment of initial guesses for self-consistent
    field calculations. Superposition of Atomic Potentials: simple yet
    efficient, J. Chem. Theory Comput. 15, 1593 (2019). DOI:
    10.1021/acs.jctc.8b01089. arXiv:1810.11659.

    This function evaluates the effective charge of a neutral atom,
    given by exchange-only LDA on top of spherically symmetric
    unrestricted Hartree-Fock calculations as described in

    S. Lehtola, L. Visscher, E. Engel, Efficient implementation of the
    superposition of atomic potentials initial guess for electronic
    structure calculations in Gaussian basis sets, J. Chem. Phys., in
    press (2020).

    The potentials have been calculated for the ground-states of
    spherically symmetric atoms at the non-relativistic level of theory
    as described in

    S. Lehtola, "Fully numerical calculations on atoms with fractional
    occupations and range-separated exchange functionals", Phys. Rev. A
    101, 012516 (2020). DOI: 10.1103/PhysRevA.101.012516

    using accurate finite-element calculations as described in

    S. Lehtola, "Fully numerical Hartree-Fock and density functional
    calculations. I. Atoms", Int. J. Quantum Chem. e25945 (2019).
    DOI: 10.1002/qua.25945

    .. note::
        This function will modify the input ks object.

    Args:
        ks : an instance of :class:`RKS`
            XC functional are controlled by ks.xc attribute.  Attribute
            ks.grids might be initialized.

    Returns:
        matrix Vsap = Vnuc + J + Vxc.
    '''
    if mol is None: mol = ks.mol
    #t0 = (time.clock(), time.time())

    if ks.grids.coords is None:
        ks.grids.build(with_non0tab=True)
        # t0 = logger.timer(ks, 'setting up grids', *t0)

    ni = ks._numint
    max_memory = ks.max_memory - lib.current_memory()[0]
    vsap = ni.nr_sap(mol, ks.grids, max_memory=max_memory)
    return vsap

# The vhfopt of standard Coulomb operator can be used here as an approximate
# opt since long-range part Coulomb is always smaller than standard Coulomb.
# It's safe to prescreen LR integrals with the integral estimation from
# standard Coulomb.
def _get_k_lr(mol, dm, omega=0, hermi=0, vhfopt=None):
    import sys
    sys.stderr.write('This function is deprecated. '
                     'It is replaced by mol.get_k(mol, dm, omege=omega)')
    dm = numpy.asarray(dm)
# Note, ks object caches the ERIs for small systems. The cached eris are
# computed with regular Coulomb operator. ks.get_jk or ks.get_k do not evalute
# the K matrix with the range separated Coulomb operator.  Here jk.get_jk
# function computes the K matrix with the modified Coulomb operator.
    nao = dm.shape[-1]
    dms = dm.reshape(-1,nao,nao)
    with mol.with_range_coulomb(omega):
        # Compute the long range part of ERIs temporarily with omega. Restore
        # the original omega when the block ends
        if vhfopt is None:
            contents = lambda: None # just a place_holder
        else:
            contents = vhfopt._this.contents
        with lib.temporary_env(contents,
                               fprescreen=_vhf._fpointer('CVHFnrs8_vk_prescreen')):
            intor = mol._add_suffix('int2e')
            vklr = jk.get_jk(mol, dms, ['ijkl,jk->il']*len(dms), intor=intor,
                             vhfopt=vhfopt)
    return numpy.asarray(vklr).reshape(dm.shape)


def energy_elec(ks, dm=None, h1e=None, vhf=None):
    r'''Electronic part of RKS energy.

    Note this function has side effects which cause mf.scf_summary updated.

    Args:
        ks : an instance of DFT class

        dm : 2D ndarray
            one-partical density matrix
        h1e : 2D ndarray
            Core hamiltonian

    Returns:
        RKS electronic energy and the 2-electron contribution
    '''
    if dm is None: dm = ks.make_rdm1()
    if h1e is None: h1e = ks.get_hcore()
    if vhf is None or getattr(vhf, 'ecoul', None) is None:
        vhf = ks.get_veff(ks.mol, dm)
    e1 = numpy.einsum('ij,ji->', h1e, dm)
    e2 = vhf.ecoul + vhf.exc
    ks.scf_summary['e1'] = e1.real
    ks.scf_summary['coul'] = vhf.ecoul.real
    ks.scf_summary['exc'] = vhf.exc.real
    logger.debug(ks, 'E1 = %s  Ecoul = %s  Exc = %s', e1, vhf.ecoul, vhf.exc)
    return (e1+e2).real, e2


NELEC_ERROR_TOL = getattr(__config__, 'dft_rks_prune_error_tol', 0.02)
def prune_small_rho_grids_(ks, mol, dm, grids):
    rho = ks._numint.get_rho(mol, dm, grids, ks.max_memory)
    n = numpy.dot(rho, grids.weights)
    if abs(n-mol.nelectron) < NELEC_ERROR_TOL*n:
        rho *= grids.weights
        idx = abs(rho) > ks.small_rho_cutoff / grids.weights.size
        logger.debug(ks, 'Drop grids %d',
                     grids.weights.size - numpy.count_nonzero(idx))
        grids.coords  = numpy.asarray(grids.coords [idx], order='C')
        grids.weights = numpy.asarray(grids.weights[idx], order='C')
        grids.non0tab = grids.make_mask(mol, grids.coords)
    return grids

def define_xc_(ks, description, xctype='LDA', hyb=0, rsh=(0,0,0)):
    libxc = ks._numint.libxc
    ks._numint = libxc.define_xc_(ks._numint, description, xctype, hyb, rsh)
    return ks


def _pdft_common_init_(mf, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis='3-21G', use_ext_basis = True):
    raise DeprecationWarning

class KohnShamPDFT(object):
    '''
    Attributes for Kohn-Sham AO-projected DFT:
        xc : str
            'X_name,C_name' for the XC functional.  Default is 'lda,vwn'
        nlc : str
            'NLC_name' for the NLC functional.  Default is '' (i.e., None)
        omega : float
            Omega of the range-separated Coulomb operator e^{-omega r_{12}^2} / r_{12}
        phyb  : float 
            Fraction of projected hybrid exchange to replace existing XC functional with. 
        paos  : str
            List, starting from zero, of the AOs in the projection. 
        grids : Grids object
            grids.level (0 - 9)  big number for large mesh grids. Default is 3

            radii_adjust
                | radi.treutler_atomic_radii_adjust (default)
                | radi.becke_atomic_radii_adjust
                | None : to switch off atomic radii adjustment

            grids.atomic_radii
                | radi.BRAGG_RADII  (default)
                | radi.COVALENT_RADII
                | None : to switch off atomic radii adjustment

            grids.radi_method  scheme for radial grids
                | radi.treutler  (default)
                | radi.delley
                | radi.mura_knowles
                | radi.gauss_chebyshev

            grids.becke_scheme  weight partition function
                | gen_grid.original_becke  (default)
                | gen_grid.stratmann

            grids.prune  scheme to reduce number of grids
                | gen_grid.nwchem_prune  (default)
                | gen_grid.sg1_prune
                | gen_grid.treutler_prune
                | None : to switch off grids pruning

            grids.symmetry  True/False  to symmetrize mesh grids (TODO)

            grids.atom_grid  Set (radial, angular) grids for particular atoms.
            Eg, grids.atom_grid = {'H': (20,110)} will generate 20 radial
            grids and 110 angular grids for H atom.

        small_rho_cutoff : float
            Drop grids if their contribution to total electrons smaller than
            this cutoff value.  Default is 1e-7.

    Examples:

    >>> mol = gto.M(atom='O 0 0 0; H 0 0 1; H 0 1 0', basis='ccpvdz', verbose=0)
    >>> mf = dft.RKS(mol)
    >>> mf.xc = 'b3lyp'
    >>> mf.kernel()
    -76.415443079840458
    '''
    _keys = {'xc', 'xcstr', 'nlc', 'grids', 'disp', 'nlcgrids', 'small_rho_cutoff', 'phyb', 'paos', 'ext_basis', 'use_ext_basis', 'SQQS'}

    def __init__(self, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis='3-21G', use_ext_basis = True):
        self.xc = xc
        self.xc_handler()
        ## Projector
        self.paos = paos
        self.phyb = phyb
        self.use_ext_basis = use_ext_basis
        self.ext_basis = ext_basis
        self._set_projector_params()

        ## Projector built as needed (dbl check)
        #self._build_proj()
        self.SQQS = None

        ## Other
        self.nlc = ''
        self.grids = gen_grid.Grids(self.mol)
        self.grids.level = getattr(__config__, 'dft_rks_RKS_grids_level',
                             self.grids.level)
        self.nlcgrids = gen_grid.Grids(self.mol)
        self.nlcgrids.level = getattr(__config__, 'dft_rks_RKS_nlcgrids_level',
                                self.nlcgrids.level)
        # Use rho to filter grids
        self.small_rho_cutoff = getattr(__config__, 'dft_rks_RKS_small_rho_cutoff', 1e-7)
        ##################################################
        # don't modify the following attributes, they are not input options
        self._numint = numint.NumInt()
        self._keys = self._keys.union(['xc', 'nlc', 'omega', 'grids', 'nlcgrids',
                               'small_rho_cutoff'])

    @property
    def omega(self):
        return self._numint.omega
    @omega.setter
    def omega(self, v):
        self._numint.omega = float(v)

    def dump_flags(self, verbose=None):
        logger.info(self, 'XC functionals = %s', self.xc)
        if self.nlc!='':
            logger.info(self, 'NLC functional = %s', self.nlc)
        # Projection.
        if self.paos is not None:
            logger.info(self, 'Projected AOs = %s', self.paos)
            logger.info(self, 'Projection fraction AOs = %7.4f', self.phyb)
        # END
        logger.info(self, 'small_rho_cutoff = %g', self.small_rho_cutoff)
        self.grids.dump_flags(verbose)
        if self.nlc!='':
            logger.info(self, '** Following is NLC Grids **')
            self.nlcgrids.dump_flags(verbose)
        return self

    define_xc_ = define_xc_

    def _set_projector_params(self):
        '''
        Set parameters related to projector (insert sanity checks here)
        and override redundant/improper inputs.

        This is called by check_sanity, so that if parameters are reset outside of self.__init
        the relevant variables are corrected. [Trying to fix a big I encountered before, related
        to compatibility with GRaCI..]

        Currently:
        a) paos is actually redunandant, should remove soon (AOs assigned by a function)
        b) converts phyb to a list if not already a list (project onto multiple edges)
        c) reset ext_basis to None if not using (otherwise, default basis is 3-21G)
        '''
        # 1. self.paos = paos
        # 2. self.phyb = phyb
        #   print("TYPE(PHYB)=",type(phyb))
        #   print("phyb=",phyb)
        phyb = self.phyb
        if type(phyb) == float:
            self.phyb = [phyb]
        else:
            assert (type(phyb) == list)
            self.phyb = phyb
        # 3. self.use_ext_basis = ...
        if self.use_ext_basis and (self.ext_basis is None):
            warnings.warn(
                f'Attribute use_external_basis toggled, but no basis set was specified.'
                 'Setting to default: 3-21G.')
            self.ext_basis = '3-21G'
        elif (not self.use_ext_basis) and (self.ext_basis is not None):
            warnings.warn(
                f'Removing external basis set: {self.ext_basis}.')
            self.ext_basis = None
        return

    def check_sanity(self):
        self._set_projector_params()
        ## Default check_sanity for dft.KohnShamDFT (below)
        out = super().check_sanity()
        #if self.do_nlc() and self.do_disp() and self._numint.libxc.is_nlc(self.xc):
        #    import warnings
        #    warnings.warn(
        #        f'nlc-type xc {self.xc} and disp {self.disp} may lead to'
        #        'double counting in NLC.')
        return out

    #assign_core_aos = project._assign_core_aos_by_label
    def get_core_aos(self):
        '''
        Function to assign core AOs.
        '''
        caos = project.assign_core_aos(self, mol = self.mol)
        return caos

    def build_proj(self, **kwargs):
        '''
        Function to build projector, SQQS = SQ,QS.

        If self.use_ext_basis == True, then SQ,QS are built in external AO basis.
        Otherwise, SQ,QS are built in the internal MO basis.

        However, if projecting onto multiple edges, this class currently overrides
        to external basis.

        '''
        ##print(type(self.phyb))
        ##print("Running pdft.rks.RKS.build_proj...")
        D = len(self.phyb)
        if (D > 1):
            warnings.warn("Building projector by edge (default to external AO basis).")
            self._build_proj_by_edge()
        else:
            self._build_proj(**kwargs)

        ## Number of projectors should match number of edges.
        #  (this really belongs inside of a 'sanity_check', but anyway..)
        assert( len(self.SQQS[0]) == len(self.SQQS[1]) )
        assert( len(self.SQQS[0]) == D )
        return

    def _build_proj(self, **kwargs):
        '''
        Private function to build projector.
        '''
        warnings.warn('''The projector builder is currently implemented for atoms. It cannot currently discriminate by element-type (O1s, N1s, etc.); however, it does discriminate by edge type (K-edge, L-edge, etc.).''')
        if self.use_ext_basis:
            sqqs = project.build_proj_in_ext_basis(self, ext_basis=self.ext_basis)
        else:
            #sqqs = project.build_proj_in_basis(self)
            sqqs = project.build_mo_proj(self, **kwargs) ## insert mo_coeff (cmos) if defined.
        SQ = [];QS = []
        SQ.append(sqqs[0])
        QS.append(sqqs[1])
        SQQS = [SQ, QS]
        self.SQQS = SQQS
        return

    def _build_proj_by_edge(self):
        '''
        Function to build projector by edge (default is external basis).
        '''
        warnings.warn('''The projector builder is currently implemented for atoms. It cannot currently discriminate by element-type (O1s, N1s, etc.);
                      however, it does discriminate by edge type (K-edge, L-edge, etc.).''')
        if not self.use_ext_basis:
            self.use_ext_basis = True
            warnings.warn("Internal (MO) basis not supported for this method. Overriding to external basis.")
        ## Iterate over atoms.
        M = self.mol.natm
        #elements = set()
        charges = set()
        for I in range(0, M):
            #elements.add( self.mol.atom_symbol(I) ) ##i.e., 'H'
            charges.add( int(self.mol.atom_charge(I)) )   ##i.e.,  Z = 1
        ## Iterate over Z (assess which edges to include..)
        edges = set()
        for Z in charges:
            if (Z > 2) and (Z < 11): # K-shell is not 'core' for H, He
                edges.add('K')
            elif (Z > 10) and (Z < 19):
                edges.add('K')
                edges.add('L')
            elif (Z > 18) and (Z < 37):
                edges.add('K')
                edges.add('L')
                edges.add('M')
            elif (Z > 36):
                edges.add('K')
                edges.add('L')
                edges.add('M')
                warnings.warn(f"Program does not currently support the assignment of core orbitals beyond the M-edge for 5th row elements: element {gto.elements.ATOMIC_NAMES[Z]} (Z={Z}).")
        ## BUILD PROJ FOR EACH EDGE
        SQ = []
        QS = []
        for X in edges:
            core_aos = project.assign_core_aos_by_edge(self, mol = self.mol, edge = X)
            sqqs = project.build_proj_in_ext_basis(self, ext_basis=self.ext_basis, caos = core_aos)
            SQ.append(sqqs[0])
            QS.append(sqqs[1])
        SQQS = [SQ, QS]
        self.SQQS = SQQS
        return

    def xc_handler(self):
        '''
        Function to handle XC functional.
        '''
        if (self.xc in lxc.XC_KEYS):
            ex, exch, corr = project._xc_handler(self.xc, include=False)
            xstr = project.format_xc(exch)
            cstr = project.format_xc(corr)
            self.xcstr = [xstr,cstr]
        elif ',' in self.xc:
            xstr,cstr = (self.xc).split(',')
            self.xcstr = [xstr,cstr]
        else:
            warnings.warn("Unknown XC type passed to RKS object.")
            self.xcstr = None
        return

    def scf(self, dm0=None, **kwargs):
        '''SCF main driver

        Kwargs:
            dm0 : ndarray
                If given, it will be used as the initial guess density matrix

        Examples:

        >>> import numpy
        >>> from pyscf import gto, scf
        >>> mol = gto.M(atom='H 0 0 0; F 0 0 1.1')
        >>> mf = scf.hf.SCF(mol)
        >>> dm_guess = numpy.eye(mol.nao_nr())
        >>> mf.kernel(dm_guess)
        converged SCF energy = -98.5521904482821
        -98.552190448282104
        '''
        cput0 = (logger.process_clock(), logger.perf_counter())

        self.dump_flags()
        self.build(self.mol)

        if dm0 is None and self.mo_coeff is not None and self.mo_occ is not None:
            # Initial guess from existing wavefunction
            dm0 = self.make_rdm1()

        if self.max_cycle > 0 or self.mo_coeff is None:
            self.converged, self.e_tot, \
                    self.mo_energy, self.mo_coeff, self.mo_occ = \
                    pscf.kernel(self, self.conv_tol, self.conv_tol_grad,
                           dm0=dm0, callback=self.callback,
                           conv_check=self.conv_check, **kwargs)
        else:
            # Avoid to update SCF orbitals in the non-SCF initialization
            # (issue #495).  But run regular SCF for initial guess if SCF was
            # not initialized.
            self.e_tot = pscf.kernel(self, self.conv_tol, self.conv_tol_grad,
                                dm0=dm0, callback=self.callback,
                                conv_check=self.conv_check, **kwargs)[1]

        logger.timer(self, 'SCF', *cput0)
        self._finalize()
        return self.e_tot
    kernel = lib.alias(scf, alias_name='kernel')

    def to_rhf(self):
        '''Convert the input mean-field object to a RHF/ROHF object.

        Note this conversion only changes the class of the mean-field object.
        The total energy and wave-function are the same as them in the input
        mean-field object.
        '''
        mf = scf.RHF(self.mol)
        mf.__dict__.update(self.to_rks().__dict__)
        mf.converged = False
        return mf

    def to_uhf(self):
        '''Convert the input mean-field object to a UHF object.

        Note this conversion only changes the class of the mean-field object.
        The total energy and wave-function are the same as them in the input
        mean-field object.
        '''
        mf = scf.UHF(self.mol)
        mf.__dict__.update(self.to_uks().__dict__)
        mf.converged = False
        return mf

    def to_ghf(self):
        '''Convert the input mean-field object to a GHF object.

        Note this conversion only changes the class of the mean-field object.
        The total energy and wave-function are the same as them in the input
        mean-field object.
        '''
        mf = scf.GHF(self.mol)
        mf.__dict__.update(self.to_gks().__dict__)
        mf.converged = False
        return mf

    def to_rks(self, xc=None):
        '''Convert the input mean-field object to a RKS/ROKS object.

        Note this conversion only changes the class of the mean-field object.
        The total energy and wave-function are the same as them in the input
        mean-field object.
        '''
        mf = scf.addons.convert_to_rhf(self)
        if xc is not None:
            mf.xc = xc
        # Projection.
        if paos is not None:
            mf.paos = paos
            mf.phyb = phyb
        # END
        if xc != self.xc or not isinstance(self, RKS):
            mf.converged = False
        return mf

    def to_uks(self, xc=None):
        '''Convert the input mean-field object to a UKS object.

        Note this conversion only changes the class of the mean-field object.
        The total energy and wave-function are the same as them in the input
        mean-field object.
        '''
        mf = scf.addons.convert_to_uhf(self)
        if xc is not None:
            mf.xc = xc
        if xc != self.xc:
            mf.converged = False
        return mf

    def to_gks(self, xc=None):
        '''Convert the input mean-field object to a GKS object.

        Note this conversion only changes the class of the mean-field object.
        The total energy and wave-function are the same as them in the input
        mean-field object.
        '''
        mf = scf.addons.convert_to_ghf(self)
        if xc is not None:
            mf.xc = xc
        if xc != self.xc:
            mf.converged = False
        return mf

    def reset(self, mol=None):
        hf.SCF.reset(self, mol)
        self.grids.reset(mol)
        self.nlcgrids.reset(mol)
        return self


def init_guess_by_vsap(mf, mol=None):
    '''Function to form the superposition of atomic potentials (SAP) guess
    density.
    '''
    if mol is None: mol = mf.mol

    vsap = mf.get_vsap()
    t = mol.intor_symmetric('int1e_kin')
    s = mf.get_ovlp(mol)
    hsap = t + vsap

    # Form guess orbitals
    mo_energy, mo_coeff = mf.eig(hsap, s)
    logger.debug(mf, 'VSAP mo energies\n{}'.format(mo_energy))

    # and guess density
    mo_occ = mf.get_occ(mo_energy, mo_coeff)
    return mf.make_rdm1(mo_coeff, mo_occ)

# Update the KohnShamDFT label in scf.hf module
hf.KohnShamDFT = KohnShamPDFT

class RKS(KohnShamPDFT, hf.RHF):
    __doc__ = '''Restricted Kohn-Sham\n''' + hf.SCF.__doc__ + KohnShamPDFT.__doc__

    def __init__(self, mol, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis = '3-21G', use_ext_basis=True):
        hf.RHF.__init__(self, mol)
        KohnShamPDFT.__init__(self, xc, phyb, paos, ext_basis, use_ext_basis)

    def dump_flags(self, verbose=None):
        hf.RHF.dump_flags(self, verbose)
        return KohnShamPDFT.dump_flags(self, verbose)

    get_veff = get_veff
    get_vsap = get_vsap
    energy_elec = energy_elec

    init_guess_by_vsap = init_guess_by_vsap

    def nuc_grad_method(self):
        from pyscf.grad import rks as rks_grad
        return rks_grad.Gradients(self)


if __name__ == '__main__':
    from pyscf import gto
    from pyscf.dft import xcfun
    mol = gto.Mole()
    #mol.verbose = 7
    #mol.output = '/dev/null' #'out_rks'

    mol.atom.extend([['He', (0.,0.,0.)], ])
    mol.basis = { 'He': 'cc-pvdz'}
    #mol.grids = { 'He': (10, 14),}
    mol.build()

    m = RKS(mol)
    m.xc = 'b88,lyp'
    print(m.scf())  # -2.8978518405

    m = RKS(mol)
    m._numint.libxc = xcfun
    m.xc = 'b88,lyp'
    print(m.scf())  # -2.8978518405
