"""
module to compute functions of the 1-TDM
"""

import os as os
import numpy as np
import functools
import operator
import importlib
from sympy import LeviCivita
import graci.interaction.interaction as interaction
import graci.utils.timing as timing
import graci.core.libs as libs
import graci.bitcitools.bitsi_init as bitsi_init
import graci.bitcitools.mrci_1tdm as mrci_1tdm
import graci.io.output as output
import graci.utils.constants as constants
import sys as sys

class Transition(interaction.Interaction):
    """Transition class for determing transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # parent attributes
        super().__init__()
        
        # name for this object
        self.label          = 'Transition'
        # user defined quanties
        self.print_orbitals = False

        # for transition it's convenient to call these variables
        # init/final instead of ket/bra. Here we just say they
        # point to the same reference as the variables in the parent
        # class
        # list of initial states
        #self.init_states      = self.ket_states = None
        # label of section to get initial states
        #self.init_label       = self.ket_label  = None
        # list of final states
        #self.final_states     = self.bra_states = None
        # label of section to get final states 
        #self.final_label      = self.bra_label  = None
        
        # often will want to include all states in bra object
        # if so convenient to have simple boolean to select that
        self.all_final_states = False

        # internal class variables -- should not be accessed
        # directly

        # list of transitions corresponding to the tdms
        self.trans_list = []
        # the tdms
        self.tdms       = None
        # various orbital quantities (NDOs, NTOs, wts, etc.)
        self.nos        = {} 
        # electric/magnetic multipole tensors -- all guages
        self.multipole  = {}
        # oscillator strengths -- all gauges
        self.oscstr     = {}

    #
    @timing.timed
    def run(self):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used."""

        self.bra_label = self.final_label
        self.ket_label = self.init_label

        mol_bra = self.bra_obj.scf.mol
        mol_ket = self.ket_obj.scf.mol

        scf_bra = self.bra_obj.scf
        scf_ket = self.ket_obj.scf

        # if bra and ket are same object, just need lower triangle
        # minus the diagonal elements
        if self.same_braket():
            list_type = 'nodiag'
        else:
            # sanity check that orbitals and geometry are the same
            if np.any(scf_bra.orbs != scf_ket.orbs):
                sys.exit('transition moments require same bra/ket orbs')

            if mol_bra.pymol().atom != mol_ket.pymol().atom:
                sys.exit('transition moments require same geometry'+
                         ' and basis set')

            # bra and ket not identical objects, require all bra-ket
            # pairs of states
            list_type = 'full'

        # initialize the bitsi library for the calculation of 1-TDMs
        bitsi_init.init(self.bra_obj, self.ket_obj, 'tdm')
        
        # if all_final_states is true, over-ride current contents of 
        # self.final_states
        if self.all_final_states:
            self.final_states = range(self.bra_obj.n_states())

        self.add_states('ket_states', self.ket_obj, self.init_states)
        self.add_states('bra_states', self.bra_obj, self.final_states)

        # this is main transition_list: stored by adiabatic label
        self.trans_list = self.build_pair_list('bra_states',
                                               'ket_states', 
                                                pairs=list_type)
        # the transitions, ordered by interacting irreps, is how
        # we call bitsi
        trans_list_sym  = self.build_pair_list('bra_states', 
                                               'ket_states', 
                                                pairs=list_type, 
                                                sym_blk=True)

        # tdms is a vector of nmo x nmo transition densities
        self.tdms = self.build_tdms(trans_list_sym)

        # build the multipole moments  -- easier to just do this once
        # for all transitions
        self.multipole = self.build_multipole(mol_ket.pymol(), scf_ket)

        # compute oscillator strengths
        self.oscstr = self.build_osc_str()

        # build natural differnence and natural transition orbitals 
        # for each pair
        self.nos = self.build_natural_orbs()

        # print orbitals if requested
        if self.print_orbitals:
            self.export_orbitals()

        # print the summary output
        self.print_log()

        # finalize the bitsi library
        bitsi_init.finalize()
 
        return

    #
    def tdm(self, b_st, k_st):
        """return the tdm for the bra and ket state indicated"""

        indx = self.trans_list.index(b_st, k_st)
        if indx is not None:
            return self.tdms[:, :, indx]
        else:
            return None
 
    #
    def transition_multipole(self, b_st, k_st, name, gauge='velocity'):
        """return the transition mulitpole tensor for pairs of
           states ordered by energy"""
       
        return self.extract_multipole(b_st, k_st, name, gauge) 

    #
    def osc_strength(self, b_st, k_st, order=0, gauge='velocity'):
        """return the transition dipole for pairs of 
           states ordered by energy"""

        return self.extract_osc(b_st, k_st, order, gauge, False)
    
    # 
    def osc_strength_iso(self, b_st, k_st, order=0, gauge='velocity'):
        """return the isotropically averaged oscillator strength in
           the specified gauge and through the specified order"""

        return self.extract_osc(b_st, k_st, order, gauge, True)

#--------------------------------------------------------------------------
# "Private" class methods
#
    # 
    @timing.timed
    def build_tdms(self, trans_list, trans_list_sym):
        """grab the TDMs from bitsi and then reshape the list of
           TDMs into a more usable format"""

        # grab the tdms
        tdm_list = mrci_1tdm.tdm(self.bra_obj, self.ket_obj, 
                                  trans_list_sym)

        # make the tdm list
        nmo      = self.bra_obj.scf.nmo
        npair    = len(trans_list)
        tdms     = np.zeros((nmo, nmo, npairs), dtype=float)

        for indx in range(len(trans_list)):
            bk_st          = trans_list[indx]
            [birr, bst]    = self.bra_obj.state_sym(bk_st[0])
            [kirr, kst]    = self.ket_obj.state_sym(bk_st[1])
            sym_indx       = trans_list_sym[birr][kirr].index([bst,kst])
            tdms[:,:,indx] = tdm_list[birr][kirr][:, :, sym_indx]

        return tdms
    
    #
    def extract_multipole(self, b_st, k_st, name, gauge='velocity'):
        """extract the transition dipole vector in the appropriate
           gauge"""

        # append a '_l' or '_v' to indicate length or velocity gauge
        keyname = name+'_'+str(gauge)[0]
        indx = self.trans_list.index(b_st, k_st)

        if keyname in self.multipole.keys() and indx is not None:
            return self.multipole[keyname][..., indx]
        else:
            return None 

    #
    def extract_osc(self, b_st, k_st, gauge, order, iso):
        """extract the transition dipole vector in the appropriate
           gauge, complete through order 'order'. If iso is True,
           return isotropically averaged value"""
    
        istr = ''
        if iso:
            istr = 'iso'

        keyname = 'f'+str(order)+istr+'_'+str(gauge)[0]
        indx = self.trans_list.index(b_st, k_st)

        if keyname in self.oscstr.keys() and indx is not None:
            return self.oscstr[keyname][..., indx]
        else:
            return None

    #
    @timing.timed
    def build_multipole(self, mol, scf):
        """builds the multipole moment transitions tensors"""

        # build the electronic tensors: dipole, quadrupole, octupole
        dipole   = self.build_dip(mol, scf)
        quadpole = self.build_quad(mol, scf)
        octpole  = self.build_oct(mol, scf)

        # build the magnetic moment tensors: dipole, quadrupole
        mdipole   = self.build_mdip(mol, scf)
        mquadpole = self.build_mquad(mol, scf)

        return {**dipole, **quadpole, **octpole, **mdipole, **mquadpole}

    # 
    @timing.timed
    def build_osc_str(self):
        """build the oscillator strengths in the length and velocity
           gauges"""

        oscstr_l = self.build_osc_str_l()
        oscstr_v = self.build_osc_str_v()

        return {**oscstr_l, **oscstr_v}


    #-----------------------------------------------------------------
    # construction of multipole moments
    #
    #

    #
    @timing.timed
    def build_dip(self, mol, scf):
        """build the electric transition dipole moments in length
           and velocity gauges"""

        nao = mol.nao
           
        # dipole moment integrals (length gauge)
        mu_ao = mol.intor('int1e_r')
        mu_mo = self.ints_ao2mo(mol, scf, mu_ao)
        # dipole moment integrals (velocity gauge) -- imaginary
        # contribution only
        p_ao  = mol.intor('int1e_ipovlp')
        p_mo  = self.ints_ao2mo(mol, scf, p_ao)

        dipole = {}
        dipole['e_dipole_l'] = self.contract_transitions(mu_mo)
        dipole['e_dipole_v'] = self.contract_transitions(p_mo)
    
        return dipole

    # 
    @timing.timed
    def build_quad(self, mol, scf):
        """builds the electric quadrupole moments in length and velocity
           gauges"""

        nao = mol.nao

        # electronic quadrupole moment tensor (length gauge)
        Q_ao  = mol.intor('int1e_rr').reshape(3, 3, nao, nao)
        Q_mo  = self.ints_ao2mo(mol, scf, Q_ao)

        # electronic quadrupole moment tensor (velocity gauge)
        # returns imaginary part only, see Eq. 28
        Qp_ao  = mol.intor('int1e_irp', comp=9, hermi=0).reshape(3,3, nao, nao)
        Qp_ao += Qp_ao.transpose(1,0,3,2)
        Qp_mo  = self.ints_ao2mo(mol, scf, Qp_ao)

        quad = {}
        quad['e_quadrupole_l'] = self.contract_transitions(Q_mo)
        quad['e_quadrupole_v'] = self.contract_transitions(Qp_mo)
    
        return quad

    #
    @timing.timed
    def build_oct(self, mol, scf):
        """bulids the electric octupole moments in length and velocity
           gauges"""

        nao = mol.nao

        # electronic octupole moment tensor (length gauge)
        O_ao  = mol.intor('int1e_rrr').reshape(3, 3, 3, nao, nao)
        O_mo  = self.ints_ao2mo(mol, scf, O_ao)

        # electronic octupole moment tensor (velocity gauge)
        # returns imaginary part only, see Eq. 42
        Op_ao  = mol.intor('int1e_irrp', comp=27, hermi=0).reshape(3,3,3, nao, nao)
        Op_ao += Op_ao.transpose(2,1,0,4,3)
        Op_ao += mol.intor('int1e_irpr', comp=27,hermi=0).reshape(3,3,3, nao, nao)
        Op_mo  = self.ints_ao2mo(mol, scf, Op_ao)

        octo = {}
        octo['e_octupole_l'] = self.contract_transitions(O_mo)
        octo['e_octupole_v'] = self.contract_transitions(Op_mo)

        return octo

    #
    @timing.timed
    def build_mdip(self, mol, scf):
        """build the magnetic transition dipole moments in the
           velocity gauge"""

        # loop over the state pairs and construct transition
        # dipole moments
        nao = mol.nao

        # Magnetic quantities
        # ------------------------------------
        # magnetic first moment -- see Eq. 37 (negletcing spin for now)
        # returns imaginary part only
        mdp_ao = mol.intor('int1e_cg_irxp', comp=3, hermi=2)
        mdp_mo = self.ints_ao2mo(mol, scf, mdp_ao)

        mdip = {}
        mdip['m_dipole_v'] = self.contract_transitions(mdp_mo)

        return mdip

    # 
    @timing.timed
    def build_mquad(self, mol, scf):
        """builds the magnetic quadrupole moments the velocity gauge"""

        # see PySCF/Libcint documentation for d-shell ordering
        # https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf
        XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ = range(9)

        # loop over the state pairs and construct transition
        # dipole moments
        nao = mol.nao

        # magnetic quadrupole tensor -- see Eq. 46
        # returns imaginary part only
        ints     = mol.intor('int1e_irrp', comp=27, hermi=0).reshape(3,9,nao,nao)
        mquad_ao = (ints[:,[YZ,ZX,XY]] - ints[:,[ZY,XZ,YX]]).transpose(1,0,2,3)
        ints     = mol.intor('int1e_irpr', comp=27, hermi=0).reshape(9,3,nao,nao)
        mquad_ao+= ints[[YZ,ZX,XY]] - ints[[ZY,XZ,YX]]
        mquad_mo = self.ints_ao2mo(mol, scf, mquad_ao)

        mquad = {}
        mquad['m_quadrupole_v'] = self.contract_transitions(mquad_mo)

        return

    #
    @timing.timed
    def build_osc_str_l(self):
        """builds the oscillator strength in the length gauge
           see J. Chem. Phys. 137, 204106 (2012) for details. This
           is _actually_ in a mixed gauge, as the electric multipoles
           are in the length gauge and magnetic multipoles are in 
           the velocity gauge. I think this is actually standard,
           but need to confirm"""

        ntrans = len(self.trans_list)

        f0_l   = np.zeros((3, ntrans), dtype=float)
        f0_iso = np.zeros((ntrans), dtype=float)
        f2_l   = np.zeros((3, ntrans), dtype=float)
        f2_iso = np.zeros((ntrans), dtype=float)

        for tpair in self.trans_list:        

            b_st = tpair[0]
            k_st = tpair[1]
            indx = self.trans_list.index([b_st, k_st])

            # set up state and energy info
            b_ener = self.bra_obj.energy(b_st)
            k_ener = self.ket_obj.energy(k_st)
            alpha  = constants.fine_str
            de     = b_ener - k_ener
            de2    = de**2
            de3    = de**3

            # electric tensors
            mu     = self.multipole['e_dipole_l'][:, indx]
            Qab    = self.multipole['e_quadrupole_l'][:, :, indx]
            Oabc   = self.multipole['e_octupole_l'][:, :, :, indx]
            # magnetic tensors
            mdp    = self.multipole['m_dipole_v'][:, indx]
            Mqp    = self.multipole['m_quadrupole_v'][:, :, indx]

            # dipole oscillator strength - length gauge
            f0_l[:, indx] = (2.*de)*mu**2
            # isotropically averaged value
            f0_iso[indx]  = np.sum(f0_l[:, indx])/3.

            # compute second-order contributions
            # electric quadrupole contribution
            fQ_iso  = np.sum(Qab**2) - (1./3.)*np.trace(Qab)**2
            fQ_iso *= (de3 * alpha**2 / 20.)

            # magnetic dipole term
            fm_iso  = (2.*de/3.) * np.sum(mdp**2) 

            # electric dipole - electric octupole term
            fuO_iso = np.sum(np.array([[mu[b]*Oabc[a,a,b]
                                   for a in range(3)]
                                   for b in range(3)]))
            fuO_iso *= -2 * de3 * alpha**2 / 45.

            # electric dipole - magnetic quadrupole term
            fuM_iso  = np.sum(np.array([[[
                                 LeviCivita(a,b,c)*mu[b]*Mqp[c,a]
                               for a in range(3)]
                               for b in range(3)]
                               for c in range(3)]))
            fuM_iso *= de2 * alpha / 3.

            # polarization dependent Osc strengths to come
            # later
            f2_l[:, indx] = np.zeros((3,), dtype=float)
            # isotropically averaged value
            f2_av = f0_iso[indx] + \
                    fQ_iso + fm_iso + fuO_iso + fuM_iso
            f2_iso[indx] = f2_av
             
        oscstr_l = {}
        oscstr_l['f0_l']    = f0_l
        oscstr_l['f0iso_l'] = f0_iso

        oscstr_l['f2_l']    = f2_l
        oscstr_l['f2iso_l'] = f2_iso

        return oscstr_l

    #
    @timing.timed
    def build_osc_str_v(self):
        """builds the oscillator strength in the velocity gauge
           see J. Chem. Phys. 143, 234103 (2015) for details."""

        ntrans = len(self.trans_list)

        f0_v   = np.zeros((3, ntrans), dtype=float)
        f0_iso = np.zeros((ntrans), dtype=float)
        f2_v   = np.zeros((3, ntrans), dtype=float)
        f2_iso = np.zeros((ntrans), dtype=float)

        for tpair in self.trans_list[0]:
    
            b_st = tpair[0]
            k_st = tpair[1]
            indx = self.trans_list.index([b_st, k_st])
 
            # set up state and energy info
            b_ener = self.bra_obj.energy(b_st)
            k_ener = self.ket_obj.energy(k_st)
            alpha  = constants.fine_str
            de     = b_ener - k_ener
            de2    = de**2
            de3    = de**3

            # electric tensors
            mu     = self.multipole['e_dipole_v'][:, indx]
            Qab    = self.multipole['e_quadrupole_v'][:, :, indx]
            Oabc   = self.multipole['e_octupole_v'][:, :, :, indx]
            # magnetic tensors
            mdp    = self.multipole['m_dipole_v'][:, indx]
            Mqp    = self.multipole['m_quadrupole_v'][:, :, indx]

            # dipole oscillator strength - length gauge
            f0_v[:, indx] = (2./de)*mu**2
            # isotropically averaged value
            f0_iso[indx]  = np.sum(f0_v[:, indx])/3.

            # compute second-order contributions
            # electric quadrupole contribution
            fQ_iso  = np.sum(Qab**2) - (1./3.)*np.trace(Qab)**2
            fQ_iso *= (de * alpha**2) / 20.

            # magnetic dipole term
            fm_iso  = (de * alpha**2 / 6.) * np.sum(mdp**2)

            # electric dipole - electric octupole term
            fuO_iso = np.sum(np.array([[mu[b]*Oabc[a,a,b]
                               for a in range(3)]
                               for b in range(3)]))
            fuO_iso *= -2 * de * alpha**2 / 45.

            # electric dipole - magnetic quadrupole term
            fuM_iso  = np.sum(np.array([[[
                             LeviCivita(a,b,c)*mu[b]*Mqp[c,a]
                               for a in range(3)]
                               for b in range(3)]
                               for c in range(3)]))
            fuM_iso *= de * alpha**2 / 9.

            # polarization dependent Osc strengths to come
            # later
            f2_v[:, indx] = np.zeros(3, dtype=float)
            # isotropically averaged value
            f2_av = f0_iso[indx] + \
                    fQ_iso + fm_iso + fuO_iso + fuM_iso
            f2_iso[bk_indx] = f2_av

        oscstr_v = {}
        oscstr_v['f0_v']    = f0_v
        oscstr_v['f0iso_v'] = f0_iso

        oscstr_v['f2_v']    = f2_v
        oscstr_v['f2iso_v'] = f2_iso

        return oscstr_v

    #
    @timing.timed
    def build_natural_orbs(self):
        """build the natural difference orbitals and natural transition
           orbitals"""

        ntrans = len(self.trans_list)
        nmo    = self.bra_obj.scf.nmo
        nao    = self.bra_obj.scf.mol.nao
        mos    = self.bra_obj.scf.orbs

        nto    = np.zeros((nao, nmo, ntrans), dtype=float)
        ndo    = np.zeros((nao, nmo, ntrans), dtype=float)
        nto_wt = np.zeros((nmo, ntrans), dtype=float)
        ndo_wt = np.zeros((nmo, ntrans), dtype=float)

        for tpair in self.trans_list:

            b_st = tpair[0]
            k_st = tpair[1]
            indx = self.trans_list.index([b_st, k_st])

            tdm = self.tdm(b_st, k_st)
            if tdm is None:
                continue

            rdm_bra = self.bra_obj.rdm1(b_st)
            rdm_ket = self.ket_obj.rdm1(k_st)

            # first perform SVD of 1TDMs to get hole and
            # particle orbitals and weights and convert
            # orbitals to AO basis
            part, s, hole = np.linalg.svd(tdm)
            part_ao       = np.matmul(mos, part)
            hole_ao       = np.matmul(mos, hole.T)

            # normalize singular values
            s  = 0.5 * np.square(s)

            # sort the NTO amplitudes by decreasing magnitude
            ordr     = np.flip(np.argsort(s))
            hole_srt = hole_ao[:,ordr]
            part_srt = part_ao[:,ordr]
            s_srt    = s[ordr] 

            # find the number of particle/hole pairs with weight
            # greater than 0.01 (this should be a parameter), 
            # and only write these to file
            npairs = sum(chk > 0.01 for chk in s_srt)

            nto_pairs = np.zeros((nao, 2*npairs), dtype=float)
            wt_pairs  = np.zeros((2*npairs), dtype=float)

            # save NTOs and weights (define hole wts as neg.)
            nto_pairs[:,0::2] =  hole_srt[:,:npairs]
            nto_pairs[:,1::2] =  part_srt[:,:npairs]
            wt_pairs[0::2]    = -s_srt[:npairs] 
            wt_pairs[1::2]    =  s_srt[:npairs]

            nto[:, :2*npairs, bk_indx] = nto_pairs
            nto_wt[:2*npairs, bk_indx] = wt_pairs

            # next construct the natural difference densities
            delta      = rdm_bra - rdm_ket
            wt, ndo_mo = np.linalg.eigh(delta)

            # transform ndo to AO basis
            ndo_ao     = np.matmul(mos, ndo_mo)

            # sort NDO wts by increasing magnitude (hole 
            # orbitals to start, then particle
            ordr    = np.argsort(wt)
            ndo_srt = ndo_ao[:,ordr]
            wt_srt  = wt[ordr]         

            npairs  = sum(chk > 0.01 for chk in wt_srt)

            ndo_pairs = np.zeros((nao, 2*npairs), dtype=float)
            wt_pairs  = np.zeros((2*npairs), dtype=float)

            ndo_pairs[:,0::2] = ndo_srt[:,:npairs]
            ndo_pairs[:,1::2] = ndo_srt[:,-1:-npairs-1:-1]
            wt_pairs[0::2]    = wt_srt[:npairs]
            wt_pairs[1::2]    = wt_srt[-1:-npairs-1:-1]

            ndo[:, :2*npairs, bk_indx] = ndo_pairs
            ndo_wt[:2*npairs, bk_indx] = wt_pairs

        ndo_nto = {}
        ndo_nto['ndo']    = ndo
        ndo_nto['ndo_wt'] = ndo_wt
        ndo_nto['nto']    = nto
        ndo_nto['nto_wt'] = nto_wt

        return ndo_nto

    #
    #@timing.timed
    def ints_ao2mo(self, mol, scf, ao_ints):
        """contract ao_tensor with the mos to convert to mo basis"""

        # mos from scf objects
        mos     = scf.orbs

        # get dimension(s) 1-particle expansion
        nao     = ao_ints.shape[-1]

        # ensure this is consistent with mo length
        if mol.nao != nao:
            print('error in ints_ao2mo: mol.nao != ao_ints.shape[-1]: '
                    +str(mol.nao)+' != '+str(nao))
            sys.exit(1)

        return np.einsum('...pq,pi,qj->...ij', ao_ints, mos, mos, 
                                                     optimize='optimal')

    #
    @timing.timed
    def contract_tdm(self, mo_ints, tdm):
        """contract mo integrals with tdm to determine 
           tensor properties"""

        if tdm.shape != mo_ints.shape[-2:]:
            print('error in contract_tdm: tdm.shape != (nmo, nmo) = '
                    +str(tdm.shape))
            sys.exit(1)

        # get the dimsensions of each rank 
        rank_dim = mo_ints.shape[:-2]
        # convert to an array
        dim_tot  = np.prod(rank_dim)

        #run through each element of the tensor
        npole = np.einsum('xij,ij->x', 
                   mo_ints.reshape((dim_tot,)+tdm.shape ), tdm)

        # reshape output back to original rank
        return npole.reshape( rank_dim )

    # 
    @timing.timed
    def contract_transitions(self, mo_ints):
        """Loop over all transitions and contract mo_ints
           with the appropriate tdm"""

        # get the dimsensions of each rank 
        rank_dim = mo_ints.shape[:-2]
        ntrans   = len(self.trans_list)
        tens     = np.zeros((rank_dim + (ntrans,)), dtype=float)

        for tpair in self.trans_list:
            indx = self.trans_list.index(tpair)
            tdm = self.tdms[:, :, indx]
            tens[..., indx] = self.contract_tdm(mo_ints, tdm)

        return tens

    #
    def print_log(self):
        """print summary output to log file"""
        # for each initial state, write out a table of 
        # oscillator strenghts and transition dipole moments
        bsym_lbl = self.bra_obj.scf.mol.irreplbl
        ksym_lbl = self.bra_obj.scf.mol.irreplbl

        # print the header
        output.print_transition_header(self.label)

        # get the list of ket states
        ket_states, ket_syms = self.get_group_states('ket_states')

        # print a 'transition table' for each initial state
        for iket in range(len(ket_states)):
            
            # shift state indices from 0..n-1 to 1..n
            init_st   = ket_states[iket]+1
            init_sym  = ksym_lbl[ket_syms[iket]]
            final_st  = []
            final_sym = []
            exc_ener  = []
            osc_str   = [[] for i in range(5)]     

            for tpair in self.trans_list:
          
                b_st = tpair[0]
                k_st = tpair[1]
                indx = self.trans_list.index(tpair)

                if k_st != ket_states[iket]:
                    continue    

                # shift state indices from 0..n-1 to 1..n
                final_st.append(b_st+1)
                final_sym.append(bsym_lbl[self.bra_obj.state_sym(b_st)[0]]) 
                exc_ener.append(self.bra_obj.energy(b_st) - 
                                self.ket_obj.energy(k_st)) 
                osc_str[0].append(self.oscstr['f0iso_l'][indx])
                osc_str[1].append(self.oscstr['f2iso_l'][indx])
                osc_str[2].append(self.oscstr['f0iso_v'][indx])
                osc_str[3].append(self.oscstr['f2iso_v'][indx])
                osc_str[4].append(self.oscstr['f0_v'][:,indx])
 
            # print a 'transition table' for each initial state
            output.print_transition_table(init_st,
                                          init_sym,
                                          final_st,
                                          final_sym,
                                          exc_ener, 
                                          osc_str[0],
                                          osc_str[1],
                                          osc_str[2],
                                          osc_str[3],
                                          osc_str[4])


        return

    #
    def export_orbitals(self, orb_format='molden', orb_dir=True):
        """export all transition/difference density orbitals that
           in the object"""

        orb_types = ['nto', 'ndo']
        for tpair in self.trans_list:
            for otype in orb_types:
                self.export_orbitals_tran(tpair[0], tpair[1], 
                                          orb_type=otype, 
                                          file_format=orb_format, 
                                          orb_dir=orb_dir)

        return

    #
    def export_orbitals_tran(self, bra, ket, orb_type='nto', 
                                    file_format='molden', orb_dir=True):
        """print the natural transition orbitals and difference
           densities to file"""
        orb_types = ['nto', 'ndo']
        otype     = str(orb_type).lower()
        indx      = self.trans_list.index([bra, ket])

        if indx is None:
            print("bra="+str(bra)+" ket="+str(ket)+
                    " not in list of transitions")
            return False

        if otype not in orb_types:
            print("orbital format: "+orb_type+" not recognized.")
            return False

        b_sym = self.bra_obj.state_sym(bra)[0]
        k_sym = self.ket_obj.state_sym(ket)[0] 

        fname = otype.lower() + \
                  '.'+ str(ket+1)+'_' +str(k_sym.lower()) + \
                  '.'+ str(bra+1)+'_' +str(b_sym.lower()) + \
                  '_'+ str(file_format).lower() + str(self.label)

        if orb_dir:
            fname = 'orbs/'+str(fname)
            if not os.path.exists('orbs'):
                os.mkdir('orbs')

        # if a file called fname exists, delete it
        if os.path.isfile(fname):
            os.remove(fname)

        wts = self.nos[otype+'_wt'][:, indx]
        ncols = sum(chk > 1.e-16 for chk in np.absolute(wts))

        # import the appropriate library for the file_format
        if file_format in output.orb_formats:
            orbtype = importlib.import_module('graci.io.'+file_format)
        else:
            print('orbital format type=' + file_format +
                                        ' not found. exiting...')
            sys.exit(1)

        orbtype.write_orbitals(fname, 
                               self.bra_obj.scf.mol, 
                               self.nos[otype][:,:ncols, indx],
                               occ=wts[:ncols])

        return
