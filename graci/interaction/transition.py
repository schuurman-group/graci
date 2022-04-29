"""
module to compute functions of the 1-TDM
"""

import os as os
import numpy as np
import sys as sys
import importlib
from sympy import LeviCivita
import graci.core.params as params
import graci.interaction.interaction as interaction
import graci.utils.timing as timing
import graci.bitcitools.bitsi_init as bitsi_init
import graci.bitcitools.mrci_1tdm as mrci_1tdm
import graci.io.output as output
import graci.utils.constants as constants

class Transition(interaction.Interaction):
    """Sotransition class for determing transition moments between spin-orbit
       copuled states. This class extends transition and will accpet either
       spinorbit objects, or, simply method objects. If the latter, will convert
       them into trivial spin orbit objects (i.e. with a single set of states and
       a single multiplicity)
    """
    def __init__(self):
        # parent attributes
        super().__init__()
    
        # input syntax should be identical to transition
        # name for this object
        self.label          = 'Transition'

        # user defined quanties
        self.print_orbitals = False
        # for transition it's convenient to call these variables
        # init/final instead of ket/bra. Here we just say they
        # point to the same reference as the variables in the parent
        # class
        # list of initial states
        self.init_states      = None
        # label of section to get initial states
        self.init_label       = None
        # list of final states
        self.final_states     = None
        # label of section to get final states 
        self.final_label      = None
        # often will want to include all states in bra object
        # if so convenient to have simple boolean to select that
        self.all_final_states = False

        # ----------------------------------------------------------
        # internal class variables -- should not be accessed
        # directly
        self.bra_obj    = None
        self.ket_obj    = None

        # list of state group sassociated with the ket object
        self.state_grps = {}
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

        # we currently require the same scf and mol objects
        # for the bra and ket -- this is them
        self.scf        = None
        self.mol        = None

    #
    @timing.timed
    def run(self, arg_list):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used.
        """
 
        # set the bra/ket objects and add the state groups associated
        # with each 
        self.bra_obj, self.ket_obj = self.add_state_grps(arg_list)

        # do some sanity checks on the bra/ket objects
        self.scf, self.mol = self.check_bra_ket()
      
        # if bra and ket are the same object, only compute the unique
        # transition densities
        list_type = 'full'
        if self.same_obj(self.bra_obj, self.ket_obj):
            list_type = 'nodiag'

        # if all_final_states is true, over-ride current contents of
        # self.final_states
        if self.all_final_states:
            self.final_states = range(self.bra_obj.n_states())

        # these are the initial and final states -- in whatever form
        # the bra and ket objects take (ci or postci objects)
        self.add_group('ket', self.ket_obj, self.init_states)
        self.add_group('bra', self.bra_obj, self.final_states)

        self.trans_list = self.build_pair_list('bra',
                                               'ket',
                                                pairs=list_type)

        # initialize the transition density matrices
        # These will be complex for spin-orbit coupled states
        self.tdms = np.zeros((self.scf.nmo, self.scf.nmo, 
                              len(self.trans_list)), dtype=np.cdouble)

        # loop over all pairs of spin free states in both the bra
        # and ket objects. In the future, may be prudent to check
        # init and final states and just which states are necessary
        for k_lbl in self.state_grps['ket']:
            ket_spin = self.get_spins(k_lbl)
            ket_ci   = self.get_obj(k_lbl) 

            for b_lbl in self.state_grps['bra']:
                bra_spin = self.get_spins(b_lbl)
                bra_ci   = self.get_obj(b_lbl)

                # if different spin manifold, tdm contribution is zero
                if ket_spin.mult != bra_spin.mult:
                    continue

                # if the same object, only build lower diagonal
                if self.same_obj(ket_ci, bra_ci): 
                    pair_type = 'nodiag'
                else:
                    pair_type = 'full'
                
                # initialize the bitsi library for the calculation of 1-TDMs
                bitsi_init.init(bra_ci, ket_ci, 'tdm')
        
                # this is main transition_list: stored by adiabatic label
                blks      = self.build_pair_list(b_lbl,
                                                 k_lbl, 
                                                 pairs=pair_type)

                # the transitions, ordered by interacting irreps, is how
                # we call bitsi
                blks_sym  = self.build_pair_list(b_lbl, 
                                                 k_lbl, 
                                                 pairs=pair_type, 
                                                 sym_blk=True)

                # tdms is a vector of nmo x nmo transition densities
                tdm_blk = self.build_tdms(bra_ci, ket_ci, blks, blks_sym)

                # rotate tdms into spin states, if need be
                self.rotate_tdms(b_lbl, k_lbl, blks, tdm_blk)

                # finalize the bitsi library
                bitsi_init.finalize()

        # build the multipole moments  -- easier to just do this once
        # for all transitions
        self.multipole = self.build_multipole()

        # compute oscillator strengths
        self.oscstr = self.build_osc_str()

        # print orbitals if requested
        if self.print_orbitals:
            # first construct the orbitals
            self.nos = {**self.build_ntos(), **self.build_ndos()}
            self.export_orbitals(orb_type='nto')
            self.export_orbitals(orb_type='ndo')

        # print the summary output
        self.print_log()

        return

    #
    def tdm(self, b_st, k_st):
        """return the tdm for the bra and ket state indicated"""

        indx = self.trans_list.index([b_st, k_st])

        if indx is not None:
            return self.tdms[:, :, indx]
        else:
            return None
 
    #
    def transition_multipole(self, b_st, k_st, name, gauge='velocity'):
        """return the transition mulitpole tensor for pairs of
           states ordered by energy"""

        # append a '_l' or '_v' to indicate length or velocity gauge
        keyname = name+'_'+str(gauge)[0]
        indx = self.trans_list.index([b_st, k_st])

        if keyname in self.multipole.keys() and indx is not None:
            return self.multipole[keyname][..., indx]
        else:
            return None

    #
    def osc_strength(self, b_st, k_st, order=0, gauge='velocity', iso=True):
        """return the transition dipole for pairs of 
           states ordered by energy in the requested gauge. If iso is
           True, return a scalar that is the isotropically averaged 
           oscillator strenth, else, return a cartesian vector for
           (x,y,z) components
        """

        istr = ''
        if iso:
            istr = 'iso'

        keyname = 'f'+str(order)+istr+'_'+str(gauge)[0]
        indx = self.trans_list.index([b_st, k_st])

        if keyname in self.oscstr.keys() and indx is not None:
            return self.oscstr[keyname][..., indx]
        else:
            return None

#--------------------------------------------------------------------------
# "Private" class methods
#
    #
    def add_state_grps(self, arg_list):
        """
        add the bra and ket state groups
        """
        if len(arg_list) != 2:
            sys.exit('transition.run() passed ' + str(len(arg_list)) +
                     ' objects. Expecting 2')

        # add both the bra and ket state groups
        lbls = ['bra', 'ket']
        for indx in range(len(arg_list)):
            arg_obj                     = arg_list[indx]
            self.state_grps[lbls[indx]] = []

            # if this is a ci method, create a new group associated
            # with the CI states
            if type(arg_obj).__name__ in params.ci_objs:
                lbl = lbls[indx]+'_grp0'
                self.state_grps[lbls[indx]].append(lbl)
                self.add_group(lbl, arg_obj, range(arg_obj.n_states()))

            # if this is a postci object, should already extend
            # interaction and have groups set up. Add all of them
            # here
            elif type(bra_ket).__name__ in params.postci_objs:
                grp_lbls = arg_obj.get_lbls()

                for ibk in range(len(grp_lbls)):
                    lbl  = lbls[indx]+'_'+grp_lbls[ibk]
                    self.state_grps[lbls[indx]].append(lbl)
                    self.add_group(lbl,
                                   arg_obj.get_obj(grp_lbls[ibk]),
                                   arg_obj.get_states(grp_lbls[ibk]))

            # if not a ci_obj and postci_obj, we don't know what to do 
            # with this
            else:
                sys.exit(' object passed to transition neither' +
                         ' ci_obj or postci_obj. Exiting...')

        return arg_list[0], arg_list[1]

    #
    def check_bra_ket(self):
        """
        do some sanity checks on the arguments passed to run
        """

        # now that we have bra and ket groups set up, do some
        # sanity checking on the orbitals/geometries involved
        if not self.same_obj(self.bra_obj, self.ket_obj):
            for bra_lbl in self.state_grps['bra']:
                scf_b = self.get_obj(bra_lbl).scf
                for ket_lbl in self.state_grps['ket']:
                    scf_k = self.get_obj(ket_lbl).scf

                    # sanity check that orbitals and geometry are
                    # the same
                    if np.any(scf_b.orbs != scf_k.orbs):
                        sys.exit('transition moments require same '+
                                 ' bra/ket orbitalss')

                    if (scf_b.mol.pymol().atom !=
                              scf_k.mol.pymol().atom):
                        sys.exit('transition moments require same '+
                                 ' geometry and basis set')

        else:
            scf_b = self.get_obj(self.state_grps['bra'][0]).scf

        return scf_b, scf_b.mol.pymol()


    @timing.timed
    def build_tdms(self, bra, ket, trans_list, trans_list_sym):
        """grab the TDMs from bitsi and then reshape the list of
           TDMs into a more usable format"""

        # grab the tdms
        tdm_list = mrci_1tdm.tdm(bra, ket, trans_list_sym)

        # make the tdm list
        nmo    = bra.scf.nmo
        npairs = len(trans_list)
        tdms   = np.zeros((nmo, nmo, npairs), dtype=float)

        for indx in range(len(trans_list)):
            bk_st          = trans_list[indx]
            [birr, bst]    = bra.state_sym(bk_st[0])
            [kirr, kst]    = ket.state_sym(bk_st[1])
            sym_indx       = trans_list_sym[birr][kirr].index([bst,kst])
            tdms[:,:,indx] = tdm_list[birr][kirr][:, :, sym_indx]

        return tdms

    @timing.timed
    def rotate_tdms(self, b_lbl, k_lbl, blk_lst, tdm_blk):
        """
        Rotate the spin-free tdms to the spin-orbit states. The 
        bra and/or ket object may be spin-free or spin-orbit
        coupled objects
        """

        # this will work for now, but feels inelegant
        if type(self.bra_obj).__name__ == 'Spinorbit':
            bra_spin = self.get_spins(b_lbl)
            grp_indx = self.state_grps['bra'].index(b_lbl)
            bra_lbl  = self.bra_obj.get_lbls()[grp_indx]
        else:
            bra_spin = None

        if type(self.ket_obj).__name__ == 'Spinorbit':
            ket_spin = self.get_spins(k_lbl)
            grp_indx = self.state_grps['ket'].index(k_lbl)
            ket_lbl  = self.ket_obj.get_lbls()[grp_indx] 
        else:
            ket_spin = None

        # run through trans_list and contribute each #
        for blk_pair in blk_lst:
            blk_ind = blk_lst.index(blk_pair)
            tdm     = tdm_blk[:, :, blk_ind]

            for pair in self.trans_list:
                ind = self.trans_list.index(pair)

                # if the bra state doesn't contribute, skip
                if bra_spin is None and blk_pair[0] != pair[0]:
                    continue
                else:
                    b_cf = [1.+0.j]

                # if the ket state doesn't contribute, skip
                if ket_spin is None and blk_pair[1] != pair[1]:
                    continue
                else:
                    k_cf = [1.+0.j]

                # pull the coefficients from the state vectors
                # we need to calculate contribution of this spin-free
                # pair to the transition density
                if bra_spin is not None:
                    b_indxs = [self.bra_obj.soc_index(bra_lbl,
                                                blk_pair[0], ms)
                                                for ms in bra_spin.M]
                    b_cf   = [self.bra_obj.soc_state(pair[0])[b_indx]
                                                for b_indx in b_indxs]

                if ket_spin is not None:
                    k_indxs = [self.ket_obj.soc_index(ket_lbl,
                                                blk_pair[1], ms)
                                                for ms in ket_spin.M]
                    k_cf   = [self.bra_obj.soc_state(pair[1])[k_indx]
                                                for k_indx in k_indxs]

                cf = sum( [bcf*kcf for kcf in k_cf for bcf in b_cf] )
                if abs(cf) > 1.e-16:
                    self.tdms[:, :, ind] += cf * tdm

        return


    #
    @timing.timed
    def build_multipole(self):
        """builds the multipole moment transitions tensors"""

        # build the electronic tensors: dipole, quadrupole, octupole
        dipole   = self.build_dip()
        quadpole = self.build_quad()
        octpole  = self.build_oct()

        # build the magnetic moment tensors: dipole, quadrupole
        mdipole   = self.build_mdip()
        mquadpole = self.build_mquad()

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
    def build_dip(self):
        """build the electric transition dipole moments in length
           and velocity gauges"""

        nao = self.mol.nao
           
        # dipole moment integrals (length gauge)
        mu_ao = self.mol.intor('int1e_r')
        mu_mo = self.ints_ao2mo(mu_ao)
        # dipole moment integrals (velocity gauge) -- imaginary
        # contribution only
        p_ao  = self.mol.intor('int1e_ipovlp')
        p_mo  = self.ints_ao2mo(p_ao)

        dipole = {}
        dipole['e_dipole_l'] = self.contract_transitions(mu_mo)
        dipole['e_dipole_v'] = self.contract_transitions(p_mo)
    
        return dipole

    # 
    @timing.timed
    def build_quad(self):
        """builds the electric quadrupole moments in length and velocity
           gauges"""

        nao = self.mol.nao

        # electronic quadrupole moment tensor (length gauge)
        Q_ao  = self.mol.intor('int1e_rr').reshape(3, 3, nao, nao)
        Q_mo  = self.ints_ao2mo(Q_ao)

        # electronic quadrupole moment tensor (velocity gauge)
        # returns imaginary part only, see Eq. 28
        Qp_ao  = self.mol.intor('int1e_irp', comp=9, hermi=0).reshape(3,3, nao, nao)
        Qp_ao += Qp_ao.transpose(1,0,3,2)
        Qp_mo  = self.ints_ao2mo(Qp_ao)

        quad = {}
        quad['e_quadrupole_l'] = self.contract_transitions(Q_mo)
        quad['e_quadrupole_v'] = self.contract_transitions(Qp_mo)
    
        return quad

    #
    @timing.timed
    def build_oct(self):
        """bulids the electric octupole moments in length and velocity
           gauges"""

        nao = self.mol.nao

        # electronic octupole moment tensor (length gauge)
        O_ao  = self.mol.intor('int1e_rrr').reshape(3, 3, 3, nao, nao)
        O_mo  = self.ints_ao2mo(O_ao)

        # electronic octupole moment tensor (velocity gauge)
        # returns imaginary part only, see Eq. 42
        Op_ao  = self.mol.intor('int1e_irrp', comp=27, hermi=0).reshape(3,3,3, nao, nao)
        Op_ao += Op_ao.transpose(2,1,0,4,3)
        Op_ao += self.mol.intor('int1e_irpr', comp=27,hermi=0).reshape(3,3,3, nao, nao)
        Op_mo  = self.ints_ao2mo(Op_ao)

        octo = {}
        octo['e_octupole_l'] = self.contract_transitions(O_mo)
        octo['e_octupole_v'] = self.contract_transitions(Op_mo)

        return octo

    #
    @timing.timed
    def build_mdip(self):
        """build the magnetic transition dipole moments in the
           velocity gauge"""

        # loop over the state pairs and construct transition
        # dipole moments
        nao = self.mol.nao

        # Magnetic quantities
        # ------------------------------------
        # magnetic first moment -- see Eq. 37 (negletcing spin for now)
        # returns imaginary part only
        mdp_ao = self.mol.intor('int1e_cg_irxp', comp=3, hermi=2)
        mdp_mo = self.ints_ao2mo(mdp_ao)

        mdip = {}
        mdip['m_dipole_v'] = self.contract_transitions(mdp_mo)

        return mdip

    # 
    @timing.timed
    def build_mquad(self):
        """builds the magnetic quadrupole moments the velocity gauge"""

        # see PySCF/Libcint documentation for d-shell ordering
        # https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf
        XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ = range(9)

        # loop over the state pairs and construct transition
        # dipole moments
        nao = self.mol.nao

        # magnetic quadrupole tensor -- see Eq. 46
        # returns imaginary part only
        ints     = self.mol.intor('int1e_irrp', comp=27, hermi=0).reshape(3,9,nao,nao)
        mquad_ao = (ints[:,[YZ,ZX,XY]] - ints[:,[ZY,XZ,YX]]).transpose(1,0,2,3)
        ints     = self.mol.intor('int1e_irpr', comp=27, hermi=0).reshape(9,3,nao,nao)
        mquad_ao+= ints[[YZ,ZX,XY]] - ints[[ZY,XZ,YX]]
        mquad_mo = self.ints_ao2mo(mquad_ao)

        mquad = {}
        mquad['m_quadrupole_v'] = self.contract_transitions(mquad_mo)

        return mquad

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

        # Levi-Civita 3x3x3 tensor
        LC     = np.asarray([[[LeviCivita(a,b,c) for a in range(3)] 
                                                 for b in range(3)] 
                                                 for c in range(3)], 
                                                 dtype=float)

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
            f0_l[:, indx] = (2.*de) * (np.conj(mu)*mu).real
            # isotropically averaged value
            f0_iso[indx]  = np.sum(f0_l[:, indx])/3.

            # compute second-order contributions
            # electric quadrupole contribution
            fQ_iso  = np.sum(np.conj(Qab)*Qab) - \
                      (1./3.)*np.conj(np.trace(Qab))*np.trace(Qab)
            fQ_iso *= (de3 * alpha**2 / 20.)

            # magnetic dipole term
            fm_iso  = (2.*de/3.) * np.sum(np.conj(mdp)*mdp) 

            # electric dipole - electric octupole term
            fuO_iso = np.sum(np.array([[mu[b]*Oabc[a,a,b]
                                   for a in range(3)]
                                   for b in range(3)]))
            fuO_iso *= -2 * de3 * alpha**2 / 45.

            # electric dipole - magnetic quadrupole term
            fuM_iso  = np.sum(np.array([[[
                               LC[a,b,c] * mu[b] * Mqp[c,a]
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
            if abs(f2_av.imag) > 0.1 * abs(f2_av):
                print('discarding imaginary component of elec. dip'+
                      ' - elec. oct term for transition '+str(tpair))
            f2_iso[indx] = f2_av.real
             
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

        # Levi-Civita 3x3x3 tensor
        LC     = np.asarray([[[LeviCivita(a,b,c) for a in range(3)]
                                                 for b in range(3)]
                                                 for c in range(3)],
                                                 dtype=float)

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
            mu     = self.multipole['e_dipole_v'][:, indx]
            Qab    = self.multipole['e_quadrupole_v'][:, :, indx]
            Oabc   = self.multipole['e_octupole_v'][:, :, :, indx]
            # magnetic tensors
            mdp    = self.multipole['m_dipole_v'][:, indx]
            Mqp    = self.multipole['m_quadrupole_v'][:, :, indx]

            # dipole oscillator strength - length gauge
            f0_v[:, indx] = (2./de) * (np.conj(mu)*mu).real
            # isotropically averaged value
            f0_iso[indx]  = np.sum(f0_v[:, indx])/3.

            # compute second-order contributions
            # electric quadrupole contribution
            fQ_iso  = np.sum(np.conj(Qab)*Qab) - \
                      (1./3.)*np.conj(np.trace(Qab))*np.trace(Qab)
            fQ_iso *= (de * alpha**2) / 20.

            # magnetic dipole term
            fm_iso  = (de * alpha**2 / 6.) * np.sum(np.conj(mdp)*mdp)

            # electric dipole - electric octupole term
            fuO_iso = np.sum(np.array([[mu[b]*Oabc[a,a,b]
                               for a in range(3)]
                               for b in range(3)]))
            fuO_iso *= -2 * de * alpha**2 / 45.

            # electric dipole - magnetic quadrupole term
            fuM_iso  = np.sum(np.array([[[
                               LC[a,b,c] * mu[b] * Mqp[c,a]
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
            if abs(f2_av.imag) > 0.1 * abs(f2_av):
                print('discarding imaginary component of elec. dip'+
                      ' - elec. oct term for transition '+str(tpair))
            f2_iso[indx] = f2_av.real

        oscstr_v = {}
        oscstr_v['f0_v']    = f0_v
        oscstr_v['f0iso_v'] = f0_iso

        oscstr_v['f2_v']    = f2_v
        oscstr_v['f2iso_v'] = f2_iso

        return oscstr_v

    #
    @timing.timed
    def build_ntos(self):
        """build the natural transition orbitals"""

        nmo    = self.scf.nmo
        nao    = self.mol.nao
        mos    = self.scf.orbs

        nto    = []
        nto_wt = []

        for tpair in self.trans_list:

            b_st = tpair[0]
            k_st = tpair[1]
            indx = self.trans_list.index([b_st, k_st])

            tdm = self.tdm(b_st, k_st)
            if tdm is None:
                continue

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
            nto_pairs[:,0::2] =  hole_srt[:,:npairs].real
            nto_pairs[:,1::2] =  part_srt[:,:npairs].real
            wt_pairs[0::2]    = -s_srt[:npairs].real 
            wt_pairs[1::2]    =  s_srt[:npairs].real

            nto.append(nto_pairs)
            nto_wt.append(wt_pairs)

        ntos           = {}
        ntos['nto']    = nto
        ntos['nto_wt'] = nto_wt

        return ntos

    #
    @timing.timed
    def build_ndos(self):
        """ build natural difference orbitals"""

        # first check if bra and ket objects define an
        # rdm method. If not, return an empty dictionary
        
        try:
            if (hasattr(self.bra_obj.__class__, 'rdm') and 
                 callable(getattr(self.bra_obj.__class__, 'rdm')) and
                  hasattr(self.ket_obj.__class__, 'rdm') and
                    callable(getattr(self.ket_obj.__class__, 'rdm'))):
                rdm_exists = True
            else:
                rdm_exists = False
        except:
            rdm_exists = False

        if not rdm_exists:
            return {}

        nmo    = self.scf.nmo
        nao    = self.mol.nao
        mos    = self.scf.orbs

        ndo    = []
        ndo_wt = []

        for tpair in self.trans_list:

            b_st = tpair[0]
            k_st = tpair[1]
            indx = self.trans_list.index([b_st, k_st])

            rdm_bra = self.bra_obj.rdm(b_st)
            rdm_ket = self.ket_obj.rdm(k_st)

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

            ndo.append(ndo_pairs)
            ndo_wt.append(wt_pairs)

        ndos = {}
        ndos['ndo']    = ndo
        ndos['ndo_wt'] = ndo_wt

        return ndos

    #
    #@timing.timed
    def ints_ao2mo(self, ao_ints):
        """contract ao_tensor with the mos to convert to mo basis"""

        # mos from scf objects
        mos     = self.scf.orbs

        # get dimension(s) 1-particle expansion
        nao     = ao_ints.shape[-1]

        # ensure this is consistent with mo length
        if self.mol.nao != nao:
            print('error in ints_ao2mo: mol.nao != ao_ints.shape[-1]: '
                    +str(self.mol.nao)+' != '+str(nao))
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
        tens     = np.zeros((rank_dim + (ntrans,)), dtype=np.cdouble)

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
        # print the header
        output.print_transition_header(self.label)

        # get the list of ket states
        ket_states = self.get_states('ket')
        bra_states = self.get_states('bra')

        # get state symmetries. If not defined, use C1 sym labels
        ket_syms = self.get_syms('ket_states')
        if ket_syms is None:
            ket_syms = [0]*len(ket_states)
            ksym_lbl = ['A']*len(ket_states)
        else:
            ksym_lbl = self.scf.mol.irreplbl

        bra_syms = self.get_syms('bra_states')
        if bra_syms is None:
            bra_syms = [0]*len(bra_states)
            bsym_lbl = ['A']*len(bra_states)
        else:
            bsym_lbl = self.scf.mol.irreplbl


        # print a 'transition table' for each initial state
        for iket in range(len(ket_states)):
            
            # shift state indices from 0..n-1 to 1..n
            init_st   = ket_states[iket]+1
            init_sym  = ksym_lbl[ket_syms[iket]]
            final_st  = []
            final_sym = []
            exc_ener  = []
            osc_str   = [[] for i in range(5)]     

            for ibra in range(len(bra_states)):

                # if [iket,ibra] not in trans_list, continue to
                # next pair
                try:
                    indx = self.trans_list.index([bra_states[ibra], 
                                                  ket_states[iket]])
                except:
                    continue

                # shift state indices from 0..n-1 to 1..n
                final_st.append(bra_states[ibra]+1)
                final_sym.append(bsym_lbl[bra_syms[ibra]])
                exc_ener.append(self.bra_obj.energy(bra_states[ibra]) - 
                                self.ket_obj.energy(ket_states[iket])) 
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
    def export_orbitals(self, orb_type='nto', orb_format='molden', 
                                              orb_dir=True):
        """export all transition/difference density orbitals that
           in the object"""

        if orb_type not in list(self.nos.keys()):
            print(' orb_type=' + str(orb_type) + ' not in ' + 
                                      str(list(self.nos.keys())))
            return False

        for tpair in self.trans_list:
            self.export_orbitals_tran(tpair[0], tpair[1], 
                                      orb_type=orb_type, 
                                      file_format=orb_format, 
                                      orb_dir=orb_dir)

        return True

    #
    def export_orbitals_tran(self, bra, ket, orb_type='nto', 
                            file_format='molden', orb_dir=True):
        """print the natural transition orbitals and difference
           densities to file"""
        orb_types = ['nto', 'ndo']

        indx      = self.trans_list.index([bra, ket])
        if indx is None:
            print("bra="+str(bra)+" ket="+str(ket)+
                    " not in list of transitions")
            return False

        if orb_type not in list(self.nos.keys()):
            print("orbital format: "+orb_type+" not recognized.")
            return False

        fname = orb_type.lower() + \
                  '_'+ str(ket+1)+'-' +str(bra+1) + \
                  '_'+ str(file_format).lower()

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
                               self.scf.mol, 
                               self.nos[orb_type][indx],
                               occ=self.nos[orb_type+'_wt'][indx])

        return
