"""
module to compute functions of the 1-TDM
"""

import numpy as np
import sys as sys
from sympy import LeviCivita
import graci.core.params as params
import graci.core.orbitals as orbitals
import graci.interaction.interaction as interaction
import graci.utils.timing as timing
import graci.interfaces.bitci.bitsi_init as bitsi_init
import graci.interfaces.bitci.mrci_1tdm as mrci_1tdm
import graci.io.output as output
import graci.utils.constants as constants

class Transition(interaction.Interaction):
    """Sotransition class for determing transition moments between spin-orbit
       coupled states. This class extends Interaction and will accept either
       spinorbit objects, or simply method objects. If the latter, will convert
       them into trivial spin orbit objects (i.e. with a single set of states
       and a single multiplicity)
    """
    def __init__(self):
        # parent attributes
        super().__init__()
    
        # user defined quanties
        self.print_orbitals = False
        # for transition it's convenient to call these variables
        # init/final instead of ket/bra. Here we just say they
        # point to the same reference as the variables in the parent
        # class
        # print threshold for NDOS and NTOs
        self.print_orb_thresh = 0.01
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
        # representaion (adiabatic or diabatic for now)
        self.representation   = 'adiabatic'

        # ----------------------------------------------------------
        # internal class variables -- should not be accessed
        # directly
        self.bra_obj    = None
        self.ket_obj    = None

        # allowed representations
        self.allowed_reps = ['adiabatic', 'diabatic']
        
        # list of state group sassociated with the ket object
        self.state_grps = {}
        # list of transitions corresponding to the tdms
        self.trans_list = []
        # the tdms
        self.tdms       = None
        # various orbital quantities (NDOs, NTOs, wts, etc.)
        self.ndos       = None
        self.ndo_wts    = None
        self.ntos       = None
        self.nto_wts    = None
        # electric/magnetic multipole tensors -- all guages
        self.multipole  = {}
        # oscillator strengths -- all gauges
        self.oscstr     = {}

        # we currently require the same orbitals and mol objects
        # for the bra and ket -- this is them
        self.mos        = None
        self.mol        = None

    def copy(self):
        """create of deepcopy of self"""
        new = self.Transition()

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))

        return new


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
        self.mos, self.mol = self.check_bra_ket()
        nmo = self.mos.shape[1]

        # sanity check on the representation
        self.check_representation()
        
        # if bra and ket are the same object, only compute the unique
        # transition densities
        list_type = 'full'
        if self.same_obj(self.bra_obj, self.ket_obj):
            list_type = 'nodiag'

        # if all_final_states is true, over-ride current contents of
        # self.final_states
        if self.all_final_states:
            self.final_states = np.asarray(range(self.bra_obj.n_states()))

        # these are the initial and final states -- in whatever form
        # the bra and ket objects take (ci or postci objects)
        self.add_group('ket', self.ket_obj, self.init_states)
        self.add_group('bra', self.bra_obj, self.final_states)

        self.trans_list = self.build_pair_list('bra',
                                               'ket',
                                                pairs=list_type)

        # initialize the transition density matrices
        # These will be complex for spin-orbit coupled states
        self.tdms = np.zeros((len(self.trans_list), nmo, nmo),
                                                       dtype=np.cdouble)

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
                bitsi_init.init(bra_ci, ket_ci, 'tdm', self.verbose)
        
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

        # build the NDOs and NTOS, print to file if print_orbitals=True
        if self.print_orbitals:
            self.nto_wts, self.ntos = self.build_ntos(basis='ao', 
                                                      print_orbs=True)
            self.ndo_wts, self.ndos = self.build_ndos(basis='ao', 
                                                      print_orbs=True)

        # print the summary output
        if self.verbose:
            ndo_wts_mo, ndo_mo = self.build_ndos(basis='mo')
            self.print_log(ndo_wts_mo, ndo_mo)

        return

    #
    def tdm(self, b_st, k_st):
        """return the tdm for the bra and ket state indicated"""

        indx = self.trans_list.index([b_st, k_st])

        if indx is not None:
            return self.tdms[indx, :, :]
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
            sys.exit('ERROR: transition.run() passed '
                     + str(len(arg_list))
                     + ' objects. Expecting 2')

        # add both the bra and ket state groups
        lbls = ['bra', 'ket']
        for indx in range(len(arg_list)):
            arg_obj                     = arg_list[indx]
            self.state_grps[lbls[indx]] = []

            # if this is a ci method, create a new group associated
            # with the CI states
            if type(arg_obj).__name__ in params.ci_objs:
                lbl = str(lbls[indx]+'_grp0')
                self.state_grps[lbls[indx]].append(lbl)
                self.add_group(lbl, arg_obj, 
                               np.asarray(range(arg_obj.n_states())))

            # if this is a postci object, should already extend
            # interaction and have groups set up. Add all of them
            # here
            elif type(arg_obj).__name__ in params.postci_objs:
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
                sys.exit('ERROR: object passed to transition neither' +
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
                mos_b = self.get_obj(bra_lbl).mos
                mol_b = self.get_obj(bra_lbl).scf.mol

                for ket_lbl in self.state_grps['ket']:
                    mos_k = self.get_obj(ket_lbl).mos
                    mol_k = self.get_obj(ket_lbl).scf.mol

                    # sanity check that orbitals and geometry are
                    # the same
                    if np.any(mos_b != mos_k):
                        sys.exit('ERROR: transition moments require same '+
                                 ' bra/ket orbitalss')

                    if (mol_b.pymol().atom != mol_k.pymol().atom):
                        sys.exit('ERROR: transition moments require same '+
                                 ' geometry and basis set')

        else:
            mol_b = self.get_obj(self.state_grps['bra'][0]).scf.mol
            mos_b = self.get_obj(self.state_grps['bra'][0]).mos

        return mos_b, mol_b


    def check_representation(self):
        """checks that the user-specified representation is
        supported"""
        
        if self.representation not in self.allowed_reps:
            sys.exit('\n ERROR: unsupported representation: '
                     +self.representation
                     +'\n Allowed values: '+str(self.allowed_reps))

        return
        
    @timing.timed
    def build_tdms(self, bra, ket, trans_list, trans_list_sym):
        """grab the TDMs from bitsi and then reshape the list of
           TDMs into a more usable format"""

        # grab the tdms
        tdm_list = mrci_1tdm.tdm(bra, ket, trans_list_sym,
                                 self.representation)

        # make the tdm list
        nmo    = bra.nmo
        npairs = len(trans_list)
        tdm    = np.zeros((npairs, nmo, nmo), dtype=float)

        for indx in range(len(trans_list)):
            bk_st          = trans_list[indx]
            [birr, bst]    = bra.state_sym(bk_st[0])
            [kirr, kst]    = ket.state_sym(bk_st[1])
            sym_indx       = trans_list_sym[birr][kirr].index([bst,kst])
            tdm[indx, :,:] = tdm_list[birr][kirr][:, :, sym_indx]

        return tdm

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
            tdm     = tdm_blk[blk_ind, :, :]

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
                    self.tdms[ind, :, :] += cf * tdm

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

    # 
    def build_ndos(self, basis='ao', print_orbs=False):
        """
        Construct the NDOS for the transitions in the trans_list

        Arguments:
            print_orbs: if true, print orbitals to file

        Returns:
            ndos:       numpy array of NDOS
            ndo_wts:    corresponding weights for the NDOS
        """
     
        if basis not in ['ao','mo']:
            sys.exit('error in build_ndos: basis must be ao or mo')

        # check to see if rdm function exists in bra/ket objects.
        # If so, also construct NDOs
        b_rdm = getattr(self.bra_obj, "rdm", None)
        k_rdm = getattr(self.ket_obj, "rdm", None)
        if b_rdm is None or k_rdm is None:
            return None, None

        # NTOs are always constructed
        ndos = []
        wts  = []

        for it in range(len(self.trans_list)):
            bst = self.trans_list[it][0]
            kst = self.trans_list[it][1]

            wt, ndo = orbitals.build_ndos(b_rdm(bst), k_rdm(kst),
                                         basis=basis, mos=self.mos)
            ndos.append(ndo)
            wts.append(wt)

            if print_orbs and basis == 'ao':
                nao   = ndo.shape[0]
                nmo   = ndo.shape[1]
 
                npair = min(sum(chk > self.print_orb_thresh 
                                           for chk in wt), int(0.5*nmo))

                orbs = np.zeros((nao, 2*npair), dtype=float)
                occ  = np.zeros((2*npair), dtype=float)

                orbs[:,0::2] = ndo[:,:npair]
                orbs[:,1::2] = ndo[:,-1:-npair-1:-1]
                occ[0::2]    = wt[:npair]
                occ[1::2]    = wt[-1:-npair-1:-1]

                fname = 'ndo_'+str(kst+1)+'_to_'+str(bst+1)+'_molden'
                orbitals.export_orbitals(fname, self.mol, orbs,
                                orb_occ=occ, 
                                orb_dir='Transition.'+str(self.label),
                                fmt='molden')

        return wts, ndos


    #
    def build_ntos(self, basis='ao', print_orbs=False):
        """ contruct the NDOs and NTOs for the transitions in the 
            trans_list list

            Arguments:
              print_orbs: if true, print orbitals to file

            Returns:
              None
        """

        if basis not in ['ao','mo']:
            sys.exit('error in build_ntos: basis must be ao or mo')

        # NTOs are always constructed
        ntos = []
        wts  = []

        for it in range(len(self.trans_list)):
            bst = self.trans_list[it][0]
            kst = self.trans_list[it][1]

            wt, nto = orbitals.build_ntos(self.tdm(bst, kst), 
                                             basis=basis, mos=self.mos)
            ntos.append(nto)
            wts.append(wt)

            if print_orbs and basis == 'ao':
                nao   = nto.shape[1]
                nmo   = nto.shape[2]

                npair = sum(chk > self.print_orb_thresh 
                                                     for chk in wt[1,:])

                orbs = np.zeros((nao, 2*npair), dtype=float)
                occ  = np.zeros((2*npair), dtype=float)

                # save NTOs and weights (define hole wts as neg.)
                orbs[:,0::2] =  nto[0,:,:npair].real
                orbs[:,1::2] =  nto[1,:,:npair].real
                occ[0::2]    =  wt[0,:npair].real
                occ[1::2]    =  wt[1,:npair].real

                fname = 'nto_'+str(kst+1)+'_to_'+str(bst+1)+'_molden'
                orbitals.export_orbitals(fname, self.mol, orbs,
                                orb_occ=occ, 
                                orb_dir='Transition.'+str(self.label),
                                fmt='molden')

        return wts, ntos 

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
        mu_ao = self.mol.pymol().intor('int1e_r')
        mu_mo = self.ints_ao2mo(mu_ao)
        # dipole moment integrals (velocity gauge) -- imaginary
        # contribution only
        p_ao  = self.mol.pymol().intor('int1e_ipovlp')
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
        Q_ao  = self.mol.pymol().intor('int1e_rr').reshape(3, 3, nao, nao)
        Q_mo  = self.ints_ao2mo(Q_ao)

        # electronic quadrupole moment tensor (velocity gauge)
        # returns imaginary part only, see Eq. 28
        Qp_ao  = self.mol.pymol().intor('int1e_irp', 
                                  comp=9, hermi=0).reshape(3,3, nao, nao)
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
        O_ao  = self.mol.pymol().intor('int1e_rrr').reshape(3, 3, 3, nao, nao)
        O_mo  = self.ints_ao2mo(O_ao)

        # electronic octupole moment tensor (velocity gauge)
        # returns imaginary part only, see Eq. 42
        Op_ao  = self.mol.pymol().intor('int1e_irrp', 
                                 comp=27, hermi=0).reshape(3,3,3, nao, nao)
        Op_ao += Op_ao.transpose(2,1,0,4,3)
        Op_ao += self.mol.pymol().intor('int1e_irpr', 
                                 comp=27,hermi=0).reshape(3,3,3, nao, nao)
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
        mdp_ao = self.mol.pymol().intor('int1e_cg_irxp', comp=3, hermi=2)
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
        ints     = self.mol.pymol().intor('int1e_irrp', 
                                 comp=27, hermi=0).reshape(3,9,nao,nao)
        mquad_ao = (ints[:,[YZ,ZX,XY]] - ints[:,[ZY,XZ,YX]]).transpose(1,0,2,3)
        ints     = self.mol.pymol().intor('int1e_irpr',
                                  comp=27, hermi=0).reshape(9,3,nao,nao)
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
                args = [s+1 for s in tpair] + [abs(f2_av), f2_av.imag]
                pstr =' Discarding imaginary component of Oscillator'+ \
                      ' Strength for transition: {0:3d}<-{1:3d}, '+    \
                      ' Abs(f2)={2:9.5f}/Im(f2)={3:9.5f}'
                print(pstr.format(*args))
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
                args = [s+1 for s in tpair] + [abs(f2_av), f2_av.imag]
                pstr =' Discarding imaginary component of Oscillator'+ \
                      ' Strength for transition: {0:3d}<-{1:3d}, '+ \
                      ' Abs(f2)={2:9.5f}/Im(f2)={3:9.5f}'
                print(pstr.format(*args))
            f2_iso[indx] = f2_av.real

        oscstr_v = {}
        oscstr_v['f0_v']    = f0_v
        oscstr_v['f0iso_v'] = f0_iso

        oscstr_v['f2_v']    = f2_v
        oscstr_v['f2iso_v'] = f2_iso

        return oscstr_v

    #
    @timing.timed
    def ints_ao2mo(self, ao_ints):
        """contract ao_tensor with the mos to convert to mo basis"""

        # get dimension(s) 1-particle expansion
        nao     = ao_ints.shape[-1]

        # ensure this is consistent with mo length
        if self.mol.nao != nao:
            print('error in ints_ao2mo: mol.nao != ao_ints.shape[-1]: '
                    +str(self.mol.nao)+' != '+str(nao))
            sys.exit(1)

        return np.einsum('...pq,pi,qj->...ij', ao_ints, 
                              self.mos, self.mos, optimize='optimal')

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
            tdm = self.tdms[indx, :, :]
            tens[..., indx] = self.contract_tdm(mo_ints, tdm)

        return tens

    #
    def print_log(self, ndo_wts, ndo_orbs):
        """print summary output to log file"""
        # for each initial state, write out a table of 
        # oscillator strenghts and transition dipole moments
        # print the header
        output.print_transition_header(self.label)

        # get the lists of bra ket states
        ket_states = self.get_states('ket')
        bra_states = self.get_states('bra')

        # get state symmetries. If not defined, use C1 sym labels
        ket_syms = self.get_syms('ket')
        if ket_syms is None:
            ket_syms = [0]*len(ket_states)
            ksym_lbl = ['A']*len(ket_states)
        else:
            ksym_lbl = self.mol.irreplbl

        bra_syms = self.get_syms('bra')
        if bra_syms is None:
            bra_syms = [0]*len(bra_states)
            bsym_lbl = ['A']*len(bra_states)
        else:
            bsym_lbl = self.mol.irreplbl

        # print a 'transition table' for each initial state
        for iket in range(len(ket_states)):
            
            # shift state indices from 0..n-1 to 1..n
            init_st   = ket_states[iket]+1
            init_sym  = ksym_lbl[ket_syms[iket]]
            final_st  = []
            final_sym = []
            exc_ener  = []
            osc_str   = [[] for i in range(5)]     
            promo_num = []

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

                # this check on the NDOs could be made stronger..
                pa, pd = orbitals.promotion_numbers(ndo_wts[indx],
                                                    ndo_orbs[indx])
                promo_num.append([pa,pd])

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
                                          osc_str[4],
                                          promo_num)


        return

