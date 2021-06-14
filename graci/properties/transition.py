"""
module for compute moments of the dipole operator
for a given electronic state
"""
import numpy as np
from sympy import LeviCivita
import graci.utils.timing as timing
import graci.core.libs as libs
import graci.citools.mrci_1tdm as mrci_1tdm
import graci.io.output as output
import graci.utils.constants as constants
import sys as sys

class Transition:
    """Transition class for determing transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # user defined quanties 
        self.init_states      = None
        self.final_states     = None
        self.init_states_sym  = None
        self.final_states_sym = None
        self.all_final_states = False
        self.init_label      = None
        self.final_label     = None

        # global variables
        # method object for the bra states
        self.bra_obj      = None
        # method object for the ket states
        self.ket_obj      = None
        # the tdms
        self.tdms         = None
        # comprehensive list of all transition pairs
        self.trans_list   = None
        # list of initial states
        self.ket_list     = None
        # symmetries of initial states
        self.ket_sym      = None
        # various orbital quantities (NDOs, NTOs, wts, etc.)
        self.nos          = {} 
        # electric/magnetic multipole tensors -- all guages
        self.multipole    = {}
        # oscillator strengths -- all gauges
        self.oscstr       = {}
       
        self.label         = 'transition'

    #-----------------------------------------------------------

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'transition'

    #
    def set_ket(self, ket_method):
        """set the method used to compute the ket state(s)"""

        try:
            self.init_label = ket_method.label
        except:
            self.init_label = None
            return

        self.ket_obj = ket_method
        return
 
    # 
    def ket_exists(self):
        """return True is self.ket_obj is not None"""

        return self.ket_obj is not None

    #
    def set_bra(self, bra_method):
        """set the method used to compute the bra state(s)"""

        try:
            self.final_label = bra_method.label
        except:
            self.final_label = None
            return

        self.bra_obj = bra_method
        return

    #
    def bra_exists(self):
        """return True if self.bra_obj is not None"""

        return self.bra_obj is not None

    #
    def run(self):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used."""

        timing.start('transition.run')

        mol_bra = self.bra_obj.mol
        mol_ket = self.ket_obj.mol

        scf_bra = self.bra_obj.scf
        scf_ket = self.ket_obj.scf
        
        # sanity check that orbitals and geometry are the same
        if np.any(scf_bra.orbs != scf_ket.orbs):
            sys.exit('transition moments require same bra/ket orbs')

        if mol_bra.pymol().atom != mol_ket.pymol().atom:
            sys.exit('transition moments require same geometry/basis')

        # initialize the bitsi library
        libs.init_bitsi(self)

        # construct the trans_list array
        # currently store them as [initial state, final state]
        self.build_trans_list()

        # grab the transition density matrices
        tdm_list = mrci_1tdm.tdm(self.bra_obj, self.ket_obj, 
                                 self.trans_list)
        self.build_tdms()

        # build the multipole moments  -- easier to just do this once
        # for all transitions
        self.build_multipole(mol_ket.pymol(), scf_ket)

        # compute oscillator strengths
        self.build_osc_str()

        # print natural differnence and natural transition orbitals 
        # for each pair
        self.build_natural_orbs()

        # print the transition and difference densities and print
        # excitations energies/osc. strengths to file
        self.print_orbitals()

        # print the summary output
        self.print_log()

        timing.stop('transition.run')

        return

    #
    def tdm1(self, b_st, k_st):
        """return the tdm for the bra and ket state indicated"""

        return self.tdms[:,:, b_st, k_st]
 
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
    def osc_strength_isotropic(self, b_st, k_st, 
                                            order=0, gauge='velocity'):
        """return the isotropically averaged oscillator strength in
           the specified gauge and through the specified order"""

        return self.extract_osc(b_st, k_st, order, gauge, True)

#--------------------------------------------------------------------------

    # 
    def build_tdms(self):
        """grab the TDMs from bitsi and then reshape the list of
           TDMs into a more usable format"""

        # grab the tdms
        tdm_list = mrci_1tdm.tdm(self.bra_obj, self.ket_obj, 
                                                     self.trans_list)

        # make the tdm list
        nbra = len(self.bra_list)
        nket = len(self.ket_list)
        nmo  = self.bra_obj.scf.nmo

        self.tdms = np.zeros((nmo, nmo, nbra, nket), dtype=float)

        for ik in range(nket):
            [kirr,kst] = self.ket_obj.state_sym(self.ket_list[ik])
            for ib in range(nbra):
                [birr,bst] = self.bra_obj.state_sym(self.bra_list[ib])
                indx       = self.trans_index(birr,kirr,bst,kst)
                self.tdms[:,:,ib,ik] = tdm_list[birr][kirr][:,:,indx]

        del(tdm_list)
        return

    # 
    def trans_index(self, b_irr, k_irr, b_st, k_st):
        """return the index in the trans_list corresponding to
           b_st,k_st"""

        # return the transition dipole vector given the irreps and
        # state indices
        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st, k_st])
        except:
            return None

        return ind

    #
    def extract_multipole(self, b_st, k_st, name, gauge='velocity'):
        """extract the transition dipole vector in the appropriate
           gauge"""

        # append a '_l' or '_v' to indicate length or velocity gauge
        keyname = name+'_'+str(gauge)[0]

        if keyname in self.multipole.keys():
            return self.multipole[keyname][..., b_st, k_st]
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

        if keyname in self.oscstr.keys():
            return self.oscstr[keyname][..., b_st, k_st]
        else:
            return None

    #
    def build_trans_list(self):
        """built an array of initial states, ordered by adiabatic energy
           the symmetries of each state will be stored in a separate 
           array"""

        # determine what initial states to include
        # -------------------------------------------------

        state_list  = []

        # iterate through the initial states, adding roots
        # to various irreps as needed -- these states are 
        # ordered by adiabatic energy
        ninit_total = sum(self.ket_obj.n_states())
        if self.init_states is not None:
            # iterate through the initial states, adding roots
            # to various irreps as nee
            for istate in self.init_states:
                if istate < ninit_total:
                    state_list.append(istate)
        
        # use the istate array (format = irrep.state) to select
        # individual states
        if self.init_states_sym is not None:
            for state in self.init_states_sym:
                irrep   = int(state)
                st      = int(10*state) - 10*irrep
                istate  = self.ket_obj.state_index(irrep,st)
                if istate is not None:
                    state_list.append(istate)

        self.ket_list = sorted(list(set(state_list)))
        self.ket_sym  = [self.ket_obj.state_sym(self.ket_list[i])[0] 
                                    for i in range(len(self.ket_list))]

        if len(self.ket_list) == 0:
            sys.exit('no initial states selected in Transition')

        # determine what final states to include
        # -------------------------------------------------

        state_list  = []

        # iterate through the final states, adding roots
        # to various irreps as needed -- these states are 
        # ordered by adiabatic energy
        nfinal_total = sum(self.bra_obj.n_states())
        if self.final_states is not None:
            # iterate through the final states, adding roots
            # to various irreps as nee
            for fstate in self.final_states:
                if fstate < nfinal_total:
                    state_list.append(fstate)

        # use the istate array (format = irrep.state) to select
        # individual states
        if self.final_states_sym is not None:
            for state in self.final_states_sym:
                irrep   = int(state)
                st      = int(10*state) - 10*irrep
                fstate  = self.bra_obj.state_index(irrep,st)
                if fstate is not None:
                    state_list.append(fstate)

        self.bra_list = sorted(list(set(state_list)))
        self.bra_sym  = [self.bra_obj.state_sym(self.bra_list[i])[0]
                                    for i in range(len(self.bra_list))]

        # create the pair list -- we build the list this way 
        # in order to facilitate call to the tdm code in bitci
        # which computes tmds for pairs of states in particular
        # irreps.

        nirr_bra = self.bra_obj.n_irrep()
        nirr_ket = self.ket_obj.n_irrep()

        self.trans_list = [[[] for irr_bra in range(nirr_bra)]
                               for irr_ket in range(nirr_ket)]

        for k_ind in range(len(self.ket_list)):
            kst  = self.ket_obj.state_sym(self.ket_list[k_ind])[1]
            ksym = self.ket_sym[k_ind]

            for b_ind in range(len(self.bra_list)):
                bst  = self.bra_obj.state_sym(self.bra_list[b_ind])[1]
                bsym = self.bra_sym[b_ind]

                self.trans_list[bsym][ksym].append([bst, kst])

        return

    #
    def build_multipole(self, mol, scf):
        """builds the multipole moment transitions tensors"""

        # build the electronic tensors: dipole, quadrupole, octupole
        self.build_dip(mol, scf)
        self.build_quad(mol, scf)
        self.build_oct(mol, scf)

        # build the magnetic moment tensors: dipole, quadrupole
        self.build_mdip(mol, scf)
        self.build_mquad(mol, scf)

        return

    # 
    def build_osc_str(self):
        """build the oscillator strengths in the length and velocity
           gauges"""

        self.build_osc_str_l()
        self.build_osc_str_v()

        return


    #
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

        self.multipole['e_dipole_l'] = self.contract_transitions(mu_mo)
        self.multipole['e_dipole_v'] = self.contract_transitions(p_mo)
    
        return

    # 
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

        self.multipole['e_quadrupole_l'] = self.contract_transitions(Q_mo)
        self.multipole['e_quadrupole_v'] = self.contract_transitions(Qp_mo)
    
        return

    #
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

        self.multipole['e_octupole_l'] = self.contract_transitions(O_mo)
        self.multipole['e_octupole_v'] = self.contract_transitions(Op_mo)

        return

    #
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

        self.multipole['m_dipole_v'] = self.contract_transitions(mdp_mo)

        return

    # 
    def build_mquad(self, mol, scf):
        """builds the electric quadrupole moments in length and velocity
           gauges"""

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

        self.multipole['m_quadrupole_v'] = self.contract_transitions(mquad_mo)

        return

    #
    def build_osc_str_l(self):
        """builds the oscillator strength in the length gauge
           see J. Chem. Phys. 137, 204106 (2012) for details. This
           is _actually_ in a mixed gauge, as the electric multipoles
           are in the length gauge and magnetic multipoles are in 
           the velocity gauge. I think this is actually standard,
           but need to confirm"""

        nbra = len(self.bra_list)
        nket = len(self.ket_list)

        f0_l   = np.zeros((3, nbra, nket), dtype=float)
        f0_iso = np.zeros((nbra, nket), dtype=float)
        f2_l   = np.zeros((3, nbra, nket), dtype=float)
        f2_iso = np.zeros((nbra, nket), dtype=float)
        
        for ib in range(nbra):
            for ik in range(nket):

                # set up state and energy info
                b_ener = self.bra_obj.energy(self.bra_list[ib])
                k_ener = self.ket_obj.energy(self.ket_list[ik])
                alpha  = constants.fine_str
                de     = b_ener - k_ener
                de2    = de**2
                de3    = de**3

                # electric tensors
                mu     = self.multipole['e_dipole_l'][:, ib, ik]
                Qab    = self.multipole['e_quadrupole_l'][:, :, ib, ik]
                Oabc   = self.multipole['e_octupole_l'][:, :, :, ib, ik]
                # magnetic tensors
                mdp    = self.multipole['m_dipole_v'][:, ib, ik]
                Mqp    = self.multipole['m_quadrupole_v'][:, :, ib, ik]

                # dipole oscillator strength - length gauge
                f0_l[:, ib, ik] = (2.*de)*mu**2
                # isotropically averaged value
                f0_iso[ib, ik]  = np.sum(f0_l[:, ib, ik])/3.

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
                f2_l[:, ib, ik] = np.zeros((3,), dtype=float)
                # isotropically averaged value
                f2_av = f0_iso[ib, ik] + \
                        fQ_iso + fm_iso + fuO_iso + fuM_iso
                f2_iso[ib, ik] = f2_av
                    
        self.oscstr['f0_l']    = f0_l
        self.oscstr['f0iso_l'] = f0_iso

        self.oscstr['f2_l']    = f2_l
        self.oscstr['f2iso_l'] = f2_iso

        return

    #
    def build_osc_str_v(self):
        """builds the oscillator strength in the velocity gauge
           see J. Chem. Phys. 143, 234103 (2015) for details."""

        nbra = len(self.bra_list)
        nket = len(self.ket_list)

        f0_v   = np.zeros((3, nbra, nket), dtype=float)
        f0_iso = np.zeros((nbra, nket), dtype=float)
        f2_v   = np.zeros((3, nbra, nket), dtype=float)
        f2_iso = np.zeros((nbra, nket), dtype=float)

        for ib in range(nbra):
            for ik in range(nket):

                # set up state and energy info
                b_ener = self.bra_obj.energy(self.bra_list[ib])
                k_ener = self.ket_obj.energy(self.ket_list[ik])
                alpha  = constants.fine_str
                de     = b_ener - k_ener
                de2    = de**2
                de3    = de**3

                # electric tensors
                mu     = self.multipole['e_dipole_v'][:, ib, ik]
                Qab    = self.multipole['e_quadrupole_v'][:, :, ib, ik]
                Oabc   = self.multipole['e_octupole_v'][:, :, :, ib, ik]
                # magnetic tensors
                mdp    = self.multipole['m_dipole_v'][:, ib, ik]
                Mqp    = self.multipole['m_quadrupole_v'][:, :, ib, ik]

                # dipole oscillator strength - length gauge
                f0_v[:, ib, ik] = (2./de)*mu**2
                # isotropically averaged value
                f0_iso[ib, ik]  = np.sum(f0_v[:, ib, ik])/3.

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
                f2_v[:, ib, ik] = np.zeros(3, dtype=float)
                # isotropically averaged value
                f2_av = f0_iso[ib, ik] + \
                        fQ_iso + fm_iso + fuO_iso + fuM_iso
                f2_iso[ib, ik] = f2_av

        self.oscstr['f0_v']    = f0_v
        self.oscstr['f0iso_v'] = f0_iso

        self.oscstr['f2_v']    = f2_v
        self.oscstr['f2iso_v'] = f2_iso

        return

    #
    def build_natural_orbs(self):
        """build the natural difference orbitals and natural transition
           orbitals"""

        nbra = len(self.bra_list)
        nket = len(self.ket_list)
        nmo  = self.bra_obj.scf.nmo
        nao  = self.bra_obj.mol.nao
        mos  = self.bra_obj.scf.orbs

        nto    = np.zeros((nao, nmo, nbra, nket), dtype=float)
        ndo    = np.zeros((nao, nmo, nbra, nket), dtype=float)
        nto_wt = np.zeros((nmo, nbra, nket), dtype=float)
        ndo_wt = np.zeros((nmo, nbra, nket), dtype=float)

        for ik in range(nket):
            for ib in range(nbra):

                tdm = self.tdm1(ib, ik)
                rdm_bra = self.bra_obj.rdm1(self.bra_list[ib])
                rdm_ket = self.ket_obj.rdm1(self.ket_list[ik])

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

                nto[:, :2*npairs, ib, ik] = nto_pairs
                nto_wt[:2*npairs, ib, ik] = wt_pairs

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

                ndo[:, :2*npairs, ib, ik] = ndo_pairs
                ndo_wt[:2*npairs, ib, ik] = wt_pairs

        self.nos['ndo']    = ndo
        self.nos['ndo_wt'] = ndo_wt
        self.nos['nto']    = nto
        self.nos['nto_wt'] = nto_wt

        return

    #
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

        return np.einsum('...pq,pi,qj->...ij', ao_ints, mos, mos)

    #
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
    def contract_transitions(self, mo_ints):
        """Loop over all transitions and contract mo_ints
           with the appropriate tdm"""

        # get the dimsensions of each rank 
        rank_dim = mo_ints.shape[:-2]
        nbra     = len(self.bra_list)
        nket     = len(self.ket_list)
        tens     = np.zeros((rank_dim + (nbra, nket)), dtype=float)

        for ib in range(nbra):
            for ik in range(nket):
                tdm = self.tdms[:, :, ib, ik]
                tens[..., ib, ik] = self.contract_tdm(mo_ints, tdm)

        return tens

    #
    def print_log(self):
        """print summary output to log file"""
        # for each initial state, write out a table of 
        # oscillator strenghts and transition dipole moments

        nbra = len(self.bra_list)
        nket = len(self.ket_list)

        k_irrlbl = self.ket_obj.mol.irreplbl
        b_irrlbl = self.bra_obj.mol.irreplbl

        # print a 'transition table' for each initial state
        for ik in range(nket):

            init_st   = self.ket_list[ik]
            init_sym  = k_irrlbl[self.ket_sym[ik]]
            final_st  = self.bra_list
            final_sym = [b_irrlbl[self.bra_sym[ib]] 
                         for ib in range(nbra)]
            exc_ener  = [self.bra_obj.energy(self.bra_list[i]) - 
                        self.ket_obj.energy(self.ket_list[ik]) 
                         for i in range(nbra)]

            # print a 'transition table' for each initial state
            output.print_transition_table(init_st,
                                          init_sym,
                                          final_st,
                                          final_sym,
                                          exc_ener, 
                                          self.oscstr['f0iso_l'][:,ik],
                                          self.oscstr['f2iso_l'][:,ik],
                                          self.oscstr['f0iso_v'][:,ik],
                                          self.oscstr['f2iso_v'][:,ik],
                                          self.oscstr['f0_v'][:,:,ik])

        return

    #
    def print_orbitals(self):
        """print the natural transition orbitals and difference
           densities to file"""

        nbra = len(self.bra_list)
        nket = len(self.ket_list)

        k_irrlbl = self.ket_obj.mol.irreplbl
        b_irrlbl = self.bra_obj.mol.irreplbl

        for ik in range(nket):
            kst  = self.ket_list[ik]
            ksym = k_irrlbl[self.ket_sym[ik]]

            for ib in range(nbra):
                bst  = self.bra_list[ib]
                bsym = b_irrlbl[self.bra_sym[ib]]

                str_suffix ='.'+ str(kst+1)+'_' +str(ksym.lower()) + \
                            '.'+ str(bst+1)+'_' +str(bsym.lower())
           
                # print out the NTOs
                # determine number of pairs
                wts = self.nos['nto_wt'][:, ib, ik]
                ncols = sum(chk > 1.e-16 for chk in np.absolute(wts))
                output.print_nos_molden(
                    'nto'+str_suffix, 
                    self.bra_obj.mol,
                    self.nos['nto'][:,:ncols, ib, ik], 
                    wts[:ncols])

                # print out the NDOs
                # determine number of pairs
                wts = self.nos['ndo_wt'][:, ib, ik]
                ncols = sum(chk > 1.e-16 for chk in np.absolute(wts))
                output.print_nos_molden(
                    'ndo'+str_suffix,
                    self.bra_obj.mol,
                    self.nos['ndo'][:,:ncols, ib, ik],
                    wts[:ncols])

        return
