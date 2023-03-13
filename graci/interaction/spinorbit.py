"""
module to compute spin-orbit coupling matrix 
elements
"""

import sys as sys
import numpy as np
import copy as copy
import graci.interaction.interaction as interaction
import graci.utils.timing as timing
import graci.interfaces.bitci.bitsi_init as bitsi_init
import graci.interfaces.bitci.mrci_soc as mrci_soc
import graci.io.output as output
import graci.utils.constants as constants

class Spinorbit(interaction.Interaction):
    """
    Spin orbit coupling class. For now will only accept calculations
    for which the 'Molecule and 'Scf' objects are the same
    """
    def __init__(self):
        # parent attributes
        super().__init__()
        
        # user defined quanties 
        #-----------------------------------------------------------
        self.print_soc_thresh = 1.
        self.mf2e             = 'atomic'
        # group of objects to couple
        self.couple_groups    = [None]
        # states from each object to couple
        self.couple_states    = [None]
        # repsentation adiabatic vs diabatic
        self.representation   = 'adiabatic'

        #
        # Private variables -- should not be directly referenced
        #                      outside the class
        #------------------------------------------------------------
        # global variables
        # store the group labels globally, because the order 
        # matters for construction of Hsoc
        self.grp_lbls     = None
        # allowed MF 2e integral schemes
        self.allowed_mf2e = ['off', 'atomic', 'full']
        # matrix to hold spin-orbit vectors
        self.so_vec       = None
        # vector to hold spin-orbit energies
        self.so_ener      = None
        # 1-RDMs for the spin-orbit coupled states
        self.dmats        = {'adiabatic' : None, 'diabatic' : None}
        # allowed representations
        self.allowed_reps = ['adiabatic'] 
        # common mol object
        self.mol          = None      
        # common mos
        self.mos          = None
 
    def copy(self):
        """create of deepcopy of self"""
        new = self.Spinorbit()

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
    def run(self, ci_objs):
        """
        returns the SOC matrix elements between the bra and ket
        states. If b_state and k_state are None, assume SOC matrix
        elements from all states in method object should be used
        """

        # section header
        if self.verbose:
            output.print_spinorbit_header(self.label)

        # sanity check on the representation
        self.check_representation()

        # sanity check that all objects in list are compatible with
        # spin-orbit computation
        self.mol, self.mos = self.check_soc_objs(ci_objs)

        # set the bra/ket objects and add the state groups associated
        # with each 
        nsoc = self.get_nsoc_states(ci_objs, self.couple_states)
        hsoc = np.zeros((nsoc, nsoc), dtype=np.cdouble)

        self.add_group('spinorbit', ci_objs,
                                    states = self.couple_states, 
                                    spin = True)

        # get the one-electron SOC matrices
        h1e = self.build_h1e()

        # loop over the ci groups in the spin orbit object
        nci_objs = self.n_ci_objs('spinorbit')
        for b in range(nci_objs):
            b_lbl = self.get_ci_lbls('spinorbit')[b]
            b_obj = self.get_ci_obj('spinorbit', b_lbl)

            # need to include the ci_grp1 == ci_grp2 (m1 == m2) block
            for k in range(b+1):
                k_lbl = self.get_ci_lbls('spinorbit')[k]
                k_obj = self.get_ci_obj('spinorbit', k_lbl)

                b_mult, b_s, b_m = self.get_ci_spins('spinorbit', b_lbl)
                k_mult, k_s, k_m = self.get_ci_spins('spinorbit', k_lbl)

                # no singlet/singlet coupling
                if b_mult == 1 and k_mult == 1:
                    continue               
 
                ptype = 'full'
                if b_lbl == k_lbl:
                    ptype = 'lower'
                ci_pairs     = self.ci_pair_list('spinorbit', b_lbl,
                                                 'spinorbit', k_lbl,
                                                  pairs=ptype)
                ci_pairs_sym = self.ci_pair_list('spinorbit', b_lbl,
                                                 'spinorbit', k_lbl,
                                                  pairs=ptype, 
                                                  sym_blk=True)

                # initialize the bitsi library for the calculation of
                # spin-orbit matrix elements
                bitsi_init.init(b_obj, k_obj, 'soc', self.verbose)

                # get the reduced matrix elements
                # <S' M'=S' psi_m||T_ij^(1)||S M=S psi'_n>
                # for the current block
                redmat = self.build_redmat(b_lbl, k_lbl, ci_pairs, 
                                                         ci_pairs_sym)

                # build this block of H_SOC
                hsoc += self.build_hsoc(nsoc, b_lbl, k_lbl, ci_pairs,
                                                          h1e, redmat)

                # finalize the bitsi library
                bitsi_init.finalize()

        # add in the on-diagonal elements
        hsoc += np.diag(self.hsoc_ondiag(nsoc))

        # diagonalise H_soc
        so_ener, so_vec = np.linalg.eigh(hsoc)
        self.set_group_states('spinorbit', so_vec, so_ener)

        # build the one-particle RDMS
        self.build_rdms()

        # print the H_SOC elements
        if self.print_soc_thresh >= 1:
            self.print_hsoc(hsoc)

        del(redmat)
        return

    #-----------------------------------------------------------------------
    #
   
    # 
    def n_states(self):
        """return the number of spin-orbit coupled states"""
 
        return len(self.get_states('spinorbit'))

    #
    def energy(self, state):
        """
        Return the spin-orbit state energy
        """

        if state > self.n_states():
            return None
        else:
            return self.get_energy('spinorbit', state)

    #
    def rdm(self, istate, rep='adiabatic'):
        """return the density matrix for the state istate"""
    
        if self.dmats[rep] is not None \
           and istate < self.n_states():
            return self.dmats[rep][istate, :, :]
        else:
            print("rdm called but density does not exist")
            return None

    #-----------------------------------------------------------------
    # 
    # methods not meant to be called from outside the class

    #
    def check_soc_objs(self, soc_objs):
        """
        Perform some sanity checks on all objs in list
        """
        for i in range(len(soc_objs)):
            for j in range(i+1):

                # sanity check that orbitals are the same
                if np.any(soc_objs[i].mos != soc_objs[j].mos):
                    sys.exit('spin-orbit coupling requires same'+
                             'bra/ket orbs')
                else:
                    mos = soc_objs[i].mos

                # sanity check that the geometry is the same
                if (soc_objs[i].scf.mol.pymol().atom !=
                                 soc_objs[j].scf.mol.pymol().atom):
                    sys.exit('spin-orbit coupling requires same '+
                             'geometry and basis set')
                else:
                    mol = soc_objs[i].scf.mol

                # sanity check that the spins are compatible
                S_i = 0.5 * (soc_objs[i].mult - 1.)
                S_j = 0.5 * (soc_objs[j].mult - 1.)

                # Delta S = -1, 0 or +1 must hold
                if S_i - S_j not in [-1., 0., 1.]:
                    sys.exit('\n ERROR: S_bra, S_ket spin combo ' \
                          ' not currently supported in spinorbit')

        # check on the MF 2e integral scheme
        if self.mf2e not in self.allowed_mf2e:
            print('\n Unrecognised mf2e value: '+self.mf2e
                  +'\n Allowed values: '+str(self.allowed_mf2e))
            sys.exit()

        return mol, mos

    #
    def get_nsoc_states(self, ci_objs, couple_states):
        """
        Determine the number of SOC states that will result
        """
        return sum([ci_objs[i].mult*len(couple_states[i])
                             for i in range(len(ci_objs))])

    #
    @timing.timed
    def build_h1e(self):
        """
        Sets up the one-electron SOC matrices h^(k), k=-1,0,+1
        """

        # PySCF Mole obeject
        pymol = self.mol.pymol()
        # No. AOs
        nao   = pymol.nao_nr()
        # No. MOs
        nmo   = self.mos.shape[1]

        # get the MF density matrix
        rho_ao = self.build_rho()

        # Initialise arrays
        h1e      = np.zeros((3, nmo, nmo), dtype=np.cdouble)
        hcart_ao = np.zeros((3, nao, nao), dtype=float)

        # One-electron contributions
        for iatm in range(pymol.natm):
            Z   = pymol.atom_charge(iatm)
            xyz = pymol.atom_coord(iatm)
            pymol.set_rinv_orig(xyz)
            hcart_ao += Z * pymol.intor('cint1e_prinvxp_sph', 3)
        hcart_ao = hcart_ao * constants.fine_str**2 / 2

        # Mean-field two-electron contributions
        if self.mf2e == 'full':
            h2e_ao    = pymol.intor("cint2e_p1vxp1_sph", comp=3, aosym="s1")
            h2e_ao   *= -0.5 * constants.fine_str**2
            hcart_ao += (np.einsum("ijklm,lm->ijk", h2e_ao, rho_ao)
                         - 1.5 * (np.einsum("ijklm, kl->ijm", h2e_ao, rho_ao)
                         + np.einsum("ijklm,mj->ilk", h2e_ao, rho_ao)))
        elif self.mf2e == 'atomic':
            hcart_ao += self.build_mf_atomic(pymol, rho_ao)

        # Transform to the MO basis
        orbs     = self.mos
        hcart_mo = orbs.T @ hcart_ao @ orbs

        # Transform to the spherical tensor representation
        ci = np.sqrt(-1+0j)
        # k = -1: x - iy
        h1e[0,:,:] = hcart_mo[0,:,:] - ci * hcart_mo[1,:,:]
        # k = +1: x + iy
        h1e[2,:,:] = hcart_mo[0,:,:] + ci * hcart_mo[1,:,:]
        # k = 0: z 
        h1e[1,:,:] = hcart_mo[2,:,:]

        return h1e

    #
    @timing.timed
    def build_redmat(self, b_lbl, k_lbl, pair_list, pair_list_sym):
        """
        grab the reduced matrix elements from bitsi and then
        reshape the list into a more usable format
        """
        bra = self.get_ci_obj('spinorbit', b_lbl)
        ket = self.get_ci_obj('spinorbit', k_lbl)

        # grab the reduced matrix elements
        redmats = mrci_soc.redmat(bra, ket, pair_list_sym)

        # Make the reduced matrix lisit
        nmo       = self.mos.shape[1]
        npairs    = len(pair_list)
        redmat_bk = np.zeros((nmo, nmo, npairs), dtype=float)

        for bkst in pair_list:
            [bir, bst] = self.get_ci_sym('spinorbit', b_lbl, bkst[0])
            [kir, kst] = self.get_ci_sym('spinorbit', k_lbl, bkst[1])
            indx       = pair_list.index(bkst)
            sindx      = pair_list_sym[bir][kir].index([bst, kst])
            redmat_bk[:,:,indx] = redmats[bir][kir][:,:,sindx].copy()

        del(redmats)

        return redmat_bk

    #
    def build_rdms(self):
        """
        construct the 1-RDMs 
        """
        g_name  = 'spinorbit'
        soc_wts = (np.conj(self.get_vectors(g_name)) * 
                           self.get_vectors(g_name)).real

        ci_lbls = self.get_ci_lbls(g_name)
        n_ci    = len(ci_lbls)
        nsf_st  = [len(self.get_ci_states(g_name, ci_lbls[i])) 
                                   for i in range(n_ci)]
       
        nst     = self.n_states()
        wts     = [np.zeros((nst, nsf_st[i]), dtype=float) 
                                   for i in range(n_ci)]

        for g_st in range(nst):
            ci_lbl, st, ci_m = self.get_ci_index(g_name, g_st)
            c_indx = ci_lbls.index(ci_lbl)
            wts[c_indx][:, st] += soc_wts[g_st, :]
          
        nmo = self.mos.shape[1] 
        self.dmats[self.representation] = np.zeros((nst, nmo, nmo), 
                                          dtype=float)

        for i in range(n_ci):
            states = self.get_ci_states(g_name, ci_lbls[i])
            ci_obj = self.get_ci_obj(g_name, ci_lbls[i])
            self.dmats[self.representation] += \
                  np.einsum('ij,jkl->ikl', wts[i], 
                           ci_obj.dmats[self.representation][states,:,:])

        return

    #
    @timing.timed
    def build_rho(self):
        """
        sets up the mean-field density matrix
        """
        nmo  = self.mos.shape[1]

        # Average density matrix across states
        rho_mo = np.zeros((nmo, nmo), dtype=float)

        nsta = 0
        for ci_lbl in self.get_ci_lbls('spinorbit'):
            for st in self.get_ci_states('spinorbit', ci_lbl):
                rho_mo += self.get_ci_obj('spinorbit', ci_lbl).rdm(st)
                nsta   += 1

        rho_mo /= nsta

        # Transform to the AO basis
        orbs   = self.mos
        rho_ao = np.matmul(np.matmul(orbs, rho_mo), orbs.T)

        return rho_ao

    #
    @timing.timed
    def build_mf_atomic(self, pymol, rho_ao):
        """
        builds the atomic one-centre approximation to the
        mean-field two-electron SOC integrals
        """

        # function result
        nao    = pymol.nao_nr()
        h2e_ao = np.zeros((3, nao, nao), dtype=float)

        # atom indices for each shell
        shell_atm = np.array([shell[0] for shell in pymol._bas])

        for iatm in range(pymol.natm):

            # shell indices for this atom
            shells     = pymol.atom_shell_ids(iatm)
            shls_slice = [shells[0], shells[-1]+1] * 4

            # fetch the integrals for this atomic centre
            ints = pymol.intor("cint2e_p1vxp1_sph", shls_slice=shls_slice,
                             comp=3, aosym="s1")
            ints *= -0.5 * constants.fine_str**2

            # density matrix for this block of AOs
            ji, jf = pymol.nao_nr_range(shells[0], shells[-1]+1)
            rho = rho_ao[ji:jf, ji:jf]

            # contraction with the density matrix
            h = (np.einsum("ijklm,lm->ijk", ints, rho)
                         - 1.5 * (np.einsum("ijklm, kl->ijm", ints, rho)
                         + np.einsum("ijklm,mj->ilk", ints, rho)))

            # mapping back up to the full set of AOs
            p = -1
            for p1 in range(ji, jf):
                p += 1
                q = -1
                for q1 in range(ji, jf):
                    q += 1
                    h2e_ao[:, p1, q1] = h[:, p, q]

        return h2e_ao

    #
    @timing.timed
    def build_hsoc(self, dim, b_lbl, k_lbl, pair_list, h1e, redmat):
        """
        Builds the iblock'th block of the SOC Hamiltonian matrix
        """

        # fill in the appropriate elements of hsoc
        hsoc = np.zeros((dim, dim), dtype=np.cdouble)

        # get the bra/ket method objects
        bra      = self.get_ci_obj('spinorbit', b_lbl)
        ket      = self.get_ci_obj('spinorbit', k_lbl)

        # get the corresponding bra/ket spin objects
        b_mult, b_s, b_M = self.get_ci_spins('spinorbit', b_lbl)
        k_mult, k_s, k_M = self.get_ci_spins('spinorbit', k_lbl)

        # Number of bra and ket multiplets in this block
        b_states = self.get_ci_states('spinorbit', b_lbl)
        k_states = self.get_ci_states('spinorbit', k_lbl)

        # get the Clebsh-Gordan coefficients needed to
        # compute the SOC matrix elements for the current block
        cg_coef = mrci_soc.clebsch_gordan(bra, ket)

        # Loop over pairs of spin-orbit coupled states in this block
        for pair in pair_list:
            indx   = pair_list.index(pair)
            bra_st = pair[0]
            ket_st = pair[1]

            for mk in k_M:
                for mb in b_M:
                    i = self.get_group_index('spinorbit', 
                                                b_lbl, bra_st, mb)
                    j = self.get_group_index('spinorbit', 
                                                k_lbl, ket_st, mk)
                    hsoc[i,j] = self.contract_redmat(h1e,
                                          b_s, k_s, mb, mk,
                                          cg_coef, redmat[:,:,indx])
                    hsoc[j,i] = np.conj(hsoc[i,j])

        return hsoc

    #
    def hsoc_ondiag(self, dim):
        """
        Adds the on-diagonal elements of H_SOC
        """

        hdiag = np.zeros(dim, dtype=np.cdouble)

        indx   = 0
        g_name = 'spinorbit'
        ci_lbls = self.get_ci_lbls(g_name)
        for ci_i in range(len(ci_lbls)):
            ci_lbl = ci_lbls[ci_i]
            ci_st  = self.get_ci_states(g_name, ci_lbl)
            for ist in range(len(ci_st)):
                mult, S, M = self.get_ci_spins(g_name, ci_lbl)
                for m_s in M:
                    hdiag[indx] = self.get_ci_energy(g_name, ci_lbl, 
                                                          ci_st[ist])
                    indx += 1

        return hdiag
        
    #
    @timing.timed
    def contract_redmat(self, h1e, S_bra, S_ket, M_bra, M_ket, 
                        cg_coef, redmat):
        """
        Computes a single element 
        < I_bra M_bra | H_SOC | I_ket M_ket >
        of the SOC Hamiltonian matrix via the contraction
        of the reduced matrices with the Clebsch-Gordan coefficient-
        scaled one-electron SOC matrices h^(k), k=-1,0,+1
        """
        
        # Sum of the scaled one-electron SOC matrices
        nmo    = self.mos.shape[1]
        hscale = np.zeros((nmo, nmo), dtype=np.cdouble)
        kval   = [-1, 0, 1]
        coe    = [1., np.sqrt(2.), -1.]
        for n in range(3):
            i, i12  = mrci_soc.clebsch_gordan_index(
                                      S_bra, M_bra, 
                                      S_ket, M_ket, kval[n])
            hscale += h1e[2-n, :, :] * cg_coef[i12, i] * coe[n]
        
        # Contraction with the reduced matrix
        hij = 0.5 * np.einsum('ij,ij', redmat, hscale)

        return hij

    #
    def create_stlbl(self):
        """
        list of all all multiplet indices and projected spin values
        """

        lbl      = []
        grp_name = 'spinorbit'
        ci_lbls  = self.get_ci_lbls(grp_name)

        for ci in range(len(ci_lbls)):
            ci_lbl     = ci_lbls[ci]
            ci_states  = self.get_ci_states(grp_name, ci_lbl)
            obj_lbl    = self.get_ci_obj(grp_name, ci_lbl).label
            mult, S, M = self.get_ci_spins(grp_name, ci_lbl)
            for ist in range(len(ci_states)):
                for iM in range(len(M)):
                    lbl.append([obj_lbl, ci_lbl, ci_states[ist], M[iM]])

        return lbl

    #
    def print_hsoc(self, hsoc):
        """
        prints the SOC values to the log file
        """

        stlbl = self.create_stlbl()

        # output the SOC matrix elements
        if self.verbose:
            output.print_spinorbit_table(hsoc, hsoc.shape[0], stlbl,
                                         self.print_soc_thresh)

        # output the eigenvectors of H_SOC
        if self.verbose:
            output.print_hsoc_eig(self.get_energy('spinorbit'), 
                                  self.get_vectors('spinorbit'),
                                  hsoc.shape[0], stlbl)
        
        return
