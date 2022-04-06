"""
module to compute spin-orbit coupling matrix 
elements
"""

import sys as sys
import numpy as np
import copy as copy
import graci.interaction.interaction as interaction
import graci.utils.timing as timing
import graci.bitcitools.bitsi_init as bitsi_init
import graci.bitcitools.mrci_soc as mrci_soc
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
        self.label         = 'spinorbit'
        self.print_thresh  = 1.
        self.mf2e          = 'atomic'
        # group of objects to couple
        self.couple_groups = [None]
        # states from each object to couple
        self.couple_states = [None]

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

    #
    @timing.timed
    def run(self, obj_list):
        """
        returns the SOC matrix elements between the bra and ket
        states. If b_state and k_state are None, assume SOC matrix
        elements from all states in method object should be used
        """

        # section header
        output.print_spinorbit_header(self.label)

        # check on the MF 2e integral scheme
        if self.mf2e not in self.allowed_mf2e:
            print('\n Unrecognised mf2e value: '+self.mf2e
                  +'\n Allowed values: '+str(self.allowed_mf2e))
            sys.exit()

        # sanity check that all objects in list are compatible with
        # spin-orbit computation
        self.check_list(obj_list)

        # add each of the groups of states
        # this is kind of a hack -- storing numpy arrays of strings
        # introduces encoding issues
        if isinstance(self.couple_groups, np.ndarray):
            self.couple_groups = self.couple_groups.tolist()
        obj_lbls = [obj.label for obj in obj_list]
        self.grp_lbls = ['grp'+str(i) for i in 
                                     range(len(self.couple_groups))]

        for igrp in range(len(self.couple_groups)):
            # this allows the same group label to appear more than once
            indx     = obj_lbls.index(self.couple_groups[igrp])
            states   = self.couple_states[igrp].tolist()
            self.add_group(self.grp_lbls[igrp], obj_list[indx], states)

        # get the one-electron SOC matrices
        h1e = self.build_h1e()
        
        # initialise the H_soc array
        hdim = 0
        for grp in self.grp_lbls:
            hdim += len(self.get_states(grp)) * self.get_spins(grp).mult
        hsoc = np.zeros((hdim, hdim), dtype=np.cdouble)

        # construct the trans_list arrays for each of the state
        # pair blocks of H_SOC
        # loop over the blocks of H_SOC
        for iblk in range(self.n_groups()):
            # loop over range(iblk+1) to include iblk==jblk blk
            for jblk in range(iblk+1):

                bra_lbl  = self.grp_lbls[iblk]
                ket_lbl  = self.grp_lbls[jblk]

                bra = self.get_obj(bra_lbl)
                ket = self.get_obj(ket_lbl)

                # skip if both manifolds of states are singlets
                if bra.mult == 1 and ket.mult == 1:
                    continue

                # do some sanity checking on the bra/ket pair
                self.check_pair(bra_lbl, ket_lbl)

                if bra_lbl == ket_lbl:
                    blk_fill = 'lower'
                else:
                    blk_fill = 'full'

                blk_list     = self.build_pair_list(bra_lbl, ket_lbl, 
                                                  pairs=blk_fill)
                blk_list_sym = self.build_pair_list(bra_lbl, ket_lbl, 
                                                  pairs=blk_fill, 
                                                  sym_blk=True)

                # initialize the bitsi library for the calculation of
                # spin-orbit matrix elements
                bitsi_init.init(bra, ket, 'soc')

                # get the reduced matrix elements
                # <S' M'=S' psi_m||T_ij^(1)||S M=S psi'_n>
                # for the current block
                redmat = self.build_redmat(bra, ket, blk_list, 
                                                     blk_list_sym)

                # build this block of H_SOC
                hsoc += self.build_hsoc(hdim, bra_lbl, ket_lbl, 
                                        blk_list, h1e, redmat)
            
                # finalize the bitsi library
                bitsi_init.finalize()

        # add in the on-diagonal elements
        hsoc += np.diag(self.hsoc_ondiag(hdim))
        
        # diagonalise H_soc
        self.so_ener, self.so_vec = np.linalg.eigh(hsoc)
        
        # print the H_SOC elements
        if self.print_thresh >= 1:
            self.print_hsoc(hsoc)

        del(redmat)
        return
   
    #
    def soc_state(self, state):
        """
        Return the spin orbit state vector in terms of spin free states
        """

        if state < self.so_vec.shape[1]:
            return self.so_vec[:, state]
        else:
            return None

    #
    def soc_basis_lbl(self, indx):
        """for a given index of the state vector, return the label of 
           of group, the state label, and value of Ms of the basis
           state
        """
 
        if indx < 0 or indx >= self.so_vec.shape[0]:
            return None

        msval = [-1, 0, 1]
        ind   = 0
        igrp  = 0
        while ind < indx:
            nstates = len(self.get_states(self.grp_lbls[igrp]))
            mult    = self.get_spins(self.grp_lbls[i]).mult
            ind    += nstates*mult
            igrp   += 1

        # decrement igrp to current grp
        igrp -= 1
        lbl   = self.grp_lbls[igrp]

        # decrement indx to start of the grp
        ind  -= nstates*mult

        # find the state index
        st    = np.floor((indx - ind) / self.get_spins(lbl).mult)

        # get Ms value (-1, 0, 1)
        ms    = indx - (ind + st*self.get_spins(lbl).mult)

        return lbl, self.get_states(lbl)[st], msval[ms]

    #
    def soc_index(self, lbl, state, m_s):
        """return the indx in the Hsoc for a given state group, state,
           and value of Ms
        """

        spin = self.get_spins(lbl)
        ist  = self.get_states(lbl).index(state)

        indx = spin.mult * ist + m_s + spin.S

        # now determine the offset for bra/ket indices
        off = 0
        for i in range(self.grp_lbls.index(lbl)):
            nstates  = len(self.get_states(self.grp_lbls[i]))
            mult     = self.get_spins(self.grp_lbls[i]).mult
            off += nstates*mult

        return int(indx + off)

    # 
    def n_states(self):
        """return the number of spin-orbit coupled states"""
 
        if self.so_ener is None:
            return 0
        else:
            return self.so_ener.shape[0]

    #
    def energy(self, state):
        """
        Return the spin-orbit state energy
        """

        if self.so_ener is None or state > self.n_states():
            return None
        else:
            return self.so_ener[state]

    #-----------------------------------------------------------------
    # 
    # methods not meant to be called from outside the class

    #
    def check_list(self, obj_lst):
        """
        Perform some sanity checks on all objs in list
        """
        for i in range(len(obj_lst)):
            for j in range(i):
                if not self.same_obj(obj_lst[i], obj_lst[j]):
                    # sanity check that orbitals and geometry are the same
                    if np.any(obj_lst[i].scf.orbs != obj_lst[j].scf.orbs):
                        sys.exit('spin-orbit coupling requires same'+
                                 'bra/ket orbs')

                    if (obj_lst[i].scf.mol.pymol().atom !=
                                     obj_lst[j].scf.mol.pymol().atom):
                        sys.exit('spin-orbit coupling requires same '+
                                 'geometry and basis set')
        return

    #
    def check_pair(self, bra_lbl, ket_lbl):
        """
        Perform some sanity checks on the bra/ket states and objects
        """

        bra_spin = self.get_spins(bra_lbl)
        ket_spin = self.get_spins(ket_lbl)

        # Delta S = -1, 0 or +1 must hold
        if bra_spin.S - ket_spin.S not in [-1., 0., 1.]:
            sys.exit('\n ERROR: S_bra, S_ket spin combination not ' \
                     +'currently supported in spinorbit')

        # In the case of Delta S = 0, S > 0 must hold
        if bra_spin.S == ket_spin.S and bra_spin.S == 0.:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination ' \
                     +'in spinorbit')
        return

    #
    @timing.timed
    def build_h1e(self):
        """
        Sets up the one-electron SOC matrices h^(k), k=-1,0,+1
        """

        # PySCF Mole obeject
        mol = self.get_obj(self.grp_lbls[0]).scf.mol.pymol().copy()
        # No. AOs
        nao = mol.nao_nr()

        # get the MF density matrix
        rho_ao = self.build_rho()

        # Initialise arrays
        h1e      = np.zeros((3, nao, nao), dtype=np.cdouble)
        hcart_ao = np.zeros((3, nao, nao), dtype=float)

        # One-electron contributions
        for iatm in range(mol.natm):
            Z   = mol.atom_charge(iatm)
            xyz = mol.atom_coord(iatm)
            mol.set_rinv_orig(xyz)
            hcart_ao += Z * mol.intor('cint1e_prinvxp_sph', 3)
        hcart_ao = hcart_ao * constants.fine_str**2 / 2

        # Mean-field two-electron contributions
        if self.mf2e == 'full':
            h2e_ao    = mol.intor("cint2e_p1vxp1_sph", comp=3, aosym="s1")
            h2e_ao   *= -0.5 * constants.fine_str**2
            hcart_ao += (np.einsum("ijklm,lm->ijk", h2e_ao, rho_ao)
                         - 1.5 * (np.einsum("ijklm, kl->ijm", h2e_ao, rho_ao)
                         + np.einsum("ijklm,mj->ilk", h2e_ao, rho_ao)))
        elif self.mf2e == 'atomic':
            hcart_ao += self.build_mf_atomic(mol, rho_ao)

        # Transform to the MO basis
        orbs     = self.get_obj(self.grp_lbls[0]).scf.orbs
        hcart_mo = np.matmul(np.matmul(orbs.T, hcart_ao), orbs)

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
    def build_rho(self):
        """
        sets up the mean-field density matrix
        """
        nsta = sum([len(self.get_states(lbl)) for lbl in self.grp_lbls])
        nmo  = self.get_obj(self.grp_lbls[0]).scf.nmo

        # Average density matrix across states
        rho_mo = np.zeros((nmo, nmo), dtype=float)

        for grp in self.grp_lbls:
            for st in self.get_states(grp):
                rho_mo += self.get_obj(grp).rdm(st)

        rho_mo /= nsta

        # Transform to the AO basis
        orbs   = self.get_obj(self.grp_lbls[0]).scf.orbs
        rho_ao = np.matmul(np.matmul(orbs, rho_mo), orbs.T)

        return rho_ao

    #
    @timing.timed
    def build_mf_atomic(self, mol, rho_ao):
        """
        builds the atomic one-centre approximation to the
        mean-field two-electron SOC integrals
        """

        # function result
        nao    = mol.nao_nr()
        h2e_ao = np.zeros((3, nao, nao), dtype=float)

        # atom indices for each shell
        shell_atm = np.array([shell[0] for shell in mol._bas])

        for iatm in range(mol.natm):

            # shell indices for this atom
            shells     = mol.atom_shell_ids(iatm)
            shls_slice = [shells[0], shells[-1]+1] * 4

            # fetch the integrals for this atomic centre
            ints = mol.intor("cint2e_p1vxp1_sph", shls_slice=shls_slice,
                             comp=3, aosym="s1")
            ints *= -0.5 * constants.fine_str**2

            # density matrix for this block of AOs
            ji, jf = mol.nao_nr_range(shells[0], shells[-1]+1)
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
    def build_redmat(self, bra, ket, pair_list, pair_list_sym):
        """
        grab the reduced matrix elements from bitsi and then
        reshape the list into a more usable format
        """

        # grab the reduced matrix elements
        redmat_list = mrci_soc.redmat(bra, ket, pair_list_sym)

        # Make the reduced matrix list
        nmo         = bra.scf.nmo
        npairs      = len(pair_list)
        redmat_blk  = np.zeros((nmo, nmo, npairs), dtype=float)

        for bkst in pair_list:
            [birr, bst] = bra.state_sym(bkst[0])
            [kirr, kst] = ket.state_sym(bkst[1])
            indx        = pair_list.index(bkst)
            indx_sym    = pair_list_sym[birr][kirr].index([bst, kst])
            redmat_blk[:,:,indx] = \
                        redmat_list[birr][kirr][:,:,indx_sym].copy()

        del(redmat_list)

        return redmat_blk

    #
    @timing.timed
    def build_hsoc(self, dim, bra_lbl, ket_lbl, pair_list, h1e, redmat):
        """
        Builds the iblock'th block of the SOC Hamiltonian matrix
        """

        # fill in the appropriate elements of hsoc
        hsoc = np.zeros((dim, dim), dtype=np.cdouble)

        # get the bra/ket method objects
        bra      = self.get_obj(bra_lbl)
        ket      = self.get_obj(ket_lbl)

        # get the corresponding bra/ket spin objects
        bra_spin = self.get_spins(bra_lbl)
        ket_spin = self.get_spins(ket_lbl)

        # Number of bra and ket multiplets in this block
        bra_states = self.get_states(bra_lbl)
        ket_states = self.get_states(ket_lbl)

        # get the Clebsh-Gordan coefficients needed to
        # compute the SOC matrix elements for the current block
        cg_coef = mrci_soc.clebsch_gordan(bra, ket)

        # Loop over pairs of spin-orbit coupled states in this block
        Mb = bra_spin.M
        Mk = ket_spin.M

        for pair in pair_list:
            indx   = pair_list.index(pair)
            bra_st = pair[0]
            ket_st = pair[1]

            for M_ket in Mk:
                for M_bra in Mb:
                    i         = self.soc_index(bra_lbl, bra_st, M_bra)
                    j         = self.soc_index(ket_lbl, ket_st, M_ket)
                    hsoc[i,j] = self.contract_redmat(
                                          h1e,
                                          bra_spin, ket_spin,
                                          M_bra, M_ket,
                                          cg_coef, redmat[:,:,indx])
                    hsoc[j,i] = np.conj(hsoc[i,j])

        return hsoc

    
    #
    def hsoc_ondiag(self, dim):
        """
        Adds the on-diagonal elements of H_SOC
        """

        hdiag = np.zeros(dim, dtype=np.cdouble)

        indx = 0
        for igrp in self.grp_lbls:
            obj = self.get_obj(igrp)
            for ist in self.get_states(igrp):
                for m in self.get_spins(igrp).M:
                    hdiag[indx] = obj.energy(ist)
                    indx += 1

        return hdiag
        
    #
    @timing.timed
    def contract_redmat(self, h1e, bra_spin, ket_spin, M_bra, M_ket, 
                        cg_coef, redmat):
        """
        Computes a single element 
        < I_bra M_bra | H_SOC | I_ket M_ket >
        of the SOC Hamiltonian matrix via the contraction
        of the reduced matrices with the Clebsch-Gordan coefficient-
        scaled one-electron SOC matrices h^(k), k=-1,0,+1
        """
        
        # Sum of the scaled one-electron SOC matrices
        nmo    = self.get_obj(self.grp_lbls[0]).scf.nmo
        hscale = np.zeros((nmo, nmo), dtype=np.cdouble)
        kval   = [-1, 0, 1]
        coe    = [1., np.sqrt(2.), -1.]
        for n in range(3):
            i, i12  = mrci_soc.clebsch_gordan_index(
                                      bra_spin.S, M_bra, 
                                      ket_spin.S, M_ket, kval[n])
            hscale += h1e[2-n, :, :] * cg_coef[i12, i] * coe[n]
        
        # Contraction with the reduced matrix
        hij = 0.5 * np.einsum('ij,ij', redmat, hscale)

        return hij

    #
    def create_stlbl(self):
        """
        list of all all multiplet indices and projected spin values
        """

        stlbl = []

        for igrp in self.grp_lbls:
            obj_lbl = self.get_obj(igrp).label
            for i in self.get_states(igrp):
                for m in self.get_spins(igrp).M:
                    stlbl.append([obj_lbl, igrp, i, m])

        return stlbl

    #
    def print_hsoc(self, hsoc):
        """
        prints the SOC values to the log file
        """

        stlbl = self.create_stlbl()

        # output the SOC matrix elements
        output.print_spinorbit_table(hsoc, hsoc.shape[0], stlbl,
                                     self.print_thresh)

        # output the eigenvectors of H_SOC
        output.print_hsoc_eig(self.so_ener, self.so_vec, hsoc.shape[0], 
                                     stlbl)
        
        return
