"""
module to compute spin-orbit coupling matrix 
elements
"""

import sys as sys
import numpy as np
import copy as copy
import graci.interaction.interaction as interaction
import graci.utils.timing as timing
import graci.core.libs as libs
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
        self.label        = 'spinorbit'
        self.print_thresh = 1.
        self.mf2e         = 'atomic'

        #
        # Private variables -- should not be directly referenced
        #                      outside the class
        #------------------------------------------------------------
        # global variables
        # One electron SOC matrices
        self.h1e          = None
        # List of multiplet indices and projected spins
        self.stlbl        = None
        # MF density matrix (AO basis)
        self.rho_ao       = None
        # allowed MF 2e integral schemes
        self.allowed_mf2e = ['off', 'atomic', 'full']
     
    #
    @timing.timed
    def run(self):
        """
        returns the SOC matrix elements between the bra and ket
        states. If b_state and k_state are None, assume SOC matrix
        elements from all states in method object should be used
        """

        # section header
        output.print_spinorbit_header(self.label)

        # bra/ket mol and scf objects
        mol_bra = self.bra_obj.scf.mol
        mol_ket = self.ket_obj.scf.mol
        scf_bra = self.bra_obj.scf
        scf_ket = self.ket_obj.scf

        if not self.same_braket():
            # sanity check that orbitals and geometry are the same
            if np.any(scf_bra.orbs != scf_ket.orbs):
                sys.exit('spin-orbit coupling requires same'+
                         'bra/ket orbs')

            if mol_bra.pymol().atom != mol_ket.pymol().atom:
                sys.exit('spin-orbit coupling requires same geometry'+
                         ' and basis set')

        # check on the MF 2e integral scheme
        if self.mf2e not in self.allowed_mf2e:
            print('\n Unrecognised mf2e value: '+self.mf2e
                  +'\n Allowed values: '+str(self.allowed_mf2e))
            sys.exit()

        bra_spin = self.SpinInfo(self.bra_obj)
        ket_spin = self.SpinInfo(self.ket_obj)

        # Delta S = -1, 0 or +1 must hold
        if bra_spin.S - ket_spin.S not in [-1., 0., 1.]:
            sys.exit('\n ERROR: S_bra, S_ket spin combination not ' \
                     +'currently supported in spinorbit')

        # In the case of Delta S = 0, S > 0 must hold
        if bra_spin.S == ket_spin.S and bra_spin.S == 0.:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination ' \
                     +'in spinorbit')

        # bra and ket multiplet indices and M-values
        self.stlbl = self.set_stlbl(bra_spin, ket_spin)

        # get the MF density matrix
        self.rho_ao = self.build_rho()
        
        # get the one-electron SOC matrices
        self.h1e = self.build_h1e()
        
        # initialise the H_soc array
        hdim = len(self.bra_states) * bra_spin.mult \
             + len(self.ket_states) * ket_spin.mult
        hsoc = np.zeros((hdim, hdim), dtype=np.cdouble)

        # construct the trans_list arrays for each of the
        # bra-bra, ket-ket and bra-ket blocks of H_SOC
        # currently store them as [initial state, final state]
        blocks = ['lower', 'lower', 'full']
        grps   = [['bra','bra'], ['ket','ket'], ['bra','ket']]

        # add the state groups
        self.add_states('bra', self.bra_obj, self.bra_states)
        self.add_states('ket', self.ket_obj, self.ket_states)
        self.add_spins('bra', bra_spin)
        self.add_spins('ket', ket_spin)

        # loop over the blocks of H_SOC
        for iblock in range(len(blocks)):
            bra_lbl  = grps[iblock][0]
            ket_lbl  = grps[iblock][1]
            blk_fill = blocks[iblock]

            b_list     = self.build_pair_list(bra_lbl, ket_lbl, 
                                              pairs=blk_fill)
            b_list_sym = self.build_pair_list(bra_lbl, ket_lbl, 
                                              pairs=blk_fill, 
                                              sym_blk=True)

            bra = self.get_bkobj(bra_lbl)
            ket = self.get_bkobj(ket_lbl)

            # skip if both manifolds of states are singlets
            if bra.mult == 1 and ket.mult == 1:
                continue
            
            # initialize the bitsi library for the calculation of
            # spin-orbit matrix elements
            bitsi_init.init(bra, ket, 'soc')

            # get the reduced matrix elements
            # <S' M'=S' psi_m||T_ij^(1)||S M=S psi'_n>
            # for the current block
            redmat = self.build_redmat(bra, ket, b_list, b_list_sym)

            # get the Clebsh-Gordan coefficients needed to
            # compute the SOC matrix elements for the current block
            cg_coef = mrci_soc.clebsch_gordan(bra, ket)

            # build this block of H_SOC
            hblk = self.build_hsoc(hdim, bra_lbl, ket_lbl, b_list, 
                                   redmat, cg_coef)
            
            # add the current block of Hsoc to the total
            hsoc += hblk

            # finalize the bitsi library
            bitsi_init.finalize()
           
        # add in the on-diagonal elements
        hsoc += np.diag(self.hsoc_ondiag(hdim))
        
        # diagonalise H_soc
        h_eig, h_vec = np.linalg.eigh(hsoc)
        
        # print the H_SOC elements
        self.print_hsoc(hsoc, h_eig, h_vec)

        del(cg_coef, redmat)

        return
    
    #
    def set_stlbl(self, bra_spin, ket_spin):
        """
        list of all all multiplet indices and projected spin values
        """
        
        stlbl = []
        
        for i in self.ket_states:
            for m in ket_spin.M:
                stlbl.append([self.ket_obj.label, i, m])
        
        for i in self.bra_states:
            for m in bra_spin.M:
                stlbl.append([self.bra_obj.label, i, m])

        return stlbl
    
    #
    def build_rho(self):
        """
        sets up the mean-field density matrix
        """
        nsta = len(self.bra_states) + len(self.ket_states)
        nmo    = self.bra_obj.scf.nmo

        # Average density matrix across states
        rho_mo = np.zeros((nmo, nmo), dtype=float)

        for i in self.ket_states:
            rho_mo += self.ket_obj.rdm(i) / nsta

        for i in self.bra_states:
            rho_mo += self.bra_obj.rdm(i) / nsta
        
        # Transform to the AO basis
        orbs   = self.bra_obj.scf.orbs
        rho_ao = np.matmul(np.matmul(orbs, rho_mo), orbs.T)

        return rho_ao
        
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
            redmat_blk[:,:,indx] = redmat_list[birr][kirr][:,:,indx_sym]

        return redmat_blk

    #
    def build_h1e(self):
        """
        Sets up the one-electron SOC matrices h^(k), k=-1,0,+1
        """

        # PySCF Mole obeject
        mol = self.bra_obj.scf.mol.pymol().copy()

        # No. AOs
        nao = mol.nao_nr()
        
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
            hcart_ao += (np.einsum("ijklm,lm->ijk", h2e_ao, self.rho_ao)
                         - 1.5 * (np.einsum("ijklm, kl->ijm", h2e_ao, self.rho_ao)
                         + np.einsum("ijklm,mj->ilk", h2e_ao, self.rho_ao)))
        elif self.mf2e == 'atomic':
            hcart_ao += self.build_mf_atomic(mol)
            
        # Transform to the MO basis
        orbs     = self.bra_obj.scf.orbs
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
    def build_mf_atomic(self, mol):
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
            rho = self.rho_ao[ji:jf, ji:jf]

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
    def build_hsoc(self, dim, bra_lbl, ket_lbl, pair_list, 
                   redmat, cg_coef):
        """
        Builds the iblock'th block of the SOC Hamiltonian matrix
        """
        
        # get the corresponding bra/ket spin objects
        bra_spin = self.get_spins(bra_lbl)
        ket_spin = self.get_spins(ket_lbl)

        # Number of bra and ket multiplets in this block
        bra_states = self.get_states(bra_lbl)[0]
        ket_states = self.get_states(ket_lbl)[0]
        nm_bra = len(bra_states)
        nm_ket = len(ket_states)

        # Loop over pairs of spin-orbit coupled states in this block
        Mb = bra_spin.M
        Mk = ket_spin.M

        # fill in the appropriate elements of hsoc
        hsoc = np.zeros((dim, dim), dtype=np.cdouble)

        for pair in pair_list:
            indx  = pair_list.index(pair)
            I_bra = bra_states.index(pair[0])
            I_ket = ket_states.index(pair[1])

            for M_ket in Mk:
                for M_bra in Mb:

                    i,j  = self.hsoc_indx(bra_lbl, ket_lbl, 
                                          I_bra, I_ket,
                                          M_bra, M_ket)
                    hsoc[i,j] = self.contract_redmat(
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

        for i in range(dim):
            indx = self.stlbl[i][1]
            
            if self.stlbl[i][0] == self.bra_obj.label:
                hdiag[i] = self.bra_obj.energies[indx]
            else:
                hdiag[i] = self.ket_obj.energies[indx]
        
        return hdiag
        
    #
    def contract_redmat(self, bra_spin, ket_spin, M_bra, M_ket, 
                        cg_coef, redmat):
        """
        Computes a single element 
        < I_bra M_bra | H_SOC | I_ket M_ket >
        of the SOC Hamiltonian matrix via the contraction
        of the reduced matrices with the Clebsch-Gordan coefficient-
        scaled one-electron SOC matrices h^(k), k=-1,0,+1
        """
        
        # Sum of the scaled one-electron SOC matrices
        nmo    = self.bra_obj.scf.nmo
        hscale = np.zeros((nmo, nmo), dtype=np.cdouble)
        kval   = [-1, 0, 1]
        coe    = [1., np.sqrt(2.), -1.]
        for n in range(3):
            i, i12  = self.cgcoe_indx(bra_spin.S, M_bra, 
                                      ket_spin.S, M_ket, kval[n])
            hscale += self.h1e[2-n, :, :] * cg_coef[i12, i] * coe[n]
        
        # Contraction with the reduced matrix
        hij = 0.5 * np.einsum('ij,ij', redmat, hscale)

        return hij

    #
    def hsoc_indx(self, bra_lbl, ket_lbl, I_bra, I_ket, M_bra, M_ket):
        """
        Given bra and ket multiplet indices, I_bra and I_ket,
        and projected spins, M_bra and M_ket, returns the
        corresponding element of the SOC Hamiltonian matrix
        """
        bra_spin = self.get_spins(bra_lbl)
        ket_spin = self.get_spins(ket_lbl)

        indx_ket = ket_spin.mult * I_ket + M_ket + ket_spin.S
        indx_bra = bra_spin.mult * I_bra + M_bra + bra_spin.S

        # this is a hacky -- should be change
        if bra_lbl == ket_lbl and bra_lbl == 'bra':
            indx_bra += len(self.ket_states)
            indx_ket += len(self.ket_states)
        elif bra_lbl != ket_lbl:
            indx_bra += len(self.ket_states)

        return int(indx_bra), int(indx_ket)

    
    #
    def cgcoe_indx(self, S, M, s1, m1, k):
        """
        returns the indices for the Clebsch-Gordan coefficient
        < s1 m1; 1 k | S M >
        """

        i = int(M + S)
        
        i1  = int(m1 + s1) + 1

        i2  = k + 2

        i12 = (i1 - 1) * 3 + i2 - 1
        
        return i, i12


    #
    def print_hsoc(self, hsoc, heig, hvec):
        """
        prints the SOC values to the log file
        """
        
        # output the SOC matrix elements
        output.print_spinorbit_table(hsoc, hsoc.shape[0], self.stlbl,
                                     self.print_thresh)

        # output the eigenvectors of H_SOC
        output.print_hsoc_eig(heig, hvec, hsoc.shape[0], self.stlbl)
        
        return
