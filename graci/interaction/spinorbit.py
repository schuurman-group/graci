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
        
        # global variables
        # Class name
        self.class_name = 'spinorbit'
        # labels of the manifolds of states in the 3 unique
        # blocks of H_SOC
        self.blkstr     = [['bra', 'bra'], ['ket', 'ket'], ['bra', 'ket']]
        self.ifdict     = {'bra' : 'final', 'ket' : 'init'}
        # Clebsch-Gordan coeficients
        self.cgcoe      = []
        # Reduced matrix elements
        self.redmat     = []
        # Bra and ket total spins
        self.S_bra      = None
        self.S_ket      = None
        # Bra and ket projected spins
        self.M_bra      = None
        self.M_ket      = None
        # Bra and ket spin multiplicities
        self.mult_bra   = None
        self.mult_ket   = None
        # One electron SOC matrices
        self.h1e        = None
        # SOC Hamiltonian matrix
        self.hdim       = None
        self.hsoc       = None
        # eigenpairs of H_SOC
        self.vec        = None
        self.eig        = None
        # List of multiplet indices and projected spins
        self.stlbl      = None
        # MF density matrix (AO basis)
        self.rho_ao     = None

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

        # check on the MF 2e integral scheme
        if self.mf2e not in self.allowed_mf2e:
            print('\n Unrecognised mf2e value: '+self.mf2e
                  +'\n Allowed values: '+str(self.allowed_mf2e))
            sys.exit()
        
        # bra and ket total spins, multiplicities, etc
        self.set_spins()

        # bra and ket multiplet indices and M-values
        self.set_stlbl()

        # get the MF density matrix
        self.build_rho()
        
        # first check to see if bra and ket are identical
        if (type(self.bra_obj).__name__ == type(self.ket_obj).__name__ 
            and self.bra_obj.label == self.ket_obj.label):
            self.braket_iden = True

        # bra/ket mol and scf objects
        mol_bra = self.bra_obj.scf.mol
        mol_ket = self.ket_obj.scf.mol
        scf_bra = self.bra_obj.scf
        scf_ket = self.ket_obj.scf

        if not self.braket_iden:
            # sanity check that orbitals and geometry are the same
            if np.any(scf_bra.orbs != scf_ket.orbs):
                sys.exit('spin-orbit coupling requires same'+
                         'bra/ket orbs')

            if mol_bra.pymol().atom != mol_ket.pymol().atom:
                sys.exit('spin-orbit coupling requires same geometry'+
                         ' and basis set')
                
        # construct the trans_list arrays for each of the
        # bra-bra, ket-ket and bra-ket blocks of H_SOC
        # currently store them as [initial state, final state]
        blocks              = ['lower', 'lower', 'full']        
        for iblock in range(3):
            bstr = self.blkstr[iblock][0]
            kstr = self.blkstr[iblock][1]
            kets, trans, trans_sym = self.build_trans_list_new(
                bstring=bstr, kstring=kstr, list_type=blocks[iblock])
            self.ket_list.append(kets)
            self.trans_list.append(trans)
            self.trans_list_sym.append(trans_sym)

        # get the one-electron SOC matrices
        self.build_h1e()

        # initialise the H_soc array
        self.hdim = len(self.final_states) * self.mult_bra \
            + len(self.init_states)  * self.mult_ket
        self.hsoc = np.zeros((self.hdim, self.hdim), dtype=np.cdouble)
        
        # loop over the blocks of H_SOC
        for iblock in range(3):

            bra = getattr(self, self.blkstr[iblock][0]+'_obj')
            ket = getattr(self, self.blkstr[iblock][1]+'_obj')

            # skip if both manifolds of states are singlets
            if bra.mult == 1 and ket.mult == 1:
                self.redmat.append([])
                self.cgcoe.append([])
                continue
            
            # initialize the bitsi library for the calculation of
            # spin-orbit matrix elements
            bitsi_init.init(bra, ket, 'soc')

            # get the reduced matrix elements
            # <S' M'=S' psi_m||T_ij^(1)||S M=S psi'_n>
            # for the current block
            self.build_redmat(bra, ket, iblock)

            # get the Clebsh-Gordan coefficients needed to
            # compute the SOC matrix elements for the current block
            coe = mrci_soc.clebsch_gordan(bra, ket)
            self.cgcoe.append(coe)

            # build this block of H_SOC
            self.build_hsoc(iblock)
            
            # finalize the bitsi library
            bitsi_init.finalize()

        # add in the on-diagonal elements
        self.hsoc_ondiag()

        # diagonalise H_soc
        self.eig, self.vec = np.linalg.eigh(self.hsoc)

        # print the H_SOC elements
        self.print_hsoc()

        del(self.cgcoe, self.redmat)
    
    #
    def set_spins(self):
        """
        sets the bra and ket S, M, etc. values
        """

        # Total spins
        self.mult_bra = self.bra_obj.mult
        self.mult_ket = self.ket_obj.mult

        # Spin multiplicities
        self.S_bra    = (self.mult_bra - 1)/2
        self.S_ket    = (self.mult_ket - 1)/2

        # All projected spins
        self.M_bra = [-self.S_bra + i for i in range(self.mult_bra)]
        self.M_ket = [-self.S_ket + i for i in range(self.mult_ket)]
        
        # Delta S = -1, 0 or +1 must hold
        if self.S_bra - self.S_ket not in [-1., 0., 1.]:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination' \
                     +'in spinorbit')
        
        # In the case of Delta S = 0, S > 0 must hold
        if self.S_bra == self.S_ket and self.S_bra == 0.:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination ' \
                     +'in spinorbit')
        
        return

    
    #
    def set_stlbl(self):
        """
        list of all all multiplet indices and projected spin values
        """
        
        self.stlbl = []
        
        for i in self.init_states:
            for m in self.M_ket:
                self.stlbl.append([self.ket_obj.label, i, m])
        
        for i in self.final_states:
            for m in self.M_bra:
                self.stlbl.append([self.bra_obj.label, i, m])

    
    #
    def build_rho(self):
        """
        sets up the mean-field density matrix
        """
        nsta = len(self.init_states) + len(self.final_states)
        nmo    = self.bra_obj.scf.nmo

        # Average density matrix across states
        rho_mo = np.zeros((nmo, nmo), dtype=float)

        for i in self.init_states:
            rho_mo += self.ket_obj.dmat[i, :, :] / nsta

        for i in self.final_states:
            rho_mo += self.bra_obj.dmat[i, :, :] / nsta
        
        # Transform to the AO basis
        orbs        = self.bra_obj.scf.orbs
        self.rho_ao = np.matmul(np.matmul(orbs, rho_mo), orbs.T)

        return
        
    #
    @timing.timed
    def build_redmat(self, bra, ket, iblock):
        """
        grab the reduced matrix elements from bitsi and then
        reshape the list into a more usable format
        """

        # grab the reduced matrix elements
        redmat_list = mrci_soc.redmat(bra, ket, self.trans_list_sym[iblock])

        # Make the reduced matrix list
        nmo         = bra.scf.nmo
        npairs      = len(self.trans_list[iblock])
        self.redmat.append(np.zeros((nmo, nmo, npairs), dtype=float))

        for bk_st in self.trans_list[iblock]:
            [birr, bst_sym] = bra.state_sym(bk_st[0])
            [kirr, kst_sym] = ket.state_sym(bk_st[1])
            indx            = self.trans_index(bk_st[0], bk_st[1], iblock)
            sym_indx        = self.trans_sym_index( birr, kirr, 
                                                    bst_sym, kst_sym, iblock)
        
            self.redmat[iblock][:, :, indx] = \
                redmat_list[birr][kirr][:, :, sym_indx]

        del(redmat_list)


    #
    def build_h1e(self):
        """
        Sets up the one-electron SOC matrices h^(k), k=-1,0,+1
        """

        # PySCF Mole obeject
        mol = self.bra_obj.scf.mol.mol_obj.copy()

        # No. AOs
        nao = mol.nao_nr()
        
        # Initialise arrays
        self.h1e = np.zeros((3, nao, nao), dtype=np.cdouble)
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
        self.h1e[0,:,:] = hcart_mo[0,:,:] - ci * hcart_mo[1,:,:]
        # k = +1: x + iy
        self.h1e[2,:,:] = hcart_mo[0,:,:] + ci * hcart_mo[1,:,:]
        # k = 0: z 
        self.h1e[1,:,:] = hcart_mo[2,:,:]
        
        return

    
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
    def build_hsoc(self, iblock):
        """
        Builds the iblock'th block of the SOC Hamiltonian matrix
        """
        
        # Number of bra and ket multiplets in this block
        fstring = self.ifdict[self.blkstr[iblock][0]]
        istring = self.ifdict[self.blkstr[iblock][1]]
        nm_bra  = len(getattr(self, fstring+'_states'))
        nm_ket  = len(getattr(self, istring+'_states'))

        # Multiplicities in this block
        mult_bra = getattr(self, 'mult_'+self.blkstr[iblock][0])
        mult_ket = getattr(self, 'mult_'+self.blkstr[iblock][1])
        
        # Loop over pairs of spin-orbit coupled states in this block
        Mb = getattr(self, 'M_'+self.blkstr[iblock][0])
        Mk = getattr(self, 'M_'+self.blkstr[iblock][1])
        for I_ket in range(nm_ket):
            for M_ket in Mk:
                for I_bra in range(nm_bra):
                    for M_bra in Mb:

                        # Skip the lower triangle of this block if
                        # the bra and ket manifolds are the same
                        if iblock in [0, 1] and I_bra < I_ket:
                            continue
                        
                        # spin-coupled state indices
                        i, j = self.hsoc_indx(mult_bra, mult_ket,
                                              I_bra, I_ket,
                                              M_bra, M_ket, iblock)

                        # block of redmat under consideration
                        rindx = I_ket * nm_bra + I_bra
                        
                        # H_ij^(SOC)
                        self.hsoc[i, j] = \
                            self.contract_redmat(mult_bra, mult_ket,
                                                 I_bra, I_ket,
                                                 M_bra, M_ket,
                                                 iblock, rindx)
                        self.hsoc[j, i] = np.conj(self.hsoc[i, j])
                        
        return

    
    #
    def hsoc_ondiag(self):
        """
        Adds the on-diagonal elements of H_SOC
        """

        for i in range(self.hdim):

            indx = self.stlbl[i][1]
            
            if self.stlbl[i][0] == self.bra_obj.label:
                self.hsoc[i ,i] = self.bra_obj.mrci_ener[indx]
            else:
                self.hsoc[i ,i] = self.ket_obj.mrci_ener[indx]
        
        return
        
    #
    def contract_redmat(self, mult_bra, mult_ket, I_bra, I_ket,
                        M_bra, M_ket, iblock, rindx):
        """
        Computes a single element 
        < I_bra M_bra | H_SOC | I_ket M_ket >
        of the SOC Hamiltonian matrix via the contraction
        of the reduced matrices with the Clebsch-Gordan coefficient-
        scaled one-electron SOC matrices h^(k), k=-1,0,+1
        """

        S_ket = (mult_ket - 1) / 2
        S_bra = (mult_bra - 1) / 2
        
        # Sum of the scaled one-electron SOC matrices
        nmo    = self.bra_obj.scf.nmo
        hscale = np.zeros((nmo, nmo), dtype=np.cdouble)
        kval   = [-1, 0, 1]
        coe    = [1., np.sqrt(2.), -1.]
        for n in range(3):
            i, i12  = self.cgcoe_indx(S_bra, M_bra, S_ket, M_ket, kval[n])
            hscale += self.h1e[2-n, :, :] * self.cgcoe[iblock][i12, i] * coe[n]
        
        # Contraction with the reduced matrix
        hij = 0.5 * np.einsum('ij,ij', self.redmat[iblock][:, :, rindx], hscale)

        return hij

    #
    def hsoc_indx(self, mult_bra, mult_ket, I_bra, I_ket, M_bra, M_ket,
                  iblock):
        """
        Given bra and ket multiplet indices, I_bra and I_ket,
        and projected spins, M_bra and M_ket, returns the
        corresponding element of the SOC Hamiltonian matrix
        """

        S_ket = (mult_ket - 1) / 2
        
        S_bra = (mult_bra - 1) / 2
        
        indx_ket = mult_ket * I_ket + M_ket + S_ket
        
        indx_bra = mult_bra * I_bra + M_bra + S_bra

        if iblock == 0:
            indx_bra += len(self.init_states)
            indx_ket += len(self.init_states)
        elif iblock == 2:
            indx_bra += len(self.init_states)

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
    def print_hsoc(self):
        """
        prints the SOC values to the log file
        """
                
        # output the SOC matrix elements
        output.print_spinorbit_table(self.hsoc, self.hdim, self.stlbl,
                                     self.print_thresh)

        # output the eigenvectors of H_SOC
        output.print_hsoc_eig(self.eig, self.vec, self.hdim, self.stlbl)
        
        return
