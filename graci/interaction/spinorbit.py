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
        self.label      = 'spinorbit'

        # global variables
        # Class name
        self.class_name = 'spinorbit'
        # Labels of the manifolds of states in the 3 unique
        # blocks of H_SOC
        self.blkstr     = [['bra', 'bra'], ['ket', 'ket'], ['bra', 'ket']]
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
        self.hsoc       = None
        
    #
    @timing.timed
    def run(self):
        """
        returns the SOC matrix elements between the bra and ket
        states. If b_state and k_state are None, assume SOC matrix
        elements from all states in method object should be used
        """

        # bra and ket total spins, multiplicities, etc
        self.set_spins()

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
            
        # loop over the blocks of H_SOC
        for iblock in range(3):

            bra = getattr(self, self.blkstr[iblock][0]+'_obj')
            ket = getattr(self, self.blkstr[iblock][1]+'_obj')

            # zero block if both manifolds of states are singlets
            if bra.mult == 1 and ket.mult == 1:
                nmo    = bra.scf.nmo
                npairs = len(self.trans_list[iblock])
                self.redmat.append(np.zeros((nmo, nmo,npairs),
                                            dtype=float))
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

        sys.exit()

    #
    def set_spins(self):
        """
        sets the bra and ket S, M, etc. values
        """

        # Total spins
        self.mult_bra = self.bra_obj.mult
        self.mult_ket = self.ket_obj.mult

        # Spin multiplicities
        self.S_bra    = (self.mult_bra-1)/2
        self.S_ket    = (self.mult_ket-1)/2

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
        hcart_mo = np.zeros((3, nao, nao), dtype=float)
        
        # One-electron contributions
        for iatm in range(mol.natm):
            Z   = mol.atom_charge(iatm)
            xyz = mol.atom_coord(iatm)
            mol.set_rinv_orig(xyz)
            hcart_ao += Z * mol.intor('cint1e_prinvxp_sph', 3)

        # Mean-field two-electron contributions
        # (We will fill this is once the rest of the code is
        # debugged and verified working)
        
        
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
    def hsoc_indx(self, I_bra, I_ket, M_bra, M_ket):
        """
        Given bra and ket multiplet indices, I_bra and I_ket,
        and projected spins, M_bra and M_ket, returns the
        corresponding element of the SOC Hamiltonian matrix
        """
    
        indx_ket = self.mult_ket * I_ket + M_ket + self.S_ket
        
        indx_bra = self.mult_bra * I_bra + M_bra + self.S_bra \
                   + len(self.init_states)
        
        return int(indx_bra), int(indx_ket)


    #
    @timing.timed
    def build_hsoc(self, iblock):
        """
        Builds the iblock'th block of the SOC Hamiltonian matrix
        """

        #-----------------------------------------------------------
        # We need to modify the following to be able to work with
        # the bra-bra, ket-ket and bra-ket blocks
        #-----------------------------------------------------------
        # i.e., things like self.init_states need to be replaced
        # with getattr(self, ...), etc.
        #-----------------------------------------------------------
        # the same goes for the hsoc_indx and cgcoe_indx functions
        #-----------------------------------------------------------
        
        ## Number of bra and ket multiplets
        #nm_bra = len(self.final_states)
        #nm_ket = len(self.init_states)
        #
        ## Number of spin-orbit coupled states
        #dim = nm_bra * self.mult_bra + nm_ket * self.mult_ket
        #
        ## Initialise hsoc (complex matrix)
        #self.hsoc = np.zeros((dim, dim), dtype=np.cdouble)
        #
        ## Loop over pairs of spin-orbit coupled states
        #for I_ket in range(nm_ket):
        #    for M_ket in self.M_ket:
        #        for I_bra in range(nm_bra):
        #            for M_bra in self.M_bra:
        #
        #                # spin-coupled state indices
        #                i, j = self.hsoc_indx(I_bra, I_ket,
        #                                      M_bra, M_ket)
        #                # H_ij^(SOC)
        #                self.hsoc[i, j] = self.contract_redmat( I_bra, I_ket,
        #                                                        M_bra, M_ket)
                        
        return

    
    #
    def contract_redmat(self, I_bra, I_ket, M_bra, M_ket):
        """
        Computes a single element 
        < I_bra M_bra | H_SOC | I_ket M_ket >
        of the SOC Hamiltonian matrix via the contraction
        of the reduced matrices with the Clebsch-Gordan coefficient-
        scaled one-electron SO matrices h^(k), k=-1,0,+1
        """

        # Sum of Clebsch-Gordan coefficient scaled one-electron
        # SO matrices
        nmo    = self.bra_obj.scf.nmo
        hscale = np.zeros((nmo, nmo), dtype=np.cdouble)
        kval = [-1, 0, 1]
        for n in range(3):
            i, i12 = self.cgcoe_indx(self.S_bra, M_bra, self.S_ket,
                                     M_ket, kval[n])
            hscale += self.h1e[n, :, :] * self.cgcoe[i12, i]
            
        print('\n How do we access the correct block of self.redmat?')
        # (Look at self.build_redmat...)
        sys.exit()

        hij = 0.
        
        return hij

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
    
