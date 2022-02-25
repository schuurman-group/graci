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
        # Clebsch-Gordan coeficients
        self.cgcoe      = None
        # Reduced matrix elements
        self.redmat     = None
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

        # Bra and ket total spins and multiplicities
        self.mult_bra = self.bra_obj.mult
        self.mult_ket = self.ket_obj.mult
        self.S_bra    = (self.mult_bra-1)/2
        self.S_ket    = (self.mult_ket-1)/2

        # Bra and ket projected spins
        self.M_bra = [-self.S_bra + i for i in range(self.mult_bra)]
        self.M_ket = [-self.S_ket + i for i in range(self.mult_ket)]
        
        # Check on spin multiplicities:
        # (1) Delta S = -1, 0 or +1 must hold
        # (2) In the case of Delta S = 0, only S > 0 must hold
        if self.S_bra - self.S_ket not in [-1., 0., 1.]:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination' \
                     +'in spinorbit')
        if self.S_bra == self.S_ket and self.S_bra == 0.:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination ' \
                     +'in spinorbit')

        # first check to see if bra and ket are identical
        if (type(self.bra_obj).__name__ == type(self.ket_obj).__name__ 
            and self.bra_obj.label == self.ket_obj.label):
            self.braket_iden = True

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
                
        # initialize the bitsi library for the calculation of
        # spin-orbit matrix elements
        bitsi_init.init(self, 'soc')

        # construct the trans_list array
        # currently store them as [initial state, final state]
        #self.build_trans_list_old()
        self.build_trans_list(list_type='full')
        
        # get the reduced matrix elements
        # <S' M'=S' psi_m||T_ij^(1)||S M=S psi'_n>
        self.build_redmat()
        
        # get the Clebsh-Gordan coefficients needed to
        # compute the SOC matrix elements
        self.cgcoe = mrci_soc.clebsch_gordan(self.bra_obj, self.ket_obj)

        # get the one-electron SOC matrices
        self.build_h1e()

        # SOC Hamiltonian matrix elements
        hsoc = self.build_hsoc()
            
    #
    @timing.timed
    def build_redmat(self):
        """
        grab the reduced matrix elements from bitsi and then
        reshape the list into a more usable format
        """

        # grab the reduced matrix elements
        redmat_list = mrci_soc.redmat(self.bra_obj, self.ket_obj,
                                      self.trans_list_sym)

        # Make the reduced matrix list
        nmo         = self.bra_obj.scf.nmo
        npairs      = len(self.trans_list)
        self.redmat = np.zeros((nmo, nmo, npairs), dtype=float)

        for bk_st in self.trans_list:
            [birr, bst_sym] = self.bra_obj.state_sym(bk_st[0])
            [kirr, kst_sym] = self.ket_obj.state_sym(bk_st[1])
            indx            = self.trans_index(bk_st[0], bk_st[1])
            sym_indx        = self.trans_sym_index( birr, kirr, 
                                                bst_sym, kst_sym)

            self.redmat[:, :, indx] = redmat_list[birr][kirr][:, :, sym_indx]

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
        self.h1e = np.zeros((3, nao, nao))
                
        # One-electron contributions
        for iatm in range(mol.natm):
            Z   = mol.atom_charge(iatm)
            xyz = mol.atom_coord(iatm)
            mol.set_rinv_orig(xyz)
            self.h1e += Z * mol.intor('cint1e_prinvxp_sph', 3)

        # Mean-field two-electron contributions
        
        
        # Transform to the MO basis

        
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
                   + len(self.init_states) - 2
        
        return int(indx_bra), int(indx_ket)


    #
    def build_hsoc(self):
        """
        Builds the entire SOC Hamiltonian matrix
        """

        # Number of bra and ket multiplets
        nm_bra = len(self.final_states)
        nm_ket = len(self.init_states)

        # Number of spin-orbit coupled states
        dim_bra = nm_bra * self.mult_bra
        dim_ket = nm_ket * self.mult_ket

        # Initialise hsoc (complex matrix)
        hsoc = np.zeros((dim_bra, dim_ket), dtype=np.cdouble)

        # Loop over pairs of spin-orbit coupled states
        for I_ket in range(nm_ket):
            for M_ket in self.M_ket:
                for I_bra in range(nm_bra):
                    for M_bra in self.M_bra:

                        # spin-coupled state indices
                        i, j = self.hsoc_indx(I_bra, I_ket,
                                              M_bra, M_ket)

                        # H_ij^(SOC)
                        hsoc[i, j] = self.contract_redmat(I_bra, I_ket,
                                                          M_bra, M_ket)
                        
        return hsoc

    
    #
    def contract_redmat(self, I_bra, I_ket, M_bra, M_ket):
        """
        Computes a single element 
        < I_bra M_bra | H_SOC | I_ket M_ket >
        of the SOC Hamiltonian matrix via the contraction
        of the reduced matrices with the Clebsch-Gordan coefficient-
        scaled one-electron SO matrices h^(k), k=-1,0,+1
        """

        hij = 0.
        
        return hij
