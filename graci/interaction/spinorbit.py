"""
module to compute spin-orbit coupling matrix 
elements
"""

import sys as sys
import numpy as np
import graci.interaction.interaction as interaction
import graci.utils.timing as timing
import graci.core.libs as libs
import graci.bitcitools.bitsi_init as bitsi_init
import graci.bitcitools.mrci_soc as mrci_soc
import graci.io.output as output
import graci.utils.constants as constants

class Spinorbit(interaction.Interaction):
    """Spin orbit coupling class. For now will only accept calculations
       for which the 'Molecule and 'Scf' objects are the same"""
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
        
    #
    @timing.timed
    def run(self):
        """returns the SOC matrix elements between the bra and ket
           states. If b_state and k_state are None, assume 
           SOC matrix elements from all states in method object should
           be used."""

        # Check on spin multiplicities:
        # (1) Delta S = -1, 0 or +1 must hold
        # (2) In the case of Delta S = 0, only S > 0 must hold
        S_bra = (self.bra_obj.mult-1)/2
        S_ket = (self.ket_obj.mult-1)/2
        if S_bra - S_ket not in [-1., 0., 1.]:
            sys.exit('\n ERROR: non-sensical S_bra, S_ket combination' \
                     +'in spinorbit')
        if S_bra == S_ket and S_bra == 0.:
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
            
    #
    @timing.timed
    def build_redmat(self):
        """grab the reduced matrix elements from bitsi and then
           reshape the list into a more usable format"""

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
