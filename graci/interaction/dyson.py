"""
Dyson orbital class
"""

import sys as sys
import numpy as np
import graci.utils.timing as timing
import graci.interaction.interaction as interaction
import graci.interfaces.bitci.bitwf_init as bitwf_init
import graci.interfaces.bitci.wf_overlap as wf_overlap
import graci.interfaces.bitci.wf_dyson as wf_dyson
import graci.io.output as output

class Dyson(interaction.Interaction):
    """
    Calculation of Dyson orbitals using the bitwf and
    overlap libraries
    """
    def __init__(self):
        super().__init__()

        # user defined quantities
        self.init_states    = None
        self.init_label     = None
        self.final_states   = None
        self.final_label    = None
        self.norm_thresh    = 0.999
        self.print_orbitals = False
        
        # ----------------------------------------------------------
        # internal class variables -- should not be accessed
        # directly
        self.bra_obj       = None
        self.ket_obj       = None
        self.bra_wfunit    = None
        self.ket_wfunit    = None
        self.dyson_orbs    = None
        self.dyson_norms   = None
        self.mo_basis      = None
        self.n_basis       = 0
                
    def copy(self):
        """create of deepcopy of self"""
        new = self.Dyson()

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
    def run(self, obj_list):
        """
        Computes the Dyson orbitals between the given
        bra and ket states
        """

        # set the bra and ket objects
        self.bra_obj, self.ket_obj = self.set_objs(obj_list)
        
        # section header
        if self.verbose:
            output.print_dyson_header(self.label)

        # check on the requested calculation
        self.check_calc()

        # mol object (has to be the same for both bra and ket)
        self.mol = self.ket_obj.scf.mol
        
        # initialise the bitwf library
        bitwf_init.init(self.bra_obj, self.ket_obj, 'dyson',
                        self.verbose)

        # set the MO basis: take the N-electron WF MOs
        # as the basis for now
        if self.bra_obj.nel < self.ket_obj.nel:
            self.mo_basis = 'ket'
            self.n_basis  = self.ket_obj.nmo
        else:
            self.mo_basis = 'bra'
            self.n_basis  = self.bra_obj.nmo
            
        # construct the list of bra-ket pairs
        self.add_group('ket', self.ket_obj, self.init_states)
        self.add_group('bra', self.bra_obj, self.final_states)
        self.trans_list = self.build_pair_list('bra', 'ket',
                                               pairs='full')
        self.trans_list_sym = self.build_pair_list('bra', 'ket',
                                               pairs='full',
                                               sym_blk=True)

        # compute the Dyson orbitals
        self.dyson_orbs, self.dyson_norms = \
            self.build_dyson_orbs(self.bra_obj,
                                  self.ket_obj,
                                  self.trans_list,
                                  self.trans_list_sym)

        # output the ionisation energies and squared Dyson orbital
        # norms
        if self.verbose:
            self.print_log()

#----------------------------------------------------------------------
# "Private" class methods
#

    #
    def set_objs(self, obj_list):
        """
        Determines the bra and ket CI objects based on their labels
        """

        bra_obj = None
        ket_obj = None

        if obj_list[0].label == self.final_label:
            bra_obj = obj_list[0]
            ket_obj = obj_list[1]
        else:
            ket_obj = obj_list[0]
            bra_obj = obj_list[1]
            
        return bra_obj, ket_obj

    #
    def check_calc(self):
        """
        Sanity check on the requested calculation
        """

        # make sure that |S_bra - S_ket| = 0.5
        if abs(self.bra_obj.mult - self.ket_obj.mult) != 1:
            sys.exit('Error: |S_bra - S_ket| != 0.5')
        
        # bitwf currently only supports equal bra and ket point groups
        if self.bra_obj.scf.mol.sym_indx != self.ket_obj.scf.mol.sym_indx:
            sys.exit('Error: unequal bra and ket point groups')
            
        return

    #
    def build_dyson_orbs(self, bra, ket, trans_list, trans_list_sym):
        """
        Grab the Dyson orbitals from bitwf and then reshape the list
        of these into a more usable format
        """

        # extract the determinant representation of the wave functions
        self.bra_wfunit, self.ket_wfunit = wf_overlap.extract(bra, ket)

        # compute the Dyson orbitals
        dyson_list = wf_dyson.dyson(bra, ket, self.mo_basis,
                                    self.n_basis, self.bra_wfunit,
                                    self.ket_wfunit, self.norm_thresh,
                                    trans_list_sym)

        # make the Dyson orbital list
        npairs     = len(trans_list)
        dyson_orbs = np.zeros((npairs, self.n_basis), dtype=float)
        for indx in range(npairs):
            bk_st               = trans_list[indx]
            [birr, bst]         = bra.state_sym(bk_st[0])
            [kirr, kst]         = ket.state_sym(bk_st[1])
            sym_indx            = trans_list_sym[birr][kirr].index([bst,kst])
            dyson_orbs[indx, :] = dyson_list[birr][kirr][:, sym_indx]

        # Dyson orbital norms
        dyson_norms = np.array([np.linalg.norm(dyson_orbs[i, :])
                                for i in range(npairs)])

        return dyson_orbs, dyson_norms

    #
    def print_log(self):
        """
        print a summary of the ionisation energies and
        squared Dyson orbital norms to the log file
        """

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

        # print an 'ionisation table' for each initial state
        for iket in range(len(ket_states)):

            # shift state indices from 0..n-1 to 1..n
            init_st   = ket_states[iket]+1
            init_sym  = ksym_lbl[ket_syms[iket]]
            final_st  = []
            final_sym = []
            exc_ener  = []
            sqnorm    = []     

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
                sqnorm.append(self.dyson_norms[indx]**2)

            # print an 'ionisation table' for each initial state
            output.print_dyson_table(init_st,
                                     init_sym,
                                     final_st,
                                     final_sym,
                                     exc_ener, 
                                     sqnorm)
        
        return
