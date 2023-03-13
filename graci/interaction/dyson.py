"""
Dyson orbital class
"""

import sys as sys
import numpy as np
import copy
import graci.core.params as params
import graci.utils.timing as timing
import graci.interaction.interaction as interaction
import graci.interfaces.bitci.bitwf_init as bitwf_init
import graci.interfaces.bitci.wf_overlap as wf_overlap
import graci.interfaces.bitci.wf_dyson as wf_dyson
import graci.core.orbitals as orbitals
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
        self.representation = 'adiabatic'

        # ----------------------------------------------------------
        # internal class variables -- should not be accessed
        # directly
        self.dyson_orbs    = None
        self.dyson_norms   = None
        self.mo_basis      = None
        self.n_basis       = 0
        self.allowed_reps = ['adiabatic']    
            
    def copy(self):
        """create of deepcopy of self"""
        new = Dyson()

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
    def run(self, bra, ket):
        """
        Computes the Dyson orbitals between the given
        bra and ket states
        """

        # section header
        if self.verbose:
            output.print_dyson_header(self.label)

        # sanity check on the representation
        self.check_representation()

        # set the bra/ket objects and add the state groups associated
        # with each 
        self.add_group('bra', [bra], states = [self.final_states])
        self.add_group('ket', [ket], states = [self.init_states])

        # check on the requested calculation
        self.mol = self.check_bra_ket()

        # can't be same object -- pairs = full
        self.trans_list = self.grp_pair_list('bra', 'ket', pairs='full')

        # initialize the dyson orbs
        ntrans = len(self.trans_list)
        self.dyson_orbs = np.zeros((ntrans, self.n_basis),dtype=complex)

        # loop over all pairs of spin free states in both the bra
        # and ket objects. In the future, may be prudent to check
        # init and final states and just which states are necessary
        for k_lbl in self.get_ci_lbls('ket'):
            mult_k, S_k, M_k = self.get_ci_spins('ket', k_lbl)
            ci_k             = self.get_ci_obj('ket', k_lbl)

            for b_lbl in self.get_ci_lbls('bra'):
                mult_b, S_b, M_b = self.get_ci_spins('bra', b_lbl)
                ci_b             = self.get_ci_obj('bra', b_lbl)

                if abs(mult_k - mult_b) != 1:
                    continue

                # initialise the bitwf library
                bitwf_init.init(ci_b, ci_k, 'dyson', self.verbose)

                # this is main transition_list: stored by adiabatic label
                ci_tran = self.ci_pair_list('bra', b_lbl, 
                                            'ket', k_lbl, 
                                             pairs='full')

                # the transitions, ordered by interacting irreps, is how
                # we call bitsi
                ci_tran_sym = self.ci_pair_list('bra', b_lbl, 
                                                'ket', k_lbl, 
                                                 pairs='full', 
                                                 sym_blk=True)

                # compute the Dyson orbitals
                dys_orbs, dys_nrms = self.build_ci_orbs(b_lbl, k_lbl, 
                                                        ci_tran, 
                                                        ci_tran_sym)

                # accumulate the dyson orbital in terms of the 
                # bra/ket states 
                self.build_grp_tensor('bra', b_lbl, 
                                      'ket', k_lbl,
                                      ci_tran, dys_orbs, 
                                      self.dyson_orbs)

                bitwf_init.finalize()

        # determine the norms of the orbitals
        self.dyson_norms = np.array([np.linalg.norm(self.dyson_orbs[i,:])
                                   for i in range(len(self.trans_list))]) 

        # print the Dyson orbitals to file
        if self.print_orbitals:
            self.export_dyson_orbs()

        # output the ionisation energies and squared Dyson orbital
        # norms
        if self.verbose:
            self.print_log()


    #
    def get_dyson_norms(self, istate=None, fstate=None):
        """
        Return requested dyson orbital norms
        """

        if istate is None and fstate is None:
            return self.dyson_norms
        else:
            if istate is None:
                chk_st  = fstate
                chk_ind = 1
            else:
                chk_st  = istate
                chk_ind = 0
            lst = [self.trans_list.index(pair)
                  for pair in self.trans_list if pair[chk_ind]==chk_st]
            return self.dyson_norms[lst]


#----------------------------------------------------------------------
# "Private" class methods
#

    #
    def check_bra_ket(self):
        """
        Sanity check on the requested calculation
        """

        valid_pair    = False
        self.mo_basis = ''
        self.n_basis  = 0

        for b_lbl in self.get_ci_lbls('bra'):
            mol_b            = self.get_ci_obj('bra',b_lbl).scf.mol
            mult_b, S_b, M_b = self.get_ci_spins('bra',b_lbl)

            for k_lbl in self.get_ci_lbls('ket'):
                mol_k            = self.get_ci_obj('ket',k_lbl).scf.mol
                mult_k, S_k, M_k = self.get_ci_spins('ket',k_lbl)

                # make sure that |S_bra - S_ket| = 0.5
                if abs(mult_b - mult_k) == 1:
                    valid_pair = True

                # right now, we're going to assume that all bra states and
                # all ket states have same number of mos and same number of 
                # electrons
                nel_b = self.get_ci_obj('bra', b_lbl).nel
                nmo_b = self.get_ci_obj('bra', b_lbl).nmo
                nel_k = self.get_ci_obj('ket', k_lbl).nel
                nmo_k = self.get_ci_obj('ket', k_lbl).nmo
                if nel_b < nel_k:
                    if self.mo_basis == 'bra':
                        msg = 'ERROR: not all bra objects have '+\
                              ' nel_bra < nel_ket'
                        sys.exit(msg)
                    self.mo_basis = 'ket'
                    self.n_basis  = nmo_k
                else:
                    if self.mo_basis == 'ket':
                        msg = 'ERROR: not all ket objects have '+\
                              ' nel_ket < nel_bra'
                        sys.exit(msg)
                    self.mo_basis = 'bra'
                    self.n_basis  = nmo_b        

                # bitwf currently only supports equal bra and ket point groups
                if mol_b.sym_indx != mol_k.sym_indx:
                    sys.exit('Error: unequal bra and ket point groups')
        
                # sanity check that orbitals and geometry are
                # the same
                if (mol_b.pymol().atom != mol_k.pymol().atom):
                    sys.exit('ERROR: transition moments require same '+
                             ' geometry and basis set')

        if not valid_pair:
            sys.exit('Error: |S_bra - S_ket| != 0.5')

        return mol_b

    #
    def build_ci_orbs(self, b_lbl, k_lbl, ci_trans, ci_trans_sym):
        """
        Grab the Dyson orbitals from bitwf and then reshape the list
        of these into a more usable format
        """

        # get ci objects
        bra = self.get_ci_obj('bra', b_lbl)
        ket = self.get_ci_obj('ket', k_lbl) 

        # extract the determinant representation of the wave functions
        bra_wfunit, ket_wfunit = \
            wf_overlap.extract(bra, ket, rep=self.representation)

        # compute the Dyson orbitals
        dys_list = wf_dyson.dyson(bra, ket, 
                                  self.mo_basis, self.n_basis,
                                  bra_wfunit, ket_wfunit, 
                                  self.norm_thresh,
                                  ci_trans_sym)

        # make the Dyson orbital list
        npairs = len(ci_trans)         
        dyson_orbs = np.zeros((npairs, self.n_basis), dtype=float)

        for indx in range(npairs):
            bk_st          = ci_trans[indx]
            [bir, bst]     = self.get_ci_sym('bra', b_lbl, bk_st[0])
            [kir, kst]     = self.get_ci_sym('ket', k_lbl, bk_st[1])
            sym_indx       = ci_trans_sym[bir][kir].index([bst, kst])
            dyson_orbs[indx, :] = dys_list[bir][kir][:, sym_indx]

        # Dyson orbital norms
        dyson_norms = np.array([np.linalg.norm(dyson_orbs[i, :])
                                for i in range(npairs)])

        return dyson_orbs, dyson_norms

    #
    #
    def export_dyson_orbs(self):
        """
        Writing of the Dyson orbitals to molden files
        """

        # AO representation of the Dyson orbitals
        nmo = self.dyson_orbs.shape[1]
        nao = self.mol.nao

        if self.mo_basis == 'ket':
            ci_lbl = self.get_ci_lbls('ket')[0]
            orbs   = self.get_ci_obj('ket', ci_lbl).scf.orbs
            dysao  = np.matmul(self.dyson_orbs, orbs.T[:nmo,:])
        else:
            ci_lbl = self.get_ci_lbls('bra')[0]
            orbs   = self.get_ci_obj('bra', ci_lbl).scf.orbs
            dysao  = np.matmul(self.dyson_orbs, orbs.T[:nmo,:])

        # export the (normalised) non-zero Dyson orbitals to Molden files
        sqnorm = np.zeros((1), dtype=float)
        for ido in range(len(self.trans_list)):
            sqnorm[0] = self.dyson_norms[ido]**2
            if sqnorm[0] < 1e-6:
                continue
            bst       = self.trans_list[ido][0]
            kst       = self.trans_list[ido][1]
            fname     = 'dyson_'+str(kst+1)+'_to_'+str(bst+1)+'_molden'
            orb       = dysao[ido,:].reshape(nao,1) / self.dyson_norms[ido]
            orbitals.export_orbitals(fname,
                                     self.mol,
                                     orb.real,
                                     orb_ener=sqnorm,
                                     orb_dir='Dyson.'+str(self.label),
                                     fmt='molden')

            if np.linalg.norm(orb.imag) > 0.2:
                orbitals.export_orbitals(fname+'_imag',
                                         self.mol,
                                         orb.imag,
                                         orb_ener=sqnorm,
                                         orb_dir='Dyson.'+str(self.label),
                                         fmt='molden')

        return


    #
    def print_log(self):
        """
        print a summary of the ionisation energies and
        squared Dyson orbital norms to the log file
        """

        # get the ket states and energies
        ket_states = self.get_states('ket')
        ket_eners  = self.get_energy('ket')

        # get the bra states and energies
        bra_states = self.get_states('bra')
        bra_eners  = self.get_energy('bra')

        # get state symmetries. If not defined, use C1 sym labels
        bsym_lbl = self.get_sym_lbl('bra')
        ksym_lbl = self.get_sym_lbl('ket')

        # print an 'ionisation table' for each initial state
        for iket in range(len(ket_states)):

            # shift state indices from 0..n-1 to 1..n
            init_st   = ket_states[iket]+1
            init_sym  = ksym_lbl[iket]
            final_st  = []
            final_sym = []
            exc_ener  = []
            sqnorm    = []     

            for ibra in range(len(bra_states)):

                indx = self.trans_list.index([bra_states[ibra], 
                                              ket_states[iket]])
                # shift state indices from 0..n-1 to 1..n
                final_st.append(bra_states[ibra]+1)
                final_sym.append(bsym_lbl[ibra])
                exc_ener.append(bra_eners[ibra] - ket_eners[iket]) 
                sqnorm.append(self.dyson_norms[indx]**2)

            # print an 'ionisation table' for each initial state
            output.print_dyson_table(init_st,
                                     init_sym,
                                     final_st,
                                     final_sym,
                                     exc_ener, 
                                     sqnorm)
        
        return
