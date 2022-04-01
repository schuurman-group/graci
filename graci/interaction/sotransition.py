"""
module to compute functions of the 1-TDM
"""

import os as os
import numpy as np
import functools
import operator
import importlib
from sympy import LeviCivita
import graci.core.params as params
import graci.interaction.interaction as interaction
import graci.interaction.spinorbit as spinorbit
import graci.utils.timing as timing
import graci.core.libs as libs
import graci.bitcitools.bitsi_init as bitsi_init
import graci.bitcitools.mrci_1tdm as mrci_1tdm
import graci.io.output as output
import graci.utils.constants as constants
import sys as sys

class Sotransition(transition.Transition):
    """Sotransition class for determing transition moments between spin-orbit
       copuled states. This class extends transition and will accpet either
       spinorbit objects, or, simply method objects. If the latter, will convert
       them into trivial spin orbit objects (i.e. with a single set of states and
       a single multiplicity)
    """
    def __init__(self):
        # parent attributes
        super().__init__()
    
        # input syntax should be identical to transition
        # name for this object
        self.label          = 'Sotransition'
      
        # ----------------------------------------------------
        # private variables
        self.init_state_vecs  = None
        self.final_state_vecs = None

    #
    @timing.timed
    def run(self, obj_list):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used.
        """

        self.bra_obj, self.ket_obj, self.scf, self.mol = \
                                           self.check_obj_list(obj_list)

        # if bra and ket are same object, just need lower triangle
        # minus the diagonal elements
        if not self.same_obj(self.bra_obj, self.ket_obj)::
            # sanity check that orbitals and geometry are the same
            if np.any(scf_bra.orbs != scf_ket.orbs):
                sys.exit('transition moments require same bra/ket orbs')

            if mol_bra.pymol().atom != mol_ket.pymol().atom:
                sys.exit('transition moments require same geometry'+
                         ' and basis set')

        # master transition list for the SO-coupled states
        if self.same_obj(self.bra_obj, self.ket_obj):
            list_type = 'nodiag'
        else:
            list_stype = 'full'
        self.add_group('ket_states', self.ket_obj, self.init_states)
        self.add_group('bra_states', self.bra_obj, self.final_states)
        self.trans_list = self.build_pair_list('bra_states',
                                               'ket_states',
                                                pairs=list_type)

        # loop over all pairs of spin free states in both the bra
        # and ket objects. In the future, may be prudent to check
        # init and final states and just which states are necessary
        ket_lbls = self.ket_obj.grp_lbls
        bra_lbls = self.bra_obj.grp_lbls

        for k_lbl in ket_lbls:
            ket_spin = sefl.get_spins(k_lbl)

            for b_lbl in ket_lbls:
                bra_spin = self.get_spins(b_lbl)
                
                # if different spins manifold, skip
                if ket_spin.S != bra_spin.S:
                    continue

                # if the same object, and 

                
                # initialize the bitsi library for the calculation of 1-TDMs
                bitsi_init.init(self.bra_obj, self.ket_obj, 'tdm')
        
                # if all_final_states is true, over-ride current contents of 
                # self.final_states
                if self.all_final_states:
                    self.final_states = range(self.bra_obj.n_states())

                self.add_group('ket_states', self.ket_obj, self.init_states)
                self.add_group('bra_states', self.bra_obj, self.final_states)

                # this is main transition_list: stored by adiabatic label
                self.trans_list = self.build_pair_list('bra_states',
                                               'ket_states', 
                                                pairs=list_type)
                # the transitions, ordered by interacting irreps, is how
                # we call bitsi
                trans_list_sym  = self.build_pair_list('bra_states', 
                                                  'ket_states', 
                                                   pairs=list_type, 
                                                    sym_blk=True)

                # tdms is a vector of nmo x nmo transition densities
                self.tdms = self.build_tdms(trans_list_sym)


        # rotate tdms into spin states
        self.rotate_tdms()

        # build the multipole moments  -- easier to just do this once
        # for all transitions
        self.multipole = self.build_multipole()

        # compute oscillator strengths
        self.oscstr = self.build_osc_str()

        # print orbitals if requested
        if self.print_orbitals:
            # construct the orbitals if they don't already exist
            if len(self.nos.keys()) == 0:
                # for now, just build both the NTOs
                self.nos = self.build_ntos()
            self.export_orbitals(orb_type='nto')

        # print the summary output
        self.print_log()

        # finalize the bitsi library
        bitsi_init.finalize()

        return

    #
    def check_obj_list(obj_list):
        """
        do some sanity checks on the objects passed to run()
        """

        if len(obj_list) > 2:
            sys.exit('transition.run() passed ' + str(len(obj_list)) + 
                     ' objects. Expecting 2')

        # we will accept method objects in lieu of spin-orbit
        # objects, but will create a spin-object on the fly with
        # a single spin group
        so_list = [None,None]
        st_list = [self.final_states, self.init_states]
        for i in range(2):
            c_name = type(obj_list[i])__name__
            if c_name != 'Spinorbit':
                if c_name in params.method_objs:
                    so_list[i] = spinorbit.Spinorbit()
                    so_list[i].label         =  obj_list[i].label
                    so_list[i].couple_groups = [obj_list[i].label]
                    so_list[i].couple_states = st_list[i]
                    so_list[i].
                else:
                    sys.exit('sotransition.run() called, but '+
                             'arguments are not spinorbit '+
                             'objs or method objs'0
            else:
                so_list[i] = obj_list[i]

        # if bra and ket are same object, just need lower triangle
        # minus the diagonal elements
        if not self.same_obj(obj_list[0], obj_list[1]):
            # sanity check that orbitals and geometry are the same
            if np.any(obj_list[0].scf.orbs != obj_list[1].scf.orbs):
                sys.exit('transition moments require same bra/ket orbs')

            if (obj_list[0].mol.pymol().atom !=
                                          obj_list[1].mol.pymol().atom):
                sys.exit('transition moments require same geometry'+
                         ' and basis set')

            # bra and ket not identical objects, require all bra-ket
            # pairs of states
            list_type = 'full'

        return obj_list[0], obj_list[1], \
               obj_list[0].scf, obj_list[0].scf.mol.pymol()



    def build_tdms(self, trans_list_sym):
        """grab the TDMs from bitsi and then reshape the list of
           TDMs into a more usable format"""

        # grab the tdms
        tdm_list = mrci_1tdm.tdm(self.bra_obj, self.ket_obj, 
                                  trans_list_sym)

        # make the tdm list
        nmo    = self.bra_obj.scf.nmo
        npairs = len(self.trans_list)
        tdms   = np.zeros((nmo, nmo, npairs), dtype=float)

        for indx in range(len(self.trans_list)):
            bk_st          = self.trans_list[indx]
            [birr, bst]    = self.bra_obj.state_sym(bk_st[0])
            [kirr, kst]    = self.ket_obj.state_sym(bk_st[1])
            sym_indx       = trans_list_sym[birr][kirr].index([bst,kst])
            tdms[:,:,indx] = tdm_list[birr][kirr][:, :, sym_indx]

        return tdms
    
