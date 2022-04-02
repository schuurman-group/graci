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
import graci.interaction.transition as transition
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

        # master transition list for the SO-coupled states
        if self.same_obj(self.bra_obj, self.ket_obj):
            list_type = 'nodiag'
        else:
            list_type = 'full'
        self.add_group('ket_states', self.ket_obj, self.init_states)
        self.add_group('bra_states', self.bra_obj, self.final_states)
        self.trans_list = self.build_pair_list('bra_states',
                                               'ket_states',
                                                pairs=list_type)

        # initialize the (spin-orbit coupled) transition density
        # matrices
        self.tdms = np.zeros((self.scf.nmo, self.scf.nmo, 
                              len(self.trans_list)), dtype=np.cdouble)

        # loop over all pairs of spin free states in both the bra
        # and ket objects. In the future, may be prudent to check
        # init and final states and just which states are necessary
        for k_lbl in self.ket_obj.grp_lbls:
            ket_so    = 'ket_'+str(k_lbl)
            ket_spin  = self.ket_obj.get_spins(k_lbl)
            ket_ci    = self.ket_obj.get_obj(k_lbl) 
            self.add_group(ket_so, 
                           ket_ci, 
                           self.ket_obj.get_states(k_lbl))

            for b_lbl in self.bra_obj.grp_lbls:
                bra_so   = 'bra_'+str(b_lbl)
                bra_spin = self.bra_obj.get_spins(b_lbl)
                bra_ci   = self.bra_obj.get_obj(b_lbl)
                self.add_group(bra_so, 
                               bra_ci, 
                               self.bra_obj.get_states(b_lbl))

                # if different spin manifold, tdm contribution is zero
                if ket_spin.mult != bra_spin.mult:
                    continue

                # if the same object, only build lower diagonal
                if self.same_obj(ket_ci, bra_ci): 
                    pair_type = 'nodiag'
                else:
                    pair_type = 'full'
                
                # initialize the bitsi library for the calculation of 1-TDMs
                bitsi_init.init(bra_ci, ket_ci, 'tdm')
        
                # this is main transition_list: stored by adiabatic label
                blk_lst = self.build_pair_list(bra_so,
                                               ket_so, 
                                               pairs=pair_type)
                # the transitions, ordered by interacting irreps, is how
                # we call bitsi
                blk_lst_sym  = self.build_pair_list(bra_so, 
                                                    ket_so, 
                                                    pairs=pair_type, 
                                                    sym_blk=True)

                # tdms is a vector of nmo x nmo transition densities
                tdms = self.build_tdms(bra_ci, ket_ci, 
                                       blk_lst, 
                                       blk_lst_sym)

                # rotate tdms into spin states
                self.rotate_tdms(b_lbl, k_lbl, blk_lst, tdms)

                # finalize the bitsi library
                bitsi_init.finalize()

        print('tdm[1]='+str(self.tdms[:,:,0]))


        # build the multipole moments  -- easier to just do this once
        # for all transitions
        self.multipole = self.build_multipole()

        # compute oscillator strengths
        self.oscstr = self.build_osc_str()

        # print orbitals if requested
        if self.print_orbitals:
            # construct the orbitals if they don't already exist
            if len(self.nos.keys()) == 0:
                self.nos = self.build_ntos()
            self.export_orbitals(orb_type='nto')

        # print the summary output
        self.print_log()

        del(tdms)
        return

    #
    def check_obj_list(self, obj_list):
        """
        do some sanity checks on the objects passed to run()
        """
        if len(obj_list) > 2:
            sys.exit('transition.run() passed ' + str(len(obj_list)) +
                     ' objects. Expecting 2')

        # if bra and ket are same spin-object, just unique bra/ket 
        # block pairs, and, the lower diagonal of the diagonal blocks
        # minus the diagonal elements
        if not self.same_obj(obj_list[0], obj_list[1]):
            for igrp in obj_list[0].grp_lbls:
                scf_b = obj_list[0].get_obj(igrp).scf
                for jgrp in obj_list[1].grp_lbls:
                    scf_k = obj_list[1].get_obj(jgrp).scf

                    # sanity check that orbitals and geometry are 
                    # the same
                    if np.any(scf_b.orbs != scf_k.orbs):
                        sys.exit('transition moments require same '+
                                 ' bra/ket orbitalss')

                    if (scf_b.mol.pymol().atom !=
                              scf_k.mol.pymol().atom):
                        sys.exit('transition moments require same '+
                                 ' geometry and basis set')

        else:
            # this is hacky and should be fixed
            grp0 = obj_list[0].grp_lbls[0]
            scf_b = obj_list[0].get_obj(grp0).scf

        return obj_list[0], obj_list[1], scf_b, scf_b.mol.pymol()

    #
    @timing.timed
    def rotate_tdms(self, bra_lbl, ket_lbl, blk_lst, tdms):
        """
        Rotate the spin-free tdms to the spin-orbit states
        """
        print('bra_lbl='+str(bra_lbl))
        print('ket_lbl='+str(ket_lbl))
        print('blk_lst='+str(blk_lst))
        print('shape tdms='+str(tdms.shape))

        bra_spin = self.bra_obj.get_spins(bra_lbl)
        ket_spin = self.ket_obj.get_spins(ket_lbl)

        # run through trans_list and contribute each #
        for pair in blk_lst:
            ind = blk_lst.index(pair)
            tdm = tdms[:, :, ind] 

            for sopair in self.trans_list:
                so_ind = self.trans_list.index(sopair)

                for ms_b in bra_spin.M:
                    b_ind = self.bra_obj.soc_index(bra_lbl, pair[0], ms_b)
                    b_cf = np.conj(self.bra_obj.soc_state(sopair[0])[b_ind])

                    for ms_k in ket_spin.M:
                        k_ind = self.ket_obj.soc_index(ket_lbl, pair[1], ms_k)
                        k_cf = self.ket_obj.soc_state(sopair[1])[k_ind]
                        self.tdms[:, :, so_ind] += tdm * b_cf * k_cf

        return

