"""
module for compute moments of the dipole operator
for a given electronic state
"""
import numpy as np
import graci.core.libs as libs
import graci.citools.mrci_1tdm as mrci_1tdm
import graci.io.output as output
import sys as sys

class Transition:
    """Transition class for determing transition moments, right now
       it will only accept calculations for which the 'Molecule' and 'Scf' 
       objects are the same"""
    def __init__(self):
        # user defined quanties 
        self.init_states      = None
        self.final_states     = None
        self.istate_array     = None
        self.fstate_array     = None
        self.all_final_states = False
        self.init_method      = None
        self.final_method     = None

        # global variables
        self.tdm          = None

        # orbital quantites stored in AO basis 
        # the natural difference orbitals
        self.ndo          = None
        # weights for each NDO
        self.ndo_wt       = None
        # the natural transition orbitals
        self.nto          = None
        # weights for each NTO
        self.nto_wt       = None
        # comprehensive list of all transition pairs
        self.trans_list   = None
        # list of initial states
        self.ket_list     = None
        # the transition dipole moments
        self.trans_mom    = None
        # the oscillator strengths
        self.osc_str      = None
        # this is the default label: the class name
        self.label        = 'transition'

    #-----------------------------------------------------------

    #
    def name(self):
        """ return the name of the class object as a string"""
        return 'transition'

    #
    def set_final_method(self, final_method):
        """set the method used to compute the ket state(s)"""
        self.final_method = final_method
        return

    #
    def set_init_method(self, init_method):
        """set the method used to compute the bra state(s)"""
        self.init_method = init_method
        return

    #
    def run(self):
        """return the transition dipole moments between the bra and
           ket states. If b_state and k_state are None, assume 
           transitions from all states in method object should be 
           used."""

        mol_bra = self.final_method.mol
        mol_ket = self.final_method.mol

        scf_bra = self.final_method.scf
        scf_ket = self.final_method.scf
        
        # sanity check that orbitals and geometry are the same
        if np.any(scf_bra.orbs != scf_ket.orbs):
            sys.exit('transition moments require same bra/ket orbs')

        if mol_bra.pymol().atom != mol_ket.pymol().atom:
            sys.exit('transition moments require same geometry/basis')

        # initialize the bitsi library
        libs.init_bitsi(self)

        # compute the dipole moment integrals
        mu_ao = mol_ket.pymol().intor('int1e_r')
        
        # contract with the MOs
        mu_mo = [np.matmul(np.matmul(scf_ket.orbs.T, mu_ao[i]), 
                                    scf_ket.orbs) for i in range(3)]

        # construct the trans_list array
        # currently store them as [initial state, final state]
        self.build_trans_list()

        # grab the transition density matrices
        self.tdms = mrci_1tdm.tdm(self)

        # compute transition dipole moments and oscillator strengths
        self.build_trans_moments(mu_mo)

        # print natural differnence and natural transition orbitals 
        # for each pair
        self.build_natural_orbs()

        # print the transition and difference densities and print
        # excitations energies/osc. strengths to file
        self.print_orbitals()

        # print the summary output
        self.print_log()

        return

    #
    def tdm1(self, b_irr, b_st, k_irr, k_st):
        """return the tdm for the bra and ket state indicated"""

        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st,k_st])
        except:
            return None

        return self.tdms[b_irr][k_irr][:,:,ind]
 
    #
    def trans_dipole(self, istate, fstate):
        """return the transition dipole for pairs of
           states ordered by energy"""
        
        [b_irr, b_st] = self.final_method.state(istate)
        [k_irr, k_st] = self.init_method.state(fstate)

        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st, k_st])
        except:
            return None

        return self.trans_dipole[b_irr][k_irr][ind][:]

    #
    def trans_dipole_sym(self, b_irr, b_st, k_irr, k_st):
        """return transition dipole matrix element"""

        # return the transition dipole vector given the irreps and
        # state indices
        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st, k_st])
        except:
            return None

        return self.trans_dipole[bra_irrep][ket_irrep][ind][:]

    #
    def osc_strength(self, istate, fstate):
        """return the transition dipole for pairs of 
           states ordered by energy"""

        [b_irr, b_st] = self.final_method.state(istate)
        [k_irr, k_st] = self.init_method.state(fstate)

        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st, k_st])
        except:
            return None

        return self.osc_str[b_irr][k_irr][ind][:]

    # 
    def osc_strength_sym(self,  b_irr, b_st, k_irr, k_st):
        """return the transition dipole moments for all pairs of
           bra_states and ket_states. If one or both are not specified,
           return all transition moments"""

        # return the transition dipole vector given the irreps and
        # state indices
        try:
            ind = self.trans_list[b_irr][k_irr].index([b_st, k_st])
        except:
            None

        return self.osc_str[bra_irrep][ket_irrep][ind][:]

#--------------------------------------------------------------------------

    def build_trans_list(self):
        """built an array of initial states, ordered by irrep"""

        nirr = self.final_method.mol.n_irrep()

        # determine what initial states to include
        # -------------------------------------------------

        self.ket_list = [[] for irr in range(nirr)]
        ket_avail = self.final_method.n_states()

        # iterate through the initial states, adding roots
        # to various irreps as needed
        if self.init_states is not None:
            # iterate through the final states, adding roots
            # to various irreps as nee
            for state in self.init_states:
                istate = self.final_method.state_n(state)
                self.ket_list[istate[0]].append(istate[1])
        
        # use the istate array (format = irrep.state) to select
        # individual states
        elif self.istate_array:
            for state in self.fstate_array:
                irrep = int(state)
                st    = int(10*state) - 10*irrep
                if st <= st_avail[irrep]:
                    self.ket_list[irrep].append(st)

        else:
            sys.exit('no initial states selected in Transition')

        # determine what final states to include
        # -------------------------------------------------

        bra_list = [[] for irr in range(nirr)]
        bra_avail = self.final_method.n_states()

        # if all_final_states is true, that is overrides anything else
        if self.all_final_states:
            bra_list = [[st for st in range(bra_avail[irr])]
                            for irr in range(nirr)]

        # if the number of final roots is given, pick the lowest n roots
        # by energy
        elif self.final_states is not None:
            # iterate through the final states, adding roots
            # to various irreps as nee
            for state in self.final_states:
                fstate = self.final_method.state_n(state)
                bra_list[fstate[0]].append(fstate[1])

        # use the fstate array (format = irrep.state) to select
        # individual states
        elif self.fstate_array is not None:
            for state in self.fstate_array:
                irrep = int(state)
                st    = int(10*state) - 10*irrep
                if st <= bra_avail[irrep]:
                    bra_list[irrep].append(st)

        else:
            sys.exit('no final states selected in Transition')

        # create the pair list
        self.trans_list = [[[] for irr_bra in range(nirr)]
                               for irr_ket in range(nirr)]
        self.state_list = []
        self.ener_list  = []
        for b_irr in range(nirr):
            for k_irr in range(nirr):
                for b_st in bra_list[b_irr]:
                    for k_st in self.ket_list[k_irr]:
                        if [k_irr, k_st] not in self.state_list:
                            self.state_list.append(
                                [k_irr, k_st])
                        if [b_irr, b_st] not in self.state_list:
                            self.state_list.append(
                                [b_irr, b_st])
                        if b_irr != k_irr or b_st != k_st:
                            self.trans_list[b_irr][k_irr].append([b_st, 
                                                                  k_st])

        return

    #
    def build_trans_moments(self, mu_mo):
        """builds the transition moment vector based on the based tdms"""

        # loop over the state pairs and construct transition 
        # dipole moments
        nirr_bra = self.final_method.n_irrep()
        nirr_ket = self.init_method.n_irrep()

        self.trans_dipole = [[[] for irr_bra in range(nirr_bra)]
                                 for irr_ket in range(nirr_ket)]

        self.osc_str      = [[[] for irr_bra in range(nirr_bra)]
                                 for irr_ket in range(nirr_ket)]

        for b_irr in range(nirr_bra):
            for k_irr in range(nirr_ket):
                for bra_ket in range(len(self.trans_list[b_irr][k_irr])):

                    states  = self.trans_list[b_irr][k_irr][bra_ket]
                    b_ener  = self.final_method.energy(b_irr, states[0])
                    k_ener  = self.init_method.energy(k_irr, states[1])

                    tdm_mat = self.tdms[b_irr][k_irr][:,:,bra_ket]
                    td_vec  = np.array([np.sum(tdm_mat*mu_mo[x])
                              for x in range(3)], dtype=float)
                    osc_vec = (2./3.)*(b_ener-k_ener)*td_vec**2

                    self.trans_dipole[b_irr][k_irr].append(td_vec)
                    self.osc_str[b_irr][k_irr].append(osc_vec)

        return

    #
    def build_natural_orbs(self):
        """build the natural difference orbitals and natural transition
           orbitals"""

        nirr_bra = self.final_method.n_irrep()
        nirr_ket = self.init_method.n_irrep()
        # should probably check bra and ket MOs are the same
        mos      = self.final_method.scf.orbs
        nmo      = self.final_method.scf.nmo

        self.nto    = [[[] for irr_bra in range(nirr_bra)]
                           for irr_ket in range(nirr_ket)]
        self.ndo    = [[[] for irr_bra in range(nirr_bra)]
                           for irr_ket in range(nirr_ket)]
        self.nto_wt = [[[] for irr_bra in range(nirr_bra)]
                           for irr_ket in range(nirr_ket)]
        self.ndo_wt = [[[] for irr_bra in range(nirr_bra)]
                           for irr_ket in range(nirr_ket)]

        for b_irr in range(nirr_bra):
            for k_irr in range(nirr_ket):
                npairs_irr = len(self.trans_list[b_irr][k_irr])
                for bra_ket in range(npairs_irr):
                    states = self.trans_list[b_irr][k_irr][bra_ket]
                    bra_st = states[0]
                    ket_st = states[1]

                    tdm = self.tdms[b_irr][k_irr][:,:,bra_ket]
                    rdm_bra = self.final_method.rdm1(b_irr, bra_st)
                    rdm_ket = self.final_method.rdm1(k_irr, ket_st)

                    # first perform SVD of 1TDMs to get hole and
                    # particle orbitals and weights
                    u_mo, s, v_mo = np.linalg.svd(tdm)
                    
                    if np.amin(s) < -1.e-3:
                        sys.exit('ERROR: build_ntos -> wt < 0')
                    else:
                        s = np.absolute(s)

                    # sort the NTO amplitudes by decreasing magnitude
                    ordr     = np.flip(np.argsort(np.sqrt(s)))
                    u_ao     = np.matmul(mos, u_mo)
                    v_ao     = np.matmul(mos, v_mo.T)
                    u_ao_srt = np.array([u_ao[:,ordr[i]] 
                                for i in range(nmo)], dtype=float).T
                    v_ao_srt = np.array([v_ao[:,ordr[i]]
                                for i in range(nmo)], dtype=float).T
                    s_srt    = np.array([np.sqrt(s[ordr[i]]) 
                                for i in range(nmo)], dtype=float)

                    self.nto[b_irr][k_irr].append([u_ao_srt, v_ao_srt])
                    self.nto_wt[b_irr][k_irr].append(s_srt)

                    # next construct the natural difference densities
                    delta       = rdm_bra - rdm_ket
                    wts, ndo_mo = np.linalg.eigh(delta)
                    ndo_ao      = np.matmul(mos, ndo_mo)

                    # sort NDO wts by decreasing magnitude
                    ordr        = np.flip(np.argsort(wts))
                    self.ndo[b_irr][k_irr].append(ndo_ao)
                    self.ndo_wt[b_irr][k_irr].append(wts)

        return

    #
    def print_log(self):
        """print summary output to log file"""
        # for each initial state, write out a table of 
        # oscillator strenghts and transition dipole moments
        nirr_bra = self.final_method.n_irrep()
        nirr_ket = self.init_method.n_irrep()
      
        k_irrlbl = self.init_method.mol.irreplbl
        b_irrlbl = self.final_method.mol.irreplbl

        for kirr in range(nirr_ket):
            for kst in range(len(self.ket_list[kirr])):
                st_table  = []
                td_table  = []
                osc_table = []
                for birr in range(nirr_bra):
                    npairs = len(self.trans_list[birr][kirr])
                    for bra_ket in range(npairs):

                        st_pair = self.trans_list[birr][kirr][bra_ket]

                        # if bra state == ket state, skip
                        if [birr, st_pair[0]] == [kirr, kst]:
                            continue
                        # if this transition pair has correct initial
                        # state, add to the list
                        if st_pair[1] == kst:
                            st_table.append([birr, st_pair[0]])
                            td_table.append(self.trans_list[birr][kirr]
                                                         [bra_ket])
                            osc_table.append(self.osc_str[birr][kirr]
                                                         [bra_ket])

                # states currently ordered by symmetry, re-order
                # by state energy
                nfinal   = len(st_table)
                st_indx  = np.array(
                            [self.final_method.state_index(*st_table[i])
                            for i in range(nfinal)])
                ordr     = np.argsort(st_indx)

                # adiabatic state labels + sym labels (state 
                # indices should start 1 (not zero)
                k_sym_st = [self.init_method.state_index(kirr, kst)+1,
                              k_irrlbl[kirr]]
                b_sym_st = [[st_indx[ordr[i]]+1,
                              b_irrlbl[st_table[ordr[i]][0]]]
                             for i in range(nfinal)]

                # sort the exc energies, tran dipoles and oscillator
                # strengths by increasing excitation energy
                exc_ordr = [
                           self.final_method.energy_n(b_sym_st[i][0]-1)-
                           self.init_method.energy_n(k_sym_st[0]-1)
                           for i in range(nfinal)]
                td_ordr  = [td_table[ordr[i]] for i in range(nfinal)]
                osc_ordr = [osc_table[ordr[i]] for i in range(nfinal)]

                # currently print table to output file. Perhaps we
                # also write to file for spectral simulations at a 
                # later date.
                output.print_transition_table(k_sym_st, b_sym_st,
                                              exc_ordr, osc_ordr)



        return

    #
    def print_orbitals(self):
        """print the natural transition orbitals and difference
           densities to file"""

        k_irrlbl = self.init_method.mol.irreplbl
        b_irrlbl = self.final_method.mol.irreplbl

        nirr_ket = len(k_irrlbl)
        nirr_bra = len(b_irrlbl)

        for k_irr in range(nirr_ket):
            for b_irr in range(nirr_bra):
                npairs = len(self.trans_list[b_irr][k_irr])
                for bra_ket in range(npairs):
                    # label the states by adiabatic label + symmetry
                    # irrep
                    b_st = self.trans_list[b_irr][k_irr][bra_ket][0]
                    k_st = self.trans_list[b_irr][k_irr][bra_ket][1]

                    bra_st_n = self.final_method.state_index(b_irr,b_st)
                    ket_st_n = self.init_method.state_index(k_irr,k_st)

                    if bra_st_n is None or ket_st_n is None:
                        sys.exit(' bra,ket=['+str(b_st)+','+str(k_st)+
                                 '] not found')
                        
                    # print out the NTOs
                    output.print_nos_wt_molden(
                        'nto', self.final_method.mol,
                        bra_st_n+1, b_irrlbl[b_irr],
                        ket_st_n+1, k_irrlbl[k_irr],
                        self.nto[b_irr][k_irr][bra_ket][0],
                        self.nto_wt[b_irr][k_irr][bra_ket])

                    # print out the NDOs
                    output.print_nos_wt_molden(
                        'ndo', self.final_method.mol,
                        bra_st_n+1, b_irrlbl[b_irr],
                        ket_st_n+1, k_irrlbl[k_irr],
                        self.ndo[b_irr][k_irr][bra_ket],
                        self.ndo_wt[b_irr][k_irr][bra_ket])

        return
