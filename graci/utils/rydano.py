"""Module for handling local basis set definitions"""

import os
import sys
import numpy as np
import operator as op
import pyscf.gto as gto
import graci.core.scf as scf
import graci.io.output as output

# the KBJ exponents are ordered by angular momentum l (l=0,1,2,3),
# and run from n=1,1.5,2,2.5,...8. The starting exponent should
# be chosen by inspecting the valence basis set, and then
# 7 additional functions (8 total) are taken
kbj_exp = [[0.245645, 0.0656456, 0.0246239, 0.0112533, 0.00585838, 
            0.00334597, 0.00204842, 0.00132364, 0.000893096, 
            0.000624313, 0.000449505, 0.000331842, 0.000250295, 
            0.000192336, 0.000150229],
           [0.430082, 0.113659, 0.0423353, 0.0192542, 0.00998821, 
            0.00568936, 0.00347568, 0.00224206, 0.00151064, 
            0.00105475, 0.000758656, 0.000559584, 0.000421752, 
            0.000323875, 0.000252822],
           [0.622557, 0.163298, 0.0605402, 0.0274457, 0.0142044, 
            0.00807659, 0.00492719, 0.00317481, 0.00213712, 
            0.00149102, 0.00107174, 0.000790066, 0.00059517, 
            0.00045685, 0.000356487],
           [0.820349, 0.214006, 0.0790699, 0.0357629, 0.0184778, 
            0.010493, 0.00639492, 0.00411721, 0.00276965, 
            0.00193124, 0.00138751, 0.00102243, 0.000769943, 
            0.000590821, 0.000460901]]

n_vals     = [0.5*(n+1) for n in range(1,16)]
ang_lbls = ['s','p','d','f','g','h','i']

class Rydano():
    """
    Generate a rydberg basis for a given molecule object
    """

    def __init__(self):
        # user variables
        self.xc         = 'hf'
        self.mult       = 2
        self.charge     = 1
        self.origin     = 'center'
        self.label      = 'rydano'
        self.verbose    = False
        self.contract   = '1s1p1d'

        # class variables


    def run(self, mol):
        """
        run an SCF calculation on the corresponding doublet
        """

        # make a copy of the molecule object we are
        cation_mol = mol.copy()

        # build the PySCF molecule object 
        cation_mol.run()

        # the valence basis will be DZP for the non-H
        # atoms and DZ for H
        cation_mol.basis = {} 
        for atm in list(set(cation_mol.asym)):
            if atm.lower().strip() == 'h':
                bas_lbl = 'dz'
            else:
                bas_lbl = 'dzp'
            cation_mol.basis[atm] = bas_lbl

        # construct basis object of uncontracted primitives
        prim_obj = \
         self.construct_prim_basis('X', 2, 
                          [0.5*(n+1) for n in range(3,11)])
        cation_mol.basis['X'] = prim_obj

        # add a ghost atom at the center of 
        # of nuclear charge and place an un-contracted
        # KBJ basis at this site
        coc_xyz = cation_mol.nuc_charge_center()
        cation_mol.asym.append('X')
        cation_mol.crds = np.append(cation_mol.crds, [coc_xyz], axis=0)

        output.print_rydano_header(self)

        # construct an SCF object
        cation_scf          = scf.Scf()
        cation_scf.xc       = self.xc
        cation_scf.mult     = self.mult
        cation_scf.charge   = self.charge
        cation_scf.verbose  = True
        cation_scf.do_ao2mo = False
        cation_scf.print_orbitals = True

        # generate the cationic density
        cation_scf.run(cation_mol, None)

        print('olap matrix='+str(cation_scf.orbs@cation_scf.orbs.T), flush=True)
        # make the AO overlap matrix
        ao_s   = cation_mol.pymol().intor('int1e_ovlp')
        ao_s_s, ao_s_c = np.linalg.eigh(ao_s)

        s_thr = 1.e-12
        s_12  = np.diag([np.sqrt(s) for s in ao_s_s])
        s_m12 = np.diag([1./np.sqrt(s) if s > s_thr else 0. 
                                                     for s in ao_s_s])

        ao_sm12 = ao_s_c @ s_m12 @ ao_s_c.T
        ao_s12  = ao_s_c @ s_12  @ ao_s_c.T

        print('X@orbs@orbs.T@X.T='+str(ao_s@cation_scf.orbs@cation_scf.orbs.T),flush=True)


        # parse the AO-basis to get relevant labels
        a_ind, a_lbl, l_val, l_shl, l_lbl = \
                                     self.parse_ao_basis(cation_scf.mol)

        print('l_shl='+str(l_shl))

        # extract virtual orbitals and set occupation
        # using the expression exp[6.9*(orb_ener[i]/min(orb_ener)-1)] 
        occ = self.set_virtual_occ(cation_scf) 

        # construct a set of ortho-normal rydberg orbitals
        #orbs_r = self.project_orbitals(ao_s12, cation_scf.orbs)

        # form the density associated with the rydberg basis 
        dao = self.form_ryd_density(cation_scf.orbs, occ, a_lbl)

        # spherically average the density
        dao_sph = self.sph_average(dao, l_val[r_indx], l_shl[r_indx])

        # form NOs via projection of different ang. mom. shells
        l_max = max(l_val)
        ano_ryd = []
        for l in range(l_max+1):
            ano_ryd.append(self.make_nos(l, dao_sph, r_orbs, r_l, r_lbl))

        return

    #
    def construct_prim_basis(self, a_lbl, l_max, n_range):
        """
        construct an uncontracted basis of KBJ functions, running 
        from 0 to l_max, and for exponents specified by the values
        of n in n_range
        """
        l_lbl = ['S', 'P', 'D', 'F', 'G', 'H', 'I']

        ostr = ''
        astr = '{:<6s}{:1s}\n'
        pstr = '{:>15.7f}'+''.join(['{:14.7f}']*len(n_range))+'\n'

        nprim = len(n_range)

        for l in range(0, l_max+1):
            ostr += astr.format(a_lbl, l_lbl[l])
            args = [0]*(nprim+1)
            for ni in range(len(n_range)):
               n_indx          = n_vals.index(n_range[ni])
               args[0]         = kbj_exp[l][n_indx]
               args[1:nprim+1] = [0.]*nprim
               args[ni+1]      = 1.
               ostr += pstr.format(*args)        

        return ostr

    #
    def parse_ao_basis(self, mol):
        """
        parse_ao_basis
        """
        ao_lbls = mol.pymol().ao_labels(fmt=True)
        print('ao_lbls='+str(ao_lbls), flush=True)

        n_ao = len(ao_lbls)

        a_indx = np.zeros(n_ao, dtype=int)
        a_lbl  = ['']*n_ao
        l_val  = np.zeros(n_ao, dtype=int)
        l_lbl  = ['']*n_ao
        l_shl  = ['']*n_ao
        
        for i in range(len(ao_lbls)):
            val = ao_lbls[i].strip().split()
            a_indx[i] = int(val[0])
            a_lbl[i]  = val[1]

            # find angular momentum
            a_fnd = False
            l     = 0
            while l < len(ang_lbls):
                if ang_lbls[l] in val[2]:
                    a_fnd = True
                    break
                else:
                    l += 1 

            if not a_fnd:
                print('error parsing ao label: '+str(vals))
                sys.exit(1)

            l_val[i]  = int(l)
            l_lbl[i]  = val[2][val[2].index(ang_lbls[l])+1:].strip()
            l_shl[i]  = int(val[2][:val[2].index(ang_lbls[l])].strip())

        return a_indx, a_lbl, l_val, l_shl, l_lbl

    # 
    def set_virtual_occ(self, cation_scf):
        """
        set the occupations on the virtual orbitals
        """

        nmo = cation_scf.nmo
        occ = np.zeros(nmo, dtype=float)

        first_vrt = None
        for imo in range(nmo):
            scf_occ  = cation_scf.orb_occ[imo]
            scf_ener = cation_scf.orb_ener[imo]
            if scf_occ == 0. and scf_ener < 0:
                if first_vrt is None:
                    first_vrt = scf_ener
                occ[imo] = np.exp(6.9*(scf_ener/first_vrt - 1.))
            if scf_ener > 0:
                break

        return occ

    #
    def ryd_density(self, orbs, occ, a_lbl):
        """
        form the rydberg density
        """
        
        nmo     = orbs.shape[1]
        nao_ryd = a_lbl.count('X')
        ryd_i   = [ind for ind,val in enumerate(a_lbl) if val=='X']
        ryd_o   = [ind for ind,val in enumerate(occ) if val > 0.] 
        dao     = np.zeros((nao_ryd, nao_ryd), dtype=float)

        orbs_r = orbs[:, ryd_o]
        occ_r  = occ[ryd_o]

        for imo in range(len(ryd_o)):
            # pull out just those AOs on the target center
            mo_ryd = orbs_r[ ryd_i, imo ] 
            dao   += occ_r[ imo ]*np.outer(mo_ryd, mo_ryd)

        return dao, ryd_i

    #
    def sph_average(self, dao, l_val, l_shl):
        """
        sph_average
        """

        nmo = dao.shape[1]
        dao_sph = np.zeros(dao.shape, dtype=float)

        for bra_l in range(0, max(l_val)+1):
            bl_ind = [ind for ind,val in enumerate(l_val) 
                                                      if val == bra_l]
            bshl  = list(set(l_shl[bl_ind]))
            for bra_n in bshl:
                bs_ind = [val for ind,val in enumerate(bl_ind) 
                                         if l_shl[val] == bra_n]

                for ket_l in range(0, max(l_val)+1):
                    kl_ind = [ind for ind,val in enumerate(l_val) 
                                                      if val == ket_l]
                    kshl  = list(set(l_shl[kl_ind]))
                    for ket_n in kshl:
                        ks_ind = [val for ind,val in enumerate(kl_ind)
                                         if l_shl[val] == ket_n]
                        indx  = np.ix_(bs_ind, ks_ind)
                        nindx = len(bs_ind) * len(ks_ind)

                        ave   = np.sum(dao[indx]) / nindx
                        dao_sph[indx] = ave

        return dao_sph

    # 
    def project_orbitals(self, ao_s, orbs, r_indx):
        return


    #
    def make_nos(self, l, dao, r_orbs, l_val, l_lbl):
        """
        make_nos
        """
        #for i in range(dao.shape[0]):
        #    print('dao['+str(i)+']='+str(dao[i,:]),flush=True)

        print('l_lbl='+str(l_lbl), flush=True)
        n_l = (l_val == l).sum()
        print('l='+str(l)+' n_l = '+str(n_l),flush=True)
        ml  = list(set([l_lbl[ind] for ind,val 
                          in enumerate(l_val) if val==l]))
        print('ml='+str(ml), flush=True)
        ml_ind = [ind for ind,val in enumerate(l_val) 
                          if val==l and l_lbl[ind]==ml[0]]
        n_ml = len(ml_ind)
        print('ml_ind='+str(ml_ind), flush=True)

        for il in ml:
            ml_ind = [ind for ind,val in enumerate(l_val)
                      if val==l and l_lbl[ind]==il]
            n_ml = len(ml_ind)
            print('ml_ind='+str(ml_ind), flush=True)

            rc_ml  = np.ix_(ml_ind, ml_ind)
            dao_ml = dao[rc_ml]
            print('dao_ml='+str(dao_ml),flush=True)
            occ, nos = np.linalg.eigh(dao_ml)        

            # put in descending order
            idx = occ.argsort()[::-1]   
            occ = occ[idx]
            nos = nos[:,idx]

            print('l='+str(l)+' ml='+str(il), flush=True)
            print('occ = '+str(occ),flush=True)
            print('nos1='+str(nos[:,0]), flush=True)
            print('nos2='+str(nos[:,1]), flush=True)
            print('nos3='+str(nos[:,2]), flush=True)

        return nos


    def tofile(self):
        """
        dump the request basis set to file
        """

        return
