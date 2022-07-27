"""Module for handling local basis set definitions"""

import os
import sys
import pyscf.gto as gto
import graci.core.scf as scf

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

n_vals = [n + 0.5 for n in range(0.5,7.5,0.5)]

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
        # to modify and use that to generate Rydberg
        # basis
        cation_mol = mol.copy()

        # build the PySCF molecule object 
        cation_mol.run()

        # the valence basis will be DZP for the non-H
        # atoms and DZ for H
        cation_mol.basis = []
        for atm in list(set(cation_mol.asym)):
            if atm.lower().strip() == 'h':
                bas_lbl = 'dz'
            else:
                bas_lbl = 'dzp'
            cation_mol.basis.extend([atm bas_lbl])

        # construct basis object of uncontracted primitives
        prim_obj = \
         self.construct_prim_basis(albl, 2, 
                                      [n+0.5 for n in range(1.5,5,0.5)])
        cation_mol.basis.extend(['X' prim_obj])

        # add a ghost atom at the center of 
        # of nuclear charge and place an un-contracted
        # KBJ basis at this site
        coc_xyz = cation_mol.nuc_charge_center()
        cation_mol.asym.append('X')
        cation_mol.cart.append(coc_xyz)

        output.print_rydano_header()

        # construct an SCF object
        ryd_scf          = scf.Scf()
        ryd_scf.xc       = self.xc
        ryd_scf.mult     = self.mult
        ryd_scf.charge   = self.charge
        ryd_scf.verbose  = True
        ryd_scf.do_ao2mo = False

        # generate the cationic density
        ryd_scf.run(mol, None)

        # 
 

        return

    #
    def construct_prim_basis(a_lbl, l_max, n_range):
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
            for ni in n_range:
               args[0]         = kbj_exp[l, n_vals.index(ni)]
               args[1:nprim+1] = [0.]*nprim
               args[ni+1]      = 1.
               ostr += pstr.format(*args)        

        return ostr

    def tofile(self):
        """
        dump the request basis set to file
        """

        return
