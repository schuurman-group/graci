"""
Module for computing HF or KS-DFT orbitals and energies
Initial pass, at least, will employ PySCF
"""

import os
import sys
import h5py
import copy
import numpy as np
import graci.core.params as params
import graci.core.orbitals as orbitals
import graci.io.output as output
import graci.utils.timing as timing
import graci.core.solvents as solvents
import graci.core.functionals as functionals
from pyscf import gto, scf, dft, symm, df
from pyscf.tools import molden
from pyscf.scf import stability

class Scf:
    """Class constructor for SCF object"""
    def __init__(self):
        # user defined input paramaters
        self.xc             = 'hf'
        self.print_orbitals = False
        self.label          = 'default'
        self.init_guess     = 'minao'
        self.diag_method    = 'cdiis'
        self.max_iter       = 50
        self.diis_start     = None
        self.diis_method    = None
        self.diis_space     = 14
        self.lvl_shift      = None
        self.verbose        = True 
        self.restart        = False
        self.damp_fac       = None
        self.mult           = 0
        self.charge         = 0
        self.grid_level     = 2
        self.x2c            = False
        self.conv_tol       = 1e-8
        self.cosmo          = False
        self.solvent        = None
        self.mol_label      = 'default'
        self.guess_label    = None
        self.direct_scf     = True
        self.chk_stable     = [False, False]
        
        # computed quantities
        self.mol          = None 
        self.nel          = None
        self.energy       = None
        self.orbs         = None
        self.orb_occ      = None
        self.orb_ener     = None
        self.orb_irrep    = []
        self.orb_sym      = []
        self.nmo          = 0
        self.nao          = 0
        self.naux         = 0
        self.rdm_ao       = None
        self.auxbasis     = None 

        # class variables
        self.valid_grids  = ['sg1_prune']
        self.valid_guess  = ['minao','1e','atom','huckel','vsap']
        self.valid_method = ['soscf', 'cdiis', 'ediis','adiis']

# Required functions #######################################################

    def copy(self):
        """create of deepcopy of self"""
        new = Scf()

        var_dict = {key:value for key,value in self.__dict__.items() 
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            if type(value).__name__ in params.valid_objs:
                setattr(new, key, value.copy())
            else:
                setattr(new, key, copy.deepcopy(value))

        return new

    @timing.timed
    def run(self, mol, guess):
        """compute the DFT energy and KS orbitals"""
        
        # set the GRaCI mol object
        self.mol = mol
       
        # the charge and multiplity of the scf sections
        # overwrites the mult/charge of the mol object
        #----------------------------------------------------

        # set the charge and multiplicity to the appropriate
        # values for this object
        self.mol.charge = self.charge

        # pyscf spin = 2*S
        self.mol.mult   = self.mult

        # build the PySCF molecule object 
        self.mol.run()

        # save the number of electrons
        self.nel = self.mol.pymol().nelectron

        # returns the PySCF molecule object needed to run SCF
        pymol = self.mol.pymol()

        # set the verbosity of the output
        if self.verbose:
            pymol.verbose = 4 
        else:
            pymol.verbose = 0

        # print standard header 
        if self.verbose:
            output.print_scf_header(self)

        # write the Cartesian coordinate to the log file
        if self.verbose:
            output.print_coords(self.mol.crds, self.mol.asym)
            
        # this is just to tell user the nature of the auxiliary basis
        if self.mol.use_df:
            # tell user what RI basis if ri_basis is not set
            self.auxbasis = self.mol.ri_basis

        # run the SCF calculation
        scf_pyscf = self.run_pyscf(pymol, guess)
        if scf_pyscf is None:
            return None       
 
        # extract orbitals, occupations and energies
        self.orbs      = scf_pyscf.mo_coeff
        
        # orb_occ are the MO occupation numbers
        self.orb_occ   = scf_pyscf.mo_occ
        # orb_ener is the array of MO energies
        self.orb_ener  = scf_pyscf.mo_energy
        # number of electrons
        self.nel       = pymol.nelectron
        # number of MOs
        self.nmo       = len(scf_pyscf.mo_occ)
        # number of auxiliary AOs
        if (self.mol.use_df and hasattr(scf_pyscf, 'with_df')):
            self.naux  = int(scf_pyscf.with_df.auxmol.nao_nr())

        # orb_sym are the symmetry labels 
        if pymol.symmetry:
            #orb_irrep are the string irrep labels - easier for check-
            #pointing if a list rather than numpy array of strings.
            self.orb_irrep = symm.label_orb_symm(pymol, pymol.irrep_name,
                                pymol.symm_orb, self.orbs).tolist()
            #orb sym are the integer indices of the irreps
            self.orb_sym   = symm.label_orb_symm(pymol, pymol.irrep_id,
                                pymol.symm_orb, self.orbs)
        else:
            self.orb_irrep = ['a'] * len(self.orbs)
            self.orb_sym   = [0] * len(self.orbs)

        # construct density matrix
        occmos = self.orbs[:,self.orb_occ>0]
        self.rdm_ao = occmos @ np.diag(self.orb_occ[self.orb_occ>0]) @ occmos.T 

        # print the summary of the output to file
        if self.verbose:
            output.print_scf_summary(self)

        # write the Molden file if requested
        if self.print_orbitals:
            orbitals.export_orbitals('mos_molden', self.mol, self.orbs,
                                     orb_sym=self.orb_irrep,
                                     orb_occ=self.orb_occ,
                                     orb_ener=self.orb_ener,
                                     orb_dir='scf.'+str(self.label), 
                                     cart=True)
        
        return self.energy

    #
    def load(self):
        """load is called on a restart to export orbitals
           and print summary/orbital information"""

        # construct the molecule object
        pymol = self.mol.pymol()

        # seit the verbosity of the output
        if self.verbose:
            pymol.verbose = 4
        else:
            pymol.verbose = 0

        # print header info even on a restart
        if self.verbose:
            output.print_scf_header(self)

        # print the summary of the output to file
        if self.verbose:
            output.print_scf_summary(self)

        # write the Molden file if requested
        if self.print_orbitals:
            orbitals.export_orbitals('mos_molden', self.mol, self.orbs, 
                                     orb_sym=self.orb_irrep, 
                                     orb_ener=self.orb_ener, 
                                     orb_dir='Scf.'+str(self.label), 
                                     cart=True)
        return

    #
    def n_states(self):
        """return number of states optmizied, i.e. 1"""
        return 1

    #
    def n_states_sym(self, irrep):
        """number of states computed"""

        if irrep == self.state_sym:
            return 1
        else:
            return 0

    #
    def energy(self):
        """return the energy of state 'state'"""

        return self.energy

    #
    def rdm1(self, basis='mo'):
        """return the density matrix for state i"""
        
        if basis == 'mo':
            return np.diag(self.orb_occ) 
        else:
            return self.rdm_ao

    #
    def natural_orbs(self, basis='ao'):
        """return hf orbitals"""

        if basis == 'ao':
            return self.orbs
        else:
            return np.identity(self.nmo)

    #
    @timing.timed
    def run_pyscf(self, pymol, guess):
        """run a PySCF HF/KSDFT computation"""

        # catch the use of XC functional aliases
        try:
            self.xc = functionals.aliases[self.xc.lower()]
        except:
            pass

        # function handle string
        if self.xc == 'hf':
            class_str  = 'scf'
            method_str = 'HF'
        else:
            class_str  = 'dft'
            method_str = 'KS'

        if self.mol.mult == 1:
            method_str = 'R'+method_str
        else:
            method_str = 'RO'+method_str

        if self.x2c:
            rel_str = '.x2c()'
        else:
            rel_str = ''

        if self.mol.use_df:
            df_str = '.density_fit(auxbasis = self.mol.ri_basis)'
        else:
            df_str=''

        if self.diag_method.lower() == 'soscf':
            diag_str = '.newton()'
        else:
            diag_str = ''

        if self.cosmo:
            func_str = 'pymol.'+method_str+'(xc=\''+self.xc+'\')' \
                + rel_str + diag_str + df_str + '.DDCOSMO()'
        else:
            func_str = class_str+'.'+method_str \
                +'(pymol)' + rel_str + diag_str + df_str

        # instantiate the scf/dft class object        
        mf = eval(func_str)

        # COSMO solvent dielectric constant
        if self.cosmo:
            mf.with_solvent.eps = solvents.dielectric[self.solvent]
            
        if self.xc != 'hf':
            # set the XC functional
            mf.xc = self.xc
        
            # set the quadrature grids
            mf.grids.level = self.grid_level
            #mf.grids.prune = dft.sg1_prune
            mf.grids.prune = dft.nwchem_prune

        # convergence threshold
        mf.conv_tol = self.conv_tol
            
        # set the integral file names
        #if hasattr(mf, 'with_df'):
        #    aoint_file = '2e_eri_'+str(self.label).strip()+'.h5'
        #    if not self.mol.use_df: 
        #        mf.with_df = None

        # set the dynamic level shift, if request
        #if self.lvl_shift is not None:
        #    scf.addons.dynamic_level_shift_(mf,
        #                              factor=float(self.lvl_shift))
        if self.lvl_shift is not None:
            #mf.level_shift = self.lvl_shift    
            mf.conv_check  = False
            scf.addons.dynamic_level_shift_(mf, factor=self.lvl_shift)

        # set the maximum number of iterations
        mf.max_cycle = self.max_iter

        # default chkfile name based on object label name
        #mf.chkfile = self.make_chkfile_name(self.label)
      
        # how to define initial guess orbitals
        if guess == None:

            dm = None

            # check if a restart file exists, if so, use same chkfile
            # name
            init_mos = self.make_chkfile_name(self.init_guess)
            
            # if restart file exists, use orbitals from that as
            # a guess
            if os.path.exists(init_mos):
                #dm = dft.roks.init_guess_by_chkfile(pymol, init_mos)
                mf.init_guess = 'chkfile'
                mf.chkfile = init_mos

            # else use user-specified init guess algorithm
            elif self.init_guess in self.valid_guess:
                mf.init_guess = self.init_guess

        # use scf object passed to run
        else:
            #dm = guess.rdm_ao 
            dm = self.guess_dm(guess)

        # SCF CONVERGENCE OPTIONS
        if self.diag_method.lower() in self.valid_method \
              and self.diag_method.lower() != 'soscf':
            mf.DIIS = eval('scf.'+self.diag_method.upper())
        if self.diis_space > 0:
            mf.diis_space = self.diis_space
        # when to start diis cycle
        if self.damp_fac is not None:
            mf.damp = self.damp_fac
        if self.diis_start is not None:
            mf.diis_start_cycle = self.diis_start

        # if this is an atom: preserve spherical symmetry
        #if self.mult != 1: 
        #    mf = scf.addons.frac_occ(mf)
        mf.direct_scf = self.direct_scf 

        # run the scf computation
        self.energy = mf.kernel(dm0=dm)
       
        # if not converged, kill things
        if not mf.converged:
            return None

        # check stability
        if any(self.chk_stable):
            mo_i, mo_e, stbl_i, stbl_e = mf.stability(
                                            return_status=True,
                                            external=self.chk_stable[1])
            # this should be cleaned up: if we're not checking for 
            # external instabilities, set stbl_e to True so mo_e is 
            # never taken
            if not self.chk_stable[1]:
                stbl_e = True

            chk_iter = 1
            while (not all([stbl_i, stbl_e]) and chk_iter <= 3):
                if self.verbose:
                    output.print_message('Orbital instability found.' + 
                              ' Re-optimizing, attempt '+str(chk_iter))

                if not stbl_i:
                    new_mo = mo_i
                else:
                    new_mo = mo_e[0]

                dm1 = mf.make_rdm1(new_mo, mf.mo_occ)
                mf  = mf.run(dm1)
                self.energy = mf.e_tot
                mo_i, mo_e, stbl_i, stbl_e = mf.stability(
                                            return_status=True, 
                                            external=self.chk_stable[1])

                if not self.chk_stable[1]:
                    stbl_e = True

                chk_iter += 1
          
            # if stable orbitals not found, return None
            if not all([stbl_i, stbl_e]):
                if self.verbose:
                    output.print_message('Stable Orbitals not found '
                     'after ' + str(chk_iter) + ' attempts. Continuing')
                #return None 

        # MO phase convention: positive dominant coefficients
        # for degenerate coefficients, pick the first occurrence
        # N.B. this is essential for diabatisation runs
        nmo  = mf.mo_coeff.shape[1]
        nao  = mf.mo_coeff.shape[0]
        imax = [np.argmax(np.abs(mf.mo_coeff[:,i])) for i in range(nmo)]
        for i in range(nmo):
            coeff_max = mf.mo_coeff[imax[i], i]
            diff = np.abs(mf.mo_coeff[:, i]) - abs(coeff_max)
            indx = [1 if abs(diff[i]) < 1e-6 else 0
                    for i in range(nao)].index(1)
            if mf.mo_coeff[indx, i] < 0.:
                mf.mo_coeff[:, i] *= -1.
        
        return mf

    #
    def guess_dm(self, guess):
        """
        returns the density matrix of the guess scf object
        projected onto the current AO basis
        
        adapted from the PySCF uhf.init_guess_by_chkfile function
        """
        
        # PySCF mol objects
        mol  = self.mol.mol_obj
        mol0 = guess.mol.mol_obj
       
        # current AO overlaps
        s = mol.intor_symmetric('int1e_ovlp')
    
        # guess MOs and occupations
        mo0  = guess.orbs
        occ0 = guess.orb_occ

        # projection
        mo       = scf.addons.project_mo_nr2nr(mol0, mo0, mol)
        norm     = np.einsum('pi,pi->i', mo.conj(), s.dot(mo))
        mo_coeff = mo / np.sqrt(norm)

        # construct the density matrix
        mo_occa = (occ0 > 1e-8).astype(np.double)
        mo_occb = occ0 - mo_occa
        dm      = scf.uhf.make_rdm1([mo_coeff, mo_coeff],
                                    [mo_occa, mo_occb])

        return dm[0] + dm[1]

        return dm   


    #
    def make_chkfile_name(self, label):
        """
        return a checkfile name based on object label
        """
        return 'Chkfile.Scf.'+str(label)
        
    #
    def mo_overlaps(self, other):
        """
        returns the overlaps with the MOs of another CI
        class object, 'other'
        """

        # mol objects
        mol0 = other.mol.pymol().copy()
        mol  = self.mol.pymol().copy()

        # AO overlaps
        smo  = gto.intor_cross('int1e_ovlp', mol0, mol)

        # AO-to-MO transformation
        smo  = np.matmul(np.matmul(other.orbs.T, smo), self.orbs)

        return smo

