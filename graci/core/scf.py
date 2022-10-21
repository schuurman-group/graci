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
from pyscf import gto, scf, dft, symm, df
from pyscf.tools import molden


class Scf:
    """Class constructor for SCF object"""
    def __init__(self):
        # user defined input paramaters
        self.xc             = 'hf'
        self.print_orbitals = False
        self.label          = 'default'
        self.verbose        = True 
        self.restart        = False
        self.mult           = 0
        self.charge         = 0
        self.grid_level     = 2
        self.conv_tol       = 1e-8
        self.cosmo          = False
        self.solvent        = None
        self.guess_label    = None
        
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
            
        # set the file names based on class label
        # save integrals -- tie them to the scf object for
        self.moint_1e     = '1e_'+str(self.label).strip()+'.h5'
        self.moint_2e_eri = '2e_eri_'+str(self.label).strip()+'.h5'

        # this is just to tell user the nature of the auxiliary basis
        if self.mol.use_df:
            # tell user what RI basis if ri_basis is not set
            self.auxbasis = self.mol.ri_basis

        # run the SCF calculation
        scf_pyscf = self.run_pyscf(pymol, guess)
        
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
        self.rdm_ao = np.zeros((self.nmo, self.nmo), dtype=float)
        for i in range(self.nmo):
            self.rdm_ao += self.orb_occ[i]*np.outer(self.orbs[:,i],\
                                                    self.orbs[:,i])

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
        
        return

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
        if self.mol.use_df:
            df_str = '.density_fit(auxbasis = self.mol.ri_basis)'
        else:
            df_str=''
        if self.cosmo:
            func_str = 'pymol.'+method_str+'(xc=\''+self.xc+'\')' \
                +df_str+ '.DDCOSMO()'
        else:
            func_str = class_str+'.'+method_str \
                +'(pymol)'+df_str

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
            mf.grids.prune = dft.sg1_prune

        # convergence threshold
        mf.conv_tol = self.conv_tol
            
        # set the name of the DF-tensor
        if hasattr(mf, 'with_df'):
            if self.mol.use_df:
                eri_name = '2e_eri_'+str(self.label).strip()+'.h5'
                mf.with_df._cderi_to_save = eri_name
            else:
                mf.with_df = None

        # optional: override the initial density matrix with
        # that of a previous SCF calculation
        if guess is None:
            dm = None
        else:
            dm = self.guess_dm(guess)

        # run the scf computation
        self.energy = mf.kernel(dm0=dm)
        
        # if not converged, kill things
        if not mf.converged:
            sys.exit('Reference SCF computation did not converge.')
    
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


