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
from pyscf import gto, scf, dft, symm, ao2mo, df
from pyscf.tools import molden


class Scf:
    """Class constructor for SCF object"""
    def __init__(self):
        # user defined input paramaters
        self.xc             = 'hf'
        self.print_orbitals = False
        self.label          = 'Scf'
        self.verbose        = True 
        self.restart        = False
        self.mult           = 0
        self.charge         = 0
        self.grid_level     = 2
        
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
        self.moint_1e     = None
        self.moint_2e_eri = None
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
    def run(self, mol):
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

        # set the file names based on class label
        # save integrals -- tie them to the scf object for
        self.moint_1e     = '1e_'+str(self.label).strip()+'.h5'
        self.moint_2e_eri = '2e_eri_'+str(self.label).strip()+'.h5'

        # this is just to tell user the nature of the auxiliary basis
        if self.mol.use_df:
            # tell user what RI basis if ri_basis is not set
            self.auxbasis = self.mol.ri_basis
            if self.auxbasis is None:
                bname = {atm :
                        gto.basis._format_basis_name(self.mol.basis[atm])
                                        for atm in self.mol.basis.keys()}
                self.auxbasis = dict()
                for atm in bname.keys():
                    if bname[atm] in df.addons.DEFAULT_AUXBASIS.keys():
                        self.auxbasis[atm] = df.addons.DEFAULT_AUXBASIS[bname[atm]][0]
                    else:
                        self.auxbasis[atm] = 'even-tempered'

        scf_pyscf = self.run_pyscf(pymol)

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
            self.naux      = int(scf_pyscf.with_df.auxmol.nao_nr())

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

        # perform AO -> MO transformation
        self.ao_to_mo(pymol, self.orbs)

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
        """load is called on a restart to generate integrals and do 
           ao2mo transformation"""

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

        # perform AO -> MO transformation
        self.ao_to_mo(pymol, self.orbs)

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
    def run_pyscf(self, pymol):
        """run a PYSCF HF/KSDFT computation"""

        # if self.xc='hf', use canonical hf orbitals
        if self.xc == 'hf':
            if self.mol.use_df:
                if self.mol.mult == 0.:
                    mf = scf.RHF(pymol).density_fit( auxbasis =
                                            self.mol.ri_basis)
                else:
                    mf = scf.ROHF(pymol).density_fit( auxbasis =
                                            self.mol.ri_basis)
            else:
                if self.mol.mult == 0.:
                    mf = scf.RHF(pymol)
                else:
                    mf = scf.ROHF(pymol)

        # this is a DFT computation
        else:
            if self.mol.use_df:
                if self.mol.mult == 0.:
                    mf = dft.RKS(pymol).density_fit( auxbasis =
                                               self.mol.ri_basis)
                else:
                    mf = dft.ROKS(pymol).density_fit( 
                                   auxbasis = self.mol.ri_basis)
            else:
                if self.mol.mult == 0.:
                    mf = dft.RKS(pymol)
                else:
                    mf = dft.ROKS(pymol)
            # set the XC functional to BHLYP
            mf.xc = self.xc

            # set the DFT quadrature grids
            mf.grids.level = self.grid_level
            mf.grids.prune = dft.sg1_prune
            
        # if using density-fitting, set the name of the DF-tensor
        if hasattr(mf, 'with_df'):
             if self.mol.use_df:
                 mf.with_df._cderi_to_save = self.moint_2e_eri
             else:
                 mf.with_df = None

        # run the dft computation
        self.energy = mf.kernel()

        # if not converged, kill things
        if not mf.converged:
            sys.exit('Reference SCF computation did not converge.')

        return mf

    #
    @timing.timed
    def ao_to_mo(self, pymol, orbs):
        """perform AO to MO integral transformation using the current
           orbitals"""

        # Do the AO -> MO transformation
        if self.mol.use_df:
            # save the auxbasis value: either user requested or the
            # PySCF default
            ij_trans = np.concatenate(([orbs], [orbs]))
            df.outcore.general(pymol, ij_trans, self.moint_2e_eri,
                        auxbasis=self.mol.ri_basis, dataname='eri_mo')
        else:
            eri_ao = pymol.intor('int2e_sph', aosym='s8')
            eri_mo = ao2mo.incore.full(eri_ao, orbs)
            with h5py.File(self.moint_2e_eri, 'w') as f:
                f['eri_mo'] = eri_mo

        # Construct the core Hamiltonian
        one_nuc_ao = pymol.intor('int1e_nuc')
        one_kin_ao = pymol.intor('int1e_kin')
        hcore_ao   = one_nuc_ao + one_kin_ao

        # transform the core Hamiltonian to MO basis explicitly 
        # and write to file
        h1_mo      = np.matmul(np.matmul(orbs.T, hcore_ao), orbs)
        with h5py.File(self.moint_1e, 'w') as f:
            f['hcore_mo'] = h1_mo

        return

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


