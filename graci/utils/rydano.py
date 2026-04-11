"""Module for handling local basis set definitions"""

import os
import sys
import numpy as np
import scipy as sp
import operator as op
import pyscf.gto as gto
import graci.core.molecule as molecule
import graci.core.scf as scf
import graci.io.output as output
import graci.utils.basis as basis
import graci.utils.constants as constants

ang_lbls = ['s','p','d','f','g','h','i']

class Rydano():
    """
    Generate a rydberg basis for a given molecule object
    """

    def __init__(self):
        # user variables
        self.xc          = 'qtp17'
        self.mult        = 2
        self.charge      = 1
        self.origin      = None
        self.label       = 'default'
        self.verbose     = True
        self.contract    = '1s1p1d'
        self.nprimitive  = 8
        self.print_ano   = False
        self.max_overlap = [1. / np.sqrt(2.)]

        # class variables
        self.anos       = []
        self.occs       = []


    def run(self, mol):
        """
        run an SCF calculation on the corresponding doublet
        """

        # make a copy of the molecule object we will use
        # to run the cation calculation
        ion_mol = mol.copy()

        # build the PySCF molecule object 
        ion_mol.run()

        # the valence basis will be DZP for the non-H
        # atoms and DZ for H
        for atm in list(set(ion_mol.asym)):
            if atm.lower().strip() == 'h':
                bas_lbl = 'dz'
            else:
                bas_lbl = 'dzp'
            ion_mol.basis[atm]    = bas_lbl
            ion_mol.ri_basis[atm] = None

        # maximum angular momentum value to include in AO basis
        n_con = basis.str_to_contract(self.contract)
        l     = [ind for ind,val in enumerate(n_con) if val>0]
        l_max = max(l)

        # add the atom to the ion_mol object. ANOs not specified, 
        # so primitives are added uncontracted. Determine which
        # primitives to include based on the original basis 
        # function -- not the dz(p) used to compute virtuals
        if self.origin is None:
            ryd_origin = np.asarray(ion_mol.nuc_charge_center())
        else:
            ryd_origin = self.origin 

        # if the rydberg basis origin is distinct from any atoms,
        # add a ghost atom, else, append basis set of atom at
        # origin
        exps = self.determine_exponents(mol, l, self.nprimitive, ryd_origin)
        nao_per_l, bas_str, bas_lst = self.make_ano_basis('X', l, exps)
        # label for the ghost atom
        xind = self.add_ano_basis(ion_mol, ryd_origin, bas_str, bas_lst)

        # print header information to the log
        if self.verbose:
            output.print_rydano_header(self)

        # construct an SCF object
        ion_scf          = scf.Scf()
        ion_scf.xc       = self.xc
        ion_scf.mult     = self.mult
        ion_scf.charge   = self.charge
        ion_scf.verbose  = False 
        ion_scf.do_ao2mo = False
        ion_scf.print_orbitals = False

        # generate the ionic density
        ion_scf.run(ion_mol, None)

        # print orbital info  / scf summary
        if self.verbose:
            output.print_scf_summary(ion_scf)

        # parse the AO-basis to get relevant labels
        a_ind, a_lbl, l_i, n_l, l_lbl = self.parse_ao_basis(ion_scf.mol)

        # extract virtual orbitals and set occupation
        # using the expression exp[6.9*(orb_ener[i]/min(orb_ener)-1)] 
        occ = self.set_virtual_occ(ion_scf) 

        # form the l-projected NOs
        S_ao = ion_mol.pymol().intor('int1e_ovlp')
        rocc, rnos = \
              self.make_nos(xind, a_ind, nao_per_l, S_ao, ion_scf.orbs, 
                            occ, l_i, l_lbl, cart=mol.ao_cart)

        # set the phase
        for li in range(len(rnos)):
            rnos[li] = self.set_phase(rnos[li])

        self.anos = []
        self.occs = []
        for li in range(len(n_con)):
           self.anos.append(rnos[li][:,:n_con[li]])
           self.occs.append(rocc[li][  :n_con[li]])      
       
        # add the Rydberg basis
        nao_per_l, bas_str, bas_lst = self.make_ano_basis(
                                                ion_mol.asym[xind], l, 
                                                exps, ano=self.anos)

        self.add_ano_basis(mol, ryd_origin, bas_str, bas_lst)
        # standard RI basis not consistent with rydberg orbitals,
        # use even-tempered
        mol.ri_basis = None

        # print basis string to file if requested
        if self.print_ano:
            with open('rydano.dat', 'w') as f:
                f.write(bas_str) 

        # print the results
        if self.verbose:
            output.print_rydano_summary(exps, rocc, rnos, self.contract)
 
        return

    #
    def kbj_exp(self, l, i):
        """generate the 'nth' KBJ expononent for specified value of l

        Args:
           i: ith value in the series: 1, 1.5, 2, 2.5, 3, 3.5,...
           l: angular momentum value

        Returns:
          ex: Rydberg exponent corresponding to n,l value
        """ 

        a = [0.584342, 0.452615, 0.382362, 0.337027, 0.304679]
        b = [0.424483, 0.309805, 0.251333, 0.215013, 0.189944]
        n = 1 + i*0.5

        return (1./( 4*n**2 )) * (1./( a[l]*n + b[l] )**2)
    
    def _contracted_overlap(self, alpha, l, shell, atom_pos_bohr,
                            origin_bohr):
        """Overlap between a normalised KBJ primitive (exponent alpha,
        angular momentum l, at the Rydberg origin) and a normalised
        contracted valence shell of the same l at atom_pos_bohr, computed
        via PySCF's integral engine.

        Args:
            alpha:         KBJ primitive exponent (bohr^-2)
            l:             angular momentum
            shell:         PySCF shell entry [l, [exp, coef], ...]
            atom_pos_bohr: atomic centre position (bohr, numpy array)
            origin_bohr:   Rydberg basis origin (bohr, numpy array)

        Returns:
            max_overlap: maximum absolute element of the overlap matrix
        """
        # KBJ primitive at the Rydberg origin (unit contraction coefficient;
        # PySCF normalises it internally)
        mol_kbj = gto.Mole()
        mol_kbj.atom  = [('H', tuple(origin_bohr))]
        mol_kbj.basis = {'H': [[l, [alpha, 1.0]]]}
        mol_kbj.unit  = 'Bohr'
        mol_kbj.verbose = 0
        mol_kbj.spin = 1
        mol_kbj.build()

        # contracted valence shell at the atomic centre
        mol_val = gto.Mole()
        mol_val.atom  = [('H', tuple(atom_pos_bohr))]
        mol_val.basis = {'H': [shell]}
        mol_val.unit  = 'Bohr'
        mol_val.verbose = 0
        mol_val.spin = 1.
        mol_val.build()

        ovlp = gto.mole.intor_cross('int1e_ovlp_sph', mol_kbj, mol_val)
        return np.max(np.abs(ovlp))

    def determine_exponents(self, mol, l, nprim, ryd_origin):
        """Determine the appropriate KBJ primitives to include for a
        given AO basis set using overlap-based criterion.

        For each angular momentum l, walks the KBJ series (from the
        most compact end) and finds the first entry whose maximum
        overlap with any same-l contracted valence shell is below
        self.max_overlap. The series of nprim primitives starts there.

        Args:
            mol:        molecule object [Molecule]
            l:          list of angular momentum values to include
            nprim:      number of primitives to include per l
            ryd_origin: Rydberg basis origin in Angstrom (numpy array)

        Return:
            exps: list of exponent lists, one per l value
        """

        # ryd_origin is in Angstrom; convert to bohr for overlap integrals
        origin_bohr = ryd_origin * constants.ang2bohr

        # atom positions in bohr from PySCF (atom_coords returns bohr)
        atom_coords_bohr = mol.pymol().atom_coords()

        # collect contracted shells per l value from the valence basis,
        # together with the atomic centre position in bohr
        shells_per_l = [[] for _ in range(len(l))]
        for iatm, atm in enumerate(mol.asym):
            if atm not in mol.basis_obj:
                continue
            pos = atom_coords_bohr[iatm]
            for shell in mol.basis_obj[atm]:
                lval = shell[0]
                if lval in l:
                    li = l.index(lval)
                    shells_per_l[li].append((shell, pos))

        # if a requested l is absent from the valence basis, borrow
        # shells from the next lower l that is present, relabelling
        # the angular momentum so the overlap is computed correctly
        for li in range(len(l)):
            if len(shells_per_l[li]) == 0 and li > 0:
                relabelled = [([l[li]] + shell[1:], pos)
                              for shell, pos in shells_per_l[li - 1]]
                shells_per_l[li] = relabelled

        # walk the KBJ series and find the starting index where the
        # maximum overlap with any same-l valence shell drops below
        # max_overlap
        kbj_i = [0] * len(l)
        max_ovr = [self.max_overlap[l] if l < len(self.max_overlap) 
                                     else self.max_overlap[-1] 
                                     for l in range(len(l))]

        for li in range(len(l)):
            while True:
                alpha = self.kbj_exp(l[li], kbj_i[li])
                max_ovlp = max(
                    (self._contracted_overlap(alpha, l[li], shell, pos,
                                             origin_bohr)
                     for shell, pos in shells_per_l[li]),
                    default=0.)
                if max_ovlp <= max_ovr[li]:
                    break
                kbj_i[li] += 1

        exps = [[self.kbj_exp(l[i], kbj_i[i] + j) for j in range(nprim)]
                for i in range(len(l))]

        return exps

    #
    def add_ano_basis(self, mol, crds, basis_str, basis_lst):
        """Add an atom to a molecule object, including the basis
           set
         
        Args:
            mol:  molecule object to amend [Molecule]
            crds: cartesian coordinates of the atom
            l:    angular momentum values to include in basis
            exps: the exponents on the primitive functions.
            nos:  list of numpy array at least max(l) long. Each 
                  column of numpy array is ANO, ordered by 'n'
    
        Returns:
            None
        """

        # numerical issues arise if we try to place a dummy atom
        # on top of an existing atom. Check minimum distance, and
        # if less than threshold, displace a bit.
        dvec = [np.linalg.norm(crd - crds) for crd in mol.crds]
        xind = None

        # if there is an existing atom at the origin, append
        # basis set to existing atom
        if min(dvec) <= 1.e-3:
            ia = np.argmin(np.array(dvec))
            self.append_ano_basis(mol, ia, basis_lst)
            xind = ia

        else:
            # add a ghost atom
            sym = 'X'

            # add atomic symbol and cartesian coordinaes
            iatm = molecule.atom_name.index(sym)
            mol.asym.append(sym)

            mol.crds   = np.append(mol.crds, [crds], axis=0)
            mol.masses = np.append(mol.masses, molecule.atom_mass[iatm])
            # construct basis object of uncontracted primitives
            mol.basis[sym]    = basis_str
            #RI basis definition limited to even-tempered for time being
            if mol.use_df and mol.ri_basis is not None:
                mol.ri_basis[sym] = None
            xind = mol.asym.index(sym)

        return xind

    #
    def append_ano_basis(self, mol, iatm, bas_lst):
        """
        append basis specified by bas_str to iatm
        """

        # get the atomic symbol of the atomic center
        isym = mol.asym[iatm]

        # if there's more than one, we'll need to create a special
        # label to distinguish the basis from the other atoms
        nsym = mol.asym.count(isym)
        basis_set = gto.basis.load(mol.basis[isym], isym)
        for ano_set in bas_lst:
            l = ano_set[0]
            ind = [bset[0] for bset in basis_set].index(l)
            basis_set.insert(ind, ano_set)
            ind += 1

        if nsym > 1:
            mol.asym[iatm] += ' 1'
        mol.basis[mol.asym[iatm]] = basis_set

        if mol.use_df and mol.ri_basis is not None:
            mol.ri_basis[mol.asym[iatm]] = None

    #
    def make_ano_basis(self, sym, l, exps, ano=None):
        """
        construct an uncontracted basis of KBJ functions, running 
        from 0 to l_max, and for exponents specified by the values
        of n in n_range
        """
        l_lbl = ['S', 'P', 'D', 'F', 'G', 'H', 'I']

        ostr = '#BASIS set generated by rydano\n'
        astr = '{:<6s}{:1s}\n'

        # determine the number of basis functions per ang. mom. 
        # shell
        nao_per_l = [0]*(max(l)+1)
        for li in l:
            nao_per_l[li] = len(exps[li])

        # if ANOs is None, or, they don't match the number of
        # primitives, use uncontracted functions
        if ano == None:
            rano = [np.identity(len(exps[li]), dtype=float)
                                             for li in range(len(l))]
        else:
            rano = ano
            for li in range(len(rano)):
                if rano[li].shape[0] != len(exps[li]):
                    print('Number of primitivies in ANO incorrect: ' +
                              str(rano[li].shape[0]) + ' != ' +
                              str(len(exps[li])))
                    sys.exit(1)

        olst = []
        # loop over angular momentum shells to build uncontracted
        # basis and NWChem file format
        for li in range(len(l)):
            nprim = rano[li].shape[0]
            ncon = rano[li].shape[1]
            pstr = '{:>15.7f}'+''.join(['{:14.7f}']*ncon)+'\n'

            # if ANOs are supplied, we're adding contracted
            # basis functions, else, if ano is 'None', we're 
            # adding uncontracted functions. This should be
            # cleaned up
            if ano != None:
                for ic in range(ncon):
                    bf = [l[li]]
                    for ip in range(nprim):
                        bf.append([exps[li][ip], rano[li][ip, ic]])
                    olst.append(bf)

            ostr += astr.format(sym, l_lbl[l[li]])
            for ip in range(nprim):
               args = []
               args.append(exps[li][ip])
               args.extend(rano[li][ip,:])
               ostr += pstr.format(*args)
               # if we're adding uncontracted functions,
               # number of ANO == nprim
               if ano == None:
                   olst.append([l[li], [exps[li][ip], 1.]])

        return nao_per_l, ostr, olst

    #
    def parse_ao_basis(self, mol):
        """parses PySCF AO basis set for ang. mom., labels, etc.

        Args:
          mol: A port GRaCI molecule class object

        Returns:
          a_indx: the int index of atomic center 0=1st, 1=2nd, etc. 
          a_lbl: the text atomic symbol for the center
          l_val: the angular momentum value, l, for each basis function
          l_shl: the quantum number 'n', or, 'shell' value of each basis
                 function of a given value of l
          l_lbl: text label for each basis function, eg. 's','px','dxy'

        """
        ao_lbls = mol.pymol().ao_labels(fmt=True)
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

            l_val[i]  = l
            l_shl[i]  = int(val[2][:val[2].index(ang_lbls[l])].strip())
            l_lbl[i]  = ao_lbls[i].split(ang_lbls[l])[1].strip()

        return a_indx, a_lbl, l_val, l_shl, l_lbl

    # 
    def set_virtual_occ(self, ion_scf):
        """
        set the occupations on the virtual orbitals
        """

        nmo = ion_scf.nmo
        occ = np.zeros(nmo, dtype=float)

        first_vrt = None
        for imo in range(nmo):
            scf_occ  = ion_scf.orb_occ[imo]
            scf_ener = ion_scf.orb_ener[imo]
            if scf_occ == 0. and scf_ener < 0:
                if first_vrt is None:
                    first_vrt = scf_ener
                occ[imo] = np.exp(6.9*(scf_ener/first_vrt - 1.))
            if scf_ener > 0:
                break

        return occ

    #
    def make_nos(self, xind, aind, nao_per_l, Smat, orbs, occ, 
                                                l_i, l_lbl, cart=False):
        """
        form the rydberg density
        """

        # determine number of AOs per angular momentum shell
        n_per_shell = [[1, 3, 5, 7, 9, 11], [1, 3, 6, 10, 15, 21]]
        if cart:
            sh_ind = 1
        else:
            sh_ind = 0

        # tabulate the number of AOs per angular momentum shell
        prim_cnt  = [0]*len(nao_per_l)
        n_prim_ao = [nao_per_l[i] * n_per_shell[sh_ind][i] 
                                  for i in range(len(nao_per_l))]

        x_ao = []
        # identify AOs in the uncontracted basis
        for ai in range(len(aind)):
            if aind[ai] == xind:
                lval = l_i[ai]

                if lval >= len(prim_cnt):
                    continue

                # this requires rydberg AOs are first in the atom
                # specification
                if prim_cnt[lval] < n_prim_ao[lval]:
                    x_ao.append(ai)
                    prim_cnt[lval] += 1

        # first extract basis properties from our selected
        # center (for now assume ghost atoms '')        
        x_mo   = [ind for ind,val in enumerate(occ) if val > 0.] 
        x_l    = [val for ind,val in enumerate(l_i) if ind in x_ao]
        x_lbl  = [val for ind,val in enumerate(l_lbl) if ind in x_ao]
        l_max  = max(x_l)
        nao_x  = len(x_ao)

        dao   = np.zeros((nao_x, nao_x),dtype=float)
        # construct the total Rydberg density
        for imo in range(len(x_mo)):
            mo   = orbs[x_ao, x_mo[imo]]
            dao += occ[x_mo[imo]] * np.outer(mo, mo)

        # compute the corresponding Rydberg orbitals 
        S_ao        = Smat[x_ao,:][:,x_ao]
        s_12, s_m12 = self.ao_s(S_ao)

        # find the orbitals that correspond to the the Rydberg density
        rocc, rmos = self.eig_h(s_12 @ dao @ s_12)
        rmos = s_m12 @ rmos
        # normalize the Rydberg density to 1
        rocc /= np.sum(rocc)

        # project out the various angular momentum components
        # and form m_l specific density matrices
        occ_all = []
        nos_all = []
        for l in range(l_max+1):

            # the various m_l components. We'll use the AO labels
            # which should be fine for sph. harmonics or cartesians
            ml_lbl = list(set([lbl for ind,lbl in enumerate(x_lbl) 
                                                     if x_l[ind] == l]))

            nao_ml = len(ml_lbl)
            if nao_ml != 2*l+1 and nao_ml != int((l+1)*(l+2)/2):
                print('Number of AOs not consistent with spherical or '+
                      'cartesian functions: l=' + str(l) +
                      ' nao=' + str(nao_ml),flush=True)
                sys.exit(1)

            # indices of the AO basis functions for each m_l value
            ml_i = [[ind for ind,val in enumerate(x_l) if
                      val==l and x_lbl[ind]==ml_lbl[ml]]
                      for ml in range(nao_ml)]

            # Ensure the number of AOs for each m_l value are the 
            # the same
            nml_all  = [len(ml_i[ml]) for ml in range(len(ml_i))]           
            if len(list(set(nml_all))) == 1:
                nml = nml_all[0]
            else:
                print('Error in rydano: number of AOs inconsistent: '+
                      'l='+str(l)+' n_ml='+str(nml_all))
                sys.exit(1)

            # project out each set of MOs corresponding to each
            # value of m_l
            dao_ml = np.zeros((nao_ml, nml, nml), dtype=float)
            S_ml   = np.zeros((nao_ml, nml, nml), dtype=float)
            for ml in range(nao_ml):

                # construct density for each value of m_l in AO basis
                for imo in range(rmos.shape[1]):
                    mo              = rmos[ml_i[ml], imo]
                    dao_ml[ml,:,:] += rocc[imo] * np.outer(mo, mo)
                    S_ml[ml,:,:]    = S_ao[ml_i[ml],:][:,ml_i[ml]] 

            # spherically average over the different m_l components
            dao_sph = self.average_matrices(dao_ml)           
            
            # determine the corresponding NOs
            S_l  = self.average_matrices(S_ml)           

            s_12, s_m12  = self.ao_s(S_ml[0])
            occ_l, nos_l = self.eig_h( s_12 @ dao_sph @ s_12)
            nos_l        = s_m12 @ nos_l          

            occ_all.append(occ_l)
            nos_all.append(nos_l)
 
        return occ_all, nos_all

    #
    def set_phase(self, vecs):
        """routine to set the phase of a set of vectors according
           to some convention (here: largest element is positive)

        Args:
           vecs: set of vectors to be phased

        Returns:
           vecs: vectors with sign convention applied to each
                 column
        """

        for i in range(vecs.shape[1]):
            mx_i = np.argmax(np.absolute(vecs[:,i]))
            if vecs[mx_i,i] < 0.:
                vecs[:,i] *= -1.

        return vecs


    #
    def eig_h(self, M):
        """diagonalize hermitian matrix H and return eigenvalues 
           in ascending order.

        Args:
            M: hermitian matrix 

        Returns:
            evals: eigenvalues in ascending order
            evecs: correspoinding eigenvectors
        """

        evals, evecs = np.linalg.eigh(M)

        # put in descending order
        idx   = evals.argsort()[::-1]
        evals = evals[idx]
        evecs = evecs[:,idx]

        return evals, evecs

    #
    def ao_s(self, Smat, thr=1.e-12):
        """construct the S^1/2 and S^-1/2 for a given overlap matrix
        
        Args:
            Smat: an AO overlap matrix
            thr:  threshold for the generalized inverse (default=1.e-12)

        Returns:
            s12:  the S^1/2 matrix
            sm12: the S^-1/2 matrix
        """

        sv, sc = np.linalg.eigh(Smat)
        s_12   = np.diag([np.sqrt(s) for s in sv])
        s_m12  = np.diag([1./np.sqrt(s) if s > thr else 0. for s in sv])

        S12    = sc @ s_12 @ sc.T
        Sm12   = sc @ s_m12 @ sc.T

        return S12, Sm12

    #
    def average_matrices(self, mats):
        """averages as list of matrices. 

        Args:
           mats: a list of nao x nao density matrices

        Returns:
           ave_mat: a single nao x nao averaged density
        """

        ave_mat = np.sum(mats, axis=0) / mats.shape[0] 
        return ave_mat 


    def tofile(self):
        """
        dump the request basis set to file
        """

        return
