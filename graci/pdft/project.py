#!/usr/bin/env python
'''
Core-Projected DFT Helper Functions
'''

import numpy as np
from pyscf import gto
from pyscf import lib
import pylibxc as libxc
from scipy import linalg

def assign_core_aos(obj):
    '''
    Function to define the subset of core atomic orbitals.

    args:    gto.mole.Mole object

    returns: list of int

    '''
    if type(obj) is gto.mole.Mole:
        core_aos, core_dict = _assign_core_aos_by_label(obj)
        return core_aos
    else:
        raise TypeError

def _assign_core_aos_by_label(obj):
    '''
    Function to define the subset of core atomic orbitals.

    Core orbitals are crudely defined as:
    '1s' for second-row elements (Li through Ne)
    '1s','2s',and '2p' for third-row elements (Na through Ar)
    '1s','2s','2p','3s',and '3p' for fourth-row elements (K through Kr)
    This **roughly** assigns core orbitals as those with an orbital energy below -100 eV.

    That is, the 2s and 2p orbitals in Carbon are not considered core orbitals.
    Note that the MO energies (MO eigenvalues) in methane are roughly:
    '1s':    290 eV
    '2s':     23 eV

    Also builds a dictionary object with entries:

    atom-id: [symbol-str, atom_charge, nl-str]
    int: [str, int, str]

    ex: 0: ['F', 9, '1s']

    __________________________________________
    args:    gto.mole.Mole object

    returns: list of int, dict

    '''
    if type(obj) is gto.mole.Mole:
        mol = obj.copy()
    else:
        try:
            mol = (obj.mol).copy()
        except:
            raise TypeError

    labs = mol.ao_labels() #returns list of str
    uflabs = mol.ao_labels(fmt=False) #returns list of tuple
    core_aos = []
    core_ao_dict = dict()

    # unformatted: [(atom-id, symbol-str, nl-str, str-of-AO-notation)], ex:(0, 'F', '1s', ''), (1, 'F', '3d', 'x2-y2'), etc.
    # formatted:   ['0 F 1s    ', '1 H 1s    ']
    for iao in range(mol.nao):
        lbl = labs[iao]
        tpl = uflabs[iao]
        elemnt = tpl[1]
        ao_type = tpl[2]
        atm_id = int(tpl[0])
        atm_Z = mol.atom_charge(atm_id)
        core_ao_defn = {'row1':('1s'), 'row2':('1s','2s','2p'), 'row3':('1s','2s','2p','3s','3p')}
        if atm_Z > 2:
            if ao_type in core_ao_defn['row1']:
                core_ao_dict[atm_id] = [elemnt, atm_Z, ao_type]
                core_aos.append(iao)
        elif atm_Z > 10 and atm_Z < 19:
            if ao_type in core_ao_defn['row2']:
                core_ao_dict[atm_id] = [elemnt, atm_Z, ao_type]
                core_aos.append(iao)
        elif atm_Z > 18 and atm_Z < 37:
            if ao_type in core_ao_defn['row3']:
                core_ao_dict[atm_id] = [elemnt, atm_Z, ao_type]
                core_aos.append(iao)
        elif atm_Z > 36:
            raise ValueError("Program does not currently support 5th row elements (Z>36).")
    core_aos = [*set(core_aos)] # Remove duplicates
    return core_aos, core_ao_dict

def build_orthogonalizer(S, default=True):
    '''
    Function to build orthogonalizer.
    '''
    s, U = linalg.eigh(S)
    if not default:
        Ndim = S.shape; N = Ndim[0]
        shalf = np.zeros(Ndim, dtype = np.complex128)
        N = self.mol.nao
        for i in range(N):
            if( abs(s[i]) > 1e-12):
                shalf[i,i] = s[i]**(-0.5)
        X = np.dot(U, np.dot(shalf, U.T))
    else:
        eig = s + 0j
        shalf = np.diag( (eig**(-0.5)) )
        ## diagonalizing matrix, X = U * s^(-1/2) * Ut
        # X = np.dot(U, np.dot(shalf, U.T) )
        X = U @ shalf @ U.T
    return X

def build_proj_in_basis(mydft):
    '''
    Build projector in the same basis.

    Args:
    mydft                         [dft.ROKS object]

    Returns:
    [QS,SQ]    Core-Projector     [np.array]
    '''
    ## Call up variables.
    m = mydft.mol
    S = mydft.get_ovlp()
    #SM = linalg.inv(S)
    N = m.nao           # total number of aos

    ## Get integral overlaps
    #s = m.intor_symmetric('int1e_ovlp')

    ## List of Indices of Core AOs
    caos = assign_core_aos(m)

    # Build overlap matrix for core aos only
    NC = len(caos)
    if(NC>0): ## if number of core aos > 0
        Q = np.zeros((N, N))
        # iterate over core aos
        for ic in range(NC):
            for jc in range(NC):
                Q[caos[ic], caos[jc]] = S[caos[ic], caos[jc]] # = 1

    ## Define the operators
    QS = np.einsum('ik,kj->ij', Q, S)
    SQ = np.einsum('ik,kj->ij', S, Q)
    ## Consolidate into single (2,N,N) dim array
    SQQS = np.array((SQ,QS))
    return SQQS

def build_proj_in_ext_basis(mydft, ext_basis = '3-21G'):
    '''
    Build projector in an external basis.

    Args:
    mydft                         [dft.ROKS object]
    ext_basis  (default:sto3g)    str


    Returns:
    [QS,SQ]    Core-Projector     [np.array]
    '''
    ## Call up variables.
    M = mydft.mol
    S = mydft.get_ovlp()
    SM = linalg.inv(S)
    N = M.nao           # total number of aos

    ## Minimal basis
    m =  M.copy()
    m.basis = ext_basis
    m.build()

    ## Get integral overlaps
    s = m.intor_symmetric('int1e_ovlp')      # 1e-integral ( | )
    Sx = gto.intor_cross('int1e_ovlp',M, m)  # cross-overlap matrix

    ## List of Indices of Core AOs
    caos = assign_core_aos(m)

    # build overlap matrix -- core aos only
    NC = len(caos)
    if(NC>0): ## if number of core aos > 0
        SC = np.zeros((NC,NC))
        SXC = np.zeros((N,NC))
        # iterate over core aos (of ext basis)
        for ic in range(NC):
            SXC[:,ic] = Sx[:, caos[ic]]
            for jc in range(NC):
                SC[ic,jc] = s[caos[ic], caos[jc]]

    ## Assume that the ith orthogonalized core AO equals the ith core AO
    # decompose S into s, U, where s = Ut * S * U
    #X = build_orthogonalizer(SC)
    #SCm = (np.dot(X.T,X) ).real
    SCm = linalg.inv(SC)

    # Core AO projection operators in current basis set
    Q = np.einsum('ia,ab,bc,dc,dj->ij',SM,SXC,SCm,SXC,SM)
    # SM: inverse of ovlap (internal basis)
    # SXC: cross-ovlap, core-aos only (external/internal basis)
    # SCm: orthogonalized ovlap, core-aos only (external basis)

    ## Define the operators
    QS = np.einsum('ik,kj->ij', Q, S)
    SQ = np.einsum('ik,kj->ij', S, Q)
    ## Consolidate into single (2,N,N) dim array
    SQQS = np.array((SQ,QS))

    return SQQS

def old_build_proj_in_basis(mydft):
    '''
    Build projector in the same basis.

    Args:
    mydft                         [dft.ROKS object]

    Returns:
    [QS,SQ]    Core-Projector     [np.array]
    '''
    ## Call up variables.
    m = mydft.mol
    S = mydft.get_ovlp()
    SM = linalg.inv(S)
    N = m.nao           # total number of aos

    ## Get integral overlaps
    #s = m.intor_symmetric('int1e_ovlp')

    ## List of Indices of Core AOs
    caos = assign_core_aos(m)

    # Build overlap matrix for core aos only
    NC = len(caos)
    if(NC>0): ## if number of core aos > 0
        SC = np.zeros((NC,NC))
        SXC = np.zeros((N,NC))
        # iterate over core aos
        for ic in range(NC):
            SXC[:,ic] = S[:, caos[ic]]
            for jc in range(NC):
                SC[ic,jc] = S[caos[ic], caos[jc]]
    ## Assume that the ith orthogonalized core AO equals the ith core AO
    # decompose S into s, U, where s = Ut * S * U
    # X = build_orthogonalizer(SC)
    # SCm = (np.dot(X.T,X)).real
    SCm = linalg.inv(SC)

    # Core AO projection operators in current basis set
    Q = np.einsum('ia,ab,bc,dc,dj->ij',SM,SXC,SCm,SXC,SM)
    # SM: inverse of overlap matrix
    # SXC: cross-overlap matrix of dim [N X NC]
    # SCm: orthogonalized core overlap of dim [NC X NC]

    ## Define the operators
    QS = np.einsum('ik,kj->ij', Q, S)
    SQ = np.einsum('ik,kj->ij', S, Q)
    ## Consolidate into single (2,N,N) dim array
    SQQS = np.array((SQ,QS))
    return SQQS

### Extra PROJ METHODS DEFINED BELOW
def _xc_handler(xc, include = False):
     '''
     Extract information about the libxc XC functional.
     Given XC Functional name, return a dictionary of X and C components.
     Must be using PyLibXC (libxc) version 6.2.2 or later. 
     The module pylibxc2 (available via pip) is not compatible.

     Given:
     xc      [str]
     include [bool]

     Returns (include = False):
     K, X, C [float, dictionary, dictionary]
     or
     Returns (include = True):
     X, C     [dictionary, dictionary]

     Example:
     >>> mf._xc_handler(xc='HYB_GGA_XC_PBEH', include = False)
        (0.25, {'gga_x_pbe': 0.75}, {'gga_c_pbe': 1.0})
     >>> mf._xc_handler(xc='HYB_GGA_XC_PBEH', include = True)
        ({'gga_x_pbe': 0.75, 'HF': 0.25}, {'gga_c_pbe': 1.0})

     Another Example (how to use pylibxc):
     >>> import pylibxc as xc
     >>> qtp17 = xc.LibXCFunctional('HYB_GGA_XC_QTP17', spin=1)
     >>> qtp.describe()
     'Functional ID: 460\n
     Functional Name: hyb_gga_xc_qtp17\nAttributes:\n
     Name: Global hybrid for vertical ionization potentials\n
     Kind: 2\n
     Family: 32\nCitations:\n
     Y. Jin and R. J. Bartlett.,  J. Chem. Phys. 149, 064111 (2018)'

     >>> qtp17.aux_funcs()
     [('lda_x', 0.38), ('lda_c_vwn_rpa', 0.19999999999999996), ('gga_c_lyp', 0.8)] 

     a_x, a_c = qtp17.get_ext_param_default_values()
              = [0.62, 0.8]
     '''
     xcfunc =  libxc.LibXCFunctional(xc, 1) #spin=self.mol.spin)
     try:
         K = xcfunc.get_hyb_exx_coef()
     except:
         K = 0.0
     aux = xcfunc.aux_funcs()
     exch = dict()
     corr = dict()
     for aux_func in aux:
         n = aux_func[0]
         coeff = round(aux_func[1],6)
         if '_x' in n:
             exch[n] = coeff
         elif '_c' in n:
             corr[n] = coeff
     if K != 0.0 and include:
         exch['HF'] = K
         return exch, corr
     elif K != 0.0 and not include:
         return K, exch, corr

def format_xc(xc_list):
    '''
    Format the dictionary containing the X and/or C components of the XC functional.

    Given dict: {'lda_c_vwn_rpa':0.200000, 'gga_c_lyp':0.800000}
    Return str: '.80*gga_c_lyp + .20*lda_c_vwn_rpa' == '.80*LYP + .20*VWN'
    '''
    part_str = "{}*{}"
    parts = []

    for aux in xc_list.keys():
        coeff = xc_list[aux]
        fmtd = part_str.format(coeff, aux)
        parts.append(fmtd)
    func_str = ' + '.join(parts)
    return func_str
