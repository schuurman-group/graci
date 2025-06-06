"""
Module for loading/finalising the bitwf library
"""

import sys as sys
import numpy as np
from pyscf import gto
import graci.core.libs as libs

#
def init(bra, ket, calctype, verbose):
    """Initialise the bitwf library"""
    
    # check on |Nel_bra - Nel_ket|
    if calctype in ['overlap'] and bra.nel != ket.nel:
        sys.exit('calctype = overlap and bra.nel != ket.nel in init_bitwf')
    if calctype == 'dyson' and np.abs(bra.nel - ket.nel) != 1:
        sys.exit('calctype = dyson and |bra.nel-ket.nel| != 1 in init_bitwf')
        
    # bra-ket MO overlaps
    molBra = bra.scf.mol.mol_obj.copy()
    molKet = ket.scf.mol.mol_obj.copy()
    smat   = gto.intor_cross('int1e_ovlp',
                             bra.scf.mol.mol_obj,
                             ket.scf.mol.mol_obj)
    smat   = np.matmul(np.matmul(bra.mos.T, smat), ket.mos)

    # bra and ket point groups
    if bra.scf.mol.sym_indx <= 0:
        pgrpBra = 1
    else:
        pgrpBra = bra.scf.mol.sym_indx + 1
    if ket.scf.mol.sym_indx <= 0:
        pgrpKet = 1
    else:
        pgrpKet = ket.scf.mol.sym_indx + 1
        
    # set all variable that have to be passed to bitwf_initialise
    multBra  = bra.mult
    multKet  = ket.mult
    nelBra   = bra.nel
    nelKet   = ket.nel
    nmoBra   = bra.nmo
    nmoKet   = ket.nmo
    smat     = np.reshape(smat, (nmoBra * nmoKet), order='F')
    mosymBra = np.array(bra.mosym)
    mosymKet = np.array(ket.mosym)

    # call to bitwf_initialise
    args = (multBra, multKet, nelBra, nelKet, nmoBra, nmoKet, smat,
            pgrpBra, pgrpKet, mosymBra, mosymKet, calctype, verbose)
    libs.lib_func('bitwf_initialise', args)
    return

def finalize():
    """
    finalize the bitwf library
    """

    args = ()
    libs.lib_func('bitwf_finalise', args)
    
    return
