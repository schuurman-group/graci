"""Module for handling local basis set definitions"""

import os
import sys
import re
import pyscf.gto as gto

def local_basis_sets(local_dir=False, return_alias=False):
    """
    Return a list of available local basis sets by alias
    name, i.e. with all special characters stripped 
    and converted to lower-case
    """
    if local_dir:
        if 'GRACI_EXT_BASIS_DIR' in os.environ:
            bdir = os.environ['GRACI_EXT_BASIS_DIR']
        else:
            return []
    else:
        bdir = os.environ['GRACI']+'/graci/utils/basis_sets'

    basis_files = [f.replace('.dat','')
                   for f in os.listdir(bdir) if
                   os.path.isfile(os.path.join(bdir, f))]
    
    if return_alias:
        return [basis_files[i].lower().replace('-','').replace('_','')
                                     for i in range(len(basis_files))]
    else:
        return basis_files

#
def load_basis(atom, name, local_dir=False):
    """
    load basis function specified by 'alias' for atom 'atom'
    """
    if local_dir:
        if 'GRACI_EXT_BASIS_DIR' in os.environ:
            bdir = os.environ['GRACI_EXT_BASIS_DIR']
        else:
            return None
    else:
        bdir = os.environ['GRACI']+'/graci/utils/basis_sets'

    alias = name.lower().replace('-','').replace('_','')
    basis_avail = local_basis_sets(local_dir=local_dir)
    basis_alias = local_basis_sets(local_dir=local_dir, return_alias=True)

    if alias not in basis_alias:
        sys.exit('basis set: ' + str(name) +
                 ' not in basis sets in directory: ' + str(bdir))

    fname = basis_alias.index(alias)
    with open(bdir+'/'+basis_avail[fname]+'.dat','r') as bf:
        bfile = bf.readlines()

    bf_atm = ''
    start  = -1
    while bf_atm != atom.strip().lower() and start<len(bfile)-1:
        start += 1
        bf_atm = bfile[start][:2].strip().lower()

    if start == len(bfile)-1:
        sys.exit('atom: ' + str(atom) +
                 ' not in basis set file: ' + str(alias+'.dat'))

    end = start
    while '#BASIS SET' not in bfile[end].strip() and end < len(bfile)-1:
        end += 1

    bstr = ''
    for i in range(start, end):
        bstr += bfile[i]

    return gto.basis.parse(bstr)

#
def str_to_contract(cstr):
    """
    attempt to convert a string of the form XsYpZd, etc. to 
    an array of contractions
    """

    angmom = ['s','p','d','f','g','h','i'] 
    ncon   = []
    constr = str(cstr)

    for ia in range(len(angmom)):

        val = constr.strip().split(angmom[ia])

        # if search string was found, try
        # to send contraction
        if len(val) > 1:
           try:
                ncon.append(int(val[0]))
           except:
                return None
           constr = val[1]
           if len(constr) == 0:
               break
           
        else:
            ncon.append(0)
 
    # if we didn't find any contractions, not valid 
    # constraction string
    if max(ncon) == 0:
        return None
    else:
        return ncon





