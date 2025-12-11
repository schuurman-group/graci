# Copyright 2014-2020 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''
Projected Hybrid Density functional theory
=========================

Simple usage::

    >>> from pyscf import gto, dft
    >>> mol = gto.M(atom='N 0 0 0; N 0 0 1', basis='def2-tzvp')
    >>> mf = dft.RKS(mol,paos='CoreAOs',phyb=1)
    >>> mf.xc = 'pbe,pbe'
    >>> mf.run()
'''
import sys
sys.path.append('../')
try:
    from pyscf.dft import libxc
    XC = libxc.XC
except (ImportError, OSError):
    pass
try:
    from pyscf.dft import xcfun
    XC = xcfun.XC
except (ImportError, OSError):
    pass
#from pyscf.dft import xc
from graci.pdft import rks
from graci.pdft import roks
from graci.pdft import uks
from graci.pdft import gks
from graci.pdft import rks_symm
from graci.pdft import uks_symm
from graci.pdft import gks_symm
from pyscf.dft import dks
from pyscf.dft import gen_grid as grid
from pyscf.dft import radi
from pyscf.df import density_fit
from pyscf.dft.gen_grid import sg1_prune, nwchem_prune, treutler_prune, \
        stratmann, original_becke, Grids
from pyscf.dft.radi import BRAGG_RADII, COVALENT_RADII, \
        delley, mura_knowles, gauss_chebyshev, treutler, treutler_ahlrichs, \
        treutler_atomic_radii_adjust, becke_atomic_radii_adjust


def KS(mol, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis = '3-21G', use_ext_basis=True):
    if mol.spin == 0:
        return RKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
    else:
        return UKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
KS.__doc__ = '''
A wrap function to create DFT object (RKS or UKS).\n
''' + rks.RKS.__doc__
DFT = KS

def RKS(mol, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis = '3-21G', use_ext_basis=True):
    if mol.nelectron == 1:
        return uks.UKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
    elif not mol.symmetry or mol.groupname == 'C1':
        if mol.spin > 0:
            return roks.ROKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
        else:
            return rks.RKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
    else:
        if mol.spin > 0:
            return rks_symm.ROKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
        else:
            return rks_symm.RKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
RKS.__doc__ = rks.RKS.__doc__

def ROKS(mol, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis = '3-21G', use_ext_basis=True):
    if mol.nelectron == 1:
        return uks.UKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
    elif not mol.symmetry or mol.groupname == 'C1':
        return roks.ROKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
    else:
        return rks_symm.ROKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
ROKS.__doc__ = roks.ROKS.__doc__

#def UKS(mol, xc='LDA,VWN'):
def UKS(mol, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis = '3-21G', use_ext_basis=True):
    if not mol.symmetry or mol.groupname == 'C1':
        return uks.UKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
    else:
        return uks_symm.UKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
UKS.__doc__ = uks.UKS.__doc__

def GKS(mol, xc='LDA,VWN', phyb=[0.0], paos=None, ext_basis = '3-21G', use_ext_basis=True):
    if not mol.symmetry or mol.groupname == 'C1':
        return gks.GKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)
    else:
        return gks_symm.GKS(mol, xc, phyb, paos, ext_basis, use_ext_basis)

GKS.__doc__ = gks.GKS.__doc__

def DKS(mol, xc='LDA,VWN'):
    from pyscf.scf import dhf
    if dhf.zquatev and mol.spin == 0:
        return dks.RDKS(mol, xc=xc)
    else:
        return dks.UDKS(mol, xc=xc)
