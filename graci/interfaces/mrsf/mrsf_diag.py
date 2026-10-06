"""
Module for the Davidson diagonalisation of the MRSF-TDDFT response
matrix
"""

import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing

@timing.timed
def diag(ci_method, irrep):
    """
    Compute the lowest ci_method.nstates[irrep] roots of the
    ci_method.mult manifold in the irrep 'irrep' (PySCF irrep id)

    Returns:
        nroots : number of roots actually computed
        ener   : excitation energies relative to the reference (Ha)
        xvec   : compressed amplitude vectors, shape (xdim, nroots)
        niter  : number of Davidson iterations
        iconv  : 1 if converged, 0 otherwise
    """

    nroots = int(ci_method.nstates[irrep])
    xdim   = int(ci_method.xdim)
    mult   = int(ci_method.mult)
    nextra = int(ci_method.nextra)
    maxvec = int(ci_method.diag_maxvec) * (nroots + nextra)
    maxit  = int(ci_method.diag_iter)
    tol    = float(ci_method.diag_tol)

    ener = np.zeros(nroots, dtype=float)
    xvec = np.zeros(xdim*nroots, dtype=float)

    args = (irrep, mult, nroots, nextra, maxvec, maxit, tol, xdim,
            ener, xvec, 0, 0)
    nroots_out, ener, xvec, niter, iconv = libs.lib_func('mrsf_diag', args)

    nroots_out = int(nroots_out)
    xvec = np.reshape(xvec, (xdim, nroots), order='F')

    return nroots_out, ener[:nroots_out], xvec[:, :nroots_out], \
           int(niter), int(iconv)
