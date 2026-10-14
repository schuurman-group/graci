"""
Module for the Davidson diagonalisation of the MRSF-TDDFT response
matrix
"""

import sys
import numpy as np
import graci.core.libs as libs
import graci.utils.timing as timing

@timing.timed
def diag(ci_method, irrep, guess=None):
    """
    Compute the lowest ci_method.nstates[irrep] roots of the
    ci_method.mult manifold in the irrep 'irrep' (PySCF irrep id).
    guess: optional (xdim, nguess) array of start vectors (compressed
    layout, already aligned to the present MOs, any irrep content: the
    library masks them and drops the empty ones)

    Returns:
        nroots : number of roots actually computed
        ener   : excitation energies relative to the reference (Ha)
        xvec   : compressed amplitude vectors, shape (xdim, nroots)
        niter  : number of Davidson iterations
        iconv  : 1 if converged, 0 otherwise
        tracked: (ener_all, xvec_all) of all nroots + nextra tracked roots
                 (start vectors of the next geometry)
    """

    nroots = int(ci_method.nstates[irrep])
    xdim   = int(ci_method.xdim)
    mult   = int(ci_method.mult)
    nextra = int(ci_method.nextra)
    maxvec = int(ci_method.diag_maxvec) * (nroots + nextra)
    maxit  = int(ci_method.diag_iter)
    tol    = float(ci_method.diag_tol)
    # correction equations: 1 = diagonal preconditioner, 2 = Jacobi-Davidson
    # with the kernel-free operator (singlets with the CV kernel only)
    precond = {'diag': 1, 'inner': 2}[str(ci_method.precond).lower()]
    inner   = int(ci_method.inner_iter)

    nguess = 0 if guess is None else int(np.shape(guess)[1])
    ncols  = max(nroots + nextra, nguess, 1)
    ener = np.zeros(ncols, dtype=float)
    xvec = np.zeros((xdim, ncols), dtype=np.float64, order='F')
    if nguess > 0:
        if np.shape(guess)[0] != xdim:
            sys.exit('mrsf_diag: the guess vectors have the wrong dimension')
        xvec[:, :nguess] = guess

    args = (irrep, mult, nroots, nextra, maxvec, maxit, tol, precond, inner,
            nguess, xdim, ener, xvec, 0, 0, 0)
    nroots_out, ener, xvec, nvec_out, niter, iconv = libs.lib_func('mrsf_diag', args)

    nroots_out = int(nroots_out)
    nvec_out   = int(nvec_out)
    xvec = np.array(xvec)

    return nroots_out, ener[:nroots_out], xvec[:, :nroots_out], \
           int(niter), int(iconv), (np.array(ener[:nvec_out]), xvec[:, :nvec_out])
