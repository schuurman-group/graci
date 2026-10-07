"""
Wrappers around the gradient kernels of libmrsf (graci/dep/mrsf/src/gradient).
Arrays cross the boundary by address in Fortran order ('dptr' type of
libs.lib_func): outputs are F-contiguous float64 arrays allocated here.
"""
import numpy as np
import graci.core.libs as libs


def fzeros(*shape):
    return np.zeros(shape, dtype=np.float64, order='F')


def farr(a):
    return np.asfortranarray(np.asarray(a, dtype=np.float64))


def grad_init(eri_file):
    """load the hole-particle DF block (second pass over the integral file);
    returns the error code (0 = ok)"""
    return libs.lib_func('mrsf_grad_init', (eri_file, 0))


def grad_free():
    libs.lib_func('mrsf_grad_free', ())


def bso(naux, nmo):
    """B^Q_{p,O_x} for all MOs p: (naux, nmo, 2)"""
    out = fzeros(naux, nmo, 2)
    libs.lib_func('mrsf_grad_bso', (out,))
    return out


def gfock(cx, dims, Ahh=None, Za=None, Zb=None, jq_add=None, want_hp=True, want_hh=False):
    """Coulomb + exchange parts of G^s[D], D^a = Ahh (HH) + sym(Za) (PH),
    D^b = sym(Zb) (PH); returns (Ghpa, Ghpb, Ghha, Ghhb), each with a
    trailing vector index"""
    nocca, nvirb, naux = dims
    nvec = 1
    for a in (Ahh, Za, Zb):
        if a is not None and np.ndim(a) == 3:
            nvec = np.shape(a)[2]
    if jq_add is not None and np.ndim(jq_add) == 2:
        nvec = np.shape(jq_add)[1]
    flags = 0

    def prep(a, shape, bit):
        nonlocal flags
        if a is None:
            return fzeros(*shape)
        flags |= bit
        return farr(np.reshape(np.asarray(a), shape, order='F'))
    Ahh_ = prep(Ahh, (nocca, nocca, nvec), 1)
    Za_ = prep(Za, (nvirb, nocca, nvec), 2)
    Zb_ = prep(Zb, (nvirb, nocca, nvec), 4)
    jq_ = prep(jq_add, (naux, nvec), 8)
    if want_hp:
        flags |= 16
    if want_hh:
        flags |= 32
    Ghpa = fzeros(nocca, nvirb, nvec); Ghpb = fzeros(nocca, nvirb, nvec)
    Ghha = fzeros(nocca, nocca, nvec); Ghhb = fzeros(nocca, nocca, nvec)
    libs.lib_func('mrsf_gfock', (nvec, float(cx), flags, Ahh_, Za_, Zb_, jq_, Ghpa, Ghpb, Ghha, Ghhb))
    return Ghpa, Ghpb, Ghha, Ghhb


def jblocks(jq, dims):
    """J_pq = sum_Q B^Q_pq jq(Q) on the HH, HP, PP blocks for a set of aux vectors"""
    nocca, nvirb, naux = dims
    jq = np.asarray(jq)
    if jq.ndim == 1:
        jq = jq[:, None]
    nvec = jq.shape[1]
    JHH = fzeros(nocca, nocca, nvec); JHP = fzeros(nocca, nvirb, nvec); JPP = fzeros(nvirb, nvirb, nvec)
    libs.lib_func('mrsf_jblocks', (nvec, farr(jq), JHH, JHP, JPP))
    return JHH, JHP, JPP


def state(cx, Xt, dims):
    """per-state exchange-channel pass; returns a dict of the outputs"""
    nocca, nvirb, naux = dims
    out = {'LaH': fzeros(nocca, nocca), 'LaP': fzeros(nvirb, nocca),
           'LbH': fzeros(nocca, nvirb), 'LbP': fzeros(nvirb, nvirb),
           'Ghh': fzeros(nocca, nocca, naux), 'Fhp': fzeros(nocca, nvirb, naux),
           'Yf': fzeros(nvirb, nocca, naux), 'Sq': fzeros(nocca, nocca, naux),
           'gpp': fzeros(naux, naux), 'KbT_HP': fzeros(nocca, nvirb),
           'KbT_HH': fzeros(nocca, nocca), 'jT': fzeros(naux), 'dq': fzeros(naux)}
    keys = ['LaH', 'LaP', 'LbH', 'LbP', 'Ghh', 'Fhp', 'Yf', 'Sq', 'gpp', 'KbT_HP', 'KbT_HH', 'jT', 'dq']
    libs.lib_func('mrsf_grad_state', (float(cx), farr(Xt)) + tuple(out[k] for k in keys))
    return out


def finish(cx, st, Ta_loc, Xt, Za, Zb, occH, dims):
    """adds the z-dependent and mean-field pieces to st['Ghh'], st['Fhp']
    (in place) and returns the metric matrix g and the Coulomb vector of
    the relaxed density"""
    nocca, nvirb, naux = dims
    g = fzeros(naux, naux); pq = fzeros(naux)
    libs.lib_func('mrsf_grad_finish',
                  (float(cx), farr(st['dq']), farr(st['jT']), farr(Ta_loc), farr(Xt), farr(Za),
                   farr(Zb), st['Sq'], farr(occH), st['Ghh'], st['Fhp'], st['gpp'], g, pq))
    return g, pq


def zvec_setup(nao, lz, fa, fb):
    """rotation space, Fock couplings and diagonal of the Z-vector operator
    in the library; returns the diagonal (lz)"""
    hdiag = fzeros(lz)
    libs.lib_func('mrsf_zvec_setup', (nao, lz, farr(fa), farr(fb), hdiag))
    return hdiag


def zvec_free():
    libs.lib_func('mrsf_zvec_free', ())


def zvec_factors(z, nao, nocca):
    """AO factors C_P Z_s (nao, nocca, 2, nvec) of the trial densities of z (lz, nvec)"""
    z = farr(z)
    nvec = z.shape[1]
    Lf = fzeros(nao, nocca, 2, nvec)
    libs.lib_func('mrsf_zvec_factors', (nvec, z, Lf))
    return Lf


def zvec_hessian(cx, z, VHP, have_xc, lz):
    """H z (lz, nvec); VHP: XC kernel HP blocks (nocca, nvirb, 2, nvec) of the trial densities"""
    z = farr(z)
    nvec = z.shape[1]
    Hz = fzeros(lz, nvec)
    libs.lib_func('mrsf_zvec_hessian', (nvec, float(cx), z, farr(VHP), bool(have_xc), Hz))
    return Hz


def aograd_free():
    libs.lib_func('mrsf_aograd_free', ())


def reffam(cx, dq, occH, dims):
    """B-space families (Ghh, g) of the ROKS reference two-electron energy"""
    nocca, nvirb, naux = dims
    Ghh = fzeros(nocca, nocca, naux); g = fzeros(naux, naux)
    libs.lib_func('mrsf_grad_reffam', (float(cx), farr(dq), farr(occH), Ghh, g))
    return Ghh, g
