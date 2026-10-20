"""
EMRSF-Pi (plan Part XVIII): the O1-O2 rotation-invariant parent of the
extended MRSF-TDDFT method. The C->V configurations are built on the
lowest eigenvector Pi = w_G G + w_D D + w_S S of the 3x3 pair block of the
singlet MRSF matrix (G = |C^2 O1^2|, D = |C^2 O2^2|, S = the open-shell
singlet pair = the MRSF OO slot), an invariant many-electron function.
Conventions: the compressed-slot phases of the library (G = slot
(pO1,hO2), D = slot (pO2,hO1), S = the singlet-folded OO slot); the pair
transition densities in that convention are gamma^GS = -sqrt2 on (O1,O2),
gamma^DS = +sqrt2 on (O2,O1), gamma^GD = 0 (pinned against the determinant
prototype ~/calculations/mrsf_dev/inv/emrsf_inv_ref.py).

Covariant one-body operators (c = response exchange fraction):
  CV block   F_Pi = F_ens + c (J - K/2)[Dg_Pi] + 1/2 W[Dg_Pi]
  couplings  F_cpl^XX = F_ens + (J - K/2)[Dg^XX],  Dg^GG = P1 - P2 = -Dg^DD, Dg^SS = 0
  (no XC term in the couplings: the explicit SOMO-pair integrals of the
  per-parent tables mix direct- and exchange-origin pieces, so a kernel
  companion of their Coulomb part has no covariant definition; checked at the
  library level, plan Part XVIII section 9)
with F_ens the closed-shell KS matrix of the ensemble density
D_ens = 2 P_C + P_O1 + P_O2, Dg_Pi the Pi pair density minus its trace part
and W the kernel-contracted (linearised) XC difference at rho_ens.
"""
import numpy as np

SQ2 = np.sqrt(2.0)


def somo_integrals(Bso, O1, O2):
    """(O_a O_b|O_c O_d), a..d in {0: O1, 1: O2}, from the SOMO columns of
    the DF tensor Bso(naux, nmo, 2) = B^Q_{p O_x}"""
    Oi = (O1, O2)

    def B(x, y):
        return Bso[:, Oi[y], x]

    ints = {}
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                for d in (0, 1):
                    ints[(a, b, c, d)] = float(B(a, b) @ B(c, d))
    return ints


def pair_block(fa, fb, O1, O2, ints, chf):
    """3x3 block (G, D, S) of the singlet MRSF matrix (compressed-slot
    convention) from the MO spin Fock matrices and the SOMO integrals"""
    def I(a, b, c, d):
        return ints[(a, b, c, d)]
    K12 = I(0, 1, 0, 1)
    AGG = fb[O1, O1] - fa[O2, O2] - chf*I(1, 1, 0, 0)
    ADD = fb[O2, O2] - fa[O1, O1] - chf*I(0, 0, 1, 1)
    AGD = -chf*K12
    AG11 = -fa[O2, O1] - chf*I(0, 1, 0, 0)
    AG22 = fb[O1, O2] - chf*I(0, 1, 1, 1)
    AD11 = fb[O2, O1] - chf*I(0, 1, 0, 0)
    AD22 = -fa[O1, O2] - chf*I(0, 1, 1, 1)
    A11 = fb[O1, O1] - fa[O1, O1] - chf*I(0, 0, 0, 0)
    A22 = fb[O2, O2] - fa[O2, O2] - chf*I(1, 1, 1, 1)
    AOO = 0.5*(A11 + A22) + chf*K12
    AGS = (AG11 - AG22)/SQ2
    ADS = (AD11 - AD22)/SQ2
    return np.array([[AGG, AGD, AGS], [AGD, ADD, ADS], [AGS, ADS, AOO]])


def pi_weights(Hp):
    """lowest eigenvector (w_G > 0) and the eigenvalues of the pair block"""
    w, v = np.linalg.eigh(0.5*(Hp + Hp.T))
    vec = np.array(v[:, 0])
    if vec[0] < 0.0:
        vec = -vec
    return vec, w


def delta_gamma_pi(w):
    """Pi pair density minus the trace part, SOMO block [[O1O1, O1O2], [O2O1, O2O2]]"""
    d = w[0]**2 - w[1]**2
    o = SQ2*w[2]*(w[1] - w[0])
    return np.array([[d, o], [o, -d]])


def somo_jk(Bso, O1, O2, gam, jblocks, dims, H, P, nmo):
    """J[gamma] and K[gamma] (MO, nmo x nmo) of a SOMO-pair density gamma
    (2x2) from the DF tensor: jq = sum_xy gamma_xy B^Q_{O_x O_y} contracted
    with the resident blocks (jblocks), K_pq = sum_xy gamma_xy sum_Q
    B^Q_{p O_x} B^Q_{q O_y}"""
    Oi = (O1, O2)
    jq = np.zeros(Bso.shape[0])
    K = np.zeros((nmo, nmo))
    for x in (0, 1):
        for y in (0, 1):
            if gam[x, y] == 0.0:
                continue
            jq += gam[x, y]*Bso[:, Oi[y], x]
            K += gam[x, y]*(Bso[:, :, x].T @ Bso[:, :, y])
    JHH, JHP, JPP = jblocks(jq, dims)
    J = np.zeros((nmo, nmo))
    J[np.ix_(H, H)] = JHH[:, :, 0]
    J[np.ix_(P, P)] = JPP[:, :, 0]
    J[np.ix_(H, P)] = JHP[:, :, 0]
    J[np.ix_(P, H)] = JHP[:, :, 0].T
    return J, 0.5*(K + K.T)


def pi_setup(sess, fa, fb, chf, linearised_xc=True):
    """the EMRSF-Pi data of a grid session: dict with w, pair (eigenvalues),
    shift, fone (CV-block operator F_Pi), fcpl (nmo, nmo, 3) and dgam;
    fa, fb: MO spin Fock matrices of the ROKS reference; chf: response
    exchange fraction"""
    import graci.interfaces.mrsf.mrsf_grad as mrsf_grad
    nmo, naux = sess.nmo, sess.naux
    H, P = sess.H, sess.P
    O1, O2 = int(sess.Oidx[0]), int(sess.Oidx[1])
    dims = (sess.nocca, sess.nvirb, naux)
    Bso = mrsf_grad.bso(naux, nmo)
    ints = somo_integrals(Bso, O1, O2)
    Hp = pair_block(fa, fb, O1, O2, ints, chf)
    w, pair = pi_weights(Hp)
    shift = float(pair[0])
    dgam = delta_gamma_pi(w)
    F_ens = sess.fock_ens()
    # F_cpl^XX: HF-type SOMO tensor of the parents G (P1 - P2), D (P2 - P1), S (0)
    Jd, Kd = somo_jk(Bso, O1, O2, np.diag([1.0, -1.0]), mrsf_grad.jblocks, dims, H, P, nmo)
    Td = Jd - 0.5*Kd
    fcpl = np.zeros((nmo, nmo, 3))
    fcpl[:, :, 0] = F_ens + Td
    fcpl[:, :, 1] = F_ens - Td
    fcpl[:, :, 2] = F_ens
    # F_Pi: c-scaled tensor part plus the linearised XC difference
    Jp, Kp = somo_jk(Bso, O1, O2, dgam, mrsf_grad.jblocks, dims, H, P, nmo)
    fone = F_ens + chf*(Jp - 0.5*Kp)
    if linearised_xc and sess.is_dft:
        C = sess.C
        dm = np.zeros((nmo, nmo))
        dm[np.ix_([O1, O2], [O1, O2])] = dgam
        fone = fone + C.T @ sess.fxc_potential_ao(C @ dm @ C.T) @ C
    return {'w': w, 'pair': pair, 'shift': shift, 'fone': fone, 'fcpl': fcpl, 'dgam': dgam, 'Hp': Hp}


def pi_fone_for_weights(sess, pi, w, chf, linearised_xc=True):
    """the CV-block operator for given parent weights (diagnostics)"""
    import graci.interfaces.mrsf.mrsf_grad as mrsf_grad
    nmo, naux = sess.nmo, sess.naux
    H, P = sess.H, sess.P
    O1, O2 = int(sess.Oidx[0]), int(sess.Oidx[1])
    dims = (sess.nocca, sess.nvirb, naux)
    Bso = mrsf_grad.bso(naux, nmo)
    dgam = delta_gamma_pi(np.asarray(w, float))
    Jp, Kp = somo_jk(Bso, O1, O2, dgam, mrsf_grad.jblocks, dims, H, P, nmo)
    fone = sess.fock_ens() + chf*(Jp - 0.5*Kp)
    if linearised_xc and sess.is_dft:
        C = sess.C
        dm = np.zeros((nmo, nmo))
        dm[np.ix_([O1, O2], [O1, O2])] = dgam
        fone = fone + C.T @ sess.fxc_potential_ao(C @ dm @ C.T) @ C
    return fone
