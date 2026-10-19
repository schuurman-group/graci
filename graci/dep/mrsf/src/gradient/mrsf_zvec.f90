!**********************************************************************
! mrsf_zvec: the Z-vector (ROKS orbital Hessian) operator of the MRSF
! gradient. The rotation space R = {C->O, C->V, O->V} is ordered as in
! the Python driver (for q in C: p in O; for q in C: p in V; for q in O:
! p in V), with the occupation weights w^s = n^s_q - n^s_p. For a batch
! of vectors z(lz, nvec):
!   (H z)_pq = 2 sum_s w^s_pq [ (f^s K - K f^s)_pq + G^s[p^z,s]_qp ],
!   K_pq = z_pq, K_qp = -z_pq,  p^{z,s} = sym(w^s z) (local Za, Zb),
! where G^s = J - x_ref K + f_xc is evaluated by gfock (two-electron
! part) and by the grid module (XC part, passed in as the HP blocks VHP);
! the beta part only contributes for q in C. zvec_factors forms the AO
! factors C_P Z_s of the trial densities for the XC kernel.
! zvec_apply is the fused product used by the solvers: weighted local
! matrices built once, the XC kernel with spin-shared factor columns
! (xc_zvec, cached AO values), the two-electron part with the particle-
! hole exchange as Q-blocked BLAS-3 products and one hole-particle
! exchange sweep for both spins (zvec_twoe), and the Fock couplings on
! the rotation blocks only (K has the blocks K_OC = z_CO, K_VC = z_CV,
! K_VO = z_OV and their negative transposes):
!   (fK - Kf)_OC = f_OO z_CO + f_OV z_CV - z_CO f_CC + z_OV^T f_VC
!   (fK - Kf)_VC = f_VO z_CO + f_VV z_CV - z_CV f_CC - z_OV f_OC
!   (fK - Kf)_VO = -f_VC z_CO^T + f_VV z_OV - z_CV f_CO - z_OV f_OO
! (beta for C->O, alpha + beta for C->V, alpha for O->V).
!**********************************************************************
module mrsf_zvec

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_gradient, only: gfock, kx_pp_hh, Bhp
  use mrsf_xcgrid, only: CPg, xc_ready, xc_zvec

  implicit none

  integer(is)              :: lz = 0, nao_z = 0
  integer(is), allocatable :: Rp(:), Rq(:), Rploc(:), Rqloc(:)
  logical, allocatable     :: RqisC(:)
  real(dp), allocatable    :: wa(:), wb(:), fa_full(:,:), fb_full(:,:)
  real(dp), allocatable    :: hdg(:)              ! the Fock diagonal of the operator (preconditioner)
  real(dp)                 :: time_zvec = 0.0_dp
  integer(is)              :: nzvec_calls = 0
  ! Fock blocks of the fused product (f^b for C->O, f^a + f^b for C->V,
  ! f^a for O->V; MO order of the classes C, O, V)
  real(dp), allocatable    :: fbOO(:,:), fbOV(:,:), fbCC(:,:), fbVC(:,:)
  real(dp), allocatable    :: fsVO(:,:), fsVV(:,:), fsCC(:,:), fsOC(:,:)
  real(dp), allocatable    :: faVC(:,:), faVV(:,:), faCO(:,:), faOO(:,:)
  ! fused-product timers: local matrices, Coulomb, particle-hole exchange,
  ! hole-particle exchange, XC kernel, Fock part + assembly; calls, vectors
  real(dp)                 :: tap(6) = 0.0_dp
  integer(is)              :: napply = 0, napply_vecs = 0

contains

!######################################################################
! zvec_setup: rotation space, weights, Fock matrices (MO, nmo x nmo)
! and the diagonal of the operator
!######################################################################
  subroutine zvec_setup(nao1, lz1, fa, fb, hdiag)

    integer(is), intent(in) :: nao1, lz1
    real(dp), intent(in)    :: fa(nmo,nmo), fb(nmo,nmo)
    real(dp), intent(out)   :: hdiag(lz1)

    integer(is) :: k, q, a, x, p, qq
    real(dp) :: na_p, na_q, nb_p, nb_q

    if (.not. init_done) call mrsf_error('zvec_setup: library not initialised')
    call zvec_free()
    nao_z = nao1
    lz = nC*2 + nC*nV + 2*nV
    if (lz /= lz1) call mrsf_error('zvec_setup: rotation-space dimension mismatch')
    allocate(Rp(lz), Rq(lz), Rploc(lz), Rqloc(lz), RqisC(lz), wa(lz), wb(lz), hdg(lz))
    allocate(fa_full(nmo,nmo), fb_full(nmo,nmo))
    fa_full = fa; fb_full = fb
    k = 0
    do q = 1, nC
       do x = 1, 2
          k = k + 1; Rploc(k) = x; Rqloc(k) = q
       enddo
    enddo
    do q = 1, nC
       do a = 3, nvirb
          k = k + 1; Rploc(k) = a; Rqloc(k) = q
       enddo
    enddo
    do x = 1, 2
       do a = 3, nvirb
          k = k + 1; Rploc(k) = a; Rqloc(k) = nC + x
       enddo
    enddo
    do k = 1, lz
       p = Pmap(Rploc(k)); qq = Hmap(Rqloc(k))
       Rp(k) = p; Rq(k) = qq
       RqisC(k) = (Rqloc(k) <= nC)
       na_p = merge(1.0_dp, 0.0_dp, occ(p) > 0.5_dp); na_q = merge(1.0_dp, 0.0_dp, occ(qq) > 0.5_dp)
       nb_p = merge(1.0_dp, 0.0_dp, occ(p) > 1.5_dp); nb_q = merge(1.0_dp, 0.0_dp, occ(qq) > 1.5_dp)
       wa(k) = na_q - na_p
       wb(k) = nb_q - nb_p
       hdiag(k) = 2.0_dp*(wa(k)*(fa(p,p) - fa(qq,qq)) + wb(k)*(fb(p,p) - fb(qq,qq)))
    enddo
    hdg = hdiag
    call fock_blocks_setup(fa, fb)

  end subroutine zvec_setup

!######################################################################
! fock_blocks_setup: the Fock blocks of the fused product, classes C =
! Hmap(1:nC), O = Pmap(1:2) (= Hmap(nC+1:nC+2)), V = Pmap(3:nvirb)
!######################################################################
  subroutine fock_blocks_setup(fa, fb)

    real(dp), intent(in) :: fa(nmo,nmo), fb(nmo,nmo)

    integer(is) :: Cm(nC), Om(2), Vm(nV)

    Cm = Hmap(1:nC); Om = Pmap(1:2); Vm = Pmap(3:nvirb)
    allocate(fbOO(2,2), fbOV(2,nV), fbCC(nC,nC), fbVC(nV,nC))
    allocate(fsVO(nV,2), fsVV(nV,nV), fsCC(nC,nC), fsOC(2,nC))
    allocate(faVC(nV,nC), faVV(nV,nV), faCO(nC,2), faOO(2,2))
    fbOO = fb(Om,Om); fbOV = fb(Om,Vm); fbCC = fb(Cm,Cm); fbVC = fb(Vm,Cm)
    fsVO = fa(Vm,Om) + fb(Vm,Om); fsVV = fa(Vm,Vm) + fb(Vm,Vm)
    fsCC = fa(Cm,Cm) + fb(Cm,Cm); fsOC = fa(Om,Cm) + fb(Om,Cm)
    faVC = fa(Vm,Cm); faVV = fa(Vm,Vm); faCO = fa(Cm,Om); faOO = fa(Om,Om)

  end subroutine fock_blocks_setup

  subroutine zvec_free()

    if (allocated(Rp)) deallocate(Rp, Rq, Rploc, Rqloc, RqisC, wa, wb, hdg)
    if (allocated(fa_full)) deallocate(fa_full, fb_full)
    if (allocated(fbOO)) deallocate(fbOO, fbOV, fbCC, fbVC, fsVO, fsVV, fsCC, fsOC, faVC, faVV, faCO, faOO)
    lz = 0

  end subroutine zvec_free

!######################################################################
! zvec_local: Za, Zb(nvirb, nocca, nvec) = w^s z scattered on the
! particle x hole blocks
!######################################################################
  subroutine zvec_local(nvec, z, Za, Zb)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: z(lz,nvec)
    real(dp), intent(out)   :: Za(nvirb,nocca,nvec), Zb(nvirb,nocca,nvec)

    integer(is) :: v, k

    Za = 0.0_dp; Zb = 0.0_dp
    do v = 1, nvec
       do k = 1, lz
          Za(Rploc(k),Rqloc(k),v) = wa(k)*z(k,v)
          Zb(Rploc(k),Rqloc(k),v) = wb(k)*z(k,v)
       enddo
    enddo

  end subroutine zvec_local

!######################################################################
! zvec_factors: AO factors L^s_v = C_P Z^s_v (nao, nocca, 2, nvec) of
! the trial densities p^{z,s} = L C_H^T + C_H L^T for the XC kernel
!######################################################################
  subroutine zvec_factors(nvec, z, Lf)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: z(lz,nvec)
    real(dp), intent(out)   :: Lf(nao_z,nocca,2,nvec)

    real(dp), allocatable :: Za(:,:,:), Zb(:,:,:)
    integer(is) :: v

    if (lz == 0) call mrsf_error('zvec_factors: call zvec_setup first')
    if (.not. xc_ready) call mrsf_error('zvec_factors: grid module not initialised')
    allocate(Za(nvirb,nocca,nvec), Zb(nvirb,nocca,nvec))
    call zvec_local(nvec, z, Za, Zb)
    do v = 1, nvec
       call dgemm('N','N', nao_z, nocca, nvirb, 1.0_dp, CPg, nao_z, Za(1,1,v), nvirb, 0.0_dp, Lf(1,1,1,v), nao_z)
       call dgemm('N','N', nao_z, nocca, nvirb, 1.0_dp, CPg, nao_z, Zb(1,1,v), nvirb, 0.0_dp, Lf(1,1,2,v), nao_z)
    enddo
    deallocate(Za, Zb)

  end subroutine zvec_factors

!######################################################################
! zvec_hessian: H z for a batch of vectors; VHP(nocca, nvirb, 2, nvec)
! holds the XC kernel HP blocks of the trial densities (used when
! have_xc), cx is the reference exchange fraction
!######################################################################
  subroutine zvec_hessian(nvec, cx, z, VHP, have_xc, Hz)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: cx, z(lz,nvec), VHP(nocca,nvirb,2,nvec)
    logical, intent(in)     :: have_xc
    real(dp), intent(out)   :: Hz(lz,nvec)

    real(dp), allocatable :: Za(:,:,:), Zb(:,:,:), Ahh(:,:,:), jq(:,:)
    real(dp), allocatable :: Ghpa(:,:,:), Ghpb(:,:,:), Ghha(:,:,:), Ghhb(:,:,:)
    real(dp), allocatable :: Kmat(:,:), dfa(:,:), dfb(:,:)
    integer(is) :: v, k, p, q
    real(dp) :: da, db, t0

    if (lz == 0) call mrsf_error('zvec_hessian: call zvec_setup first')
    t0 = wall_time()
    allocate(Za(nvirb,nocca,nvec), Zb(nvirb,nocca,nvec), Ahh(nocca,nocca,nvec), jq(naux,nvec))
    allocate(Ghpa(nocca,nvirb,nvec), Ghpb(nocca,nvirb,nvec), Ghha(nocca,nocca,nvec), Ghhb(nocca,nocca,nvec))
    call zvec_local(nvec, z, Za, Zb)
    ! two-electron response: J - cx K of the trial densities on the HP blocks
    Ahh = 0.0_dp; jq = 0.0_dp
    call gfock(nvec, cx, 2+4+16, Ahh, Za, Zb, jq, Ghpa, Ghpb, Ghha, Ghhb)
    if (have_xc) then
       Ghpa = Ghpa + VHP(:,:,1,:)
       Ghpb = Ghpb + VHP(:,:,2,:)
    endif
    ! Fock couplings and assembly
    allocate(Kmat(nmo,nmo), dfa(nmo,nmo), dfb(nmo,nmo))
    do v = 1, nvec
       Kmat = 0.0_dp
       do k = 1, lz
          Kmat(Rp(k),Rq(k)) = z(k,v)
          Kmat(Rq(k),Rp(k)) = -z(k,v)
       enddo
       call dgemm('N','N', nmo, nmo, nmo, 1.0_dp, fa_full, nmo, Kmat, nmo, 0.0_dp, dfa, nmo)
       call dgemm('N','N', nmo, nmo, nmo, -1.0_dp, Kmat, nmo, fa_full, nmo, 1.0_dp, dfa, nmo)
       call dgemm('N','N', nmo, nmo, nmo, 1.0_dp, fb_full, nmo, Kmat, nmo, 0.0_dp, dfb, nmo)
       call dgemm('N','N', nmo, nmo, nmo, -1.0_dp, Kmat, nmo, fb_full, nmo, 1.0_dp, dfb, nmo)
       do k = 1, lz
          p = Rp(k); q = Rq(k)
          da = dfa(p,q) + Ghpa(Rqloc(k),Rploc(k),v)
          db = dfb(p,q)
          if (RqisC(k)) db = db + Ghpb(Rqloc(k),Rploc(k),v)
          Hz(k,v) = 2.0_dp*(wa(k)*da + wb(k)*db)
       enddo
    enddo
    deallocate(Za, Zb, Ahh, jq, Ghpa, Ghpb, Ghha, Ghhb, Kmat, dfa, dfb)
    time_zvec = time_zvec + wall_time() - t0
    nzvec_calls = nzvec_calls + 1

  end subroutine zvec_hessian

!######################################################################
! zab_from_z: weighted local matrices Zab(nvirb, nocca, 2, nvec) of a
! batch of rotation vectors (block copies; rotation order C->O, C->V,
! O->V with z_CO(2,nC), z_CV(nV,nC), z_OV(nV,2) column-major)
!######################################################################
  subroutine zab_from_z(nvec, z, Zab)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: z(lz,nvec)
    real(dp), intent(out)   :: Zab(nvirb,nocca,2,nvec)

    integer(is) :: v, q, x, oCV, oOV

    oCV = 2*nC; oOV = 2*nC + nC*nV
    Zab = 0.0_dp
    do v = 1, nvec
       do q = 1, nC
          Zab(1:2,q,2,v)      = z(2*(q-1)+1:2*q,v)                 ! z_CO (beta)
          Zab(3:nvirb,q,1,v)  = z(oCV+nV*(q-1)+1:oCV+nV*q,v)       ! z_CV (alpha)
          Zab(3:nvirb,q,2,v)  = z(oCV+nV*(q-1)+1:oCV+nV*q,v)       ! z_CV (beta)
       enddo
       do x = 1, 2
          Zab(3:nvirb,nC+x,1,v) = z(oOV+nV*(x-1)+1:oOV+nV*x,v)     ! z_OV (alpha)
       enddo
    enddo

  end subroutine zab_from_z

!######################################################################
! zvec_twoe: G^s_HP = J[D^a + D^b] - cx K[D^s] on the HP blocks for the
! trial densities of Zab (as gfock with flags 2+4+16): Coulomb by two
! dgemms; particle-hole exchange sum_Q (B^Q_HP Z^s) B^Q_HP for all
! vectors and spins in Q blocks (Bp((i,q),a) = B^{Q0+q}_{ia}: W = Bp Z,
! then per (s,v) K += Wt Bp with K = nocca*nq); hole-particle exchange
! [sum_Q B^Q_PP Z B^Q_HH]^T in one kx_pp_hh call with 2 nvec columns
!######################################################################
  subroutine zvec_twoe(nvec, cx, Zab, Ghp)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: cx, Zab(nvirb,nocca,2,nvec)
    real(dp), intent(out)   :: Ghp(nocca,nvirb,2,nvec)

    integer(is) :: v, s, i, a, j, q, Q0, nQb, nq, sv, ns
    real(dp)    :: t0
    real(dp), allocatable :: ZT(:,:,:), jq(:,:), JHP(:,:,:), Bp(:,:), W(:,:), Wt(:,:), Sz(:,:,:,:)

    ns = 2*nvec
    ! Coulomb
    t0 = wall_time()
    allocate(ZT(nocca,nvirb,nvec), jq(naux,nvec), JHP(nocca,nvirb,nvec))
    do v = 1, nvec
       do a = 1, nvirb
          do i = 1, nocca
             ZT(i,a,v) = Zab(a,i,1,v) + Zab(a,i,2,v)
          enddo
       enddo
    enddo
    call dgemm('T','N', naux, nvec, nocca*nvirb, 2.0_dp, Bhp, nocca*nvirb, ZT, nocca*nvirb, 0.0_dp, jq, naux)
    call dgemm('N','N', nocca*nvirb, nvec, naux, 1.0_dp, Bhp, nocca*nvirb, jq, naux, 0.0_dp, JHP, nocca*nvirb)
    deallocate(ZT, jq)
    tap(2) = tap(2) + wall_time() - t0

    ! particle-hole exchange (accumulated in Ghp)
    t0 = wall_time()
    Ghp = 0.0_dp
    nQb = max(1, min(naux, (3072 + nocca - 1)/nocca))
    allocate(Bp(nocca*nQb,nvirb), W(nocca*nQb,nocca*ns), Wt(nocca,nocca*nQb))
    do Q0 = 0, naux-1, nQb
       nq = min(nQb, naux - Q0)
       do a = 1, nvirb
          do q = 1, nq
             Bp(nocca*(q-1)+1:nocca*q,a) = Bhp(1:nocca,a,Q0+q)
          enddo
       enddo
       call dgemm('N','N', nocca*nq, nocca*ns, nvirb, 1.0_dp, Bp, nocca*nQb, Zab, nvirb, 0.0_dp, &
            W, nocca*nQb)
       do sv = 1, ns
          do q = 1, nq
             do j = 1, nocca
                Wt(1:nocca,j+nocca*(q-1)) = W(nocca*(q-1)+1:nocca*q,j+nocca*(sv-1))
             enddo
          enddo
          s = mod(sv-1, 2) + 1; v = (sv-1)/2 + 1
          call dgemm('N','N', nocca, nvirb, nocca*nq, 1.0_dp, Wt, nocca, Bp, nocca*nQb, 1.0_dp, &
               Ghp(1,1,s,v), nocca)
       enddo
    enddo
    deallocate(Bp, W, Wt)
    tap(3) = tap(3) + wall_time() - t0

    ! hole-particle exchange, both spins in one sweep
    t0 = wall_time()
    allocate(Sz(nvirb,nocca,2,nvec))
    Sz = 0.0_dp
    call kx_pp_hh(ns, 1.0_dp, Zab, Sz)
    do v = 1, nvec
       do s = 1, 2
          Ghp(:,:,s,v) = JHP(:,:,v) - cx*(Ghp(:,:,s,v) + transpose(Sz(:,:,s,v)))
       enddo
    enddo
    deallocate(Sz, JHP)
    tap(4) = tap(4) + wall_time() - t0

  end subroutine zvec_twoe

!######################################################################
! fock_rot: the Fock couplings of one rotation vector on the rotation
! blocks: bOC = (f^bK - Kf^b)_OC, bVC = (f^sK - Kf^s)_VC with f^s = f^a +
! f^b, bVO = (f^aK - Kf^a)_VO
!######################################################################
  subroutine fock_rot(zCO, zCV, zOV, bOC, bVC, bVO)

    real(dp), intent(in)  :: zCO(2,nC), zCV(nV,nC), zOV(nV,2)
    real(dp), intent(out) :: bOC(2,nC), bVC(nV,nC), bVO(nV,2)

    ! (f^b K - K f^b)_OC = f_OO z_CO + f_OV z_CV - z_CO f_CC + z_OV^T f_VC
    call dgemm('N','N', 2, nC, 2, 1.0_dp, fbOO, 2, zCO, 2, 0.0_dp, bOC, 2)
    call dgemm('N','N', 2, nC, nV, 1.0_dp, fbOV, 2, zCV, nV, 1.0_dp, bOC, 2)
    call dgemm('N','N', 2, nC, nC, -1.0_dp, zCO, 2, fbCC, nC, 1.0_dp, bOC, 2)
    call dgemm('T','N', 2, nC, nV, 1.0_dp, zOV, nV, fbVC, nV, 1.0_dp, bOC, 2)
    ! (f K - K f)_VC = f_VO z_CO + f_VV z_CV - z_CV f_CC - z_OV f_OC
    call dgemm('N','N', nV, nC, 2, 1.0_dp, fsVO, nV, zCO, 2, 0.0_dp, bVC, nV)
    call dgemm('N','N', nV, nC, nV, 1.0_dp, fsVV, nV, zCV, nV, 1.0_dp, bVC, nV)
    call dgemm('N','N', nV, nC, nC, -1.0_dp, zCV, nV, fsCC, nC, 1.0_dp, bVC, nV)
    call dgemm('N','N', nV, nC, 2, -1.0_dp, zOV, nV, fsOC, 2, 1.0_dp, bVC, nV)
    ! (f^a K - K f^a)_VO = -f_VC z_CO^T + f_VV z_OV - z_CV f_CO - z_OV f_OO
    call dgemm('N','T', nV, 2, nC, -1.0_dp, faVC, nV, zCO, 2, 0.0_dp, bVO, nV)
    call dgemm('N','N', nV, 2, nV, 1.0_dp, faVV, nV, zOV, nV, 1.0_dp, bVO, nV)
    call dgemm('N','N', nV, 2, nC, -1.0_dp, zCV, nV, faCO, nC, 1.0_dp, bVO, nV)
    call dgemm('N','N', nV, 2, 2, -1.0_dp, zOV, nV, faOO, 2, 1.0_dp, bVO, nV)

  end subroutine fock_rot

!######################################################################
! zvec_apply: the fused operator product Hz(lz, nvec) (same result as
! zvec_factors + XC kernel + zvec_hessian); the XC part requires the
! cached AO values of the grid module (have_xc)
!######################################################################
  subroutine zvec_apply(nvec, cx, have_xc, z, Hz)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: cx, z(lz,nvec)
    logical, intent(in)     :: have_xc
    real(dp), intent(out)   :: Hz(lz,nvec)

    integer(is) :: v, q, x, a, oCV, oOV
    real(dp)    :: t0, t1
    real(dp), allocatable :: Zab(:,:,:,:), Ghp(:,:,:,:), VHP(:,:,:,:), bOC(:,:), bVC(:,:), bVO(:,:)

    if (lz == 0) call mrsf_error('zvec_apply: call zvec_setup first')
    t1 = wall_time()
    oCV = 2*nC; oOV = 2*nC + nC*nV
    t0 = wall_time()
    allocate(Zab(nvirb,nocca,2,nvec), Ghp(nocca,nvirb,2,nvec))
    call zab_from_z(nvec, z, Zab)
    tap(1) = tap(1) + wall_time() - t0
    call zvec_twoe(nvec, cx, Zab, Ghp)
    if (have_xc) then
       t0 = wall_time()
       allocate(VHP(nocca,nvirb,2,nvec))
       call xc_zvec(nvec, Zab, VHP)
       Ghp = Ghp + VHP
       deallocate(VHP)
       tap(5) = tap(5) + wall_time() - t0
    endif
    t0 = wall_time()
    allocate(bOC(2,nC), bVC(nV,nC), bVO(nV,2))
    do v = 1, nvec
       call fock_rot(z(1,v), z(oCV+1,v), z(oOV+1,v), bOC, bVC, bVO)
       do q = 1, nC
          do x = 1, 2
             Hz(2*(q-1)+x,v) = 2.0_dp*(bOC(x,q) + Ghp(q,x,2,v))
          enddo
          do a = 1, nV
             Hz(oCV+nV*(q-1)+a,v) = 2.0_dp*(bVC(a,q) + Ghp(q,2+a,1,v) + Ghp(q,2+a,2,v))
          enddo
       enddo
       do x = 1, 2
          do a = 1, nV
             Hz(oOV+nV*(x-1)+a,v) = 2.0_dp*(bVO(a,x) + Ghp(nC+x,2+a,1,v))
          enddo
       enddo
    enddo
    deallocate(Zab, Ghp, bOC, bVC, bVO)
    tap(6) = tap(6) + wall_time() - t0
    time_zvec = time_zvec + wall_time() - t1
    napply = napply + 1
    napply_vecs = napply_vecs + nvec

  end subroutine zvec_apply

end module mrsf_zvec

!**********************************************************************
! C-bound entry points
!**********************************************************************
module mrsf_zvec_interface

  use iso_c_binding
  use mrsf_constants
  use mrsf_global
  use mrsf_zvec

  implicit none

contains

  subroutine mrsf_zvec_setup_c(nao1, lz1, fa, fb, hdiag) bind(c, name='mrsf_zvec_setup')
    integer(is), intent(in) :: nao1, lz1
    real(dp), intent(in)    :: fa(nmo,nmo), fb(nmo,nmo)
    real(dp), intent(out)   :: hdiag(lz1)
    call zvec_setup(nao1, lz1, fa, fb, hdiag)
  end subroutine mrsf_zvec_setup_c

  subroutine mrsf_zvec_free_c() bind(c, name='mrsf_zvec_free')
    call zvec_free()
  end subroutine mrsf_zvec_free_c

  subroutine mrsf_zvec_factors_c(nvec, z, Lf) bind(c, name='mrsf_zvec_factors')
    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: z(lz,nvec)
    real(dp), intent(out)   :: Lf(nao_z,nocca,2,nvec)
    call zvec_factors(nvec, z, Lf)
  end subroutine mrsf_zvec_factors_c

  subroutine mrsf_zvec_hessian_c(nvec, cx, z, VHP, have_xc, Hz) bind(c, name='mrsf_zvec_hessian')
    integer(is), intent(in)     :: nvec
    real(dp), intent(in)        :: cx, z(lz,nvec), VHP(nocca,nvirb,2,nvec)
    logical(c_bool), intent(in) :: have_xc
    real(dp), intent(out)       :: Hz(lz,nvec)
    call zvec_hessian(nvec, cx, z, VHP, logical(have_xc), Hz)
  end subroutine mrsf_zvec_hessian_c

  subroutine mrsf_zvec_apply_c(nvec, cx, have_xc, z, Hz) bind(c, name='mrsf_zvec_apply')
    integer(is), intent(in)     :: nvec
    real(dp), intent(in)        :: cx, z(lz,nvec)
    logical(c_bool), intent(in) :: have_xc
    real(dp), intent(out)       :: Hz(lz,nvec)
    call zvec_apply(nvec, cx, logical(have_xc), z, Hz)
  end subroutine mrsf_zvec_apply_c

  ! timers of the fused product: local matrices, Coulomb, particle-hole
  ! exchange, hole-particle exchange, XC kernel, Fock + assembly, calls,
  ! vectors; XC kernel sub-timers (thread seconds) in t(9:13)
  subroutine mrsf_zvec_timings_c(t) bind(c, name='mrsf_zvec_timings')
    use mrsf_xcgrid, only: tzv
    real(dp), intent(out) :: t(16)
    t = 0.0_dp
    t(1:6) = tap
    t(7) = real(napply, dp)
    t(8) = real(napply_vecs, dp)
    t(9:13) = tzv
  end subroutine mrsf_zvec_timings_c

  subroutine mrsf_zvec_reset_timings_c() bind(c, name='mrsf_zvec_reset_timings')
    use mrsf_xcgrid, only: tzv
    tap = 0.0_dp; napply = 0; napply_vecs = 0; tzv = 0.0_dp
  end subroutine mrsf_zvec_reset_timings_c

end module mrsf_zvec_interface
