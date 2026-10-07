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
!**********************************************************************
module mrsf_zvec

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_gradient, only: gfock
  use mrsf_xcgrid, only: CPg, xc_ready

  implicit none

  integer(is)              :: lz = 0, nao_z = 0
  integer(is), allocatable :: Rp(:), Rq(:), Rploc(:), Rqloc(:)
  logical, allocatable     :: RqisC(:)
  real(dp), allocatable    :: wa(:), wb(:), fa_full(:,:), fb_full(:,:)
  real(dp)                 :: time_zvec = 0.0_dp
  integer(is)              :: nzvec_calls = 0

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
    allocate(Rp(lz), Rq(lz), Rploc(lz), Rqloc(lz), RqisC(lz), wa(lz), wb(lz))
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

  end subroutine zvec_setup

  subroutine zvec_free()

    if (allocated(Rp)) deallocate(Rp, Rq, Rploc, Rqloc, RqisC, wa, wb)
    if (allocated(fa_full)) deallocate(fa_full, fb_full)
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

end module mrsf_zvec_interface
