!**********************************************************************
! mrsf_gradient: MO-space density-fitted kernels for the analytic
! RO-MRSF-TDDFT nuclear gradient
!
! All quantities are in the local layouts of the energy code:
!   holes     i in H = C + [O1,O2]   (nocca),  particles a in P = [O1,O2] + V (nvirb)
!   B^Q_{HH} = BooQ(nocca,nocca,naux), B^Q_{HP} = Bhp(nocca,nvirb,naux),
!   B^Q_{PP} = paired planes Bvv (+ Dall diagonals) of the energy code.
! "Families" are the B-space three-index densities
!   Gamma~^Q = [HH: Ghh(:,:,Q)] + [HP+PH: F(:,:,Q), both orientations
!   summed] + [PP: Yf^Q X~^T + X~ Yf^Q^T],
! from which the gradient driver forms Gamma^P = L^-T Gamma~ and
! gamma = L^-T g L^-1 with g = <B^Q, Gamma~^Q'>.
! The gradient path requires double-precision integral storage.
!**********************************************************************
module mrsf_gradient

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_space
  use mrsf_sigma, only: ensure_work, plane_apply
  use mrsf_etensor, only: etensor_kx, etensor_ready, etensor_nfc

  implicit none

  real(dp), allocatable :: Bhp(:,:,:)    ! (nocca,nvirb,naux)  B^Q_{Hmap(i),Pmap(a)}
  real(dp), allocatable :: BooQ(:,:,:)   ! (nocca,nocca,naux)  plane-contiguous B^Q_{HH}
  logical               :: grad_ready = .false.

contains

!######################################################################
! grad_init: second streaming pass over the Ao2mo file collecting the
! hole-particle block, plus the plane-contiguous hole-hole copy
!######################################################################
  subroutine grad_init(fname, ierr)

    character(len=*), intent(in) :: fname
    integer(is), intent(out)     :: ierr

    integer(is)           :: unit, dims(2), nrec, cpr, irec, c0, c1, ncol
    integer(is)           :: n_ij, ip, iq, ij, k, i, a, i2, a2, Q, blk, Ql
    integer(is), allocatable :: ptab(:), qtab(:)
    real(dp), allocatable :: buf(:,:)
    logical               :: exists
    real(dp)              :: t0

    ierr = 0
    t0 = wall_time()
    if (.not. ints_loaded) then
       ierr = 1; return
    endif
    if (store_sp) then
       ierr = 2; return
    endif
    inquire(file=trim(fname), exist=exists)
    if (.not. exists) then
       ierr = 3; return
    endif

    call grad_free()
    allocate(Bhp(nocca,nvirb,naux), source=0.0_dp)
    allocate(BooQ(nocca,nocca,naux), source=0.0_dp)

    ! plane-contiguous copy of the Q-blocked hole-hole integrals
    do Q = 1, naux
       blk = (Q-1)/nQ + 1
       Ql  = Q - (blk-1)*nQ
       BooQ(:,:,Q) = Boo(Ql,:,:,blk)
    enddo

    ! pair tables (1-based, p >= q, ij = p(p-1)/2 + q)
    n_ij = nmo*(nmo+1)/2
    allocate(ptab(n_ij), qtab(n_ij))
    ij = 0
    do ip = 1, nmo
       do iq = 1, ip
          ij = ij + 1
          ptab(ij) = ip
          qtab(ij) = iq
       enddo
    enddo

    call freeunit(unit)
    open(unit, file=trim(fname), form='unformatted', status='old')
    read(unit) dims(1)
    read(unit) dims(2)
    if (dims(1) /= naux .or. dims(2) /= n_ij) then
       close(unit); ierr = 4; return
    endif
    read(unit) nrec
    read(unit) cpr
    allocate(buf(naux,cpr))
    do irec = 1, nrec
       c0   = (irec-1)*cpr + 1
       c1   = min(irec*cpr, n_ij)
       ncol = c1 - c0 + 1
       read(unit) buf(1:naux, 1:ncol)
       !$omp parallel do private(ij, ip, iq, i, a, i2, a2)
       do k = 1, ncol
          ij = c0 + k - 1
          ip = ptab(ij); iq = qtab(ij)
          i = Hinv(ip); a = Pinv(iq)
          if (i > 0 .and. a > 0) Bhp(i,a,:) = buf(:,k)
          i2 = Hinv(iq); a2 = Pinv(ip)
          if (i2 > 0 .and. a2 > 0) Bhp(i2,a2,:) = buf(:,k)
       enddo
       !$omp end parallel do
    enddo
    close(unit)
    deallocate(buf, ptab, qtab)

    grad_ready = .true.
    if (verbose) write(6,'(/,2x,a,f10.2,a,f8.3,a)') 'MRSF gradient blocks loaded in ', &
         wall_time() - t0, ' s (Bhp ', 8.0_dp*nocca*nvirb*naux/1.0e9_dp, ' GB)'

  end subroutine grad_init

  subroutine grad_free()
    if (allocated(Bhp)) deallocate(Bhp)
    if (allocated(BooQ)) deallocate(BooQ)
    grad_ready = .false.
  end subroutine grad_free

!######################################################################
! grad_bso: B^Q_{p,O_x} for all MOs p (for the pairing terms)
!######################################################################
  subroutine grad_bso(Bso)

    real(dp), intent(out) :: Bso(naux,nmo,2)
    integer(is) :: x, p, hO

    do x = 1, 2
       hO = nC + x
       do p = 1, nmo
          if (Hinv(p) > 0) then
             Bso(:,p,x) = BooQ(Hinv(p),hO,:)
          else
             Bso(:,p,x) = Bhp(hO,Pinv(p),:)
          endif
       enddo
    enddo

  end subroutine grad_bso

!######################################################################
! kx_pp_hh: S(:,:,v) += scale * sum_Q B^Q_PP Z(:,:,v) B^Q_HH
! (the sigma exchange kernel with a general particle-hole matrix)
!######################################################################
  subroutine kx_pp_hh(nvec, scale, Z, S)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: scale
    real(dp), intent(in)    :: Z(nvirb,nocca,nvec)
    real(dp), intent(inout) :: S(nvirb,nocca,nvec)

    integer(is)           :: blk, Q0, Ql, Q, kk, k, v, ldT, c0, j0, w0
    real(dp), allocatable :: dd(:)

    if (use_etensor .and. etensor_ready .and. etensor_nfc == 0) then
       call etensor_kx(nvec, scale, Z, S)
       return
    endif
    call ensure_work(nvec)
    allocate(dd(nvirb))
    ldT = nvirb * nQ
    do blk = 1, nblk
       Q0 = (blk-1)*nQ
       if (vv_full) then
          do Ql = 1, nQ
             Q = Q0 + Ql
             if (Q > naux) exit
             call dgemm('N','N', nvirb, nocca*nvec, nvirb, 1.0_dp, Bvv(1,1,Q), nvirb, &
                  Z, nvirb, 0.0_dp, Twork(1,Ql,1,1), ldT)
          enddo
       else
          do kk = 1, npblk
             k = (blk-1)*npblk + kk
             if (k > nplane) exit
             Q  = 2*k - 1
             Ql = 2*kk - 1
             call plane_apply(nocca*nvec, Bvv(1,1,k), Z, Q, Ql, ldT, dd)
          enddo
       endif
       ! step 1 wrote the nocca*nvec columns contiguously (flat column index
       ! over the last two dimensions of Twork, whose third extent is ncol,
       ! not nocca, in an extended-method session): vector v starts at the
       ! flat column nocca*(v-1)+1
       do v = 1, nvec
          c0 = nocca*(v-1)
          j0 = mod(c0, ncol) + 1
          w0 = c0/ncol + 1
          call dgemm('N','N', nvirb, nocca, nQ*nocca, scale, Twork(1,1,j0,w0), nvirb, &
               Boo(1,1,1,blk), nQ*nocca, 1.0_dp, S(1,1,v), nvirb)
       enddo
    enddo
    deallocate(dd)

  end subroutine kx_pp_hh

!######################################################################
! gfock: Coulomb and exchange parts of the spin Fock response
!   G^s[D] = J[D^a + D^b] - cx K[D^s]
! for the symmetric MO densities
!   D^a = Ahh (on HxH) + Za (on PxH) + Za^T (on HxP)
!   D^b =                Zb (on PxH) + Zb^T (on HxP)
! plus an externally supplied Coulomb vector jq_add (e.g. of a PP
! density).  Outputs on the HP blocks (i in H, a in P) and, if wanted,
! the HH blocks.  flags bits: 1 Ahh, 2 Za, 4 Zb, 8 jq_add, 16 want HP,
! 32 want HH.
!######################################################################
  subroutine gfock(nvec, cx, flags, Ahh, Za, Zb, jq_add, Ghpa, Ghpb, Ghha, Ghhb)

    integer(is), intent(in) :: nvec, flags
    real(dp), intent(in)    :: cx
    real(dp), intent(in)    :: Ahh(nocca,nocca,nvec), Za(nvirb,nocca,nvec), &
                               Zb(nvirb,nocca,nvec), jq_add(naux,nvec)
    real(dp), intent(out)   :: Ghpa(nocca,nvirb,nvec), Ghpb(nocca,nvirb,nvec), &
                               Ghha(nocca,nocca,nvec), Ghhb(nocca,nocca,nvec)

    logical :: hasA, hasZa, hasZb, hasJ, wantHP, wantHH
    integer(is) :: v, Q, i, a
    real(dp), allocatable :: jq(:,:), ZT(:,:,:), JHP(:,:,:), JHH(:,:,:)
    real(dp), allocatable :: KaHP(:,:,:), KbHP(:,:,:), KaHH(:,:,:), KbHH(:,:,:)
    real(dp), allocatable :: W1(:,:), W2(:,:), Sz(:,:,:)

    hasA   = iand(flags, 1_is) /= 0
    hasZa  = iand(flags, 2_is) /= 0
    hasZb  = iand(flags, 4_is) /= 0
    hasJ   = iand(flags, 8_is) /= 0
    wantHP = iand(flags, 16_is) /= 0
    wantHH = iand(flags, 32_is) /= 0

    allocate(jq(naux,nvec), source=0.0_dp)
    allocate(ZT(nocca,nvirb,nvec), source=0.0_dp)

    ! Coulomb vector j^Q = sum_pq B^Q_pq D_pq
    if (hasA) call dgemm('T','N', naux, nvec, nocca*nocca, 1.0_dp, BooQ, nocca*nocca, &
         Ahh, nocca*nocca, 1.0_dp, jq, naux)
    if (hasZa .or. hasZb) then
       do v = 1, nvec
          do a = 1, nvirb
             do i = 1, nocca
                ZT(i,a,v) = 0.0_dp
                if (hasZa) ZT(i,a,v) = ZT(i,a,v) + Za(a,i,v)
                if (hasZb) ZT(i,a,v) = ZT(i,a,v) + Zb(a,i,v)
             enddo
          enddo
       enddo
       call dgemm('T','N', naux, nvec, nocca*nvirb, 2.0_dp, Bhp, nocca*nvirb, &
            ZT, nocca*nvirb, 1.0_dp, jq, naux)
    endif
    if (hasJ) jq = jq + jq_add

    ! Coulomb matrices on the HP / HH blocks
    allocate(JHP(nocca,nvirb,nvec), source=0.0_dp)
    allocate(JHH(nocca,nocca,nvec), source=0.0_dp)
    if (wantHP) call dgemm('N','N', nocca*nvirb, nvec, naux, 1.0_dp, Bhp, nocca*nvirb, &
         jq, naux, 0.0_dp, JHP, nocca*nvirb)
    if (wantHH) call dgemm('N','N', nocca*nocca, nvec, naux, 1.0_dp, BooQ, nocca*nocca, &
         jq, naux, 0.0_dp, JHH, nocca*nocca)

    ! exchange matrices
    allocate(KaHP(nocca,nvirb,nvec), KbHP(nocca,nvirb,nvec), source=0.0_dp)
    allocate(KaHH(nocca,nocca,nvec), KbHH(nocca,nocca,nvec), source=0.0_dp)
    allocate(W1(nocca,nocca), W2(nocca,nocca))

    if (hasA) then
       do v = 1, nvec
          do Q = 1, naux
             ! W1 = B^Q_HH Ahh
             call dgemm('N','N', nocca, nocca, nocca, 1.0_dp, BooQ(1,1,Q), nocca, &
                  Ahh(1,1,v), nocca, 0.0_dp, W1, nocca)
             if (wantHP) call dgemm('N','N', nocca, nvirb, nocca, 1.0_dp, W1, nocca, &
                  Bhp(1,1,Q), nocca, 1.0_dp, KaHP(1,1,v), nocca)
             if (wantHH) call dgemm('N','N', nocca, nocca, nocca, 1.0_dp, W1, nocca, &
                  BooQ(1,1,Q), nocca, 1.0_dp, KaHH(1,1,v), nocca)
          enddo
       enddo
    endif

    ! PH part of Z (D_{a i} = Z(a,i)):  sum_Q (B^Q_HP Z) B^Q_HP   /  (B^Q_HP Z) B^Q_HH (+ h.c.)
    if (hasZa .or. hasZb) then
       do v = 1, nvec
          do Q = 1, naux
             if (hasZa) then
                call dgemm('N','N', nocca, nocca, nvirb, 1.0_dp, Bhp(1,1,Q), nocca, &
                     Za(1,1,v), nvirb, 0.0_dp, W2, nocca)
                if (wantHP) call dgemm('N','N', nocca, nvirb, nocca, 1.0_dp, W2, nocca, &
                     Bhp(1,1,Q), nocca, 1.0_dp, KaHP(1,1,v), nocca)
                if (wantHH) then
                   call dgemm('N','N', nocca, nocca, nocca, 1.0_dp, W2, nocca, &
                        BooQ(1,1,Q), nocca, 1.0_dp, KaHH(1,1,v), nocca)
                   call dgemm('T','T', nocca, nocca, nocca, 1.0_dp, BooQ(1,1,Q), nocca, &
                        W2, nocca, 1.0_dp, KaHH(1,1,v), nocca)
                endif
             endif
             if (hasZb) then
                call dgemm('N','N', nocca, nocca, nvirb, 1.0_dp, Bhp(1,1,Q), nocca, &
                     Zb(1,1,v), nvirb, 0.0_dp, W2, nocca)
                if (wantHP) call dgemm('N','N', nocca, nvirb, nocca, 1.0_dp, W2, nocca, &
                     Bhp(1,1,Q), nocca, 1.0_dp, KbHP(1,1,v), nocca)
                if (wantHH) then
                   call dgemm('N','N', nocca, nocca, nocca, 1.0_dp, W2, nocca, &
                        BooQ(1,1,Q), nocca, 1.0_dp, KbHH(1,1,v), nocca)
                   call dgemm('T','T', nocca, nocca, nocca, 1.0_dp, BooQ(1,1,Q), nocca, &
                        W2, nocca, 1.0_dp, KbHH(1,1,v), nocca)
                endif
             endif
          enddo
       enddo
       ! HP part of Z (D_{i a} = Z(a,i)):  [sum_Q B^Q_PP Z B^Q_HH]^T on HP
       if (wantHP) then
          allocate(Sz(nvirb,nocca,nvec))
          if (hasZa) then
             Sz = 0.0_dp
             call kx_pp_hh(nvec, 1.0_dp, Za, Sz)
             do v = 1, nvec
                KaHP(:,:,v) = KaHP(:,:,v) + transpose(Sz(:,:,v))
             enddo
          endif
          if (hasZb) then
             Sz = 0.0_dp
             call kx_pp_hh(nvec, 1.0_dp, Zb, Sz)
             do v = 1, nvec
                KbHP(:,:,v) = KbHP(:,:,v) + transpose(Sz(:,:,v))
             enddo
          endif
          deallocate(Sz)
       endif
    endif

    if (wantHP) then
       Ghpa = JHP - cx*KaHP
       Ghpb = JHP - cx*KbHP
    endif
    if (wantHH) then
       Ghha = JHH - cx*KaHH
       Ghhb = JHH - cx*KbHH
    endif

    deallocate(jq, ZT, JHP, JHH, KaHP, KbHP, KaHH, KbHH, W1, W2)

  end subroutine gfock

!######################################################################
! jblocks: J_pq = sum_Q B^Q_pq jq(Q) on the HH, HP and PP blocks for a
! set of auxiliary vectors
!######################################################################
  subroutine jblocks(nvec, jq, JHH, JHP, JPP)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: jq(naux,nvec)
    real(dp), intent(out)   :: JHH(nocca,nocca,nvec), JHP(nocca,nvirb,nvec), &
                               JPP(nvirb,nvirb,nvec)
    integer(is) :: v, k, a, b, Q

    call dgemm('N','N', nocca*nocca, nvec, naux, 1.0_dp, BooQ, nocca*nocca, &
         jq, naux, 0.0_dp, JHH, nocca*nocca)
    call dgemm('N','N', nocca*nvirb, nvec, naux, 1.0_dp, Bhp, nocca*nvirb, &
         jq, naux, 0.0_dp, JHP, nocca*nvirb)
    if (vv_full) then
       call dgemm('N','N', nvirb*nvirb, nvec, naux, 1.0_dp, Bvv, nvirb*nvirb, &
            jq, naux, 0.0_dp, JPP, nvirb*nvirb)
    else
       ! paired planes: plane k holds B^{2k-1} (lower triangle, incl. the
       ! diagonal) and B^{2k} (strict upper triangle). Two dgemms over the
       ! planes give T1(a,b) = sum_k Bvv(a,b,k) jq(2k-1) and
       ! T2(a,b) = sum_k Bvv(a,b,k) jq(2k); for a > b the lower element of
       ! J is T1(a,b) + T2(b,a).
       block
         real(dp), allocatable :: jodd(:,:), jeven(:,:), T1(:,:,:), T2(:,:,:)
         allocate(jodd(nplane,nvec), jeven(nplane,nvec), source=0.0_dp)
         allocate(T1(nvirb,nvirb,nvec), T2(nvirb,nvirb,nvec))
         do k = 1, nplane
            jodd(k,:) = jq(2*k-1,:)
            if (2*k <= naux) jeven(k,:) = jq(2*k,:)
         enddo
         call dgemm('N','N', nvirb*nvirb, nvec, nplane, 1.0_dp, Bvv, nvirb*nvirb, &
              jodd, nplane, 0.0_dp, T1, nvirb*nvirb)
         call dgemm('N','N', nvirb*nvirb, nvec, nplane, 1.0_dp, Bvv, nvirb*nvirb, &
              jeven, nplane, 0.0_dp, T2, nvirb*nvirb)
         do v = 1, nvec
            !$omp parallel do private(a)
            do b = 1, nvirb
               do a = b+1, nvirb
                  JPP(a,b,v) = T1(a,b,v) + T2(b,a,v)
                  JPP(b,a,v) = JPP(a,b,v)
               enddo
            enddo
            !$omp end parallel do
            do a = 1, nvirb
               JPP(a,a,v) = 0.0_dp
               do Q = 1, naux
                  JPP(a,a,v) = JPP(a,a,v) + Dall(Q,Pmap(a))*jq(Q,v)
               enddo
            enddo
         enddo
         deallocate(jodd, jeven, T1, T2)
       end block
    endif

  end subroutine jblocks

!######################################################################
! grad_state: per-state exchange-channel quantities in one pass over
! the Q-blocks:
!   T^Q = B^Q_PP X~,  U^Q = X~^T T^Q,  Y^Q = X~ B^Q_HH,  S^Q = B^Q_HP X~
!   LaH = -2cx sum_Q B^Q_HH U^Q,   LaP = -2cx sum_Q (B^Q_HP)^T U^Q
!   LbH = -2cx sum_Q S^Q (Y^Q)^T,  LbP = -2cx sum_Q T^Q (Y^Q)^T
!   Ghh^Q = -cx/2 U^Q,  Yf^Q = -cx/4 Y^Q + dqf(Q)/4 X~
!   (dqf: Coulomb vector of the density whose J pair with X~X~^T enters the
!   families, the reference d^Q for MRSF; dq returns the reference d^Q)
!   gpp(Q,Q') = 2 sum_aj T^Q_aj Yf^Q'_aj
!   KbT_HP = sum_Q S^Q (T^Q)^T,  KbT_HH = sum_Q S^Q (S^Q)^T   (K[X~X~^T] blocks)
!   jT(Q) = sum_aj T^Q_aj X~_aj                                (Coulomb vector of X~X~^T)
!######################################################################
  subroutine grad_state(cx, dqf, Xt, LaH, LaP, LbH, LbP, Ghh, Fhp, Yf, Sq, gpp, &
       KbT_HP, KbT_HH, jT, dq)

    real(dp), intent(in)  :: cx, dqf(naux), Xt(nvirb,nocca)
    real(dp), intent(out) :: dq(naux)
    real(dp), intent(out) :: LaH(nocca,nocca), LaP(nvirb,nocca), LbH(nocca,nvirb), &
                             LbP(nvirb,nvirb)
    real(dp), intent(out) :: Ghh(nocca,nocca,naux), Fhp(nocca,nvirb,naux), &
                             Yf(nvirb,nocca,naux), Sq(nocca,nocca,naux), gpp(naux,naux)
    real(dp), intent(out) :: KbT_HP(nocca,nvirb), KbT_HH(nocca,nocca), jT(naux)

    integer(is)           :: blk, Q0, Ql, Q, kk, k, j, ldT, nQb, p
    real(dp), allocatable :: dd(:), Yw(:,:,:), U(:,:), Tc(:,:,:)
    real(dp), external    :: ddot

    ! Coulomb vector of the reference density: d^Q = sum_p B^Q_pp occ_p
    do Q = 1, naux
       dq(Q) = 0.0_dp
       do p = 1, nmo
          dq(Q) = dq(Q) + Dall(Q,p)*occ(p)
       enddo
    enddo

    call ensure_work(1_is)
    ldT = nvirb*nQ
    allocate(dd(nvirb), U(nocca,nocca), Tc(nvirb,nocca,nQ))
    allocate(Yw(nvirb,nocca,naux))

    ! Y^Q = X~ B^Q_HH for all Q in one dgemm
    call dgemm('N','N', nvirb, nocca*naux, nocca, 1.0_dp, Xt, nvirb, BooQ, nocca, &
         0.0_dp, Yw, nvirb)
    do Q = 1, naux
       Yf(:,:,Q) = -0.25_dp*cx*Yw(:,:,Q) + 0.25_dp*dqf(Q)*Xt
    enddo
    ! S^Q = B^Q_HP X~
    do Q = 1, naux
       call dgemm('N','N', nocca, nocca, nvirb, 1.0_dp, Bhp(1,1,Q), nocca, Xt, nvirb, &
            0.0_dp, Sq(1,1,Q), nocca)
    enddo
    LbH = 0.0_dp; KbT_HH = 0.0_dp
    do Q = 1, naux
       call dgemm('N','T', nocca, nvirb, nocca, -2.0_dp*cx, Sq(1,1,Q), nocca, Yw(1,1,Q), nvirb, &
            1.0_dp, LbH, nocca)
       call dgemm('N','T', nocca, nocca, nocca, 1.0_dp, Sq(1,1,Q), nocca, Sq(1,1,Q), nocca, &
            1.0_dp, KbT_HH, nocca)
    enddo

    LaH = 0.0_dp; LaP = 0.0_dp; LbP = 0.0_dp; KbT_HP = 0.0_dp; gpp = 0.0_dp
    Fhp = 0.0_dp

    do blk = 1, nblk
       Q0  = (blk-1)*nQ
       nQb = min(nQ, naux - Q0)
       ! T^Q = B^Q_PP X~ into Twork(:,Ql,:,1)
       if (vv_full) then
          do Ql = 1, nQb
             call dgemm('N','N', nvirb, nocca, nvirb, 1.0_dp, Bvv(1,1,Q0+Ql), nvirb, &
                  Xt, nvirb, 0.0_dp, Twork(1,Ql,1,1), ldT)
          enddo
       else
          do kk = 1, npblk
             k = (blk-1)*npblk + kk
             if (k > nplane) exit
             Q  = 2*k - 1
             Ql = 2*kk - 1
             call plane_apply(nocca, Bvv(1,1,k), Xt, Q, Ql, ldT, dd)
          enddo
       endif
       do Ql = 1, nQb
          Q = Q0 + Ql
          ! U^Q = X~^T T^Q
          call dgemm('T','N', nocca, nocca, nvirb, 1.0_dp, Xt, nvirb, Twork(1,Ql,1,1), ldT, &
               0.0_dp, U, nocca)
          Ghh(:,:,Q) = -0.5_dp*cx*U
          call dgemm('N','N', nocca, nocca, nocca, -2.0_dp*cx, BooQ(1,1,Q), nocca, U, nocca, &
               1.0_dp, LaH, nocca)
          call dgemm('T','N', nvirb, nocca, nocca, -2.0_dp*cx, Bhp(1,1,Q), nocca, U, nocca, &
               1.0_dp, LaP, nvirb)
          call dgemm('N','T', nvirb, nvirb, nocca, -2.0_dp*cx, Twork(1,Ql,1,1), ldT, Yw(1,1,Q), nvirb, &
               1.0_dp, LbP, nvirb)
          call dgemm('N','T', nocca, nvirb, nocca, 1.0_dp, Sq(1,1,Q), nocca, Twork(1,Ql,1,1), ldT, &
               1.0_dp, KbT_HP, nocca)
          jT(Q) = 0.0_dp
          do j = 1, nocca
             jT(Q) = jT(Q) + ddot(nvirb, Twork(1,Ql,j,1), 1_is, Xt(1,j), 1_is)
             Tc(:,j,Ql) = Twork(:,Ql,j,1)
          enddo
       enddo
       ! gpp(Q0+1:Q0+nQb, :) = 2 sum_aj T^Q_aj Yf^Q'_aj
       call dgemm('T','N', nQb, naux, nvirb*nocca, 2.0_dp, Tc, nvirb*nocca, Yf, nvirb*nocca, &
            0.0_dp, gpp(Q0+1,1), naux)
    enddo

    deallocate(dd, U, Tc, Yw)

  end subroutine grad_state

!######################################################################
! grad_finish: z-dependent and mean-field contributions to the families
! and the metric matrix g.
!   pq(Q) = sum_HH B^Q.Ta + jT(Q) + 2 sum_HP B^Q.(Za+Zb)       (Coulomb vector of P_t)
!   J-type:  Ghh += dq/2 Ta + pq/2 diag(occ_H),  Fhp += dq (Za+Zb)^T
!   K-type (alpha, D^a = 1_H):
!        Ghh += -cx/2 [Ta B^Q_HH + B^Q_HH Ta + Za^T B^Q_HP^T + B^Q_HP Za]
!        Fhp += -cx   [B^Q_HH Za^T]
!   K-type (beta, D^b = 1_C; rows/cols restricted to C):
!        Ghh|_CC += -cx/2 [Zb^T B^Q_HP^T + B^Q_HP Zb]
!        Fhp|_C  += -cx   [S^Q X~^T + B^Q_HH Zb^T]
!   g(Q,Q') = sum_HH B^Q Ghh^Q' + sum_HP B^Q Fhp^Q' + gpp(Q,Q')
!######################################################################
  subroutine grad_finish(cx, dq, jT, Ta, Xt, Za, Zb, Sq, occH, Ghh, Fhp, gpp, g, pq)

    real(dp), intent(in)    :: cx, dq(naux), jT(naux), Ta(nocca,nocca), Xt(nvirb,nocca)
    real(dp), intent(in)    :: Za(nvirb,nocca), Zb(nvirb,nocca), Sq(nocca,nocca,naux)
    real(dp), intent(in)    :: occH(nocca)
    real(dp), intent(inout) :: Ghh(nocca,nocca,naux), Fhp(nocca,nvirb,naux)
    real(dp), intent(in)    :: gpp(naux,naux)
    real(dp), intent(out)   :: g(naux,naux), pq(naux)

    integer(is)           :: Q, i, j, a
    real(dp), allocatable :: ZT(:,:), ZaT(:,:), ZbT(:,:), W(:,:), Wb(:,:), WF(:,:)
    real(dp), external    :: ddot

    allocate(ZT(nocca,nvirb), ZaT(nocca,nvirb), ZbT(nocca,nvirb))
    allocate(W(nocca,nocca), Wb(nocca,nocca), WF(nocca,nvirb))
    ZaT = transpose(Za); ZbT = transpose(Zb); ZT = ZaT + ZbT

    ! Coulomb vector of the relaxed density
    call dgemv('T', nocca*nocca, naux, 1.0_dp, BooQ, nocca*nocca, Ta, 1_is, 0.0_dp, pq, 1_is)
    call dgemv('T', nocca*nvirb, naux, 2.0_dp, Bhp, nocca*nvirb, ZT, 1_is, 1.0_dp, pq, 1_is)
    pq = pq + jT

    do Q = 1, naux
       ! J-type pieces
       Ghh(:,:,Q) = Ghh(:,:,Q) + 0.5_dp*dq(Q)*Ta
       do i = 1, nocca
          Ghh(i,i,Q) = Ghh(i,i,Q) + 0.5_dp*pq(Q)*occH(i)
       enddo
       Fhp(:,:,Q) = Fhp(:,:,Q) + dq(Q)*ZT
       ! alpha mean-field exchange: Ta B_HH + B_HH Ta
       call dgemm('N','N', nocca, nocca, nocca, 1.0_dp, Ta, nocca, BooQ(1,1,Q), nocca, 0.0_dp, W, nocca)
       Ghh(:,:,Q) = Ghh(:,:,Q) - 0.5_dp*cx*(W + transpose(W))
       ! B_HP Za (+ transpose)
       call dgemm('N','N', nocca, nocca, nvirb, 1.0_dp, Bhp(1,1,Q), nocca, Za, nvirb, 0.0_dp, W, nocca)
       Ghh(:,:,Q) = Ghh(:,:,Q) - 0.5_dp*cx*(W + transpose(W))
       ! Fhp += -cx B_HH Za^T
       call dgemm('N','N', nocca, nvirb, nocca, -cx, BooQ(1,1,Q), nocca, ZaT, nocca, 1.0_dp, Fhp(1,1,Q), nocca)
       ! beta mean-field exchange, restricted to the core rows/columns
       call dgemm('N','N', nocca, nocca, nvirb, 1.0_dp, Bhp(1,1,Q), nocca, Zb, nvirb, 0.0_dp, Wb, nocca)
       Wb = Wb + transpose(Wb)
       do j = 1, nC
          do i = 1, nC
             Ghh(i,j,Q) = Ghh(i,j,Q) - 0.5_dp*cx*Wb(i,j)
          enddo
       enddo
       call dgemm('N','T', nocca, nvirb, nocca, 1.0_dp, Sq(1,1,Q), nocca, Xt, nvirb, 0.0_dp, WF, nocca)
       call dgemm('N','N', nocca, nvirb, nocca, 1.0_dp, BooQ(1,1,Q), nocca, ZbT, nocca, 1.0_dp, WF, nocca)
       do a = 1, nvirb
          do i = 1, nC
             Fhp(i,a,Q) = Fhp(i,a,Q) - cx*WF(i,a)
          enddo
       enddo
    enddo

    ! metric matrix
    call dgemm('T','N', naux, naux, nocca*nocca, 1.0_dp, BooQ, nocca*nocca, Ghh, nocca*nocca, &
         0.0_dp, g, naux)
    call dgemm('T','N', naux, naux, nocca*nvirb, 1.0_dp, Bhp, nocca*nvirb, Fhp, nocca*nvirb, &
         1.0_dp, g, naux)
    g = g + gpp

    deallocate(ZT, ZaT, ZbT, W, Wb, WF)

  end subroutine grad_finish

!######################################################################
! grad_reffam: B-space families of the ROKS reference two-electron
! energy  E = 1/2 (D_t J[D_t]) - cx/2 sum_s (D^s K[D^s]):
!   Ghh^Q = dq(Q)/2 diag(occ_H) - cx/2 [B^Q_HH + B^Q_CC]   (alpha: 1_H, beta: 1_C)
!   g(Q,Q') = sum_HH B^Q Ghh^Q'
!######################################################################
  subroutine grad_reffam(cx, dq, occH, Ghh, g)

    real(dp), intent(in)  :: cx, dq(naux), occH(nocca)
    real(dp), intent(out) :: Ghh(nocca,nocca,naux), g(naux,naux)
    integer(is) :: Q, i, j

    do Q = 1, naux
       Ghh(:,:,Q) = -0.5_dp*cx*BooQ(:,:,Q)
       do j = 1, nC
          do i = 1, nC
             Ghh(i,j,Q) = Ghh(i,j,Q) - 0.5_dp*cx*BooQ(i,j,Q)
          enddo
       enddo
       do i = 1, nocca
          Ghh(i,i,Q) = Ghh(i,i,Q) + 0.5_dp*dq(Q)*occH(i)
       enddo
    enddo
    call dgemm('T','N', naux, naux, nocca*nocca, 1.0_dp, BooQ, nocca*nocca, Ghh, nocca*nocca, &
         0.0_dp, g, naux)

  end subroutine grad_reffam

!######################################################################
! grad_dq: Coulomb vector of a diagonal MO density, dq(Q) = sum_p B^Q_pp occ_p
!######################################################################
  subroutine grad_dq(occ1, dq)

    real(dp), intent(in)  :: occ1(nmo)
    real(dp), intent(out) :: dq(naux)

    call dgemv('N', naux, nmo, 1.0_dp, Dall, naux, occ1, 1_is, 0.0_dp, dq, 1_is)

  end subroutine grad_dq

!######################################################################
! grad_booq: copy of the hole-hole block B^Q_hh' (nocca, nocca, naux)
!######################################################################
  subroutine grad_booq(B)

    real(dp), intent(out) :: B(nocca,nocca,naux)

    if (.not. grad_ready) call mrsf_error('grad_booq: gradient blocks not initialised')
    B = BooQ

  end subroutine grad_booq

!######################################################################
! grad_jq: Coulomb vectors jq(Q,v) = sum_hh' B^Q_hh' Ahh(h,h',v)
!                                   + sum_hp B^Q_hp Ahp(p,h,v)
! (hole-hole and particle-hole MO matrices in the local orders; no
! symmetrisation: a symmetric density with HP + PH blocks needs 2 Ahp)
!######################################################################
  subroutine grad_jq(nvec, Ahh, Ahp, jq)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: Ahh(nocca,nocca,nvec), Ahp(nvirb,nocca,nvec)
    real(dp), intent(out)   :: jq(naux,nvec)

    real(dp), allocatable   :: AhpT(:,:,:)
    integer(is)             :: v

    if (.not. grad_ready) call mrsf_error('grad_jq: gradient blocks not initialised')
    allocate(AhpT(nocca,nvirb,nvec))
    do v = 1, nvec
       AhpT(:,:,v) = transpose(Ahp(:,:,v))
    enddo
    call dgemm('T','N', naux, nvec, nocca*nocca, 1.0_dp, BooQ, nocca*nocca, Ahh, nocca*nocca, &
         0.0_dp, jq, naux)
    call dgemm('T','N', naux, nvec, nocca*nvirb, 1.0_dp, Bhp, nocca*nvirb, AhpT, nocca*nvirb, &
         1.0_dp, jq, naux)
    deallocate(AhpT)

  end subroutine grad_jq

!######################################################################
! grad_bvec: out(Q, t, v) = sum_{q in cls} B^Q_{t q} U(q, v) for all MOs t
!   cls = 1: q runs over the doubly occupied MOs (local hole order 1..nC)
!   cls = 2: q runs over the virtuals (local particle order 3..nvirb)
! one pass over BooQ / Bhp, and for cls = 2 one sweep over the vir-vir
! planes (dsymm with nv columns per plane)
!######################################################################
  subroutine grad_bvec(nvu, cls, U, out)

    integer(is), intent(in) :: nvu, cls
    real(dp), intent(in)    :: U(*)
    real(dp), intent(out)   :: out(naux,nmo,nvu)

    integer(is)           :: Q, h, a, v, k, ncl
    real(dp), allocatable :: Uemb(:,:), W(:,:), T1(:,:), T2(:,:)

    if (.not. grad_ready) call mrsf_error('grad_bvec: gradient blocks not initialised')
    out = 0.0_dp
    if (cls == 1) then
       ncl = nC
       if (ncl == 0) return
       allocate(W(nocca,nvu), T1(nvirb,nvu))
       do Q = 1, naux
          ! t in H: sum_q BooQ(h, q, Q) U(q, v)
          call dgemm('N','N', nocca, nvu, ncl, 1.0_dp, BooQ(1,1,Q), nocca, U, ncl, 0.0_dp, W, nocca)
          ! t in V: sum_q Bhp(q, a, Q) U(q, v)
          call dgemm('T','N', nvirb, nvu, ncl, 1.0_dp, Bhp(1,1,Q), nocca, U, ncl, 0.0_dp, T1, nvirb)
          do v = 1, nvu
             do h = 1, nocca
                out(Q,Hmap(h),v) = W(h,v)
             enddo
             do a = 3, nvirb
                out(Q,Pmap(a),v) = T1(a,v)
             enddo
          enddo
       enddo
       deallocate(W, T1)
    else if (cls == 2) then
       ncl = nV
       if (ncl == 0) return
       allocate(Uemb(nvirb,nvu), source=0.0_dp)
       call copy_rows(ncl, nvu, U, Uemb)
       allocate(W(nocca,nvu), T1(nvirb,nvu), T2(nvirb,nvu))
       ! t in H: sum_a Bhp(h, a, Q) Uemb(a, v)
       do Q = 1, naux
          call dgemm('N','N', nocca, nvu, nvirb, 1.0_dp, Bhp(1,1,Q), nocca, Uemb, nvirb, 0.0_dp, W, nocca)
          do v = 1, nvu
             do h = 1, nocca
                out(Q,Hmap(h),v) = W(h,v)
             enddo
          enddo
       enddo
       ! t in V: sum_a' B^Q_{a a'} Uemb(a', v) from the planes
       if (vv_full) then
          do Q = 1, naux
             call dgemm('N','N', nvirb, nvu, nvirb, 1.0_dp, Bvv(1,1,Q), nvirb, Uemb, nvirb, &
                  0.0_dp, T1, nvirb)
             do v = 1, nvu
                do a = 3, nvirb
                   out(Q,Pmap(a),v) = T1(a,v)
                enddo
             enddo
          enddo
       else
          do k = 1, nplane
             Q = 2*k - 1
             call dsymm('L','L', nvirb, nvu, 1.0_dp, Bvv(1,1,k), nvirb, Uemb, nvirb, 0.0_dp, T1, nvirb)
             do v = 1, nvu
                do a = 3, nvirb
                   out(Q,Pmap(a),v) = T1(a,v)
                enddo
             enddo
             if (Q + 1 <= naux) then
                call dsymm('L','U', nvirb, nvu, 1.0_dp, Bvv(1,1,k), nvirb, Uemb, nvirb, 0.0_dp, T2, nvirb)
                do v = 1, nvu
                   do a = 3, nvirb
                      out(Q+1,Pmap(a),v) = T2(a,v) + (Dall(Q+1,Pmap(a)) - Dall(Q,Pmap(a)))*Uemb(a,v)
                   enddo
                enddo
             endif
          enddo
       endif
       deallocate(Uemb, W, T1, T2)
    else
       call mrsf_error('grad_bvec: cls must be 1 (core) or 2 (virtual)')
    endif

  end subroutine grad_bvec

!######################################################################
! grad_bdot: out(t, w) = sum_Q sum_{s in cls} B^Q_{t s} W(s, Q, w) for all
! MOs t (cls as in grad_bvec); for cls = 2 one sweep over the planes
!######################################################################
  subroutine grad_bdot(nw, cls, W, out)

    integer(is), intent(in) :: nw, cls
    real(dp), intent(in)    :: W(*)
    real(dp), intent(out)   :: out(nmo,nw)

    integer(is)           :: Q, h, a, v, k, ncl
    real(dp), allocatable :: Wq(:,:), OH(:,:), OP(:,:), T1(:,:), T2(:,:)

    if (.not. grad_ready) call mrsf_error('grad_bdot: gradient blocks not initialised')
    out = 0.0_dp
    ncl = merge(nC, nV, cls == 1)
    if (cls /= 1 .and. cls /= 2) call mrsf_error('grad_bdot: cls must be 1 (core) or 2 (virtual)')
    if (ncl == 0) return
    allocate(OH(nocca,nw), OP(nvirb,nw), source=0.0_dp)
    if (cls == 1) then
       allocate(Wq(ncl,nw))
       do Q = 1, naux
          call gather_q(ncl, nw, Q, W, Wq)
          call dgemm('N','N', nocca, nw, ncl, 1.0_dp, BooQ(1,1,Q), nocca, Wq, ncl, 1.0_dp, OH, nocca)
          call dgemm('T','N', nvirb, nw, ncl, 1.0_dp, Bhp(1,1,Q), nocca, Wq, ncl, 1.0_dp, OP, nvirb)
       enddo
       deallocate(Wq)
    else
       allocate(Wq(nvirb,nw), source=0.0_dp)
       allocate(T1(nvirb,nw), T2(nvirb,nw))
       do Q = 1, naux
          call gather_q_emb(ncl, nw, Q, W, Wq)
          call dgemm('N','N', nocca, nw, nvirb, 1.0_dp, Bhp(1,1,Q), nocca, Wq, nvirb, 1.0_dp, OH, nocca)
       enddo
       if (vv_full) then
          do Q = 1, naux
             call gather_q_emb(ncl, nw, Q, W, Wq)
             call dgemm('N','N', nvirb, nw, nvirb, 1.0_dp, Bvv(1,1,Q), nvirb, Wq, nvirb, 1.0_dp, OP, nvirb)
          enddo
       else
          do k = 1, nplane
             Q = 2*k - 1
             call gather_q_emb(ncl, nw, Q, W, Wq)
             call dsymm('L','L', nvirb, nw, 1.0_dp, Bvv(1,1,k), nvirb, Wq, nvirb, 1.0_dp, OP, nvirb)
             if (Q + 1 <= naux) then
                call gather_q_emb(ncl, nw, Q+1, W, Wq)
                call dsymm('L','U', nvirb, nw, 1.0_dp, Bvv(1,1,k), nvirb, Wq, nvirb, 1.0_dp, OP, nvirb)
                do v = 1, nw
                   do a = 3, nvirb
                      OP(a,v) = OP(a,v) + (Dall(Q+1,Pmap(a)) - Dall(Q,Pmap(a)))*Wq(a,v)
                   enddo
                enddo
             endif
          enddo
       endif
       deallocate(Wq, T1, T2)
    endif
    do v = 1, nw
       do h = 1, nocca
          out(Hmap(h),v) = OH(h,v)
       enddo
       do a = 3, nvirb
          out(Pmap(a),v) = OP(a,v)
       enddo
    enddo
    deallocate(OH, OP)

  end subroutine grad_bdot

  ! U(ncl, nv) -> rows 3..nvirb of Uemb(nvirb, nv)
  subroutine copy_rows(ncl, nv, U, Uemb)
    integer(is), intent(in) :: ncl, nv
    real(dp), intent(in)    :: U(ncl,nv)
    real(dp), intent(inout) :: Uemb(nvirb,nv)
    Uemb(3:2+ncl,:) = U
  end subroutine copy_rows

  ! Wq(:, w) = W(:, Q, w)
  subroutine gather_q(ncl, nw, Q, W, Wq)
    integer(is), intent(in) :: ncl, nw, Q
    real(dp), intent(in)    :: W(ncl,naux,nw)
    real(dp), intent(out)   :: Wq(ncl,nw)
    Wq = W(:,Q,:)
  end subroutine gather_q

  ! Wq(3:, w) = W(:, Q, w) (rows 1:2 stay zero)
  subroutine gather_q_emb(ncl, nw, Q, W, Wq)
    integer(is), intent(in) :: ncl, nw, Q
    real(dp), intent(in)    :: W(ncl,naux,nw)
    real(dp), intent(inout) :: Wq(nvirb,nw)
    Wq(3:2+ncl,:) = W(:,Q,:)
  end subroutine gather_q_emb

end module mrsf_gradient
