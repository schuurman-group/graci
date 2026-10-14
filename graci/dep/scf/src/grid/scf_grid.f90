!**********************************************************************
! scf_grid: the quadrature grid (PySCF's points and weights), blocks of
! nbk = 112 points (two of PySCF's 56-point screening blocks), the
! significant shells and AOs of every block, the compacted AO cache
! (values and gradients of the significant AOs only, layout
! (np, nsig, ncomp) per block) filled once per geometry when it fits the
! memory budget (else the values are evaluated block by block in every
! pass), and the XC passes of the SCF:
!   ROKS exact:    E_xc, V_a, V_b
!   ROKS shortcut: E_xc, V_c = (V_a + V_b)/2 and (V_a - V_b) C_O, where
!                  C_O = the last ka - kb columns of the factor X
!   RKS:           E_xc, V
! One pass reads every AO block once: densities from the factor X
! (D_a = X X^T, D_b = X(:,1:kb) X(:,1:kb)^T; RKS: D = 2 X X^T), libxc,
! weighted AO products and V += phi^T aow (V = T + T^T at the end).
! OpenMP over blocks with sequential MKL and per-thread accumulators.
!**********************************************************************
module scf_grid

  use scf_constants
  use scf_io
  use scf_gto
  use scf_xc
  !$ use omp_lib

  implicit none

  integer(is), parameter   :: nbk = 112
  integer(is)              :: ngrid = 0, nblk = 0, ncomp_g = 4, nsigmax = 0
  real(dp), allocatable    :: gx(:), gy(:), gz(:), gw(:)
  integer(is), allocatable :: bbeg(:), bend(:)
  integer(is), allocatable :: shp(:), shl(:)       ! CSR: significant shells per block
  integer(is), allocatable :: aop(:), aol(:)       ! CSR: significant AOs per block (1-based)
  integer(is), allocatable :: rnp(:), rpos(:), rlen(:) ! CSR: runs of consecutive AOs per block
  logical                  :: cached = .false.
  real(dp), allocatable    :: cache(:)
  integer(ib), allocatable :: coff(:)
  character(len=255)       :: gstamp = ''
  real(dp)                 :: time_aofill = 0.0_dp, time_xcpass = 0.0_dp
  real(dp)                 :: flop_xc = 0.0_dp   ! flops of one phi^T aow product over the grid
  real(dp)                 :: frac_sig = 0.0_dp, cost_sig = 0.0_dp
  ! thread-summed sub-timers of the XC passes: AO evaluation, MO values,
  ! density + libxc, weighted AOs, V products, scatter of V
  real(dp)                 :: tsub(6) = 0.0_dp

  external :: dgemm

contains

!######################################################################
! grid_set: grid (xyz(3,ng): PySCF's C-ordered coordinates), blocks,
! screening, AO lists and (within budget bytes) the AO cache; ncomp = 4
! (GGA) or 1 (LDA); ierr 0 ok
!######################################################################
  subroutine grid_set(ng, xyz, w, ncomp, budget, stamp, ierr)

    integer(is), intent(in)      :: ng, ncomp
    real(dp), intent(in)         :: xyz(3,ng), w(ng), budget
    character(len=*), intent(in) :: stamp
    integer(is), intent(out)     :: ierr

    logical, allocatable     :: keep(:,:)
    integer(is), allocatable :: nsh(:), nsg(:)
    integer(is) :: iblk, g0, g1, np, ish, k, n, m, deriv
    integer(ib) :: tot
    real(dp)    :: t0

    ierr = 0
    call grid_free()
    t0 = wall_time()
    ngrid = ng
    ncomp_g = ncomp
    allocate(gx(ng), gy(ng), gz(ng), gw(ng))
    gx = xyz(1,:); gy = xyz(2,:); gz = xyz(3,:); gw = w
    nblk = (ng + nbk - 1)/nbk
    allocate(bbeg(nblk), bend(nblk), keep(nbas_b,nblk), nsh(nblk), nsg(nblk))
    do iblk = 1, nblk
       bbeg(iblk) = (iblk-1)*nbk + 1
       bend(iblk) = min(iblk*nbk, ng)
    enddo

    ! significant shells per block
    !$omp parallel do private(iblk, g0, g1, np, ish, n, m) schedule(dynamic)
    do iblk = 1, nblk
       g0 = bbeg(iblk); g1 = bend(iblk); np = g1 - g0 + 1
       call gto_screen(np, gx(g0), gy(g0), gz(g0), keep(1,iblk))
       n = 0; m = 0
       do ish = 1, nbas_b
          if (keep(ish,iblk)) then
             n = n + 1
             m = m + aoloc_b(ish+1) - aoloc_b(ish)
          endif
       enddo
       nsh(iblk) = n
       nsg(iblk) = m
    enddo
    !$omp end parallel do

    ! CSR lists
    allocate(shp(nblk+1), aop(nblk+1))
    shp(1) = 1; aop(1) = 1
    do iblk = 1, nblk
       shp(iblk+1) = shp(iblk) + nsh(iblk)
       aop(iblk+1) = aop(iblk) + nsg(iblk)
    enddo
    allocate(shl(max(shp(nblk+1)-1,1)), aol(max(aop(nblk+1)-1,1)))
    do iblk = 1, nblk
       n = shp(iblk) - 1
       m = aop(iblk) - 1
       do ish = 1, nbas_b
          if (keep(ish,iblk)) then
             n = n + 1
             shl(n) = ish
             do k = aoloc_b(ish)+1, aoloc_b(ish+1)
                m = m + 1
                aol(m) = k
             enddo
          endif
       enddo
    enddo
    ! runs of consecutive AO indices in every block's list (contiguous
    ! scatter of the block matrices)
    allocate(rnp(nblk+1))
    rnp(1) = 1
    do iblk = 1, nblk
       n = 0
       do k = aop(iblk), aop(iblk+1)-1
          if (k == aop(iblk)) then
             n = n + 1
          else if (aol(k) /= aol(k-1) + 1) then
             n = n + 1
          endif
       enddo
       rnp(iblk+1) = rnp(iblk) + n
    enddo
    allocate(rpos(max(rnp(nblk+1)-1,1)), rlen(max(rnp(nblk+1)-1,1)))
    do iblk = 1, nblk
       n = rnp(iblk) - 1
       do k = aop(iblk), aop(iblk+1)-1
          if (k == aop(iblk)) then
             n = n + 1
             rpos(n) = k - aop(iblk) + 1
             rlen(n) = 1
          else if (aol(k) /= aol(k-1) + 1) then
             n = n + 1
             rpos(n) = k - aop(iblk) + 1
             rlen(n) = 1
          else
             rlen(n) = rlen(n) + 1
          endif
       enddo
    enddo
    nsigmax = maxval(nsg)
    flop_xc = 0.0_dp
    do iblk = 1, nblk
       flop_xc = flop_xc + 2.0_dp*real(bend(iblk)-bbeg(iblk)+1,dp)*real(nsg(iblk),dp)**2
    enddo
    frac_sig = real(sum(int(nsg,ib)),dp)/(real(nblk,dp)*real(nao_b,dp))
    cost_sig = sum((real(nsg,dp)/real(nao_b,dp))**2)/real(nblk,dp)

    ! AO cache
    tot = 0
    allocate(coff(nblk))
    do iblk = 1, nblk
       coff(iblk) = tot
       tot = tot + int(bend(iblk)-bbeg(iblk)+1,ib)*int(nsg(iblk),ib)*int(ncomp,ib)
    enddo
    cached = (8.0_dp*real(tot,dp) <= budget)
    if (cached) then
       allocate(cache(tot))
       deriv = 0
       if (ncomp == 4) deriv = 1
       !$omp parallel do private(iblk, g0, g1, np) schedule(dynamic)
       do iblk = 1, nblk
          g0 = bbeg(iblk); g1 = bend(iblk); np = g1 - g0 + 1
          call gto_eval(np, gx(g0), gy(g0), gz(g0), nsh(iblk), shl(shp(iblk)), deriv, np, nsg(iblk), &
               cache(coff(iblk)+1))
       enddo
       !$omp end parallel do
    endif
    gstamp = stamp
    time_aofill = wall_time() - t0
    deallocate(keep, nsh, nsg)

  end subroutine grid_set

  subroutine grid_free()
    if (allocated(gx)) deallocate(gx, gy, gz, gw)
    if (allocated(bbeg)) deallocate(bbeg, bend)
    if (allocated(shp)) deallocate(shp, shl, aop, aol)
    if (allocated(rnp)) deallocate(rnp, rpos, rlen)
    if (allocated(cache)) deallocate(cache)
    if (allocated(coff)) deallocate(coff)
    cached = .false.
    ngrid = 0
    nblk = 0
    gstamp = ''
  end subroutine grid_free

!######################################################################
! grid_pass_roks: mode 0 (exact): V1 = V_a, V2 = V_b; mode 1
! (shortcut): V1 = V_c, dCO = (V_a - V_b) X(:,kb+1:ka); exc = E_xc
!######################################################################
  subroutine grid_pass_roks(nao, ka, kb, X, mode, exc, V1, V2, dCO)

    integer(is), intent(in) :: nao, ka, kb, mode
    real(dp), intent(in)    :: X(nao,ka)
    real(dp), intent(out)   :: exc, V1(nao,nao), V2(nao,nao), dCO(nao,*)

    real(dp), allocatable :: Vt1(:,:,:), Vt2(:,:,:), dt(:,:,:), ext(:), tt(:,:)
    real(dp), allocatable :: aobuf(:), Xs(:), psi(:), aow(:), Vl(:), t1(:), t2(:)
    integer(is) :: nt, tid, iblk, g0, np, nsg, no, j, t, nr
    real(dp)    :: t0, ta
    logical     :: same

    t0 = wall_time()
    no = ka - kb
    ! equal spin densities (the guess): V_b = V_a
    same = (mode == 0 .and. ka == kb)
    nt = 1
    !$ nt = omp_get_max_threads()
    allocate(Vt1(nao,nao,nt), ext(nt))
    if (mode == 0) then
       allocate(Vt2(nao,nao,nt))
    else
       allocate(Vt2(1,1,1))
    endif
    allocate(dt(nao,max(no,1),nt), tt(6,nt))
    Vt1 = 0.0_dp; Vt2 = 0.0_dp; dt = 0.0_dp; ext = 0.0_dp; tt = 0.0_dp

    !$omp parallel private(tid, iblk, g0, np, nsg, nr, aobuf, Xs, psi, aow, Vl, t1, t2, ta)
    tid = 1
    !$ tid = omp_get_thread_num() + 1
    call blas_sequential()
    if (.not. cached) allocate(aobuf(nbk*nsigmax*ncomp_g))
    allocate(Xs(nsigmax*ka), psi(nbk*ka*ncomp_g), aow(nbk*nsigmax), Vl(nsigmax*nsigmax), &
         t1(nbk*max(no,1)), t2(nsigmax*max(no,1)))
    !$omp do schedule(dynamic)
    do iblk = 1, nblk
       g0  = bbeg(iblk)
       np  = bend(iblk) - g0 + 1
       nsg = aop(iblk+1) - aop(iblk)
       if (nsg == 0) cycle
       nr  = rnp(iblk+1) - rnp(iblk)
       if (cached) then
          call kern(np, nsg, aol(aop(iblk)), nr, rpos(rnp(iblk)), rlen(rnp(iblk)), cache(coff(iblk)+1), &
               gw(g0), Xs, psi, aow, Vl, t1, t2, Vt1(1,1,tid), Vt2(1,1,min(tid,size(Vt2,3))), dt(1,1,tid), &
               ext(tid), tt(1,tid))
       else
          ta = wall_time()
          call gto_eval(np, gx(g0), gy(g0), gz(g0), shp(iblk+1)-shp(iblk), shl(shp(iblk)), &
               merge(1_is, 0_is, ncomp_g == 4), np, nsg, aobuf)
          tt(1,tid) = tt(1,tid) + wall_time() - ta
          call kern(np, nsg, aol(aop(iblk)), nr, rpos(rnp(iblk)), rlen(rnp(iblk)), aobuf, gw(g0), Xs, psi, &
               aow, Vl, t1, t2, Vt1(1,1,tid), Vt2(1,1,min(tid,size(Vt2,3))), dt(1,1,tid), ext(tid), tt(1,tid))
       endif
    enddo
    !$omp end do
    if (allocated(aobuf)) deallocate(aobuf)
    deallocate(Xs, psi, aow, Vl, t1, t2)
    call blas_restore()
    !$omp end parallel
    tsub = tsub + sum(tt, dim=2)

    ! reductions and V = T + T^T
    !$omp parallel do private(j, t)
    do j = 1, nao
       V1(:,j) = Vt1(:,j,1)
       do t = 2, nt
          V1(:,j) = V1(:,j) + Vt1(:,j,t)
       enddo
       if (mode == 0 .and. .not. same) then
          V2(:,j) = Vt2(:,j,1)
          do t = 2, nt
             V2(:,j) = V2(:,j) + Vt2(:,j,t)
          enddo
       endif
    enddo
    !$omp end parallel do
    call symmetrise_sum(nao, V1)
    if (mode == 0) then
       if (same) then
          V2 = V1
       else
          call symmetrise_sum(nao, V2)
       endif
    endif
    if (mode == 1 .and. no > 0) then
       dCO(1:nao,1:no) = dt(:,1:no,1)
       do t = 2, nt
          dCO(1:nao,1:no) = dCO(1:nao,1:no) + dt(:,1:no,t)
       enddo
    endif
    exc = sum(ext)
    deallocate(Vt1, Vt2, dt, ext, tt)
    time_xcpass = time_xcpass + wall_time() - t0

  contains

    subroutine kern(np1, nsg1, idx, nr1, rp1, rl1, ao, w, Xs1, psi1, aow1, Vl1, t11, t21, Va, Vb, dd, ex, tk)
      integer(is), intent(in) :: np1, nsg1, idx(nsg1), nr1, rp1(nr1), rl1(nr1)
      real(dp), intent(in)    :: ao(np1,nsg1,ncomp_g), w(np1)
      real(dp), intent(inout) :: Xs1(nsg1,ka), psi1(np1,ka,ncomp_g), aow1(np1,nsg1), Vl1(nsg1,nsg1)
      real(dp), intent(inout) :: t11(np1,*), t21(nsg1,*)
      real(dp), intent(inout) :: Va(nao,nao), Vb(nao,*), dd(nao,*), ex, tk(6)
      real(dp) :: ra(np1,4), rb(np1,4), wa(np1,4), wb(np1,4), wc(np1,4)
      integer(is) :: i, jj, c, g
      real(dp) :: tb

      tb = wall_time()
      do i = 1, ka
         do jj = 1, nsg1
            Xs1(jj,i) = X(idx(jj),i)
         enddo
      enddo
      do c = 1, ncomp_g
         call dgemm('N', 'N', np1, ka, nsg1, 1.0_dp, ao(1,1,c), np1, Xs1, nsg1, 0.0_dp, psi1(1,1,c), np1)
      enddo
      tk(2) = tk(2) + wall_time() - tb
      tb = wall_time()
      ra = 0.0_dp; rb = 0.0_dp
      do i = 1, ka
         do g = 1, np1
            ra(g,1) = ra(g,1) + psi1(g,i,1)*psi1(g,i,1)
         enddo
         do c = 2, ncomp_g
            do g = 1, np1
               ra(g,c) = ra(g,c) + 2.0_dp*psi1(g,i,1)*psi1(g,i,c)
            enddo
         enddo
         if (i == kb) rb = ra
      enddo
      if (kb == ka) rb = ra
      if (kb == 0) rb = 0.0_dp
      call xc_eval_pol(np1, ra, rb, w, wa, wb, ex)
      tk(3) = tk(3) + wall_time() - tb
      if (mode == 0) then
         call accumulate_v(np1, nsg1, idx, nr1, rp1, rl1, ao, wa, aow1, Vl1, Va, tk)
         if (.not. same) call accumulate_v(np1, nsg1, idx, nr1, rp1, rl1, ao, wb, aow1, Vl1, Vb, tk)
      else
         wc = 0.5_dp*(wa + wb)
         call accumulate_v(np1, nsg1, idx, nr1, rp1, rl1, ao, wc, aow1, Vl1, Va, tk)
         if (ka > kb) then
            wc = wa - wb
            call make_aow(np1, nsg1, ao, wc, aow1)
            ! (V_a - V_b) C_O = phi^T (aow_d C_O) + aow_d^T psi_O
            call dgemm('N', 'N', np1, ka-kb, nsg1, 1.0_dp, aow1, np1, Xs1(1,kb+1), nsg1, 0.0_dp, t11, np1)
            call dgemm('T', 'N', nsg1, ka-kb, np1, 1.0_dp, ao(1,1,1), np1, t11, np1, 0.0_dp, t21, nsg1)
            call dgemm('T', 'N', nsg1, ka-kb, np1, 1.0_dp, aow1, np1, psi1(1,kb+1,1), np1, 1.0_dp, t21, nsg1)
            do c = 1, ka-kb
               do jj = 1, nsg1
                  dd(idx(jj),c) = dd(idx(jj),c) + t21(jj,c)
               enddo
            enddo
         endif
      endif

    end subroutine kern

  end subroutine grid_pass_roks

!######################################################################
! grid_pass_rks: closed shell, D = 2 X X^T: exc = E_xc, V
!######################################################################
  subroutine grid_pass_rks(nao, k, X, exc, V)

    integer(is), intent(in) :: nao, k
    real(dp), intent(in)    :: X(nao,k)
    real(dp), intent(out)   :: exc, V(nao,nao)

    real(dp), allocatable :: Vt(:,:,:), ext(:), tt(:,:)
    real(dp), allocatable :: aobuf(:), Xs(:), psi(:), aow(:), Vl(:)
    integer(is) :: nt, tid, iblk, g0, np, nsg, j, t, nr
    real(dp)    :: t0, ta

    t0 = wall_time()
    nt = 1
    !$ nt = omp_get_max_threads()
    allocate(Vt(nao,nao,nt), ext(nt), tt(6,nt))
    Vt = 0.0_dp; ext = 0.0_dp; tt = 0.0_dp

    !$omp parallel private(tid, iblk, g0, np, nsg, nr, aobuf, Xs, psi, aow, Vl, ta)
    tid = 1
    !$ tid = omp_get_thread_num() + 1
    call blas_sequential()
    if (.not. cached) allocate(aobuf(nbk*nsigmax*ncomp_g))
    allocate(Xs(nsigmax*k), psi(nbk*k*ncomp_g), aow(nbk*nsigmax), Vl(nsigmax*nsigmax))
    !$omp do schedule(dynamic)
    do iblk = 1, nblk
       g0  = bbeg(iblk)
       np  = bend(iblk) - g0 + 1
       nsg = aop(iblk+1) - aop(iblk)
       if (nsg == 0) cycle
       nr  = rnp(iblk+1) - rnp(iblk)
       if (cached) then
          call kern(np, nsg, aol(aop(iblk)), nr, rpos(rnp(iblk)), rlen(rnp(iblk)), cache(coff(iblk)+1), &
               gw(g0), Xs, psi, aow, Vl, Vt(1,1,tid), ext(tid), tt(1,tid))
       else
          ta = wall_time()
          call gto_eval(np, gx(g0), gy(g0), gz(g0), shp(iblk+1)-shp(iblk), shl(shp(iblk)), &
               merge(1_is, 0_is, ncomp_g == 4), np, nsg, aobuf)
          tt(1,tid) = tt(1,tid) + wall_time() - ta
          call kern(np, nsg, aol(aop(iblk)), nr, rpos(rnp(iblk)), rlen(rnp(iblk)), aobuf, gw(g0), Xs, psi, &
               aow, Vl, Vt(1,1,tid), ext(tid), tt(1,tid))
       endif
    enddo
    !$omp end do
    if (allocated(aobuf)) deallocate(aobuf)
    deallocate(Xs, psi, aow, Vl)
    call blas_restore()
    !$omp end parallel
    tsub = tsub + sum(tt, dim=2)

    !$omp parallel do private(j, t)
    do j = 1, nao
       V(:,j) = Vt(:,j,1)
       do t = 2, nt
          V(:,j) = V(:,j) + Vt(:,j,t)
       enddo
    enddo
    !$omp end parallel do
    call symmetrise_sum(nao, V)
    exc = sum(ext)
    deallocate(Vt, ext, tt)
    time_xcpass = time_xcpass + wall_time() - t0

  contains

    subroutine kern(np1, nsg1, idx, nr1, rp1, rl1, ao, w, Xs1, psi1, aow1, Vl1, Vacc, ex, tk)
      integer(is), intent(in) :: np1, nsg1, idx(nsg1), nr1, rp1(nr1), rl1(nr1)
      real(dp), intent(in)    :: ao(np1,nsg1,ncomp_g), w(np1)
      real(dp), intent(inout) :: Xs1(nsg1,k), psi1(np1,k,ncomp_g), aow1(np1,nsg1), Vl1(nsg1,nsg1)
      real(dp), intent(inout) :: Vacc(nao,nao), ex, tk(6)
      real(dp) :: r(np1,4), wv(np1,4)
      integer(is) :: i, jj, c, g
      real(dp) :: tb

      tb = wall_time()
      do i = 1, k
         do jj = 1, nsg1
            Xs1(jj,i) = X(idx(jj),i)
         enddo
      enddo
      do c = 1, ncomp_g
         call dgemm('N', 'N', np1, k, nsg1, 1.0_dp, ao(1,1,c), np1, Xs1, nsg1, 0.0_dp, psi1(1,1,c), np1)
      enddo
      tk(2) = tk(2) + wall_time() - tb
      tb = wall_time()
      r = 0.0_dp
      do i = 1, k
         do g = 1, np1
            r(g,1) = r(g,1) + 2.0_dp*psi1(g,i,1)*psi1(g,i,1)
         enddo
         do c = 2, ncomp_g
            do g = 1, np1
               r(g,c) = r(g,c) + 4.0_dp*psi1(g,i,1)*psi1(g,i,c)
            enddo
         enddo
      enddo
      call xc_eval_unpol(np1, r, w, wv, ex)
      tk(3) = tk(3) + wall_time() - tb
      call accumulate_v(np1, nsg1, idx, nr1, rp1, rl1, ao, wv, aow1, Vl1, Vacc, tk)
    end subroutine kern

  end subroutine grid_pass_rks

!######################################################################
! accumulate_v: aow = sum_c wv(:,c) phi_c; Vacc(idx,idx) += phi_0^T aow
!######################################################################
  subroutine accumulate_v(np, nsg, idx, nrun, rp, rl, ao, wv, aow, Vl, Vacc, tk)

    integer(is), intent(in) :: np, nsg, idx(nsg), nrun, rp(nrun), rl(nrun)
    real(dp), intent(in)    :: ao(np,nsg,ncomp_g), wv(np,4)
    real(dp), intent(inout) :: aow(np,nsg), Vl(nsg,nsg), Vacc(nao_b,nao_b), tk(6)
    integer(is) :: i, j, r, a0, p0, jj
    real(dp) :: tb

    tb = wall_time()
    call make_aow(np, nsg, ao, wv, aow)
    tk(4) = tk(4) + wall_time() - tb
    tb = wall_time()
    if (nsg == nao_b) then
       call dgemm('T', 'N', nsg, nsg, np, 1.0_dp, ao, np, aow, np, 1.0_dp, Vacc, nao_b)
       tk(5) = tk(5) + wall_time() - tb
    else
       call dgemm('T', 'N', nsg, nsg, np, 1.0_dp, ao, np, aow, np, 0.0_dp, Vl, nsg)
       tk(5) = tk(5) + wall_time() - tb
       tb = wall_time()
       do j = 1, nsg
          jj = idx(j)
          do r = 1, nrun
             p0 = rp(r)
             a0 = idx(p0)
             do i = 0, rl(r)-1
                Vacc(a0+i,jj) = Vacc(a0+i,jj) + Vl(p0+i,j)
             enddo
          enddo
       enddo
       tk(6) = tk(6) + wall_time() - tb
    endif

  end subroutine accumulate_v

!######################################################################
! make_aow: aow = sum_c wv(:,c) phi_c (c = value and gradient components)
!######################################################################
  subroutine make_aow(np, nsg, ao, wv, aow)

    integer(is), intent(in) :: np, nsg
    real(dp), intent(in)    :: ao(np,nsg,ncomp_g), wv(np,4)
    real(dp), intent(out)   :: aow(np,nsg)
    integer(is) :: j, c, g

    do j = 1, nsg
       do g = 1, np
          aow(g,j) = wv(g,1)*ao(g,j,1)
       enddo
       do c = 2, ncomp_g
          do g = 1, np
             aow(g,j) = aow(g,j) + wv(g,c)*ao(g,j,c)
          enddo
       enddo
    enddo

  end subroutine make_aow

  subroutine symmetrise_sum(n, A)
    integer(is), intent(in) :: n
    real(dp), intent(inout) :: A(n,n)
    integer(is) :: i, j
    real(dp)    :: t
    do j = 1, n
       do i = 1, j
          t = A(i,j) + A(j,i)
          A(i,j) = t
          A(j,i) = t
       enddo
    enddo
  end subroutine symmetrise_sum

end module scf_grid
