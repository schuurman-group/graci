!**********************************************************************
! mrsf_xcgrid: exchange-correlation grid contractions of the MRSF
! gradient. PySCF supplies, once, the quadrature weights, the libxc
! functional derivatives and the AO values on the grid (optionally
! cached here); the kernel action on factor-pair densities and the XC
! probe (gradient) terms are evaluated here without forming any
! nao x nao matrix.
!
! Grid data (PySCF conventions): ncomp = 1 (LDA: AO values) or 4 (GGA:
! value, d/dx, d/dy, d/dz); ndc = ncomp density components (rho, drho);
!   fxc(g, y, t, x, s) = d2 f / d rho_{s,x} d rho_{t,y}   (PySCF fxc[s,x,t,y,g])
!   vxc(g, x, s)       = d f / d rho_{s,x}                (PySCF vxc[s,x,g])
! AO blocks arrive in PySCF's native buffer layout, ao(g, mu, c) with
! the grid point fastest (zero copy); all elementwise work runs on
! unit-stride grid vectors and all contractions are dgemm calls.
!
! Symmetric AO densities are given as factor pairs D = L R^T + R L^T with
! hole-width factors L, R (nao x k); R may be the hole MO block C_H, whose
! grid values are cached (psiH). The kernel potential of a density,
!   V^t_{mu nu} = sum_g [ wv0 phi_mu phi_nu + sum_c wv_c (d_c phi_mu phi_nu + phi_mu d_c phi_nu) ],
!   wv(g, y, t) = w(g) sum_{s,x} rho1(g, x, s) fxc(g, y, t, x, s),
! is projected onto the HP and HH blocks through the (nocca x nao)
! matrices M = sum_g [ psi_H^T aow + (C_H^T aow^T) phi ],
! aow = wv0/2 phi + sum_c wv_c d_c phi, so that V_HP = M C_P, V_HH = M C_H.
!
! Storage: cached AO values block-contiguous in the flat array aoc, block
! ib starting at aoff(ib)+1 with layout (nb, nao, ncomp); psiH likewise
! (nb, nocca, ncomp) at poff(ib)+1. A kernel evaluation is bracketed by
! xc_begin / [xc_cached | xc_block per block] / xc_end, so that the AO
! values can also be streamed in from PySCF when the cache does not fit
! the memory budget; the probe terms (xc_probe_*) always stream the AO
! values with the derivatives they need.
!**********************************************************************
module mrsf_xcgrid

  use iso_c_binding, only: c_int
  use mrsf_constants
  use mrsf_global
  use mrsf_io
!$ use omp_lib

  implicit none

#ifdef USE_MKL
  interface
     function mkl_threads_local(nt) bind(c, name='MKL_Set_Num_Threads_Local') result(prev)
       import :: c_int
       integer(c_int), value :: nt
       integer(c_int)        :: prev
     end function mkl_threads_local
  end interface
#endif

  integer, parameter       :: i8 = selected_int_kind(18)
  integer(is)              :: ngrid = 0, nblocks = 0, ncomp = 0, ndc = 0, nao_g = 0
  integer(is), allocatable :: gbeg(:), gend(:)
  integer(i8), allocatable :: aoff(:), poff(:)
  real(dp), allocatable    :: aoc(:)               ! cached AO values (see above)
  real(dp), allocatable    :: psiH(:)              ! hole MOs on the grid
  real(dp), allocatable    :: wgt(:)               ! (ngrid)
  real(dp), allocatable    :: fxcw(:,:,:,:,:)      ! (ndc,2,ndc,2,ngrid): w * fxc at (y,t,x,s,g)
  real(dp), allocatable    :: fxcw_cs(:,:,:)       ! (ndc,ndc,ngrid): w * (f_aa + f_ab) (closed-shell mode)
  real(dp), allocatable    :: vxc(:,:,:)           ! (ngrid, ndc, 2)
  ! closed-shell kernel set at the density of the configuration G of the
  ! extended method (EMRSF gradients): spin-0 derivatives of the total
  ! density in PySCF's convention, w f^(0) and w k^(0), and the unweighted
  ! v_xc; installed with xc_g_set after xc_init (same grid)
  real(dp), allocatable    :: vxcG(:,:)             ! (ngrid, ndc)
  real(dp), allocatable    :: fxcwG(:,:,:)          ! (ndc, ndc, ngrid) w f (y,x,g)
  real(dp), allocatable    :: kxcwG(:,:,:,:)        ! (ndc, ndc, ndc, ngrid) w k (z,y,x,g)
  logical                  :: g_ready = .false., g_have_kxc = .false.
  ! G-probe factors registered for a probe pass (xc_gprobe_set)
  integer(is)              :: kMg = 0, kSg = 0, nstg = 0
  real(dp), allocatable    :: LMg(:,:,:), RMg(:,:,:), LSg(:,:,:)
  real(dp), allocatable    :: taccG(:,:,:)          ! (nao, 3, nstg)
  real(dp)                 :: time_xcg = 0.0_dp
  integer(is)              :: nxcg_calls = 0
  real(dp), allocatable    :: CHg(:,:), CPg(:,:)   ! (nao, nocca), (nao, nvirb)
  logical                  :: xc_ready = .false., ao_cached = .false.
  logical                  :: ref_ready = .false.     ! reference-density kernel installed
  logical                  :: cs_ready = .false.      ! closed-shell kernel of the CV block installed
  ! quadratic-form probe (Part XIII): densities P = B S B^T, B = [C_H | Bexp]
  logical                  :: probe_q = .false., gq = .false., gq_kxc = .false.
  integer(is)              :: nexp_q = 0, maxnb_q = 0, nst_q = 0
  real(dp), allocatable    :: Bexp_q(:,:), Sq(:,:,:,:), Bq(:,:,:,:), Wq(:,:,:,:)
  integer(is), allocatable :: nbq(:,:), idxq(:,:,:)
  ! evaluation in progress
  integer(is)              :: kcur = 0, nvcur = 0
  real(dp), allocatable    :: Lcur(:,:,:,:), Rcur(:,:,:,:)   ! (nao, k, 2, nvec)
  integer(is), allocatable :: rch(:,:)                        ! (2, nvec): 1 -> R = C_H
  real(dp), allocatable    :: Macc(:,:,:,:)                   ! (nocca, nao, 2, nvec)
  real(dp), allocatable    :: tacc(:,:,:)                     ! (nao, 3, nst+1)
  real(dp)                 :: tstart = 0.0_dp
  ! second-derivative AO component of d_x d_c (1-based PySCF order:
  ! xx xy xz yy yz zz = 5..10)
  integer(is), parameter   :: d2idx(3,3) = reshape((/5,6,7, 6,8,9, 7,9,10/), (/3,3/))
  ! timings
  real(dp)                 :: time_xc = 0.0_dp, time_probe = 0.0_dp
  real(dp)                 :: tcs(5) = 0.0_dp   ! closed-shell kernel: psiL, rho+wv, aow, projections, final
  integer(is)              :: nxc_calls = 0, nxc_vecs = 0, nprobe_calls = 0
  real(dp)                 :: tpart(6) = 0.0_dp   ! kernel: psi, rho, wv, aow, projection, misc

contains

!######################################################################
! xc_init: grid dimensions and blocks, weights, kernel arrays, MO blocks
!######################################################################
  subroutine xc_init(nao1, ngrid1, ncomp1, nv1, nblocks1, gbeg1, gend1, wgt1, fxc1, vxc1, &
       CH1, CP1, cache_ao)

    integer(is), intent(in) :: nao1, ngrid1, ncomp1, nv1, nblocks1
    integer(is), intent(in) :: gbeg1(nblocks1), gend1(nblocks1)
    real(dp), intent(in)    :: wgt1(ngrid1), fxc1(ngrid1,nv1,2,nv1,2), vxc1(ngrid1,nv1,2)
    real(dp), intent(in)    :: CH1(nao1,nocca), CP1(nao1,nvirb)
    logical, intent(in)     :: cache_ao

    call xc_init_grid(nao1, ngrid1, ncomp1, nv1, nblocks1, gbeg1, gend1, wgt1, CH1, CP1, cache_ao)
    call xc_set_ref_kernel(ngrid1, nv1, vxc1, fxc1)

  end subroutine xc_init

!######################################################################
! xc_init_grid: grid dimensions and blocks, weights, MO blocks and the
! AO cache (no kernel: the reference kernel is installed separately by
! xc_set_ref_kernel, the G-density set by xc_g_set)
!######################################################################
  subroutine xc_init_grid(nao1, ngrid1, ncomp1, nv1, nblocks1, gbeg1, gend1, wgt1, CH1, CP1, cache_ao)

    integer(is), intent(in) :: nao1, ngrid1, ncomp1, nv1, nblocks1
    integer(is), intent(in) :: gbeg1(nblocks1), gend1(nblocks1)
    real(dp), intent(in)    :: wgt1(ngrid1)
    real(dp), intent(in)    :: CH1(nao1,nocca), CP1(nao1,nvirb)
    logical, intent(in)     :: cache_ao

    integer(is) :: ib

    call xc_free()
    nao_g = nao1; ngrid = ngrid1; ncomp = ncomp1; ndc = nv1; nblocks = nblocks1
    if (ndc /= ncomp) call mrsf_error('xc_init: nv must equal ncomp')
    allocate(gbeg(nblocks), gend(nblocks), aoff(nblocks), poff(nblocks))
    gbeg = gbeg1; gend = gend1
    do ib = 1, nblocks
       aoff(ib) = int(ncomp,i8)*int(nao_g,i8)*int(gbeg(ib)-1,i8)
       poff(ib) = int(ncomp,i8)*int(nocca,i8)*int(gbeg(ib)-1,i8)
    enddo
    allocate(wgt(ngrid)); wgt = wgt1
    allocate(CHg(nao_g,nocca), CPg(nao_g,nvirb))
    CHg = CH1; CPg = CP1
    ao_cached = cache_ao
    if (ao_cached) then
       allocate(aoc(int(ncomp,i8)*int(nao_g,i8)*int(ngrid,i8)))
       allocate(psiH(int(ncomp,i8)*int(nocca,i8)*int(ngrid,i8)))
    endif
    xc_ready = .true.

  end subroutine xc_init_grid

!######################################################################
! xc_set_ref_kernel: potentials and weighted kernel of the reference
! (spin-resolved, PySCF layout vxc(g,x,s), fxc(g,y,t,x,s))
!######################################################################
  subroutine xc_set_ref_kernel(ngrid1, nv1, vxc1, fxc1)

    integer(is), intent(in) :: ngrid1, nv1
    real(dp), intent(in)    :: vxc1(ngrid1,nv1,2), fxc1(ngrid1,nv1,2,nv1,2)

    integer(is) :: g, s, x, t, y

    if (.not. xc_ready) call mrsf_error('xc_set_ref_kernel: grid not initialised')
    if (ngrid1 /= ngrid .or. nv1 /= ndc) call mrsf_error('xc_set_ref_kernel: dimension mismatch')
    if (allocated(vxc)) deallocate(vxc)
    if (allocated(fxcw)) deallocate(fxcw)
    allocate(vxc(ngrid,ndc,2)); vxc = vxc1
    allocate(fxcw(ndc,2,ndc,2,ngrid))
    !$omp parallel do private(s,x,t,y)
    do g = 1, ngrid
       do s = 1, 2
          do x = 1, ndc
             do t = 1, 2
                do y = 1, ndc
                   fxcw(y,t,x,s,g) = wgt(g)*fxc1(g,y,t,x,s)
                enddo
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
    ref_ready = .true.

  end subroutine xc_set_ref_kernel

!######################################################################
! xc_add_block: cache the AO values of block ib, ao(nb, nao, ncomp), and
! the hole MOs on its points
!######################################################################
  subroutine xc_add_block(ib, ao)

    integer(is), intent(in) :: ib
    real(dp), intent(in)    :: ao(*)

    integer(is) :: nb

    if (.not. ao_cached) return
    if (ib < 1 .or. ib > nblocks) call mrsf_error('xc_add_block: block index out of range')
    nb = gend(ib) - gbeg(ib) + 1
    call add_block(nb, ao, aoc(aoff(ib)+1), psiH(poff(ib)+1))

  end subroutine xc_add_block

  subroutine add_block(nb, ao, aob, psib)

    integer(is), intent(in) :: nb
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp)
    real(dp), intent(out)   :: aob(nb,nao_g,ncomp), psib(nb,nocca,ncomp)

    integer(is) :: c, mu

    !$omp parallel do private(mu)
    do c = 1, ncomp
       do mu = 1, nao_g
          aob(:,mu,c) = ao(:,mu,c)
       enddo
    enddo
    !$omp end parallel do
    do c = 1, ncomp
       call dgemm('N','N', nb, nocca, nao_g, 1.0_dp, ao(1,1,c), nb, CHg, nao_g, 0.0_dp, psib(1,1,c), nb)
    enddo

  end subroutine add_block

!######################################################################
! xc_free / xc_cleanup
!######################################################################
  subroutine xc_free()

    call xc_cleanup()
    if (allocated(gbeg)) deallocate(gbeg, gend, aoff, poff)
    if (allocated(wgt)) deallocate(wgt)
    if (allocated(fxcw)) deallocate(fxcw)
    if (allocated(fxcw_cs)) deallocate(fxcw_cs)
    if (allocated(vxc)) deallocate(vxc)
    if (allocated(aoc)) deallocate(aoc)
    if (allocated(psiH)) deallocate(psiH)
    if (allocated(CHg)) deallocate(CHg, CPg)
    call xc_g_free()
    xc_ready = .false.; ao_cached = .false.; ref_ready = .false.; cs_ready = .false.

  end subroutine xc_free

  subroutine xc_cleanup()

    if (allocated(Lcur)) deallocate(Lcur, Rcur, rch)
    if (allocated(Macc)) deallocate(Macc)
    if (allocated(tacc)) deallocate(tacc)
    if (allocated(LMg)) deallocate(LMg, RMg, LSg)
    if (allocated(taccG)) deallocate(taccG)
    if (allocated(Bexp_q)) deallocate(Bexp_q, nbq, idxq, Sq, Bq, Wq)
    kcur = 0; nvcur = 0; nstg = 0; probe_q = .false.; gq = .false.; gq_kxc = .false.

  end subroutine xc_cleanup

!######################################################################
! set_factors: store the factor pairs of an evaluation, Lfac/Rfac(nao,
! k, 2, nvec); the R factors flagged by rch1(s,v) = 1 are replaced by C_H
!######################################################################
  subroutine set_factors(nvec, k, Lfac, Rfac, rch1)

    integer(is), intent(in) :: nvec, k, rch1(2,nvec)
    real(dp), intent(in)    :: Lfac(nao_g,k,2,nvec), Rfac(nao_g,k,2,nvec)

    integer(is) :: v, s

    if (.not. xc_ready) call mrsf_error('set_factors: grid not initialised')
    if (.not. ref_ready) call mrsf_error('set_factors: reference kernel not installed')
    if (any(rch1 == 1) .and. k /= nocca) call mrsf_error('set_factors: R = C_H requires k = nocca')
    call xc_cleanup()
    tstart = wall_time()
    kcur = k; nvcur = nvec
    allocate(Lcur(nao_g,k,2,nvec), Rcur(nao_g,k,2,nvec), rch(2,nvec))
    Lcur = Lfac; Rcur = Rfac; rch = rch1
    do v = 1, nvec
       do s = 1, 2
          if (rch(s,v) == 1) Rcur(:,:,s,v) = CHg
       enddo
    enddo

  end subroutine set_factors

!######################################################################
! xc_begin: start a kernel evaluation for nvec factor-pair densities
!######################################################################
  subroutine xc_begin(nvec, k, Lfac, Rfac, rch1)

    integer(is), intent(in) :: nvec, k, rch1(2,nvec)
    real(dp), intent(in)    :: Lfac(nao_g,k,2,nvec), Rfac(nao_g,k,2,nvec)

    call set_factors(nvec, k, Lfac, Rfac, rch1)
    allocate(Macc(nocca,nao_g,2,nvec))
    Macc = 0.0_dp

  end subroutine xc_begin

!######################################################################
! xc_cached: all blocks from the cache
!######################################################################
  subroutine xc_cached()

    integer(is) :: ib, nb

    if (.not. ao_cached) call mrsf_error('xc_cached: AO values are not cached')
    if (.not. allocated(Macc)) call mrsf_error('xc_cached: no kernel evaluation in progress')
    call cached_parallel()

  end subroutine xc_cached

!######################################################################
! cached_parallel: all cached blocks, OpenMP over blocks with sequential
! MKL inside; per-thread flat work buffers (mapped onto explicit-shape
! dummies of kernel_block_thread with the block's own nb) and per-thread
! accumulators reduced into Macc once at the end
!######################################################################
  subroutine cached_parallel()

    integer(is)           :: ib, nb, nbmax, nk, nthr
    integer(c_int)        :: prev
    logical               :: par
    real(dp)              :: tl(5)
    real(dp), allocatable :: psiL(:), psiR(:), rho1(:), wv(:), Bst(:), Ast(:), Mth(:,:,:,:)

    nk = kcur*2*nvcur
    nbmax = 0
    do ib = 1, nblocks
       nbmax = max(nbmax, gend(ib) - gbeg(ib) + 1)
    enddo
    nthr = 1
    !$ nthr = omp_get_max_threads()
    par = (nthr > 1) .and. (nblocks >= nthr)
    tl  = 0.0_dp

    !$omp parallel if(par) default(shared) private(ib, nb, prev, psiL, psiR, rho1, wv, Bst, Ast, Mth) &
    !$omp reduction(+:tl)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(1_c_int)
#endif
    allocate(psiL(nbmax*nk*ncomp), psiR(max(nbmax*kcur*ncomp,1_is)), rho1(nbmax*ndc*2), &
         wv(nbmax*ndc*2), Bst(2*nbmax*nao_g*2), Ast(nocca*2*nbmax*2), Mth(nocca,nao_g,2,nvcur))
    Mth = 0.0_dp
    !$omp do schedule(dynamic)
    do ib = 1, nblocks
       nb = gend(ib) - gbeg(ib) + 1
       call kernel_block_thread(nb, gbeg(ib), aoc(aoff(ib)+1), psiH(poff(ib)+1), psiL, psiR, &
            rho1, wv, Bst, Ast, Mth, tl)
    enddo
    !$omp end do
    !$omp critical
    Macc = Macc + Mth
    !$omp end critical
    deallocate(psiL, psiR, rho1, wv, Bst, Ast, Mth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel
    tpart(1:5) = tpart(1:5) + tl(1:5)

  end subroutine cached_parallel

!######################################################################
! xc_block: one streamed block, ao(nb, nao, ncomp)
!######################################################################
  subroutine xc_block(ib, ao)

    integer(is), intent(in) :: ib
    real(dp), intent(in)    :: ao(*)

    integer(is) :: nb, c
    real(dp), allocatable :: psiHw(:,:,:)

    if (ib < 1 .or. ib > nblocks) call mrsf_error('xc_block: block index out of range')
    if (.not. allocated(Macc)) call mrsf_error('xc_block: no kernel evaluation in progress')
    nb = gend(ib) - gbeg(ib) + 1
    allocate(psiHw(nb,nocca,ncomp))
    do c = 1, ncomp
       call dgemm('N','N', nb, nocca, nao_g, 1.0_dp, ao(1+nb*nao_g*(c-1)), nb, CHg, nao_g, &
            0.0_dp, psiHw(1,1,c), nb)
    enddo
    call kernel_block(ib, nb, ao, psiHw)
    deallocate(psiHw)

  end subroutine xc_block

!######################################################################
! xc_end: project the accumulated matrices onto the HP / HH blocks
!######################################################################
  subroutine xc_end(want_hp, want_hh, VHP, VHH)

    logical, intent(in)   :: want_hp, want_hh
    real(dp), intent(out) :: VHP(nocca,nvirb,2,nvcur), VHH(nocca,nocca,2,nvcur)

    integer(is) :: v, t

    if (.not. allocated(Macc)) call mrsf_error('xc_end: no kernel evaluation in progress')
    do v = 1, nvcur
       do t = 1, 2
          if (want_hp) call dgemm('N','N', nocca, nvirb, nao_g, 1.0_dp, Macc(1,1,t,v), nocca, &
               CPg, nao_g, 0.0_dp, VHP(1,1,t,v), nocca)
          if (want_hh) call dgemm('N','N', nocca, nocca, nao_g, 1.0_dp, Macc(1,1,t,v), nocca, &
               CHg, nao_g, 0.0_dp, VHH(1,1,t,v), nocca)
       enddo
    enddo
    time_xc = time_xc + wall_time() - tstart
    nxc_calls = nxc_calls + 1
    nxc_vecs = nxc_vecs + nvcur
    call xc_cleanup()

  end subroutine xc_end

!######################################################################
! kernel_block: densities, kernel potentials and projections of one
! grid block, ao(nb, nao, ncomp), psiHb(nb, nocca, ncomp)
!######################################################################
  subroutine kernel_block(ib, nb, ao, psiHb)

    integer(is), intent(in) :: ib, nb
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp), psiHb(nb,nocca,ncomp)

    integer(is) :: g0, nk, v, s, t, c, g, gg, x, y, mu, jc
    real(dp), allocatable :: psiL(:,:,:), psiR(:,:,:), rho1(:,:,:), wv(:,:,:), aow(:,:,:), aowCH(:,:)
    real(dp) :: acc, t0, t1

    g0 = gbeg(ib)
    nk = kcur*2*nvcur
    t0 = wall_time()

    ! MO-factor values on the block for all densities at once
    allocate(psiL(nb,nk,ncomp))
    do c = 1, ncomp
       call dgemm('N','N', nb, nk, nao_g, 1.0_dp, ao(1,1,c), nb, Lcur, nao_g, 0.0_dp, psiL(1,1,c), nb)
    enddo
    if (any(rch == 0)) allocate(psiR(nb,kcur,ncomp))
    allocate(rho1(nb,ndc,2), wv(nb,ndc,2), aow(nb,nao_g,2), aowCH(nocca,nb))
    t1 = wall_time(); tpart(1) = tpart(1) + t1 - t0; t0 = t1

    do v = 1, nvcur
       ! trial densities of both spins
       do s = 1, 2
          jc = kcur*(s-1) + 2*kcur*(v-1) + 1
          if (rch(s,v) == 1) then
             call block_rho(nb, kcur, psiL(1,jc,1), nb*nk, psiHb, nb*nocca, rho1(1,1,s))
          else
             do c = 1, ncomp
                call dgemm('N','N', nb, kcur, nao_g, 1.0_dp, ao(1,1,c), nb, Rcur(1,1,s,v), nao_g, &
                     0.0_dp, psiR(1,1,c), nb)
             enddo
             call block_rho(nb, kcur, psiL(1,jc,1), nb*nk, psiR, nb*kcur, rho1(1,1,s))
          endif
       enddo
       t1 = wall_time(); tpart(2) = tpart(2) + t1 - t0; t0 = t1
       ! weighted kernel potentials wv(g,y,t) = w sum_{s,x} rho1(g,x,s) fxc(g,y,t,x,s)
       !$omp parallel do private(gg,t,y,s,x,acc)
       do g = 1, nb
          gg = g0 + g - 1
          do t = 1, 2
             do y = 1, ndc
                acc = 0.0_dp
                do s = 1, 2
                   do x = 1, ndc
                      acc = acc + rho1(g,x,s)*fxcw(y,t,x,s,gg)
                   enddo
                enddo
                wv(g,y,t) = acc
             enddo
          enddo
       enddo
       !$omp end parallel do
       t1 = wall_time(); tpart(3) = tpart(3) + t1 - t0; t0 = t1
       ! aow(:,mu,t) = wv0/2 phi + sum_c wv_c d_c phi
       !$omp parallel do private(t,c)
       do mu = 1, nao_g
          do t = 1, 2
             aow(:,mu,t) = 0.5_dp*wv(:,1,t)*ao(:,mu,1)
             do c = 2, ncomp
                aow(:,mu,t) = aow(:,mu,t) + wv(:,c,t)*ao(:,mu,c)
             enddo
          enddo
       enddo
       !$omp end parallel do
       t1 = wall_time(); tpart(4) = tpart(4) + t1 - t0; t0 = t1
       ! M += psiH^T aow + (C_H^T aow^T) phi
       do t = 1, 2
          call dgemm('T','N', nocca, nao_g, nb, 1.0_dp, psiHb, nb, aow(1,1,t), nb, &
               1.0_dp, Macc(1,1,t,v), nocca)
          call dgemm('T','T', nocca, nb, nao_g, 1.0_dp, CHg, nao_g, aow(1,1,t), nb, &
               0.0_dp, aowCH, nocca)
          call dgemm('N','N', nocca, nao_g, nb, 1.0_dp, aowCH, nocca, ao(1,1,1), nb, &
               1.0_dp, Macc(1,1,t,v), nocca)
       enddo
       t1 = wall_time(); tpart(5) = tpart(5) + t1 - t0; t0 = t1
    enddo

    deallocate(psiL, rho1, wv, aow, aowCH)
    if (allocated(psiR)) deallocate(psiR)
    tpart(6) = tpart(6) + wall_time() - t0

  end subroutine kernel_block

!######################################################################
! kernel_block_thread: one grid block for all nvcur densities, run by
! one thread with sequential BLAS: psiL = ao L for all densities (one
! dgemm per component), then per vector and spin the trial density, the
! kernel potential, aow (first half of Bst) and the weighted hole-MO
! products aowH from the cached psiH (second half of Ast), and ONE
! projection dgemm per spin, M(i,nu,t) += [psi_H^T | aowH_t] [aow_t ; phi].
! tl(1:5): factor values, densities, potentials, weighted products,
! projections (thread seconds)
!######################################################################
  subroutine kernel_block_thread(nb, g0, ao, psiHb, psiL, psiR, rho1, wv, Bst, Ast, Mth, tl)

    integer(is), intent(in) :: nb, g0
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp), psiHb(nb,nocca,ncomp)
    real(dp), intent(inout) :: psiL(nb,kcur*2*nvcur,ncomp), psiR(nb,kcur,ncomp)
    real(dp), intent(inout) :: rho1(nb,ndc,2), wv(nb,ndc,2)
    real(dp), intent(inout) :: Bst(2*nb,nao_g,2), Ast(nocca,2*nb,2)
    real(dp), intent(inout) :: Mth(nocca,nao_g,2,nvcur), tl(5)

    integer(is) :: nk, v, s, t, c, g, gg, x, y, mu, jc, i
    real(dp)    :: acc, t1

    nk = kcur*2*nvcur

    ! block-constant halves of the stacked operands
    do t = 1, 2
       do g = 1, nb
          Ast(1:nocca,g,t) = psiHb(g,1:nocca,1)
       enddo
       do mu = 1, nao_g
          Bst(nb+1:2*nb,mu,t) = ao(1:nb,mu,1)
       enddo
    enddo

    ! MO-factor values on the block for all densities at once
    t1 = wall_time()
    do c = 1, ncomp
       call dgemm('N','N', nb, nk, nao_g, 1.0_dp, ao(1,1,c), nb, Lcur, nao_g, 0.0_dp, psiL(1,1,c), nb)
    enddo
    tl(1) = tl(1) + wall_time() - t1

    do v = 1, nvcur
       ! trial densities of both spins
       t1 = wall_time()
       do s = 1, 2
          jc = kcur*(s-1) + 2*kcur*(v-1) + 1
          if (rch(s,v) == 1) then
             call block_rho(nb, kcur, psiL(1,jc,1), nb*nk, psiHb, nb*nocca, rho1(1,1,s))
          else
             do c = 1, ncomp
                call dgemm('N','N', nb, kcur, nao_g, 1.0_dp, ao(1,1,c), nb, Rcur(1,1,s,v), nao_g, &
                     0.0_dp, psiR(1,1,c), nb)
             enddo
             call block_rho(nb, kcur, psiL(1,jc,1), nb*nk, psiR, nb*kcur, rho1(1,1,s))
          endif
       enddo
       tl(2) = tl(2) + wall_time() - t1
       ! weighted kernel potentials wv(g,y,t) = w sum_{s,x} rho1(g,x,s) fxc(g,y,t,x,s)
       t1 = wall_time()
       do g = 1, nb
          gg = g0 + g - 1
          do t = 1, 2
             do y = 1, ndc
                acc = 0.0_dp
                do s = 1, 2
                   do x = 1, ndc
                      acc = acc + rho1(g,x,s)*fxcw(y,t,x,s,gg)
                   enddo
                enddo
                wv(g,y,t) = acc
             enddo
          enddo
       enddo
       tl(3) = tl(3) + wall_time() - t1
       ! aow_t (first half of Bst) and the weighted hole-MO products aowH_t
       ! (second half of Ast) from the cached psiH components
       t1 = wall_time()
       do t = 1, 2
          do mu = 1, nao_g
             Bst(1:nb,mu,t) = 0.5_dp*wv(1:nb,1,t)*ao(1:nb,mu,1)
             do c = 2, ncomp
                Bst(1:nb,mu,t) = Bst(1:nb,mu,t) + wv(1:nb,c,t)*ao(1:nb,mu,c)
             enddo
          enddo
          do g = 1, nb
             do i = 1, nocca
                acc = 0.5_dp*wv(g,1,t)*psiHb(g,i,1)
                do c = 2, ncomp
                   acc = acc + wv(g,c,t)*psiHb(g,i,c)
                enddo
                Ast(i,nb+g,t) = acc
             enddo
          enddo
       enddo
       tl(4) = tl(4) + wall_time() - t1
       ! projections, one dgemm per spin with K = 2 nb
       t1 = wall_time()
       do t = 1, 2
          call dgemm('N','N', nocca, nao_g, 2*nb, 1.0_dp, Ast(1,1,t), nocca, Bst(1,1,t), 2*nb, &
               1.0_dp, Mth(1,1,t,v), nocca)
       enddo
       tl(5) = tl(5) + wall_time() - t1
    enddo

  end subroutine kernel_block_thread

!######################################################################
! block_rho: rho(g, x) of the density L R^T + R L^T from the factor
! values psiL(g, kk, c), psiR(g, kk, c) on the block (grid point fastest,
! component stride csL / csR)
!######################################################################
  subroutine block_rho(nb, k, psiL, csL, psiR, csR, rho)

    integer(is), intent(in) :: nb, k, csL, csR
    real(dp), intent(in)    :: psiL(*), psiR(*)
    real(dp), intent(out)   :: rho(nb,ndc)

    integer(is) :: kk, c, oL, oR, oLc, oRc

    rho = 0.0_dp
    do kk = 1, k
       oL = nb*(kk-1); oR = nb*(kk-1)
       rho(:,1) = rho(:,1) + 2.0_dp*psiL(oL+1:oL+nb)*psiR(oR+1:oR+nb)
       do c = 2, ndc
          oLc = oL + csL*(c-1); oRc = oR + csR*(c-1)
          rho(:,c) = rho(:,c) + 2.0_dp*(psiL(oLc+1:oLc+nb)*psiR(oR+1:oR+nb) &
               + psiL(oL+1:oL+nb)*psiR(oRc+1:oRc+nb))
       enddo
    enddo

  end subroutine block_rho

!######################################################################
! XC probe term of the gradient and the reference XC gradient on the
! fixed grid. For a weighted potential wv(g, 0:3) and a symmetric AO
! density X, PySCF's AO-derivative matrices (rks_grad._gga_grad_sum_)
! contracted with X give the per-AO vectors
!   t_x(mu) = sum_g [ d_x phi_mu A_mu + aow2_x,mu Phi_mu ],
!   Phi_mu = sum_nu X_mu,nu phi_nu,  A_mu = sum_nu X_mu,nu aow_nu,
!   aow = wv0/2 phi + sum_c wv_c d_c phi,  aow2_x = wv0/2 d_x phi + sum_c wv_c d_x d_c phi,
! and dE/dR_{A,x} = -2 sum_{mu in A} t_x(mu). Output entry 1: the
! reference term (vxc, D^s); entry v+1: the state term (vxc, P^s_v) +
! (f_xc[P_v], D^s); both summed over the spins s. The state densities
! P^s_v = L R^T + R L^T are given as factor pairs (xc_probe_begin); the
! AO values with second derivatives (GGA, 10 components) or first
! derivatives (LDA, 4) are streamed block by block (xc_probe_block).
!######################################################################
  subroutine xc_probe_begin(nst, k, Lfac, Rfac, rch1)

    integer(is), intent(in) :: nst, k, rch1(2,nst)
    real(dp), intent(in)    :: Lfac(nao_g,k,2,nst), Rfac(nao_g,k,2,nst)

    call set_factors(nst, k, Lfac, Rfac, rch1)
    allocate(tacc(nao_g,3,nst+1))
    tacc = 0.0_dp

  end subroutine xc_probe_begin

  subroutine xc_probe_block(ib, nc2, ao)

    integer(is), intent(in) :: ib, nc2
    real(dp), intent(in)    :: ao(*)

    integer(is) :: nb

    if (.not. allocated(tacc)) call mrsf_error('xc_probe_block: no probe evaluation in progress')
    if (ib < 1 .or. ib > nblocks) call mrsf_error('xc_probe_block: block index out of range')
    if ((ncomp == 1 .and. nc2 < 4) .or. (ncomp == 4 .and. nc2 < 10)) &
         call mrsf_error('xc_probe_block: AO derivative components missing')
    nb = gend(ib) - gbeg(ib) + 1
    call probe_block(ib, nb, nc2, ao)

  end subroutine xc_probe_block

  subroutine probe_block(ib, nb, nc2, ao)

    integer(is), intent(in) :: ib, nb, nc2
    real(dp), intent(in)    :: ao(nb,nao_g,nc2)

    integer(is) :: g0, nk, v, s, t, c, g, gg, x, y, jc, nocc(2)
    real(dp), allocatable :: psiHb(:,:,:), psiL(:,:,:), psiR(:,:,:), rho1(:,:,:), wvR(:,:,:), wv1(:,:,:)
    real(dp), allocatable :: aowR(:,:,:), aow2R(:,:,:,:), aow1(:,:,:), aow2_1(:,:,:,:)
    real(dp), allocatable :: PhiD(:,:,:), Phi(:,:), Amat(:,:), tmpL(:,:), tmpR(:,:)
    real(dp) :: acc

    g0 = gbeg(ib)
    nocc(1) = nocca; nocc(2) = nocca - 2

    ! hole MOs on the block (value + gradient components)
    allocate(psiHb(nb,nocca,ncomp))
    if (ao_cached) then
       call copy_block(nb*nocca*ncomp, psiH(poff(ib)+1), psiHb)
    else
       do c = 1, ncomp
          call dgemm('N','N', nb, nocca, nao_g, 1.0_dp, ao(1,1,c), nb, CHg, nao_g, 0.0_dp, psiHb(1,1,c), nb)
       enddo
    endif

    ! reference potentials, their AO products and the reference densities
    allocate(wvR(nb,ndc,2), aowR(nb,nao_g,2), aow2R(nb,nao_g,3,2))
    do s = 1, 2
       do c = 1, ndc
          wvR(:,c,s) = wgt(g0:g0+nb-1)*vxc(g0:g0+nb-1,c,s)
       enddo
    enddo
    call make_aow(nb, nc2, ao, wvR, aowR, aow2R)
    allocate(PhiD(nb,nao_g,2), Phi(nb,nao_g), Amat(nb,nao_g))
    allocate(tmpL(nb,max(kcur,nocca)), tmpR(nb,max(kcur,nocca)))
    do s = 1, 2
       call dgemm('N','T', nb, nao_g, nocc(s), 1.0_dp, psiHb, nb, CHg, nao_g, 0.0_dp, PhiD(1,1,s), nb)
       call dgemm('N','N', nb, nocc(s), nao_g, 1.0_dp, aowR(1,1,s), nb, CHg, nao_g, 0.0_dp, tmpL, nb)
       call dgemm('N','T', nb, nao_g, nocc(s), 1.0_dp, tmpL, nb, CHg, nao_g, 0.0_dp, Amat, nb)
       call probe_reduce(nb, nc2, ao, Amat, aow2R(1,1,1,s), PhiD(1,1,s), tacc(1,1,1))
    enddo

    ! state densities
    if (nvcur > 0) then
       nk = kcur*2*nvcur
       allocate(psiL(nb,nk,ncomp), psiR(nb,kcur,ncomp))
       do c = 1, ncomp
          call dgemm('N','N', nb, nk, nao_g, 1.0_dp, ao(1,1,c), nb, Lcur, nao_g, 0.0_dp, psiL(1,1,c), nb)
       enddo
       allocate(rho1(nb,ndc,2), wv1(nb,ndc,2), aow1(nb,nao_g,2), aow2_1(nb,nao_g,3,2))
       do v = 1, nvcur
          do s = 1, 2
             jc = kcur*(s-1) + 2*kcur*(v-1) + 1
             if (rch(s,v) == 1) then
                psiR = psiHb
             else
                do c = 1, ncomp
                   call dgemm('N','N', nb, kcur, nao_g, 1.0_dp, ao(1,1,c), nb, Rcur(1,1,s,v), nao_g, &
                        0.0_dp, psiR(1,1,c), nb)
                enddo
             endif
             call block_rho(nb, kcur, psiL(1,jc,1), nb*nk, psiR, nb*kcur, rho1(1,1,s))
             ! (vxc, P^s_v): Phi = psiR L^T + psiL R^T, A = (aow R) L^T + (aow L) R^T
             call dgemm('N','T', nb, nao_g, kcur, 1.0_dp, psiR, nb, Lcur(1,1,s,v), nao_g, 0.0_dp, Phi, nb)
             call dgemm('N','T', nb, nao_g, kcur, 1.0_dp, psiL(1,jc,1), nb, Rcur(1,1,s,v), nao_g, &
                  1.0_dp, Phi, nb)
             call dgemm('N','N', nb, kcur, nao_g, 1.0_dp, aowR(1,1,s), nb, Rcur(1,1,s,v), nao_g, &
                  0.0_dp, tmpR, nb)
             call dgemm('N','N', nb, kcur, nao_g, 1.0_dp, aowR(1,1,s), nb, Lcur(1,1,s,v), nao_g, &
                  0.0_dp, tmpL, nb)
             call dgemm('N','T', nb, nao_g, kcur, 1.0_dp, tmpR, nb, Lcur(1,1,s,v), nao_g, 0.0_dp, Amat, nb)
             call dgemm('N','T', nb, nao_g, kcur, 1.0_dp, tmpL, nb, Rcur(1,1,s,v), nao_g, 1.0_dp, Amat, nb)
             call probe_reduce(nb, nc2, ao, Amat, aow2R(1,1,1,s), Phi, tacc(1,1,v+1))
          enddo
          ! (f_xc[P_v], D^t): weighted kernel potential of P_v
          !$omp parallel do private(gg,t,y,s,x,acc)
          do g = 1, nb
             gg = g0 + g - 1
             do t = 1, 2
                do y = 1, ndc
                   acc = 0.0_dp
                   do s = 1, 2
                      do x = 1, ndc
                         acc = acc + rho1(g,x,s)*fxcw(y,t,x,s,gg)
                      enddo
                   enddo
                   wv1(g,y,t) = acc
                enddo
             enddo
          enddo
          !$omp end parallel do
          call make_aow(nb, nc2, ao, wv1, aow1, aow2_1)
          do t = 1, 2
             call dgemm('N','N', nb, nocc(t), nao_g, 1.0_dp, aow1(1,1,t), nb, CHg, nao_g, 0.0_dp, tmpL, nb)
             call dgemm('N','T', nb, nao_g, nocc(t), 1.0_dp, tmpL, nb, CHg, nao_g, 0.0_dp, Amat, nb)
             call probe_reduce(nb, nc2, ao, Amat, aow2_1(1,1,1,t), PhiD(1,1,t), tacc(1,1,v+1))
          enddo
       enddo
       deallocate(psiL, psiR, rho1, wv1, aow1, aow2_1)
    endif
    if (allocated(taccG)) call gprobe_block(nb, g0, nc2, ao, psiHb)
    deallocate(psiHb, wvR, aowR, aow2R, PhiD, Phi, Amat, tmpL, tmpR)

  end subroutine probe_block

  subroutine copy_block(n, src, dst)
    integer(is), intent(in) :: n
    real(dp), intent(in)    :: src(n)
    real(dp), intent(out)   :: dst(n)
    dst = src
  end subroutine copy_block

!######################################################################
! xc_probe_end: return the accumulated per-AO vectors
!######################################################################
  subroutine xc_probe_end(t)

    real(dp), intent(out) :: t(nao_g,3,nvcur+1)

    if (.not. allocated(tacc)) call mrsf_error('xc_probe_end: no probe evaluation in progress')
    t = tacc
    time_probe = time_probe + wall_time() - tstart
    nprobe_calls = nprobe_calls + 1
    call xc_cleanup()

  end subroutine xc_probe_end

!######################################################################
! make_aow: aow(:,mu,s) = wv0/2 phi + sum_c wv_c d_c phi and
! aow2(:,mu,x,s) = wv0/2 d_x phi + sum_c wv_c d_x d_c phi, both spins
!######################################################################
  subroutine make_aow(nb, nc2, ao, wv, aow, aow2)

    integer(is), intent(in) :: nb, nc2
    real(dp), intent(in)    :: ao(nb,nao_g,nc2), wv(nb,ndc,2)
    real(dp), intent(out)   :: aow(nb,nao_g,2), aow2(nb,nao_g,3,2)

    integer(is) :: mu, s, c, x

    !$omp parallel do private(s,c,x)
    do mu = 1, nao_g
       do s = 1, 2
          aow(:,mu,s) = 0.5_dp*wv(:,1,s)*ao(:,mu,1)
          do c = 2, ncomp
             aow(:,mu,s) = aow(:,mu,s) + wv(:,c,s)*ao(:,mu,c)
          enddo
          do x = 1, 3
             aow2(:,mu,x,s) = 0.5_dp*wv(:,1,s)*ao(:,mu,1+x)
             if (ncomp > 1) then
                do c = 1, 3
                   aow2(:,mu,x,s) = aow2(:,mu,x,s) + wv(:,c+1,s)*ao(:,mu,d2idx(x,c))
                enddo
             endif
          enddo
       enddo
    enddo
    !$omp end parallel do

  end subroutine make_aow

!######################################################################
! probe_reduce: t(mu,x) += sum_g [ d_x phi(g,mu) A(g,mu) + aow2(g,mu,x) Phi(g,mu) ]
!######################################################################
  subroutine probe_reduce(nb, nc2, ao, Amat, aow2, Phi, t)

    integer(is), intent(in)  :: nb, nc2
    real(dp), intent(in)     :: ao(nb,nao_g,nc2), Amat(nb,nao_g), aow2(nb,nao_g,3), Phi(nb,nao_g)
    real(dp), intent(inout)  :: t(nao_g,3)

    integer(is) :: mu, x

    !$omp parallel do private(x)
    do mu = 1, nao_g
       do x = 1, 3
          t(mu,x) = t(mu,x) + dot_product(ao(:,mu,1+x), Amat(:,mu)) + dot_product(aow2(:,mu,x), Phi(:,mu))
       enddo
    enddo
    !$omp end parallel do

  end subroutine probe_reduce

!######################################################################
! xc_closed_shell_init: pre-summed spin-free kernel array of the
! closed-shell mode, fxcw_cs(y,x,g) = w (f_aa + f_ab)(y,x)
!######################################################################
  subroutine xc_closed_shell_init()

    integer(is) :: g, x, y

    if (.not. xc_ready) call mrsf_error('xc_closed_shell_init: grid not initialised')
    if (.not. ref_ready) call mrsf_error('xc_closed_shell_init: reference kernel not installed')
    if (allocated(fxcw_cs)) deallocate(fxcw_cs)
    allocate(fxcw_cs(ndc,ndc,ngrid))
    !$omp parallel do private(x,y)
    do g = 1, ngrid
       do x = 1, ndc
          do y = 1, ndc
             fxcw_cs(y,x,g) = fxcw(y,1,x,1,g) + fxcw(y,1,x,2,g)
          enddo
       enddo
    enddo
    !$omp end parallel do
    cs_ready = .true.

  end subroutine xc_closed_shell_init

!######################################################################
! xc_kernel_cv: closed-shell CV kernel of the extended method on the
! CV amplitudes Y_v of nvec vectors (GGA form),
!   Vcv(j,b,v) = sum_g [ psi_j(g) aow_v(g,nu) + aowC_v(j,g) phi_nu(g) ] C_nu b
! with aow = w [v1 phi/2 + sum_c v_c d_c phi] and aowC the same weighting
! of the cached active-core MO values. Block-parallel: OpenMP over the
! cached grid blocks with sequential MKL inside, per-thread work arrays
! allocated once and per-thread accumulators reduced at the end. Per
! block: psiL = ao L for all vectors (one dgemm per component), then
! per vector the trial density, the kernel potential, aow (first half
! of the stacked right operand Bst), aowC (second half of the stacked
! left operand Ast) and ONE projection dgemm with K = 2 nb:
!   M(j,nu) += [psi_C^T | aowC] [aow ; phi].
!######################################################################
  subroutine xc_kernel_cv(nvec, X, Vcv)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: X(nvirb,ncol,nvec)
    real(dp), intent(out)   :: Vcv(nC,nV,nvec)

    real(dp), allocatable :: Lall(:,:,:), Macc(:,:,:), CCact(:,:), Yact(:,:), Vact(:,:,:)
    real(dp), allocatable :: psiL(:,:,:), psiC(:,:,:), rho1(:,:), wv(:,:), Bst(:,:), Ast(:,:), Mth(:,:,:)
    integer(is)           :: v, ib, nb, nbmax, nk, jj, nca, nthr
    integer(c_int)        :: prev
    real(dp)              :: t0, t1, tl(4)
    logical               :: par

    if (.not. xc_ready) call mrsf_error('xc_kernel_cv: grid not initialised')
    if (.not. ao_cached) call mrsf_error('xc_kernel_cv: AO values not cached')
    if (.not. allocated(fxcw_cs)) call mrsf_error('xc_kernel_cv: closed-shell kernel not initialised')
    t0 = wall_time()

    ! active core columns only (frozen-core holes carry no CV amplitudes)
    nca = nC_act
    nk  = nca*nvec
    allocate(Lall(nao_g,nca,nvec), Macc(nca,nao_g,nvec), CCact(nao_g,nca), Yact(nV,nca))
    do jj = 1, nca
       CCact(:,jj) = CHg(:,cact(jj))
    enddo
    ! factors L_v = C_V Y_v / 2 for all vectors, (nao, nC_act, nvec)
    do v = 1, nvec
       do jj = 1, nca
          Yact(:,jj) = X(3:nvirb,nocca+cact(jj),v)
       enddo
       call dgemm('N','N', nao_g, nca, nV, 0.5_dp, CPg(1,3), nao_g, Yact, nV, &
            0.0_dp, Lall(1,1,v), nao_g)
    enddo
    Macc = 0.0_dp

    nbmax = 0
    do ib = 1, nblocks
       nbmax = max(nbmax, gend(ib) - gbeg(ib) + 1)
    enddo
    nthr = 1
    !$ nthr = omp_get_max_threads()
    par = (nthr > 1) .and. (nblocks >= nthr)
    tl  = 0.0_dp

    !$omp parallel if(par) default(shared) private(ib, nb, prev, psiL, psiC, rho1, wv, Bst, Ast, Mth) &
    !$omp reduction(+:tl)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(1_c_int)
#endif
    allocate(psiL(nbmax,nk,ncomp), psiC(nbmax,nca,ncomp), rho1(nbmax,ndc), wv(nbmax,ndc), &
         Bst(2*nbmax,nao_g), Ast(nca,2*nbmax), Mth(nca,nao_g,nvec))
    Mth = 0.0_dp
    !$omp do schedule(dynamic)
    do ib = 1, nblocks
       nb = gend(ib) - gbeg(ib) + 1
       call cv_block(nb, nbmax, nk, nca, nvec, gbeg(ib), aoc(aoff(ib)+1), psiH(poff(ib)+1), &
            Lall, psiL, psiC, rho1, wv, Bst, Ast, Mth, tl)
    enddo
    !$omp end do
    !$omp critical
    Macc = Macc + Mth
    !$omp end critical
    deallocate(psiL, psiC, rho1, wv, Bst, Ast, Mth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel
    tcs(1:4) = tcs(1:4) + tl(1:4)

    t1 = wall_time()
    allocate(Vact(nca,nV,nvec))
    Vcv = 0.0_dp
    do v = 1, nvec
       call dgemm('N','N', nca, nV, nao_g, 1.0_dp, Macc(1,1,v), nca, CPg(1,3), nao_g, &
            0.0_dp, Vact(1,1,v), nca)
       do jj = 1, nca
          Vcv(cact(jj),:,v) = Vact(jj,:,v)
       enddo
    enddo
    tcs(5) = tcs(5) + wall_time() - t1
    deallocate(Lall, Macc, CCact, Yact, Vact)
    time_xc = time_xc + wall_time() - t0
    nxc_calls = nxc_calls + 1
    nxc_vecs = nxc_vecs + nvec

  end subroutine xc_kernel_cv

!######################################################################
! cv_block: one grid block of the CV kernel for all vectors, run by one
! thread with sequential BLAS. tl(1:4): factor values, density and
! potential, weighted products, projection (thread seconds).
!######################################################################
  subroutine cv_block(nb, nbmax, nk, nca, nvec, g0, ao, psiHb, Lall, psiL, psiC, rho1, wv, &
       Bst, Ast, Mth, tl)

    integer(is), intent(in) :: nb, nbmax, nk, nca, nvec, g0
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp), psiHb(nb,nocca,ncomp)
    real(dp), intent(in)    :: Lall(nao_g,nca,nvec)
    real(dp), intent(inout) :: psiL(nbmax,nk,ncomp), psiC(nbmax,nca,ncomp)
    real(dp), intent(inout) :: rho1(nbmax,ndc), wv(nbmax,ndc)
    real(dp), intent(inout) :: Bst(2*nbmax,nao_g), Ast(nca,2*nbmax), Mth(nca,nao_g,nvec)
    real(dp), intent(inout) :: tl(4)

    integer(is) :: c, jj, mu, g, gg, v, jc, ix, iy
    real(dp)    :: t1, acc

    ! active core MO values of the block and the block-constant halves of
    ! the stacked operands: Ast(:,1:nb) = psi_C^T, Bst(nb+1:2nb,:) = phi
    do c = 1, ncomp
       do jj = 1, nca
          psiC(1:nb,jj,c) = psiHb(1:nb,cact(jj),c)
       enddo
    enddo
    do g = 1, nb
       Ast(1:nca,g) = psiC(g,1:nca,1)
    enddo
    do mu = 1, nao_g
       Bst(nb+1:2*nb,mu) = ao(1:nb,mu,1)
    enddo

    ! factor values psiL = ao L for all vectors and components
    t1 = wall_time()
    do c = 1, ncomp
       call dgemm('N','N', nb, nk, nao_g, 1.0_dp, ao(1,1,c), nb, Lall, nao_g, &
            0.0_dp, psiL(1,1,c), nbmax)
    enddo
    tl(1) = tl(1) + wall_time() - t1

    do v = 1, nvec
       jc = nca*(v-1)
       ! trial density and gradient, rho = 2 psi_L psi_C
       t1 = wall_time()
       rho1(1:nb,:) = 0.0_dp
       do jj = 1, nca
          rho1(1:nb,1) = rho1(1:nb,1) + 2.0_dp*psiL(1:nb,jc+jj,1)*psiC(1:nb,jj,1)
          do c = 2, ncomp
             rho1(1:nb,c) = rho1(1:nb,c) + 2.0_dp*(psiL(1:nb,jc+jj,c)*psiC(1:nb,jj,1) &
                  + psiL(1:nb,jc+jj,1)*psiC(1:nb,jj,c))
          enddo
       enddo
       ! kernel potential wv(g,y) = w sum_x rho(g,x) (f_aa+f_ab)(y,x)
       do g = 1, nb
          gg = g0 + g - 1
          do iy = 1, ndc
             acc = 0.0_dp
             do ix = 1, ndc
                acc = acc + rho1(g,ix)*fxcw_cs(iy,ix,gg)
             enddo
             wv(g,iy) = acc
          enddo
       enddo
       tl(2) = tl(2) + wall_time() - t1
       ! weighted AO products (first half of Bst) and weighted active-core
       ! MO products (second half of Ast)
       t1 = wall_time()
       do mu = 1, nao_g
          Bst(1:nb,mu) = 0.5_dp*wv(1:nb,1)*ao(1:nb,mu,1)
          do c = 2, ncomp
             Bst(1:nb,mu) = Bst(1:nb,mu) + wv(1:nb,c)*ao(1:nb,mu,c)
          enddo
       enddo
       do g = 1, nb
          do jj = 1, nca
             acc = 0.5_dp*wv(g,1)*psiC(g,jj,1)
             do c = 2, ncomp
                acc = acc + wv(g,c)*psiC(g,jj,c)
             enddo
             Ast(jj,nb+g) = acc
          enddo
       enddo
       tl(3) = tl(3) + wall_time() - t1
       ! projection: M(j,nu) += sum_g psi_j(g) aow(g,nu) + sum_g aowC(j,g) phi_nu(g)
       t1 = wall_time()
       call dgemm('N','N', nca, nao_g, 2*nb, 1.0_dp, Ast, nca, Bst, 2*nbmax, &
            1.0_dp, Mth(1,1,v), nca)
       tl(4) = tl(4) + wall_time() - t1
    enddo

  end subroutine cv_block

!######################################################################
! make_aow1: aow(:,mu) = wv0/2 phi + sum_c wv_c d_c phi (one spin)
!######################################################################
  subroutine make_aow1(nb, ao, wv, aow)

    integer(is), intent(in) :: nb
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp), wv(nb,ndc)
    real(dp), intent(out)   :: aow(nb,nao_g)

    integer(is) :: mu, c

    !$omp parallel do private(c)
    do mu = 1, nao_g
       aow(:,mu) = 0.5_dp*wv(:,1)*ao(:,mu,1)
       do c = 2, ncomp
          aow(:,mu) = aow(:,mu) + wv(:,c)*ao(:,mu,c)
       enddo
    enddo
    !$omp end parallel do

  end subroutine make_aow1

!######################################################################
! xc_timings: print the kernel sub-step timers (diagnostics)
!######################################################################
  subroutine xc_timings()

    write(6,'(/,2x,a)') 'XC grid kernel timings (s)'
    write(6,'(2x,a,f10.3)') 'factor values (dgemm)  : ', tpart(1)
    write(6,'(2x,a,f10.3)') 'densities              : ', tpart(2)
    write(6,'(2x,a,f10.3)') 'kernel potentials      : ', tpart(3)
    write(6,'(2x,a,f10.3)') 'weighted AO products   : ', tpart(4)
    write(6,'(2x,a,f10.3)') 'projections (dgemm)    : ', tpart(5)
    write(6,'(2x,a,f10.3)') 'other                  : ', tpart(6)
    write(6,'(2x,a,f10.3,a,i0,a)') 'kernel total           : ', time_xc, ' (', nxc_calls, ' calls)'
    write(6,'(2x,a,f10.3,a,i0,a)') 'probe total            : ', time_probe, ' (', nprobe_calls, ' calls)'
    if (nxcg_calls > 0) write(6,'(2x,a,f10.3,a,i0,a)') 'G-density potentials   : ', time_xcg, ' (', nxcg_calls, ' calls)'
    if (sum(tcs) > 0.0_dp) then
       write(6,'(2x,a)') 'closed-shell CV kernel (s): factor values, densities+potentials, AO products, projections, final'
       write(6,'(2x,5f10.3)') tcs
    endif
    flush(6)

  end subroutine xc_timings

!######################################################################
! xc_g_set: closed-shell kernel set at the G density (PySCF spin-0
! derivatives of the total density; transposed PySCF arrays
! vxc1(ngrid,nv), fxc1(ngrid,nv,nv), kxc1(ngrid,nv,nv,nv)); the grid
! weights are folded into the kernels. have_kxc = .false. (triplets,
! HF) leaves the third derivatives out.
!######################################################################
  subroutine xc_g_set(ngrid1, nv1, vxc1, fxc1, kxc1, have_kxc)

    integer(is), intent(in) :: ngrid1, nv1
    real(dp), intent(in)    :: vxc1(ngrid1,nv1), fxc1(ngrid1,nv1,nv1), kxc1(*)
    logical, intent(in)     :: have_kxc

    integer(is) :: g, x, y

    if (.not. xc_ready) call mrsf_error('xc_g_set: grid not initialised')
    if (ngrid1 /= ngrid .or. nv1 /= ndc) call mrsf_error('xc_g_set: dimension mismatch')
    call xc_g_free()
    allocate(vxcG(ngrid,ndc), fxcwG(ndc,ndc,ngrid))
    vxcG = vxc1
    !$omp parallel do private(x,y)
    do g = 1, ngrid
       do x = 1, ndc
          do y = 1, ndc
             fxcwG(y,x,g) = wgt(g)*fxc1(g,y,x)
          enddo
       enddo
    enddo
    !$omp end parallel do
    g_have_kxc = have_kxc
    if (have_kxc) then
       allocate(kxcwG(ndc,ndc,ndc,ngrid))
       call fold_kxc(ngrid, ndc, kxc1, wgt, kxcwG)
    endif
    g_ready = .true.

  end subroutine xc_g_set

  subroutine fold_kxc(ng1, nv1, kxc1, w, out)

    integer(is), intent(in) :: ng1, nv1
    real(dp), intent(in)    :: kxc1(ng1,nv1,nv1,nv1), w(ng1)
    real(dp), intent(out)   :: out(nv1,nv1,nv1,ng1)

    integer(is) :: g, x, y, z

    !$omp parallel do private(x,y,z)
    do g = 1, ng1
       do x = 1, nv1
          do y = 1, nv1
             do z = 1, nv1
                out(z,y,x,g) = w(g)*kxc1(g,z,y,x)
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

  end subroutine fold_kxc

  subroutine xc_g_free()

    if (allocated(vxcG)) deallocate(vxcG)
    if (allocated(fxcwG)) deallocate(fxcwG)
    if (allocated(kxcwG)) deallocate(kxcwG)
    g_ready = .false.; g_have_kxc = .false.

  end subroutine xc_g_free

!######################################################################
! g_wv_f: weighted kernel potential wv(g,y) = sum_x w f(y,x) rho(g,x)
! g_wv_k: weighted third-derivative potential
!         wv(g,z) = sum_{x,y} w k(z,y,x) rho(g,x) rho(g,y)
! (one spin; rho are the components of the trial density on the block)
!######################################################################
  subroutine g_wv_f(nb, g0, rho, wv)

    integer(is), intent(in) :: nb, g0
    real(dp), intent(in)    :: rho(nb,ndc)
    real(dp), intent(out)   :: wv(nb,ndc)

    integer(is) :: g, gg, x, y
    real(dp)    :: acc

    do g = 1, nb
       gg = g0 + g - 1
       do y = 1, ndc
          acc = 0.0_dp
          do x = 1, ndc
             acc = acc + fxcwG(y,x,gg)*rho(g,x)
          enddo
          wv(g,y) = acc
       enddo
    enddo

  end subroutine g_wv_f

  subroutine g_wv_k(nb, g0, rho, wv)

    integer(is), intent(in) :: nb, g0
    real(dp), intent(in)    :: rho(nb,ndc)
    real(dp), intent(out)   :: wv(nb,ndc)

    integer(is) :: g, gg, x, y, z
    real(dp)    :: acc

    if (.not. g_have_kxc) call mrsf_error('g_wv_k: third derivatives not installed')
    do g = 1, nb
       gg = g0 + g - 1
       do z = 1, ndc
          acc = 0.0_dp
          do x = 1, ndc
             do y = 1, ndc
                acc = acc + kxcwG(z,y,x,gg)*rho(g,x)*rho(g,y)
             enddo
          enddo
          wv(g,z) = acc
       enddo
    enddo

  end subroutine g_wv_k

!######################################################################
! g_project: M(i,nu) += sum_g [ psi_i(g) aow(g,nu) + aowH_i(g) phi_nu(g) ]
! for a weighted potential wv (one spin): aow into the first half of the
! stacked right operand Bst (second half = phi, block constant), aowH
! into the second half of the stacked left operand Ast (first half =
! psi_H^T, block constant), one dgemm with K = 2 nb
!######################################################################
  subroutine g_project(nb, ao, psiHb, wv, Bst, Ast, Mth)

    integer(is), intent(in) :: nb
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp), psiHb(nb,nocca,ncomp), wv(nb,ndc)
    real(dp), intent(inout) :: Bst(2*nb,nao_g), Ast(nocca,2*nb), Mth(nocca,nao_g)

    integer(is) :: mu, c, g, i
    real(dp)    :: acc

    do mu = 1, nao_g
       Bst(1:nb,mu) = 0.5_dp*wv(1:nb,1)*ao(1:nb,mu,1)
       do c = 2, ncomp
          Bst(1:nb,mu) = Bst(1:nb,mu) + wv(1:nb,c)*ao(1:nb,mu,c)
       enddo
    enddo
    do g = 1, nb
       do i = 1, nocca
          acc = 0.5_dp*wv(g,1)*psiHb(g,i,1)
          do c = 2, ncomp
             acc = acc + wv(g,c)*psiHb(g,i,c)
          enddo
          Ast(i,nb+g) = acc
       enddo
    enddo
    call dgemm('N','N', nocca, nao_g, 2*nb, 1.0_dp, Ast, nocca, Bst, 2*nb, 1.0_dp, Mth, nocca)

  end subroutine g_project

!######################################################################
! xc_g_potential: kernel potentials at the G density for the EMRSF
! gradient in one block-parallel pass over the cached AO blocks
! (Phase-5 pattern: OpenMP over blocks, sequential MKL inside,
! per-thread buffers and accumulators reduced once):
!   MM(i,nu) = sum_g [psi_i aowM(nu) + aowH_i phi_nu],  wv_M = w f rho(M)
!            (= the projection of v[f rho(M)] on the hole MOs: V_M(H, all) = MM C)
!   MK(i,nu) = the same for wv_K = w k rho(D_s)^2      (v[k rho(D_s)^2] on H x all)
!   V1(mu,nu) = sum_g [aowS(mu) phi_nu + phi_mu aowS(nu)],  wv_S = w f rho(D_s)
!            (v[f rho(D_s)] as a full AO matrix)
! with M = LM RM^T + RM LM^T (kM columns) and D_s = LS C_C^T + C_C LS^T
! (kS = nocca - 2 columns; C_C = the doubly occupied MOs, cached on the
! grid as the first columns of psiH).
!######################################################################
  subroutine xc_g_potential(kM, LM, RM, kS, LS, MM, MK, V1)

    integer(is), intent(in) :: kM, kS
    real(dp), intent(in)    :: LM(nao_g,kM), RM(nao_g,kM), LS(nao_g,kS)
    real(dp), intent(out)   :: MM(nocca,nao_g), MK(nocca,nao_g), V1(nao_g,nao_g)

    integer(is)           :: ib, nb, nbmax, nthr
    integer(c_int)        :: prev
    logical               :: par
    real(dp)              :: t0
    real(dp), allocatable :: psiLM(:), psiRM(:), psiLS(:), rhoM(:), rhoS(:), wv(:), Bst(:), Ast(:), aowS(:)
    real(dp), allocatable :: MthM(:,:), MthK(:,:), Vth(:,:)

    if (.not. g_ready) call mrsf_error('xc_g_potential: G kernel set not installed')
    if (.not. ao_cached) call mrsf_error('xc_g_potential: the AO values must be cached (raise mem_budget)')
    if (kS /= nocca - 2) call mrsf_error('xc_g_potential: kS must equal the number of doubly occupied MOs')
    t0 = wall_time()
    nbmax = 0
    do ib = 1, nblocks
       nbmax = max(nbmax, gend(ib) - gbeg(ib) + 1)
    enddo
    nthr = 1
    !$ nthr = omp_get_max_threads()
    par = (nthr > 1) .and. (nblocks >= nthr)
    MM = 0.0_dp; MK = 0.0_dp; V1 = 0.0_dp

    !$omp parallel if(par) default(shared) private(ib, nb, prev, psiLM, psiRM, psiLS, rhoM, rhoS, wv, &
    !$omp Bst, Ast, aowS, MthM, MthK, Vth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(1_c_int)
#endif
    allocate(psiLM(nbmax*kM*ncomp), psiRM(nbmax*kM*ncomp), psiLS(nbmax*kS*ncomp), rhoM(nbmax*ndc), &
         rhoS(nbmax*ndc), wv(nbmax*ndc), Bst(2*nbmax*nao_g), Ast(nocca*2*nbmax), aowS(nbmax*nao_g), &
         MthM(nocca,nao_g), MthK(nocca,nao_g), Vth(nao_g,nao_g))
    MthM = 0.0_dp; MthK = 0.0_dp; Vth = 0.0_dp
    !$omp do schedule(dynamic)
    do ib = 1, nblocks
       nb = gend(ib) - gbeg(ib) + 1
       call gpot_block(nb, gbeg(ib), aoc(aoff(ib)+1), psiH(poff(ib)+1), kM, LM, RM, kS, LS, &
            psiLM, psiRM, psiLS, rhoM, rhoS, wv, Bst, Ast, aowS, MthM, MthK, Vth)
    enddo
    !$omp end do
    !$omp critical
    MM = MM + MthM
    MK = MK + MthK
    V1 = V1 + Vth
    !$omp end critical
    deallocate(psiLM, psiRM, psiLS, rhoM, rhoS, wv, Bst, Ast, aowS, MthM, MthK, Vth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel
    V1 = V1 + transpose(V1)
    time_xcg = time_xcg + wall_time() - t0
    nxcg_calls = nxcg_calls + 1

  end subroutine xc_g_potential

  subroutine gpot_block(nb, g0, ao, psiHb, kM, LM, RM, kS, LS, psiLM, psiRM, psiLS, rhoM, rhoS, &
       wv, Bst, Ast, aowS, MthM, MthK, Vth)

    integer(is), intent(in) :: nb, g0, kM, kS
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp), psiHb(nb,nocca,ncomp)
    real(dp), intent(in)    :: LM(nao_g,kM), RM(nao_g,kM), LS(nao_g,kS)
    real(dp), intent(inout) :: psiLM(nb,kM,ncomp), psiRM(nb,kM,ncomp), psiLS(nb,kS,ncomp)
    real(dp), intent(inout) :: rhoM(nb,ndc), rhoS(nb,ndc), wv(nb,ndc)
    real(dp), intent(inout) :: Bst(2*nb,nao_g), Ast(nocca,2*nb), aowS(nb,nao_g)
    real(dp), intent(inout) :: MthM(nocca,nao_g), MthK(nocca,nao_g), Vth(nao_g,nao_g)

    integer(is) :: c, g, mu

    ! block-constant halves of the stacked operands
    do g = 1, nb
       Ast(1:nocca,g) = psiHb(g,1:nocca,1)
    enddo
    do mu = 1, nao_g
       Bst(nb+1:2*nb,mu) = ao(1:nb,mu,1)
    enddo
    ! factor values on the block
    do c = 1, ncomp
       call dgemm('N','N', nb, kM, nao_g, 1.0_dp, ao(1,1,c), nb, LM, nao_g, 0.0_dp, psiLM(1,1,c), nb)
       call dgemm('N','N', nb, kM, nao_g, 1.0_dp, ao(1,1,c), nb, RM, nao_g, 0.0_dp, psiRM(1,1,c), nb)
       call dgemm('N','N', nb, kS, nao_g, 1.0_dp, ao(1,1,c), nb, LS, nao_g, 0.0_dp, psiLS(1,1,c), nb)
    enddo
    call block_rho(nb, kM, psiLM, nb*kM, psiRM, nb*kM, rhoM)
    call block_rho(nb, kS, psiLS, nb*kS, psiHb, nb*nocca, rhoS)
    ! v[f rho(M)] on the hole MOs
    call g_wv_f(nb, g0, rhoM, wv)
    call g_project(nb, ao, psiHb, wv, Bst, Ast, MthM)
    ! v[k rho(D_s)^2] on the hole MOs
    if (g_have_kxc) then
       call g_wv_k(nb, g0, rhoS, wv)
       call g_project(nb, ao, psiHb, wv, Bst, Ast, MthK)
    endif
    ! v[f rho(D_s)] as a full AO matrix (symmetrised by the caller)
    call g_wv_f(nb, g0, rhoS, wv)
    do mu = 1, nao_g
       aowS(1:nb,mu) = 0.5_dp*wv(1:nb,1)*ao(1:nb,mu,1)
       do c = 2, ncomp
          aowS(1:nb,mu) = aowS(1:nb,mu) + wv(1:nb,c)*ao(1:nb,mu,c)
       enddo
    enddo
    call dgemm('T','N', nao_g, nao_g, nb, 1.0_dp, aowS, nb, ao(1,1,1), nb, 1.0_dp, Vth, nao_g)

  end subroutine gpot_block

!######################################################################
! xc_gprobe_set: register, after xc_probe_begin, the EMRSF densities of
! the probe pass: per state M = LM RM^T + RM LM^T (kM columns) and
! D_s = LS C_C^T + C_C LS^T (kS = nocca - 2). The probe then also
! accumulates, per state,
!   taccG(:,:,v) = (vxc_G, M) + (f_G[M], D_G) + (f_G[D_s], D_s) + (k_G[D_s, D_s]/2, D_G)
! with D_G = 2 C_G C_G^T (G = the first nocca - 1 hole MOs), in the same
! per-AO form as the reference probe (dE/dR_{A,x} = -2 sum_{mu in A} t).
!######################################################################
  subroutine xc_gprobe_set(nst, kM, LM, RM, kS, LS, use_kxc)

    integer(is), intent(in) :: nst, kM, kS
    logical, intent(in)     :: use_kxc
    real(dp), intent(in)    :: LM(nao_g,kM,nst), RM(nao_g,kM,nst), LS(nao_g,kS,nst)

    if (.not. g_ready) call mrsf_error('xc_gprobe_set: G kernel set not installed')
    if (.not. allocated(tacc)) call mrsf_error('xc_gprobe_set: call xc_probe_begin first')
    if (nst /= nvcur) call mrsf_error('xc_gprobe_set: state count differs from xc_probe_begin')
    if (kS /= nocca - 2) call mrsf_error('xc_gprobe_set: kS must equal the number of doubly occupied MOs')
    if (use_kxc .and. .not. g_have_kxc) call mrsf_error('xc_gprobe_set: third functional derivatives not installed')
    gq_kxc = use_kxc
    if (allocated(LMg)) deallocate(LMg, RMg, LSg)
    if (allocated(taccG)) deallocate(taccG)
    nstg = nst; kMg = kM; kSg = kS
    allocate(LMg(nao_g,kM,nst), RMg(nao_g,kM,nst), LSg(nao_g,kS,nst), taccG(nao_g,3,nst))
    LMg = LM; RMg = RM; LSg = LS
    taccG = 0.0_dp

  end subroutine xc_gprobe_set

!######################################################################
! xc_gprobe_get: the accumulated G-probe vectors (call before xc_probe_end)
!######################################################################
  subroutine xc_gprobe_get(tG)

    real(dp), intent(out) :: tG(nao_g,3,nstg)

    if (.not. allocated(taccG)) call mrsf_error('xc_gprobe_get: no G probe in progress')
    tG = taccG

  end subroutine xc_gprobe_get

  subroutine gprobe_block(nb, g0, nc2, ao, psiHb)

    integer(is), intent(in) :: nb, g0, nc2
    real(dp), intent(in)    :: ao(nb,nao_g,nc2), psiHb(nb,nocca,ncomp)

    integer(is) :: v, c, nGo
    real(dp), allocatable :: wvG(:,:), aowG(:,:), aow2G(:,:,:), PhiG(:,:), Phi(:,:), Amat(:,:), tmp(:,:), tmp2(:,:)
    real(dp), allocatable :: psiLM(:,:,:), psiRM(:,:,:), psiLS(:,:,:), rho(:,:), wv1(:,:), aow1(:,:), aow2_1(:,:,:)

    nGo = nocca - 1
    allocate(wvG(nb,ndc), aowG(nb,nao_g), aow2G(nb,nao_g,3))
    do c = 1, ndc
       wvG(:,c) = wgt(g0:g0+nb-1)*vxcG(g0:g0+nb-1,c)
    enddo
    call make_aow_1s(nb, nc2, ao, wvG, aowG, aow2G)
    allocate(PhiG(nb,nao_g), Phi(nb,nao_g), Amat(nb,nao_g))
    allocate(tmp(nb,max(kMg,nocca)), tmp2(nb,max(kMg,nocca)))
    ! Phi_G = psi_G C_G^T
    call dgemm('N','T', nb, nao_g, nGo, 1.0_dp, psiHb, nb, CHg, nao_g, 0.0_dp, PhiG, nb)
    allocate(psiLM(nb,kMg,ncomp), psiRM(nb,kMg,ncomp), psiLS(nb,kSg,ncomp), rho(nb,ndc), wv1(nb,ndc), &
         aow1(nb,nao_g), aow2_1(nb,nao_g,3))
    do v = 1, nstg
       do c = 1, ncomp
          call dgemm('N','N', nb, kMg, nao_g, 1.0_dp, ao(1,1,c), nb, LMg(1,1,v), nao_g, 0.0_dp, psiLM(1,1,c), nb)
          call dgemm('N','N', nb, kMg, nao_g, 1.0_dp, ao(1,1,c), nb, RMg(1,1,v), nao_g, 0.0_dp, psiRM(1,1,c), nb)
          call dgemm('N','N', nb, kSg, nao_g, 1.0_dp, ao(1,1,c), nb, LSg(1,1,v), nao_g, 0.0_dp, psiLS(1,1,c), nb)
       enddo
       ! (vxc_G, M): Phi = psiRM LM^T + psiLM RM^T, A = (aowG RM) LM^T + (aowG LM) RM^T
       call dgemm('N','T', nb, nao_g, kMg, 1.0_dp, psiRM, nb, LMg(1,1,v), nao_g, 0.0_dp, Phi, nb)
       call dgemm('N','T', nb, nao_g, kMg, 1.0_dp, psiLM, nb, RMg(1,1,v), nao_g, 1.0_dp, Phi, nb)
       call dgemm('N','N', nb, kMg, nao_g, 1.0_dp, aowG, nb, RMg(1,1,v), nao_g, 0.0_dp, tmp, nb)
       call dgemm('N','N', nb, kMg, nao_g, 1.0_dp, aowG, nb, LMg(1,1,v), nao_g, 0.0_dp, tmp2, nb)
       call dgemm('N','T', nb, nao_g, kMg, 1.0_dp, tmp, nb, LMg(1,1,v), nao_g, 0.0_dp, Amat, nb)
       call dgemm('N','T', nb, nao_g, kMg, 1.0_dp, tmp2, nb, RMg(1,1,v), nao_g, 1.0_dp, Amat, nb)
       call probe_reduce(nb, nc2, ao, Amat, aow2G, Phi, taccG(1,1,v))
       ! (f_G[M], D_G): D_G = 2 C_G C_G^T, Phi_G carries C_G C_G^T (factor 2)
       call block_rho(nb, kMg, psiLM, nb*kMg, psiRM, nb*kMg, rho)
       call g_wv_f(nb, g0, rho, wv1)
       wv1 = 2.0_dp*wv1
       call make_aow_1s(nb, nc2, ao, wv1, aow1, aow2_1)
       call dgemm('N','N', nb, nGo, nao_g, 1.0_dp, aow1, nb, CHg, nao_g, 0.0_dp, tmp, nb)
       call dgemm('N','T', nb, nao_g, nGo, 1.0_dp, tmp, nb, CHg, nao_g, 0.0_dp, Amat, nb)
       call probe_reduce(nb, nc2, ao, Amat, aow2_1, PhiG, taccG(1,1,v))
       ! kernel terms of the singlet CV block only
       if (.not. gq_kxc) cycle
       ! (f_G[D_s], D_s): Phi = psi_C LS^T + psiLS C_C^T, A = (aow C_C) LS^T + (aow LS) C_C^T
       call block_rho(nb, kSg, psiLS, nb*kSg, psiHb, nb*nocca, rho)
       call g_wv_f(nb, g0, rho, wv1)
       call make_aow_1s(nb, nc2, ao, wv1, aow1, aow2_1)
       call dgemm('N','T', nb, nao_g, kSg, 1.0_dp, psiHb, nb, LSg(1,1,v), nao_g, 0.0_dp, Phi, nb)
       call dgemm('N','T', nb, nao_g, kSg, 1.0_dp, psiLS, nb, CHg, nao_g, 1.0_dp, Phi, nb)
       call dgemm('N','N', nb, kSg, nao_g, 1.0_dp, aow1, nb, CHg, nao_g, 0.0_dp, tmp, nb)
       call dgemm('N','N', nb, kSg, nao_g, 1.0_dp, aow1, nb, LSg(1,1,v), nao_g, 0.0_dp, tmp2, nb)
       call dgemm('N','T', nb, nao_g, kSg, 1.0_dp, tmp, nb, LSg(1,1,v), nao_g, 0.0_dp, Amat, nb)
       call dgemm('N','T', nb, nao_g, kSg, 1.0_dp, tmp2, nb, CHg, nao_g, 1.0_dp, Amat, nb)
       call probe_reduce(nb, nc2, ao, Amat, aow2_1, Phi, taccG(1,1,v))
       ! (k_G[D_s, D_s]/2, D_G) = (k_G[D_s, D_s], C_G C_G^T)
       call g_wv_k(nb, g0, rho, wv1)
       call make_aow_1s(nb, nc2, ao, wv1, aow1, aow2_1)
       call dgemm('N','N', nb, nGo, nao_g, 1.0_dp, aow1, nb, CHg, nao_g, 0.0_dp, tmp, nb)
       call dgemm('N','T', nb, nao_g, nGo, 1.0_dp, tmp, nb, CHg, nao_g, 0.0_dp, Amat, nb)
       call probe_reduce(nb, nc2, ao, Amat, aow2_1, PhiG, taccG(1,1,v))
    enddo
    deallocate(wvG, aowG, aow2G, PhiG, Phi, Amat, tmp, tmp2, psiLM, psiRM, psiLS, rho, wv1, aow1, aow2_1)

  end subroutine gprobe_block

!######################################################################
! make_aow_1s: one-spin version of make_aow
!######################################################################
  subroutine make_aow_1s(nb, nc2, ao, wv, aow, aow2)

    integer(is), intent(in) :: nb, nc2
    real(dp), intent(in)    :: ao(nb,nao_g,nc2), wv(nb,ndc)
    real(dp), intent(out)   :: aow(nb,nao_g), aow2(nb,nao_g,3)

    integer(is) :: mu, c, x

    !$omp parallel do private(c,x)
    do mu = 1, nao_g
       aow(:,mu) = 0.5_dp*wv(:,1)*ao(:,mu,1)
       do c = 2, ncomp
          aow(:,mu) = aow(:,mu) + wv(:,c)*ao(:,mu,c)
       enddo
       do x = 1, 3
          aow2(:,mu,x) = 0.5_dp*wv(:,1)*ao(:,mu,1+x)
          if (ncomp > 1) then
             do c = 1, 3
                aow2(:,mu,x) = aow2(:,mu,x) + wv(:,c+1)*ao(:,mu,d2idx(x,c))
             enddo
          endif
       enddo
    enddo
    !$omp end parallel do

  end subroutine make_aow_1s

!######################################################################
! xc_g_density: density (value and gradient) of the closed-shell
! configuration G = the first nocca-1 hole MOs doubly occupied, from the
! cached hole-MO values: rho(g,1) = 2 sum_i psi_i^2, rho(g,1+x) = 4 sum_i psi_i d_x psi_i
!######################################################################
  subroutine xc_g_density(rho)

    real(dp), intent(out) :: rho(ngrid,ndc)

    integer(is) :: ib, nb

    if (.not. xc_ready) call mrsf_error('xc_g_density: grid not initialised')
    if (.not. ao_cached) call mrsf_error('xc_g_density: the AO values must be cached (raise mem_budget)')
    !$omp parallel do private(nb)
    do ib = 1, nblocks
       nb = gend(ib) - gbeg(ib) + 1
       call gdens_block(nb, gbeg(ib), nocca-1, psiH(poff(ib)+1), rho)
    enddo
    !$omp end parallel do

  end subroutine xc_g_density

  subroutine gdens_block(nb, g0, nGo, psib, rho)

    integer(is), intent(in)  :: nb, g0, nGo
    real(dp), intent(in)     :: psib(nb,nocca,ncomp)
    real(dp), intent(inout)  :: rho(ngrid,ndc)

    integer(is) :: g, x

    do g = 1, nb
       rho(g0+g-1,1) = 2.0_dp*sum(psib(g,1:nGo,1)**2)
       do x = 2, ndc
          rho(g0+g-1,x) = 4.0_dp*sum(psib(g,1:nGo,1)*psib(g,1:nGo,x))
       enddo
    enddo

  end subroutine gdens_block

!######################################################################
! xc_g_vxc_ao: AO matrix of the XC potential of the G density,
! V = T + T^T, T = sum_g aow_G(g) phi(g)^T (PySCF nr_rks convention),
! from the cached AO values and the installed vxc_G; block-parallel
!######################################################################
  subroutine xc_g_vxc_ao(V)

    real(dp), intent(out) :: V(nao_g,nao_g)

    integer(is)           :: ib, nb, nbmax, nthr
    integer(c_int)        :: prev
    logical               :: par
    real(dp)              :: t0
    real(dp), allocatable :: wv(:), aow(:), Vth(:,:)

    if (.not. g_ready) call mrsf_error('xc_g_vxc_ao: G kernel set not installed')
    if (.not. ao_cached) call mrsf_error('xc_g_vxc_ao: the AO values must be cached (raise mem_budget)')
    t0 = wall_time()
    nbmax = 0
    do ib = 1, nblocks
       nbmax = max(nbmax, gend(ib) - gbeg(ib) + 1)
    enddo
    nthr = 1
    !$ nthr = omp_get_max_threads()
    par = (nthr > 1) .and. (nblocks >= nthr)
    V = 0.0_dp

    !$omp parallel if(par) default(shared) private(ib, nb, prev, wv, aow, Vth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(1_c_int)
#endif
    allocate(wv(nbmax*ndc), aow(nbmax*nao_g), Vth(nao_g,nao_g))
    Vth = 0.0_dp
    !$omp do schedule(dynamic)
    do ib = 1, nblocks
       nb = gend(ib) - gbeg(ib) + 1
       call gvxc_block(nb, gbeg(ib), aoc(aoff(ib)+1), wv, aow, Vth)
    enddo
    !$omp end do
    !$omp critical
    V = V + Vth
    !$omp end critical
    deallocate(wv, aow, Vth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel
    V = V + transpose(V)
    time_xcg = time_xcg + wall_time() - t0
    nxcg_calls = nxcg_calls + 1

  end subroutine xc_g_vxc_ao

  subroutine gvxc_block(nb, g0, ao, wv, aow, Vth)

    integer(is), intent(in) :: nb, g0
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp)
    real(dp), intent(inout) :: wv(nb,ndc), aow(nb,nao_g), Vth(nao_g,nao_g)

    integer(is) :: c, mu

    do c = 1, ndc
       wv(:,c) = wgt(g0:g0+nb-1)*vxcG(g0:g0+nb-1,c)
    enddo
    do mu = 1, nao_g
       aow(:,mu) = 0.5_dp*wv(:,1)*ao(:,mu,1)
       do c = 2, ncomp
          aow(:,mu) = aow(:,mu) + wv(:,c)*ao(:,mu,c)
       enddo
    enddo
    call dgemm('T','N', nao_g, nao_g, nb, 1.0_dp, aow, nb, ao(1,1,1), nb, 1.0_dp, Vth, nao_g)

  end subroutine gvxc_block

!######################################################################
! xc_cs_from_g: the closed-shell kernel array of the extended method's
! CV block from the installed G set, fxcw_cs = fac * fxcwG
!######################################################################
  subroutine xc_cs_from_g(fac)

    real(dp), intent(in) :: fac

    integer(is) :: g

    if (.not. g_ready) call mrsf_error('xc_cs_from_g: G kernel set not installed')
    if (allocated(fxcw_cs)) deallocate(fxcw_cs)
    allocate(fxcw_cs(ndc,ndc,ngrid))
    !$omp parallel do
    do g = 1, ngrid
       fxcw_cs(:,:,g) = fac*fxcwG(:,:,g)
    enddo
    !$omp end parallel do
    cs_ready = .true.

  end subroutine xc_cs_from_g

!######################################################################
! Quadratic-form probe (Part XIII): every probe density is P = B S B^T in
! the basis B = [C_H (hole MOs, cached on the grid) | Bexp (explicit AO
! columns)]: per state the alpha and beta relaxed densities (kinds 1, 2)
! and, for the extended method, the F^G density M (kind 3) and the CV
! trial density D_s (kind 4). The AO values with second derivatives are
! streamed in superblocks of consecutive cache blocks (xc_probe_super),
! processed in parallel over the cache blocks with sequential BLAS.
! Per term: rho = sum_ab psi_a S_ab psi_b, Phi = psi_B W^T, A = (aow B) W^T,
! W = B S; the kernel terms (f[P], D) use the reference densities as today.
!######################################################################
  subroutine xc_probe_begin_q(nst, nexp, Bexp, maxnb, nb, idx, S, with_g, with_kxc)

    integer(is), intent(in) :: nst, nexp, maxnb
    real(dp), intent(in)    :: Bexp(nao_g,max(nexp,1))
    integer(is), intent(in) :: nb(4,nst), idx(maxnb,4,nst)
    real(dp), intent(in)    :: S(maxnb,maxnb,4,nst)
    logical, intent(in)     :: with_g, with_kxc

    integer(is) :: v, t, a, j, nbt

    if (.not. xc_ready) call mrsf_error('xc_probe_begin_q: grid not initialised')
    if (.not. ref_ready) call mrsf_error('xc_probe_begin_q: reference kernel not installed')
    if (.not. ao_cached) call mrsf_error('xc_probe_begin_q: the AO values must be cached (raise mem_budget)')
    if (with_g .and. .not. g_ready) call mrsf_error('xc_probe_begin_q: G kernel set not installed')
    if (with_g .and. with_kxc .and. .not. g_have_kxc) &
         call mrsf_error('xc_probe_begin_q: third functional derivatives not installed')
    call xc_cleanup()
    nst_q = nst; nexp_q = nexp; maxnb_q = maxnb; gq = with_g; gq_kxc = (with_g .and. with_kxc)
    nvcur = nst; nstg = merge(nst, 0_is, with_g); kcur = 0
    allocate(Bexp_q(nao_g,max(nexp,1)), nbq(4,nst), idxq(maxnb,4,nst), Sq(maxnb,maxnb,4,nst))
    allocate(Bq(nao_g,maxnb,4,nst), Wq(nao_g,maxnb,4,nst))
    Bexp_q = Bexp; nbq = nb; idxq = idx; Sq = S
    Bq = 0.0_dp; Wq = 0.0_dp
    do v = 1, nst
       do t = 1, 4
          nbt = nb(t,v)
          if (nbt > maxnb) call mrsf_error('xc_probe_begin_q: term too long')
          do a = 1, nbt
             j = idx(a,t,v)
             if (j < 1 .or. j > nocca + nexp) call mrsf_error('xc_probe_begin_q: bad basis index')
             if (j <= nocca) then
                Bq(:,a,t,v) = CHg(:,j)
             else
                Bq(:,a,t,v) = Bexp(:,j-nocca)
             endif
          enddo
          if (nbt > 0) call dgemm('N','N', nao_g, nbt, nbt, 1.0_dp, Bq(1,1,t,v), nao_g, S(1,1,t,v), maxnb, &
               0.0_dp, Wq(1,1,t,v), nao_g)
       enddo
    enddo
    allocate(tacc(nao_g,3,nst+1)); tacc = 0.0_dp
    if (with_g) then
       allocate(taccG(nao_g,3,nst)); taccG = 0.0_dp
    endif
    probe_q = .true.
    tstart = wall_time()

  end subroutine xc_probe_begin_q

!######################################################################
! xc_probe_super: the AO values (nc2 components) of the ng points of the
! cache blocks ib0 .. ib0+nsub-1, processed in parallel over the blocks
!######################################################################
  subroutine xc_probe_super(ib0, nsub, ng, nc2, ao)

    integer(is), intent(in) :: ib0, nsub, ng, nc2
    real(dp), intent(in)    :: ao(ng,nao_g,nc2)

    integer(is)           :: ib, nb, nbmax, nthr, off, nst
    integer(c_int)        :: prev
    logical               :: par
    real(dp), allocatable :: aob(:), psiHb(:), psiE(:), psiB(:), Tq(:), rho(:), wv(:), aowR(:), aow2R(:), &
                             aow1(:), aow21(:), PhiD(:), Phi(:), Amat(:), tmp(:), wvG(:), aowG(:), aow2G(:), &
                             PhiG(:), tth(:,:,:), tGth(:,:,:)

    if (.not. probe_q) call mrsf_error('xc_probe_super: no quadratic-form probe in progress')
    if (ib0 < 1 .or. ib0 + nsub - 1 > nblocks) call mrsf_error('xc_probe_super: block range out of range')
    if (gend(ib0+nsub-1) - gbeg(ib0) + 1 /= ng) call mrsf_error('xc_probe_super: point count mismatch')
    if ((ncomp == 1 .and. nc2 < 4) .or. (ncomp == 4 .and. nc2 < 10)) &
         call mrsf_error('xc_probe_super: AO derivative components missing')
    nst = nst_q
    nbmax = 0
    do ib = ib0, ib0+nsub-1
       nbmax = max(nbmax, gend(ib) - gbeg(ib) + 1)
    enddo
    nthr = 1
    !$ nthr = omp_get_max_threads()
    par = (nthr > 1) .and. (nsub > 1)

    !$omp parallel if(par) default(shared) private(ib, nb, off, prev, aob, psiHb, psiE, psiB, Tq, rho, wv, &
    !$omp aowR, aow2R, aow1, aow21, PhiD, Phi, Amat, tmp, wvG, aowG, aow2G, PhiG, tth, tGth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(1_c_int)
#endif
    allocate(aob(nbmax*nao_g*nc2), psiHb(nbmax*nocca*ncomp), psiE(nbmax*max(nexp_q,1)*ncomp), &
         psiB(nbmax*maxnb_q*ncomp), Tq(nbmax*maxnb_q), rho(nbmax*ndc*2), wv(nbmax*ndc*2), &
         aowR(nbmax*nao_g*2), aow2R(nbmax*nao_g*3*2), aow1(nbmax*nao_g*2), aow21(nbmax*nao_g*3*2), &
         PhiD(nbmax*nao_g*2), Phi(nbmax*nao_g), Amat(nbmax*nao_g), tmp(nbmax*max(maxnb_q,nocca)), &
         wvG(nbmax*ndc), aowG(nbmax*nao_g), aow2G(nbmax*nao_g*3), PhiG(nbmax*nao_g), &
         tth(nao_g,3,nst+1), tGth(nao_g,3,max(nst,1)))
    tth = 0.0_dp; tGth = 0.0_dp
    !$omp do schedule(dynamic)
    do ib = ib0, ib0+nsub-1
       nb  = gend(ib) - gbeg(ib) + 1
       off = gbeg(ib) - gbeg(ib0)
       call copy_sub(ng, nb, off, nao_g*nc2, ao, aob)
       call copy_block(nb*nocca*ncomp, psiH(poff(ib)+1), psiHb)
       call probe_q_block(nb, gbeg(ib), nc2, aob, psiHb, psiE, psiB, Tq, rho, wv, aowR, aow2R, aow1, aow21, &
            PhiD, Phi, Amat, tmp, wvG, aowG, aow2G, PhiG, tth, tGth)
    enddo
    !$omp end do
    !$omp critical
    tacc = tacc + tth
    if (gq) taccG = taccG + tGth
    !$omp end critical
    deallocate(aob, psiHb, psiE, psiB, Tq, rho, wv, aowR, aow2R, aow1, aow21, PhiD, Phi, Amat, tmp, &
         wvG, aowG, aow2G, PhiG, tth, tGth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel

  end subroutine xc_probe_super

  subroutine copy_sub(ng, nb, off, ncol, ao, aob)
    integer(is), intent(in) :: ng, nb, off, ncol
    real(dp), intent(in)    :: ao(ng,ncol)
    real(dp), intent(out)   :: aob(nb,ncol)
    integer(is) :: j
    do j = 1, ncol
       aob(:,j) = ao(off+1:off+nb,j)
    enddo
  end subroutine copy_sub

  subroutine probe_q_block(nb, g0, nc2, ao, psiHb, psiE, psiB, Tq, rho, wv, aowR, aow2R, aow1, aow21, &
       PhiD, Phi, Amat, tmp, wvG, aowG, aow2G, PhiG, tth, tGth)

    integer(is), intent(in) :: nb, g0, nc2
    real(dp), intent(in)    :: ao(nb,nao_g,nc2), psiHb(nb,nocca,ncomp)
    real(dp), intent(inout) :: psiE(nb,max(nexp_q,1),ncomp), psiB(nb,maxnb_q,ncomp), Tq(nb,maxnb_q)
    real(dp), intent(inout) :: rho(nb,ndc,2), wv(nb,ndc,2), aowR(nb,nao_g,2), aow2R(nb,nao_g,3,2)
    real(dp), intent(inout) :: aow1(nb,nao_g,2), aow21(nb,nao_g,3,2), PhiD(nb,nao_g,2), Phi(nb,nao_g)
    real(dp), intent(inout) :: Amat(nb,nao_g), tmp(nb,max(maxnb_q,nocca))
    real(dp), intent(inout) :: wvG(nb,ndc), aowG(nb,nao_g), aow2G(nb,nao_g,3), PhiG(nb,nao_g)
    real(dp), intent(inout) :: tth(nao_g,3,nst_q+1), tGth(nao_g,3,max(nst_q,1))

    integer(is) :: v, s, t, c, g, gg, x, y, nocc(2), nGo, nbt
    real(dp)    :: acc

    nocc(1) = nocca; nocc(2) = nocca - 2
    nGo = nocca - 1

    ! reference potentials, their AO products and the reference densities
    do s = 1, 2
       do c = 1, ndc
          wv(:,c,s) = wgt(g0:g0+nb-1)*vxc(g0:g0+nb-1,c,s)
       enddo
    enddo
    call make_aow_s(nb, nc2, ao, wv, aowR, aow2R)
    do s = 1, 2
       call dgemm('N','T', nb, nao_g, nocc(s), 1.0_dp, psiHb, nb, CHg, nao_g, 0.0_dp, PhiD(1,1,s), nb)
       call dgemm('N','N', nb, nocc(s), nao_g, 1.0_dp, aowR(1,1,s), nb, CHg, nao_g, 0.0_dp, tmp, nb)
       call dgemm('N','T', nb, nao_g, nocc(s), 1.0_dp, tmp, nb, CHg, nao_g, 0.0_dp, Amat, nb)
       call probe_reduce_s(nb, nc2, ao, Amat, aow2R(1,1,1,s), PhiD(1,1,s), tth(1,1,1))
    enddo
    if (nst_q == 0) return

    ! explicit basis columns on the block
    if (nexp_q > 0) then
       do c = 1, ncomp
          call dgemm('N','N', nb, nexp_q, nao_g, 1.0_dp, ao(1,1,c), nb, Bexp_q, nao_g, 0.0_dp, psiE(1,1,c), nb)
       enddo
    endif
    if (gq) then
       do c = 1, ndc
          wvG(:,c) = wgt(g0:g0+nb-1)*vxcG(g0:g0+nb-1,c)
       enddo
       call make_aow_1s_s(nb, nc2, ao, wvG, aowG, aow2G)
       call dgemm('N','T', nb, nao_g, nGo, 1.0_dp, psiHb, nb, CHg, nao_g, 0.0_dp, PhiG, nb)
    endif

    do v = 1, nst_q
       ! (vxc^s, P^s_v) and the densities of P^s_v
       do s = 1, 2
          nbt = nbq(s,v)
          call term_density(nb, nbt, idxq(1,s,v), Sq(1,1,s,v), psiHb, psiE, psiB, Tq, rho(1,1,s))
          call term_phi_a(nb, nbt, psiB, Bq(1,1,s,v), Wq(1,1,s,v), aowR(1,1,s), tmp, Phi, Amat)
          call probe_reduce_s(nb, nc2, ao, Amat, aow2R(1,1,1,s), Phi, tth(1,1,v+1))
       enddo
       ! (f_xc[P_v], D^t): weighted kernel potential of P_v
       do g = 1, nb
          gg = g0 + g - 1
          do t = 1, 2
             do y = 1, ndc
                acc = 0.0_dp
                do s = 1, 2
                   do x = 1, ndc
                      acc = acc + rho(g,x,s)*fxcw(y,t,x,s,gg)
                   enddo
                enddo
                wv(g,y,t) = acc
             enddo
          enddo
       enddo
       call make_aow_s(nb, nc2, ao, wv, aow1, aow21)
       do t = 1, 2
          call dgemm('N','N', nb, nocc(t), nao_g, 1.0_dp, aow1(1,1,t), nb, CHg, nao_g, 0.0_dp, tmp, nb)
          call dgemm('N','T', nb, nao_g, nocc(t), 1.0_dp, tmp, nb, CHg, nao_g, 0.0_dp, Amat, nb)
          call probe_reduce_s(nb, nc2, ao, Amat, aow21(1,1,1,t), PhiD(1,1,t), tth(1,1,v+1))
       enddo
       if (.not. gq) cycle
       ! --- G-density terms of the extended method ---
       ! (vxc_G, M)
       nbt = nbq(3,v)
       call term_density(nb, nbt, idxq(1,3,v), Sq(1,1,3,v), psiHb, psiE, psiB, Tq, rho(1,1,1))
       call term_phi_a(nb, nbt, psiB, Bq(1,1,3,v), Wq(1,1,3,v), aowG, tmp, Phi, Amat)
       call probe_reduce_s(nb, nc2, ao, Amat, aow2G, Phi, tGth(1,1,v))
       ! (f_G[M], D_G), D_G = 2 C_G C_G^T (Phi_G carries C_G C_G^T: factor 2)
       call g_wv_f(nb, g0, rho(1,1,1), wv(1,1,1))
       wv(:,:,1) = 2.0_dp*wv(:,:,1)
       call make_aow_1s_s(nb, nc2, ao, wv(1,1,1), aow1(1,1,1), aow21(1,1,1,1))
       call dgemm('N','N', nb, nGo, nao_g, 1.0_dp, aow1(1,1,1), nb, CHg, nao_g, 0.0_dp, tmp, nb)
       call dgemm('N','T', nb, nao_g, nGo, 1.0_dp, tmp, nb, CHg, nao_g, 0.0_dp, Amat, nb)
       call probe_reduce_s(nb, nc2, ao, Amat, aow21(1,1,1,1), PhiG, tGth(1,1,v))
       if (.not. gq_kxc) cycle
       ! (f_G[D_s], D_s)
       nbt = nbq(4,v)
       call term_density(nb, nbt, idxq(1,4,v), Sq(1,1,4,v), psiHb, psiE, psiB, Tq, rho(1,1,2))
       call g_wv_f(nb, g0, rho(1,1,2), wv(1,1,1))
       call make_aow_1s_s(nb, nc2, ao, wv(1,1,1), aow1(1,1,1), aow21(1,1,1,1))
       call term_phi_a(nb, nbt, psiB, Bq(1,1,4,v), Wq(1,1,4,v), aow1(1,1,1), tmp, Phi, Amat)
       call probe_reduce_s(nb, nc2, ao, Amat, aow21(1,1,1,1), Phi, tGth(1,1,v))
       ! (k_G[D_s, D_s]/2, D_G) = (k_G[D_s, D_s], C_G C_G^T)
       call g_wv_k(nb, g0, rho(1,1,2), wv(1,1,1))
       call make_aow_1s_s(nb, nc2, ao, wv(1,1,1), aow1(1,1,1), aow21(1,1,1,1))
       call dgemm('N','N', nb, nGo, nao_g, 1.0_dp, aow1(1,1,1), nb, CHg, nao_g, 0.0_dp, tmp, nb)
       call dgemm('N','T', nb, nao_g, nGo, 1.0_dp, tmp, nb, CHg, nao_g, 0.0_dp, Amat, nb)
       call probe_reduce_s(nb, nc2, ao, Amat, aow21(1,1,1,1), PhiG, tGth(1,1,v))
    enddo

  end subroutine probe_q_block

  ! psi_B (values + gradients) of a term and its density rho = psi_B S psi_B^T
  subroutine term_density(nb, nbt, idx, S, psiHb, psiE, psiB, Tq, rho)

    integer(is), intent(in) :: nb, nbt, idx(maxnb_q)
    real(dp), intent(in)    :: S(maxnb_q,maxnb_q), psiHb(nb,nocca,ncomp), psiE(nb,max(nexp_q,1),ncomp)
    real(dp), intent(inout) :: psiB(nb,maxnb_q,ncomp), Tq(nb,maxnb_q)
    real(dp), intent(out)   :: rho(nb,ndc)

    call qf_density(nb, nbt, maxnb_q, nexp_q, idx, S, psiHb, psiE, psiB, Tq, rho)

  end subroutine term_density

  ! Phi = psi_B W^T and A = (aow B) W^T of a term
  subroutine term_phi_a(nb, nbt, psiB, B, W, aow, tmp, Phi, Amat)

    integer(is), intent(in) :: nb, nbt
    real(dp), intent(in)    :: psiB(nb,maxnb_q,ncomp), B(nao_g,maxnb_q), W(nao_g,maxnb_q), aow(nb,nao_g)
    real(dp), intent(inout) :: tmp(nb,max(maxnb_q,nocca))
    real(dp), intent(out)   :: Phi(nb,nao_g), Amat(nb,nao_g)

    if (nbt == 0) then
       Phi = 0.0_dp; Amat = 0.0_dp
       return
    endif
    call dgemm('N','T', nb, nao_g, nbt, 1.0_dp, psiB, nb, W, nao_g, 0.0_dp, Phi, nb)
    call dgemm('N','N', nb, nbt, nao_g, 1.0_dp, aow, nb, B, nao_g, 0.0_dp, tmp, nb)
    call dgemm('N','T', nb, nao_g, nbt, 1.0_dp, tmp, nb, W, nao_g, 0.0_dp, Amat, nb)

  end subroutine term_phi_a

  ! serial variants of make_aow / make_aow_1s / probe_reduce (called inside
  ! the parallel region of xc_probe_super)
  subroutine make_aow_s(nb, nc2, ao, wv, aow, aow2)

    integer(is), intent(in) :: nb, nc2
    real(dp), intent(in)    :: ao(nb,nao_g,nc2), wv(nb,ndc,2)
    real(dp), intent(out)   :: aow(nb,nao_g,2), aow2(nb,nao_g,3,2)

    integer(is) :: mu, s, c, x

    do mu = 1, nao_g
       do s = 1, 2
          aow(:,mu,s) = 0.5_dp*wv(:,1,s)*ao(:,mu,1)
          do c = 2, ncomp
             aow(:,mu,s) = aow(:,mu,s) + wv(:,c,s)*ao(:,mu,c)
          enddo
          do x = 1, 3
             aow2(:,mu,x,s) = 0.5_dp*wv(:,1,s)*ao(:,mu,1+x)
             if (ncomp > 1) then
                do c = 1, 3
                   aow2(:,mu,x,s) = aow2(:,mu,x,s) + wv(:,c+1,s)*ao(:,mu,d2idx(x,c))
                enddo
             endif
          enddo
       enddo
    enddo

  end subroutine make_aow_s

  subroutine make_aow_1s_s(nb, nc2, ao, wv, aow, aow2)

    integer(is), intent(in) :: nb, nc2
    real(dp), intent(in)    :: ao(nb,nao_g,nc2), wv(nb,ndc)
    real(dp), intent(out)   :: aow(nb,nao_g), aow2(nb,nao_g,3)

    integer(is) :: mu, c, x

    do mu = 1, nao_g
       aow(:,mu) = 0.5_dp*wv(:,1)*ao(:,mu,1)
       do c = 2, ncomp
          aow(:,mu) = aow(:,mu) + wv(:,c)*ao(:,mu,c)
       enddo
       do x = 1, 3
          aow2(:,mu,x) = 0.5_dp*wv(:,1)*ao(:,mu,1+x)
          if (ncomp > 1) then
             do c = 1, 3
                aow2(:,mu,x) = aow2(:,mu,x) + wv(:,c+1)*ao(:,mu,d2idx(x,c))
             enddo
          endif
       enddo
    enddo

  end subroutine make_aow_1s_s

  subroutine probe_reduce_s(nb, nc2, ao, Amat, aow2, Phi, t)

    integer(is), intent(in)  :: nb, nc2
    real(dp), intent(in)     :: ao(nb,nao_g,nc2), Amat(nb,nao_g), aow2(nb,nao_g,3), Phi(nb,nao_g)
    real(dp), intent(inout)  :: t(nao_g,3)

    integer(is) :: mu, x

    do mu = 1, nao_g
       do x = 1, 3
          t(mu,x) = t(mu,x) + dot_product(ao(:,mu,1+x), Amat(:,mu)) + dot_product(aow2(:,mu,x), Phi(:,mu))
       enddo
    enddo

  end subroutine probe_reduce_s

!######################################################################
! xc_g_potential_q: the G-density kernel potentials of the EMRSF
! gradient (as xc_g_potential) for the densities given as quadratic
! forms in the basis [C_H (cached) | Bexp]: M (term nbM, idxM, SM) and
! D_s (nbS, idxS, SS). Instead of the full AO matrix of v[f rho(D_s)]
! its products with the explicit columns and the hole MOs are returned:
!   MM(i,nu) = projection of v[f rho(M)] on the hole MOs (V(H, all) = MM C)
!   MK(i,nu) = the same for v[k rho(D_s)^2]           (when want_k)
!   VE(mu,j) = (V1 Bexp)(mu,j),  VH(mu,i) = (V1 C_H)(mu,i),  V1 = v[f rho(D_s)] (AO)
! block-parallel over the cached AO blocks, sequential MKL per thread.
!######################################################################
  subroutine xc_g_potential_q(nexp, Bexp, maxnb, nbM, idxM, SM, nbS, idxS, SS, want_k, MM, MK, VE, VH)

    integer(is), intent(in) :: nexp, maxnb, nbM, nbS, idxM(maxnb), idxS(maxnb)
    real(dp), intent(in)    :: Bexp(nao_g,max(nexp,1)), SM(maxnb,maxnb), SS(maxnb,maxnb)
    logical, intent(in)     :: want_k
    real(dp), intent(out)   :: MM(nocca,nao_g), MK(nocca,nao_g), VE(nao_g,max(nexp,1)), VH(nao_g,nocca)

    integer(is)           :: ib, nb, nbmax, nthr
    integer(c_int)        :: prev
    logical               :: par
    real(dp)              :: t0
    real(dp), allocatable :: psiE(:), psiB(:), Tq(:), rhoM(:), rhoS(:), wv(:), Bst(:), Ast(:), aowS(:), &
                             aowE(:), aowH(:), MthM(:,:), MthK(:,:), VEth(:,:), VHth(:,:)

    if (.not. g_ready) call mrsf_error('xc_g_potential_q: G kernel set not installed')
    if (.not. ao_cached) call mrsf_error('xc_g_potential_q: the AO values must be cached (raise mem_budget)')
    if (want_k .and. .not. g_have_kxc) call mrsf_error('xc_g_potential_q: third derivatives not installed')
    t0 = wall_time()
    nbmax = 0
    do ib = 1, nblocks
       nbmax = max(nbmax, gend(ib) - gbeg(ib) + 1)
    enddo
    nthr = 1
    !$ nthr = omp_get_max_threads()
    par = (nthr > 1) .and. (nblocks >= nthr)
    MM = 0.0_dp; MK = 0.0_dp; VE = 0.0_dp; VH = 0.0_dp

    !$omp parallel if(par) default(shared) private(ib, nb, prev, psiE, psiB, Tq, rhoM, rhoS, wv, Bst, Ast, &
    !$omp aowS, aowE, aowH, MthM, MthK, VEth, VHth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(1_c_int)
#endif
    allocate(psiE(nbmax*max(nexp,1)*ncomp), psiB(nbmax*maxnb*ncomp), Tq(nbmax*maxnb), rhoM(nbmax*ndc), &
         rhoS(nbmax*ndc), wv(nbmax*ndc), Bst(2*nbmax*nao_g), Ast(nocca*2*nbmax), aowS(nbmax*nao_g), &
         aowE(nbmax*max(nexp,1)), aowH(nbmax*nocca), MthM(nocca,nao_g), MthK(nocca,nao_g), &
         VEth(nao_g,max(nexp,1)), VHth(nao_g,nocca))
    MthM = 0.0_dp; MthK = 0.0_dp; VEth = 0.0_dp; VHth = 0.0_dp
    !$omp do schedule(dynamic)
    do ib = 1, nblocks
       nb = gend(ib) - gbeg(ib) + 1
       call gpot_q_block(nb, gbeg(ib), aoc(aoff(ib)+1), psiH(poff(ib)+1), nexp, Bexp, maxnb, nbM, idxM, SM, &
            nbS, idxS, SS, want_k, psiE, psiB, Tq, rhoM, rhoS, wv, Bst, Ast, aowS, aowE, aowH, &
            MthM, MthK, VEth, VHth)
    enddo
    !$omp end do
    !$omp critical
    MM = MM + MthM
    MK = MK + MthK
    VE = VE + VEth
    VH = VH + VHth
    !$omp end critical
    deallocate(psiE, psiB, Tq, rhoM, rhoS, wv, Bst, Ast, aowS, aowE, aowH, MthM, MthK, VEth, VHth)
#ifdef USE_MKL
    if (par) prev = mkl_threads_local(0_c_int)
#endif
    !$omp end parallel
    time_xcg = time_xcg + wall_time() - t0
    nxcg_calls = nxcg_calls + 1

  end subroutine xc_g_potential_q

  subroutine gpot_q_block(nb, g0, ao, psiHb, nexp, Bexp, maxnb, nbM, idxM, SM, nbS, idxS, SS, want_k, &
       psiE, psiB, Tq, rhoM, rhoS, wv, Bst, Ast, aowS, aowE, aowH, MthM, MthK, VEth, VHth)

    integer(is), intent(in) :: nb, g0, nexp, maxnb, nbM, nbS, idxM(maxnb), idxS(maxnb)
    real(dp), intent(in)    :: ao(nb,nao_g,ncomp), psiHb(nb,nocca,ncomp)
    real(dp), intent(in)    :: Bexp(nao_g,max(nexp,1)), SM(maxnb,maxnb), SS(maxnb,maxnb)
    logical, intent(in)     :: want_k
    real(dp), intent(inout) :: psiE(nb,max(nexp,1),ncomp), psiB(nb,maxnb,ncomp), Tq(nb,maxnb)
    real(dp), intent(inout) :: rhoM(nb,ndc), rhoS(nb,ndc), wv(nb,ndc)
    real(dp), intent(inout) :: Bst(2*nb,nao_g), Ast(nocca,2*nb), aowS(nb,nao_g)
    real(dp), intent(inout) :: aowE(nb,max(nexp,1)), aowH(nb,nocca)
    real(dp), intent(inout) :: MthM(nocca,nao_g), MthK(nocca,nao_g), VEth(nao_g,max(nexp,1)), VHth(nao_g,nocca)

    integer(is) :: c, g, mu, j

    ! block-constant halves of the stacked operands of g_project
    do g = 1, nb
       Ast(1:nocca,g) = psiHb(g,1:nocca,1)
    enddo
    do mu = 1, nao_g
       Bst(nb+1:2*nb,mu) = ao(1:nb,mu,1)
    enddo
    ! explicit basis columns on the block
    if (nexp > 0) then
       do c = 1, ncomp
          call dgemm('N','N', nb, nexp, nao_g, 1.0_dp, ao(1,1,c), nb, Bexp, nao_g, 0.0_dp, psiE(1,1,c), nb)
       enddo
    endif
    ! densities of M and D_s
    call qf_density(nb, nbM, maxnb, nexp, idxM, SM, psiHb, psiE, psiB, Tq, rhoM)
    call qf_density(nb, nbS, maxnb, nexp, idxS, SS, psiHb, psiE, psiB, Tq, rhoS)
    ! v[f rho(M)] on the hole MOs
    call g_wv_f(nb, g0, rhoM, wv)
    call g_project(nb, ao, psiHb, wv, Bst, Ast, MthM)
    ! v[k rho(D_s)^2] on the hole MOs
    if (want_k) then
       call g_wv_k(nb, g0, rhoS, wv)
       call g_project(nb, ao, psiHb, wv, Bst, Ast, MthK)
    endif
    ! v[f rho(D_s)] times the explicit columns and the hole MOs:
    ! (V1 B)(mu,j) = sum_g [aowS(g,mu) psi_j(g) + phi_mu(g) aowpsi_j(g)]
    call g_wv_f(nb, g0, rhoS, wv)
    do mu = 1, nao_g
       aowS(1:nb,mu) = 0.5_dp*wv(1:nb,1)*ao(1:nb,mu,1)
       do c = 2, ncomp
          aowS(1:nb,mu) = aowS(1:nb,mu) + wv(1:nb,c)*ao(1:nb,mu,c)
       enddo
    enddo
    do j = 1, nexp
       aowE(1:nb,j) = 0.5_dp*wv(1:nb,1)*psiE(1:nb,j,1)
       do c = 2, ncomp
          aowE(1:nb,j) = aowE(1:nb,j) + wv(1:nb,c)*psiE(1:nb,j,c)
       enddo
    enddo
    do j = 1, nocca
       aowH(1:nb,j) = 0.5_dp*wv(1:nb,1)*psiHb(1:nb,j,1)
       do c = 2, ncomp
          aowH(1:nb,j) = aowH(1:nb,j) + wv(1:nb,c)*psiHb(1:nb,j,c)
       enddo
    enddo
    if (nexp > 0) then
       call dgemm('T','N', nao_g, nexp, nb, 1.0_dp, aowS, nb, psiE(1,1,1), nb, 1.0_dp, VEth, nao_g)
       call dgemm('T','N', nao_g, nexp, nb, 1.0_dp, ao(1,1,1), nb, aowE, nb, 1.0_dp, VEth, nao_g)
    endif
    call dgemm('T','N', nao_g, nocca, nb, 1.0_dp, aowS, nb, psiHb(1,1,1), nb, 1.0_dp, VHth, nao_g)
    call dgemm('T','N', nao_g, nocca, nb, 1.0_dp, ao(1,1,1), nb, aowH, nb, 1.0_dp, VHth, nao_g)

  end subroutine gpot_q_block

  ! psi_B (values + gradients) of a quadratic-form term and its density
  ! rho = sum_ab psi_a S_ab psi_b (value and gradient components)
  subroutine qf_density(nb, nbt, maxnb, nexp, idx, S, psiHb, psiE, psiB, Tq, rho)

    integer(is), intent(in) :: nb, nbt, maxnb, nexp, idx(maxnb)
    real(dp), intent(in)    :: S(maxnb,maxnb), psiHb(nb,nocca,ncomp), psiE(nb,max(nexp,1),ncomp)
    real(dp), intent(inout) :: psiB(nb,maxnb,ncomp), Tq(nb,maxnb)
    real(dp), intent(out)   :: rho(nb,ndc)

    integer(is) :: a, c, j

    do c = 1, ncomp
       do a = 1, nbt
          j = idx(a)
          if (j <= nocca) then
             psiB(:,a,c) = psiHb(:,j,c)
          else
             psiB(:,a,c) = psiE(:,j-nocca,c)
          endif
       enddo
    enddo
    rho = 0.0_dp
    if (nbt == 0) return
    call dgemm('N','N', nb, nbt, nbt, 1.0_dp, psiB, nb, S, maxnb, 0.0_dp, Tq, nb)
    do a = 1, nbt
       rho(:,1) = rho(:,1) + Tq(:,a)*psiB(:,a,1)
       do c = 2, ndc
          rho(:,c) = rho(:,c) + 2.0_dp*Tq(:,a)*psiB(:,a,c)
       enddo
    enddo

  end subroutine qf_density

end module mrsf_xcgrid
