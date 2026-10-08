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

  use mrsf_constants
  use mrsf_global
  use mrsf_io

  implicit none

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
  real(dp), allocatable    :: CHg(:,:), CPg(:,:)   ! (nao, nocca), (nao, nvirb)
  logical                  :: xc_ready = .false., ao_cached = .false.
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

    integer(is) :: g, s, x, t, y, ib

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
    allocate(CHg(nao_g,nocca), CPg(nao_g,nvirb))
    CHg = CH1; CPg = CP1
    ao_cached = cache_ao
    if (ao_cached) then
       allocate(aoc(int(ncomp,i8)*int(nao_g,i8)*int(ngrid,i8)))
       allocate(psiH(int(ncomp,i8)*int(nocca,i8)*int(ngrid,i8)))
    endif
    xc_ready = .true.

  end subroutine xc_init

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
    xc_ready = .false.; ao_cached = .false.

  end subroutine xc_free

  subroutine xc_cleanup()

    if (allocated(Lcur)) deallocate(Lcur, Rcur, rch)
    if (allocated(Macc)) deallocate(Macc)
    if (allocated(tacc)) deallocate(tacc)
    kcur = 0; nvcur = 0

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
    do ib = 1, nblocks
       nb = gend(ib) - gbeg(ib) + 1
       call kernel_block(ib, nb, aoc(aoff(ib)+1), psiH(poff(ib)+1))
    enddo

  end subroutine xc_cached

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

  end subroutine xc_closed_shell_init

!######################################################################
! xc_kernel_cv: closed-shell kernel of the extended method on a batch
! of CV trial amplitudes.  X(nvirb,ncol,nvec) holds the vectors; the
! CV columns nocca+1..ncol, rows 3..nvirb, are Y(b,j) (b in V, j in C).
! Trial density rho_t = 2 sum_{jb} Y_bj phi_j phi_b = L R^T + R L^T with
! L = C_V Y, R = C_C (the first nC cached hole MOs); with both spin
! densities equal to rho_t the alpha potential is 2 (f_aa + f_ab) *
! (.|jb) Y, so L is taken as C_V Y / 2 and
!   Vcv(j,b,v) = (jb| f_aa + f_ab |ld) Y_ld.
! Loop order block-outer, vector-inner: each cached AO block is read
! once per batch; all elementwise work is on unit-stride grid vectors.
!######################################################################
  subroutine xc_kernel_cv(nvec, X, Vcv)

    integer(is), intent(in) :: nvec
    real(dp), intent(in)    :: X(nvirb,ncol,nvec)
    real(dp), intent(out)   :: Vcv(nC,nV,nvec)

    real(dp), allocatable :: Lall(:,:,:), Macc(:,:,:)
    real(dp), allocatable :: psiL(:,:,:), rho1(:,:), wv(:,:), aow(:,:), aowC(:,:)
    integer(is)           :: v, ib, nb, g0, c, g, gg, ix, iy, mu, jc, nk
    integer(i8)           :: ao0, ps0
    real(dp)              :: acc, t0, t1

    if (.not. xc_ready) call mrsf_error('xc_kernel_cv: grid not initialised')
    if (.not. ao_cached) call mrsf_error('xc_kernel_cv: AO values not cached')
    if (.not. allocated(fxcw_cs)) call mrsf_error('xc_kernel_cv: closed-shell kernel not initialised')
    t0 = wall_time()

    ! factors L_v = C_V Y_v / 2 for all vectors, (nao, nC, nvec)
    nk = nC*nvec
    allocate(Lall(nao_g,nC,nvec), Macc(nC,nao_g,nvec))
    do v = 1, nvec
       call dgemm('N','N', nao_g, nC, nV, 0.5_dp, CPg(1,3), nao_g, X(3,nocca+1,v), nvirb, &
            0.0_dp, Lall(1,1,v), nao_g)
    enddo
    Macc = 0.0_dp

    do ib = 1, nblocks
       nb  = gend(ib) - gbeg(ib) + 1
       g0  = gbeg(ib)
       ao0 = aoff(ib)
       ps0 = poff(ib)
       allocate(psiL(nb,nk,ncomp), rho1(nb,ndc), wv(nb,ndc), aow(nb,nao_g), aowC(nC,nb))
       t1 = wall_time()
       ! factor values on the block, all vectors at once
       do c = 1, ncomp
          call dgemm('N','N', nb, nk, nao_g, 1.0_dp, aoc(ao0 + 1 + int(nb,i8)*int(nao_g,i8)*int(c-1,i8)), nb, &
               Lall, nao_g, 0.0_dp, psiL(1,1,c), nb)
       enddo
       tcs(1) = tcs(1) + wall_time() - t1
       do v = 1, nvec
          t1 = wall_time()
          jc = nC*(v-1) + 1
          ! trial density from the factors and the cached core MOs
          call block_rho(nb, nC, psiL(1,jc,1), nb*nk, psiH(ps0+1), nb*nocca, rho1)
          ! weighted kernel potential wv(g,y) = w sum_x rho(g,x) (f_aa+f_ab)(y,x)
          !$omp parallel do private(gg,iy,ix,acc)
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
          !$omp end parallel do
          tcs(2) = tcs(2) + wall_time() - t1
          t1 = wall_time()
          call make_aow1(nb, aoc(ao0+1), wv, aow)
          tcs(3) = tcs(3) + wall_time() - t1
          t1 = wall_time()
          ! M_C += psi_C^T aow + (C_C^T aow^T) phi
          call dgemm('T','N', nC, nao_g, nb, 1.0_dp, psiH(ps0+1), nb, aow, nb, &
               1.0_dp, Macc(1,1,v), nC)
          call dgemm('T','T', nC, nb, nao_g, 1.0_dp, CHg, nao_g, aow, nb, 0.0_dp, aowC, nC)
          call dgemm('N','N', nC, nao_g, nb, 1.0_dp, aowC, nC, aoc(ao0+1), nb, &
               1.0_dp, Macc(1,1,v), nC)
          tcs(4) = tcs(4) + wall_time() - t1
       enddo
       deallocate(psiL, rho1, wv, aow, aowC)
    enddo

    t1 = wall_time()
    do v = 1, nvec
       call dgemm('N','N', nC, nV, nao_g, 1.0_dp, Macc(1,1,v), nC, CPg(1,3), nao_g, &
            0.0_dp, Vcv(1,1,v), nC)
    enddo
    tcs(5) = tcs(5) + wall_time() - t1
    deallocate(Lall, Macc)
    time_xc = time_xc + wall_time() - t0
    nxc_calls = nxc_calls + 1
    nxc_vecs = nxc_vecs + nvec

  end subroutine xc_kernel_cv

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
    if (sum(tcs) > 0.0_dp) then
       write(6,'(2x,a)') 'closed-shell CV kernel (s): factor values, densities+potentials, AO products, projections, final'
       write(6,'(2x,5f10.3)') tcs
    endif
    flush(6)

  end subroutine xc_timings

end module mrsf_xcgrid
