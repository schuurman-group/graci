!**********************************************************************
! mrsf_davidson: in-core block Davidson solver for the lowest roots of
! the (symmetric) MRSF-TDDFT A matrix of one (multiplicity, irrep)
! block. Exact-diagonal preconditioner, locking of converged roots,
! re-orthogonalised Gram-Schmidt expansion, Ritz-vector collapse.
! Vectors have xdim_tot entries (xdim for the standard method, xdim +
! nvirb*nC with the extended CV columns).
!**********************************************************************
module mrsf_davidson

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_space
  use mrsf_sigma
  use mrsf_extended, only: kernel_off

  implicit none

contains

!######################################################################
! davidson_solve
!   nroots : requested roots (reduced if the active space is smaller)
!   nextra : additional tracked roots
!   maxvec : maximum subspace dimension (raised to >= 2*nsolve)
!   tol    : convergence threshold on the residual norm
!   nguess : start vectors supplied in the first nguess columns of vecs
!            (warm start: the vectors of a previous geometry aligned to
!            the present MOs, in the compressed layout); they are masked,
!            orthonormalised, wrong-irrep/dependent ones dropped, at most
!            nsolve kept, and the block is completed with unit vectors on
!            the lowest active diagonal entries not yet represented
!   ener   : eigenvalues (omega, relative to the reference energy) of the
!            nvec_out = nsolve tracked roots (the first nroots converged)
!   vecs   : compressed Ritz vectors of the nvec_out tracked roots
!   iconv  : 1 if converged, 0 otherwise
!######################################################################
  subroutine davidson_solve(mult, irrep, nroots, nextra, maxvec_in, maxiter, tol, &
       precond, inner_iter, nguess, ener, vecs, nvec_out, niter, iconv)

    integer(is), intent(in)    :: mult, irrep, nextra, maxvec_in, maxiter
    integer(is), intent(in)    :: precond, inner_iter, nguess
    integer(is), intent(inout) :: nroots
    real(dp), intent(in)       :: tol
    real(dp), intent(out)      :: ener(*)
    real(dp), intent(inout)    :: vecs(xdim_tot,*)
    integer(is), intent(out)   :: nvec_out, niter, iconv

    integer(is), allocatable :: act(:), cols(:)
    real(dp), allocatable    :: d(:), V(:,:), AV(:,:), G(:,:), Gc(:,:), theta(:), work(:)
    real(dp), allocatable    :: R(:,:), AR(:,:), res(:,:), rnorm(:), dnew(:,:), ovl(:,:)
    real(dp), allocatable    :: dtmp(:), gvec(:), govl(:,:)
    integer(is)              :: nact, nsolve, maxvec, currdim, ist, iend, nnew, nkept, ngacc
    integer(is)              :: ia, a, i, k, l, iter, lwork, info, kmin, ndim, ninner_tot
    real(dp)                 :: denom, rmax, dnorm, proj, dmin, gnorm
    logical                  :: collapsed, use_inner

    iconv = 0
    niter = 0
    nvec_out = 0
    ndim  = xdim_tot
    ! the Jacobi-Davidson correction with the kernel-free operator pays
    ! off only where the exact operator carries the CV kernel
    use_inner  = (precond == 2) .and. extended .and. use_kernel .and. (mult == 1)
    ninner_tot = 0

    ! active slots
    allocate(act(ndim))
    nact = 0
    do i = 1, ncol
       do a = 1, nvirb
          if (slot_active(mult, irrep, a, i)) then
             nact = nact + 1
             act(nact) = a + nvirb*(i-1)
          endif
       enddo
    enddo
    if (nact == 0) then
       nroots = 0
       return
    endif

    nroots = min(nroots, nact)
    nsolve = min(nroots + nextra, nact)
    maxvec = max(maxvec_in, 2*nsolve)
    maxvec = min(maxvec, nact)

    allocate(d(ndim), V(ndim,maxvec), AV(ndim,maxvec), G(maxvec,maxvec), Gc(maxvec,maxvec))
    allocate(theta(maxvec), R(ndim,nsolve), AR(ndim,nsolve), res(ndim,nsolve))
    allocate(rnorm(nsolve), dnew(ndim,nsolve), ovl(maxvec,nsolve), dtmp(nact), cols(nsolve))
    V = 0.0_dp; AV = 0.0_dp; G = 0.0_dp

    ! diagonal and preconditioner
    call diagonal(mult, d)
    do ia = 1, ndim
       if (slot_irrep(ia) /= irrep .and. irrep >= 0) d(ia) = 1.0e20_dp
    enddo

    ! start vectors: (1) the supplied guess vectors (warm start), masked
    ! for this multiplicity/irrep, orthonormalised, the empty or dependent
    ! ones dropped, at most nsolve kept
    ngacc = 0
    allocate(gvec(ndim))
    do k = 1, nguess
       if (ngacc == nsolve) exit
       gvec = vecs(:,k)
       call apply_mask(mult, irrep, 1_is, gvec)
       gnorm = sqrt(dot_product(gvec, gvec))
       if (gnorm < 1.0e-3_dp) cycle
       gvec = gvec / gnorm
       do l = 1, 2
          do i = 1, ngacc
             proj = dot_product(V(:,i), gvec)
             gvec = gvec - proj * V(:,i)
          enddo
       enddo
       gnorm = sqrt(dot_product(gvec, gvec))
       if (gnorm < 1.0e-6_dp) cycle
       ngacc = ngacc + 1
       V(:,ngacc) = gvec / gnorm
    enddo
    ! (2) unit vectors on the lowest active diagonal elements, skipping
    ! the slots already represented by the guess block
    do k = 1, nact
       dtmp(k) = d(act(k))
    enddo
    k = ngacc
    do while (k < nsolve)
       kmin = 0
       dmin = huge(1.0_dp)
       do l = 1, nact
          if (dtmp(l) < dmin) then
             dmin = dtmp(l)
             kmin = l
          endif
       enddo
       if (kmin == 0) exit
       dtmp(kmin) = huge(1.0_dp)
       gvec = 0.0_dp
       gvec(act(kmin)) = 1.0_dp
       do i = 1, k
          proj = dot_product(V(:,i), gvec)
          gvec = gvec - proj * V(:,i)
       enddo
       gnorm = sqrt(dot_product(gvec, gvec))
       if (gnorm < 0.5_dp) cycle
       k = k + 1
       V(:,k) = gvec / gnorm
    enddo
    if (k < nsolve) then
       ! the active space is too small to hold nsolve orthonormal vectors
       nsolve = k
       nroots = min(nroots, nsolve)
    endif
    if (ngacc > 0) then
       allocate(govl(ngacc,nsolve))
       govl = 0.0_dp
    endif

    ! LAPACK workspace
    lwork = max(1_is, 3*maxvec + 64)
    allocate(work(lwork))

    if (verbose) then
       write(6,'(/,2x,a,i0,a,i0,a,i0,a,i0)') 'Davidson: mult = ', mult, &
            ', irrep = ', irrep, ', active dim = ', nact, ', roots = ', nroots
       if (nguess > 0) write(6,'(2x,a,i0,a,i0,a)') 'warm start: ', ngacc, ' of ', nguess, &
            ' supplied vectors accepted as start vectors'
       write(6,'(2x,a)') ' iter   dim   max|r|         lowest eigenvalues'
    endif

    currdim = nsolve
    ist = 1
    iend = nsolve
    ! keep the accepted guess block for the overlap diagnostic
    if (ngacc > 0) then
       deallocate(gvec)
       allocate(gvec(ndim*ngacc))
       gvec = reshape(V(:,1:ngacc), (/ ndim*ngacc /))
    endif

    do iter = 1, maxiter
       niter = iter
       nnew = iend - ist + 1

       ! sigma vectors for the new subspace vectors
       call sigma_batch(nnew, mult, irrep, V(:,ist:iend), AV(:,ist:iend))

       ! subspace matrix update
       call dgemm('T','N', currdim, nnew, ndim, 1.0_dp, V, ndim, AV(1,ist), ndim, &
            0.0_dp, G(1,ist), maxvec)
       do k = ist, iend
          do l = 1, currdim
             if (l < ist) then
                G(k,l) = G(l,k)
             else if (l /= k) then
                G(k,l) = 0.5_dp * (G(k,l) + G(l,k))
                G(l,k) = G(k,l)
             endif
          enddo
       enddo

       ! diagonalise
       Gc(1:currdim,1:currdim) = G(1:currdim,1:currdim)
       call dsyev('V','U', currdim, Gc, maxvec, theta, work, lwork, info)
       if (info /= 0) call mrsf_error('dsyev failed in the Davidson solver')

       ! Ritz vectors and residuals
       call dgemm('N','N', ndim, nsolve, currdim, 1.0_dp, V, ndim, Gc, maxvec, 0.0_dp, R, ndim)
       call dgemm('N','N', ndim, nsolve, currdim, 1.0_dp, AV, ndim, Gc, maxvec, 0.0_dp, AR, ndim)
       do k = 1, nsolve
          res(:,k) = AR(:,k) - theta(k) * R(:,k)
          rnorm(k) = sqrt(dot_product(res(:,k), res(:,k)))
       enddo
       rmax = maxval(rnorm(1:nroots))

       if (verbose) then
          write(6,'(2x,i4,2x,i5,2x,es10.2,2x,6f12.6)') iter, currdim, rmax, &
               theta(1:min(6_is,nroots))
          flush(6)
       endif

       if (rmax < tol) then
          iconv = 1
          exit
       endif

       ! corrections for the unconverged tracked roots: diagonal
       ! preconditioner, or the Jacobi-Davidson correction equation solved
       ! with the kernel-free operator (singlets with the CV kernel)
       nnew = 0
       do k = 1, nsolve
          if (rnorm(k) < tol) cycle
          nnew = nnew + 1
          cols(nnew) = k
       enddo
       if (nnew == 0) exit
       if (use_inner) then
          call jd_correction(mult, irrep, ndim, nsolve, nnew, cols, theta, R, res, d, &
               inner_iter, dnew)
          ninner_tot = ninner_tot + inner_iter*nnew
       else
          do k = 1, nnew
             dnew(:,k) = 0.0_dp
             do l = 1, nact
                ia = act(l)
                denom = d(ia) - theta(cols(k))
                if (abs(denom) < 1.0e-4_dp) denom = sign(1.0e-4_dp, denom)
                dnew(ia,k) = -res(ia,cols(k)) / denom
             enddo
          enddo
       endif

       ! collapse if the subspace would overflow
       collapsed = .false.
       if (currdim + nnew > maxvec) then
          V(:,1:nsolve)  = R
          AV(:,1:nsolve) = AR
          G = 0.0_dp
          do k = 1, nsolve
             G(k,k) = theta(k)
          enddo
          currdim = nsolve
          collapsed = .true.
       endif

       ! orthogonalise the corrections against the subspace (twice)
       do k = 1, 2
          call dgemm('T','N', currdim, nnew, ndim, 1.0_dp, V, ndim, dnew, ndim, 0.0_dp, ovl, maxvec)
          call dgemm('N','N', ndim, nnew, currdim, -1.0_dp, V, ndim, ovl, maxvec, 1.0_dp, dnew, ndim)
       enddo

       ! modified Gram-Schmidt among the corrections, drop small vectors
       nkept = 0
       do k = 1, nnew
          do l = 1, nkept
             proj = dot_product(dnew(:,l), dnew(:,k))
             dnew(:,k) = dnew(:,k) - proj * dnew(:,l)
          enddo
          dnorm = sqrt(dot_product(dnew(:,k), dnew(:,k)))
          if (dnorm < 1.0e-8_dp) cycle
          nkept = nkept + 1
          if (nkept /= k) dnew(:,nkept) = dnew(:,k)
          dnew(:,nkept) = dnew(:,nkept) / dnorm
       enddo
       if (nkept == 0) exit

       V(:,currdim+1:currdim+nkept) = dnew(:,1:nkept)
       ist  = currdim + 1
       iend = currdim + nkept
       currdim = iend
    enddo

    nvec_out = nsolve
    ener(1:nsolve) = theta(1:nsolve)
    vecs(1:ndim,1:nsolve) = R(:,1:nsolve)

    if (verbose) then
       if (iconv == 1) then
          write(6,'(2x,a,i0,a)') 'Davidson converged in ', niter, ' iterations'
       else
          write(6,'(2x,a,i0,a)') 'Davidson NOT converged after ', niter, ' iterations'
       endif
       if (use_inner) write(6,'(2x,a,i0,a)') 'inner correction equations: ', ninner_tot, &
            ' kernel-free operator applications'
       if (ngacc > 0) then
          ! norm of the projection of each converged root on the guess span
          call dgemm('T','N', ngacc, nsolve, ndim, 1.0_dp, gvec, ndim, R, ndim, 0.0_dp, govl, ngacc)
          write(6,'(2x,a,8f8.4)') 'warm start: |P_guess R_j| per root: ', &
               (sqrt(dot_product(govl(:,k), govl(:,k))), k = 1, min(nroots, 8_is))
       endif
    endif

    deallocate(act, d, V, AV, G, Gc, theta, work, R, AR, res, rnorm, dnew, ovl, dtmp, cols, gvec)
    if (allocated(govl)) deallocate(govl)

  end subroutine davidson_solve

!######################################################################
! jd_correction: Jacobi-Davidson-type corrections of the tracked roots
! cols(1:ncol),
!   P (A0 - theta_k) P t_k = -P r_k,   P = 1 - R R^T,   t_k _|_ R,
! A0 = the sigma operator without the CV kernel, solved with m steps of
! MINRES (Paige-Saunders) preconditioned with P |d - theta_k|^-1 P. All
! columns advance in lockstep so that every step is one kernel-free
! sigma batch; the inner work is booked separately from the exact sigma
! vectors (time_inner, ninner_vecs).
!######################################################################
  subroutine jd_correction(mult, irrep, ndim, nsolve, ncol, cols, theta, R, res, d, m, T)

    integer(is), intent(in) :: mult, irrep, ndim, nsolve, ncol, cols(ncol), m
    real(dp), intent(in)    :: theta(*), R(ndim,nsolve), res(ndim,*), d(ndim)
    real(dp), intent(out)   :: T(ndim,ncol)

    real(dp), allocatable :: r1(:,:), r2(:,:), y(:,:), v(:,:), w(:,:), w1(:,:), w2(:,:)
    real(dp), allocatable :: b(:,:), kd(:,:), az(:,:)
    real(dp), allocatable :: th(:), oldb(:), beta(:), dbar(:), epsln(:), phibar(:), cs(:), sn(:), alfa(:)
    logical, allocatable  :: done(:)
    real(dp)              :: oldeps, delta, gbar, gam, phi
    integer(is)           :: c, k, it

    allocate(r1(ndim,ncol), r2(ndim,ncol), y(ndim,ncol), v(ndim,ncol), w(ndim,ncol), &
         w1(ndim,ncol), w2(ndim,ncol), b(ndim,ncol), kd(ndim,ncol), az(ndim,ncol))
    allocate(th(ncol), oldb(ncol), beta(ncol), dbar(ncol), epsln(ncol), phibar(ncol), &
         cs(ncol), sn(ncol), alfa(ncol), done(ncol))

    do c = 1, ncol
       k = cols(c)
       th(c)   = theta(k)
       kd(:,c) = max(abs(d - th(c)), 1.0e-4_dp)
       b(:,c)  = -res(:,k)
    enddo
    call project(b)

    ! start vector: the projected diagonal-preconditioned residual
    T = b / kd
    call project(T)
    call apply_m(T, az)
    r1 = b - az
    y  = r1 / kd
    call project(y)
    done = .false.
    do c = 1, ncol
       beta(c) = sqrt(max(dot_product(r1(:,c), y(:,c)), 0.0_dp))
       if (beta(c) <= 0.0_dp) done(c) = .true.
    enddo
    oldb = 0.0_dp; dbar = 0.0_dp; epsln = 0.0_dp; phibar = beta
    cs = -1.0_dp; sn = 0.0_dp
    w = 0.0_dp; w2 = 0.0_dp; r2 = r1

    do it = 1, m
       if (all(done)) exit
       do c = 1, ncol
          if (done(c)) then
             v(:,c) = 0.0_dp
          else
             v(:,c) = y(:,c) / beta(c)
          endif
       enddo
       call apply_m(v, y)
       do c = 1, ncol
          if (done(c)) cycle
          if (it >= 2) y(:,c) = y(:,c) - (beta(c)/oldb(c)) * r1(:,c)
          alfa(c) = dot_product(v(:,c), y(:,c))
          y(:,c)  = y(:,c) - (alfa(c)/beta(c)) * r2(:,c)
          r1(:,c) = r2(:,c)
          r2(:,c) = y(:,c)
          y(:,c)  = r2(:,c) / kd(:,c)
       enddo
       call project(y)
       do c = 1, ncol
          if (done(c)) cycle
          oldb(c) = beta(c)
          beta(c) = sqrt(max(dot_product(r2(:,c), y(:,c)), 0.0_dp))
          oldeps   = epsln(c)
          delta    = cs(c)*dbar(c) + sn(c)*alfa(c)
          gbar     = sn(c)*dbar(c) - cs(c)*alfa(c)
          epsln(c) = sn(c)*beta(c)
          dbar(c)  = -cs(c)*beta(c)
          gam      = max(sqrt(gbar**2 + beta(c)**2), epsilon(1.0_dp))
          cs(c)    = gbar/gam
          sn(c)    = beta(c)/gam
          phi      = cs(c)*phibar(c)
          phibar(c) = sn(c)*phibar(c)
          w1(:,c) = w2(:,c)
          w2(:,c) = w(:,c)
          w(:,c)  = (v(:,c) - oldeps*w1(:,c) - delta*w2(:,c)) / gam
          T(:,c)  = T(:,c) + phi*w(:,c)
          if (beta(c) <= 0.0_dp .or. phibar(c) <= 1.0e-14_dp) done(c) = .true.
       enddo
    enddo
    call project(T)

    deallocate(r1, r2, y, v, w, w1, w2, b, kd, az, th, oldb, beta, dbar, epsln, phibar, &
         cs, sn, alfa, done)

  contains

    subroutine project(Z)
      ! Z <- (1 - R R^T) Z
      real(dp), intent(inout) :: Z(ndim,ncol)
      real(dp), allocatable   :: o(:,:)
      allocate(o(nsolve,ncol))
      call dgemm('T','N', nsolve, ncol, ndim, 1.0_dp, R, ndim, Z, ndim, 0.0_dp, o, nsolve)
      call dgemm('N','N', ndim, ncol, nsolve, -1.0_dp, R, ndim, o, nsolve, 1.0_dp, Z, ndim)
      deallocate(o)
    end subroutine project

    subroutine apply_m(Z, MZ)
      ! MZ = P (A0 Z) - theta Z   for Z in the range of P
      real(dp), intent(in)  :: Z(ndim,ncol)
      real(dp), intent(out) :: MZ(ndim,ncol)
      real(dp)              :: t0
      integer(is)           :: cc
      t0 = wall_time()
      kernel_off = .true.
      call sigma_batch(ncol, mult, irrep, Z, MZ)
      kernel_off = .false.
      t0 = wall_time() - t0
      nsigma_vecs  = nsigma_vecs - ncol
      nsigma_calls = nsigma_calls - 1
      time_sigma   = time_sigma - t0
      ninner_vecs  = ninner_vecs + ncol
      time_inner   = time_inner + t0
      call project(MZ)
      do cc = 1, ncol
         MZ(:,cc) = MZ(:,cc) - th(cc)*Z(:,cc)
      enddo
    end subroutine apply_m

  end subroutine jd_correction

end module mrsf_davidson
