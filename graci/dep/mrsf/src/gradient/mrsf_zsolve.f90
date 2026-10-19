!**********************************************************************
! mrsf_zsolve: the Z-vector equations H x = b of the MRSF gradients and
! derivative couplings, all right-hand sides of a gradient section at
! once (H = the ROKS orbital Hessian, applied by zvec_apply).
!
! Method: one shared subspace for all columns, minimal residual. The
! basis V is built from the preconditioned residuals of the columns that
! are not yet converged (orthogonalised twice against V and among
! themselves; dependent directions dropped); AV is kept as an incremental
! QR factorisation AV = Q Rm (classical Gram-Schmidt twice), so that for
! every column y = Rm^-1 Q^T r0 minimises |r0 - AV y|_2 over the whole
! subspace, x = x0 + V y and r = r0 - Q Q^T r0. One operator product per
! new basis vector; converged columns are locked. Unlike conjugate
! gradients the method does not assume a definite operator: at
! displaced geometries the orbital Hessian of the reference can have a
! negative eigenvalue (the ROKS solution is then a saddle point), and it
! is symmetric only up to the SCF orbital gradient.
!
! Preconditioner: the Fock diagonal of the operator (|d| >= 1e-4), or
! (precond = 2, with the XC kernel) ninner steps of Jacobi-preconditioned
! conjugate gradients on the kernel-free operator (Fock + J - cx K, no
! XC kernel; several times cheaper), started from zero, per direction.
! When the basis would exceed maxsub vectors the method restarts from
! the current solutions.
!
! Diagnostics: the Ritz values of D^-1 H in the final subspace (lowest
! two and the number of negative ones) are returned in info.
!**********************************************************************
module mrsf_zsolve

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_zvec

  implicit none

  real(dp)    :: time_zsolve = 0.0_dp, time_zprod = 0.0_dp, time_zinner = 0.0_dp
  integer(is) :: nzs_full = 0, nzs_inner = 0

  real(dp), external :: dnrm2, ddot
  external :: dgemm, dgemv, dtrsv, dsygv, daxpy, dscal

contains

!######################################################################
! zvec_solve: B(lz, ncol) right-hand sides; X: start vectors if nx0 = 1
! (zero columns: none), solutions on return; tol: |r|_2 per column;
! maxit: expansion steps; precond: 1 diagonal, 2 kernel-free inner
! steps (ninner); maxsub: largest basis before a restart.
! Outputs: niter (steps while active), resid (final |r|), hist(it+1, j)
! (residual after step it, 0 once converged), info: 1 full products,
! 2 kernel-free products, 3 restarts, 4 basis size, 5-6 lowest Ritz
! values of D^-1 H, 7 number of negative Ritz values, 8 seconds in the
! operator, 9 seconds in the inner steps, 10 total seconds
!######################################################################
  subroutine zvec_solve(ncol, cx, have_xc, B, X, nx0, tol, maxit, precond, ninner, maxsub, &
       niter, resid, hist, info, ierr)

    integer(is), intent(in)  :: ncol, nx0, maxit, precond, ninner, maxsub
    real(dp), intent(in)     :: cx, tol, B(lz,ncol)
    logical, intent(in)      :: have_xc
    real(dp), intent(inout)  :: X(lz,ncol)
    integer(is), intent(out) :: niter(ncol), ierr
    real(dp), intent(out)    :: resid(ncol), hist(maxit+1,ncol), info(16)

    real(dp), allocatable    :: d(:), X0(:,:), R0(:,:), R(:,:), V(:,:), Q(:,:), Rm(:,:)
    real(dp), allocatable    :: W(:,:), AW(:,:), G(:,:), c(:), y(:)
    integer(is), allocatable :: act(:), idx(:)
    logical, allocatable     :: conv(:)
    integer(is) :: n, m, it, j, k, kk, na, nnew, mmax, nrestart, pass, nz, col
    real(dp)    :: t0, t1, nrm0, nrm, h
    logical     :: inner

    t0 = wall_time()
    ierr = 0
    info = 0.0_dp
    if (lz == 0) then
       ierr = 1
       return
    endif
    n = lz
    mmax = max(maxsub, 2*ncol)
    inner = (precond == 2) .and. have_xc .and. (ninner > 0)
    allocate(d(n), X0(n,ncol), R0(n,ncol), R(n,ncol), V(n,mmax), Q(n,mmax), Rm(mmax,mmax), &
         act(ncol), idx(ncol), conv(ncol), c(mmax), y(mmax))
    d = hdg
    where (abs(d) < 1.0e-4_dp) d = sign(1.0e-4_dp, d)
    nzs_full = 0; nzs_inner = 0; time_zprod = 0.0_dp; time_zinner = 0.0_dp

    ! start vectors and initial residuals
    X0 = 0.0_dp
    R0 = B
    if (nx0 > 0) then
       nz = 0
       do j = 1, ncol
          if (dnrm2(n, X(1,j), 1) > 0.0_dp) then
             nz = nz + 1
             idx(nz) = j
          endif
       enddo
       if (nz > 0) then
          allocate(W(n,nz), AW(n,nz))
          do k = 1, nz
             W(:,k) = X(:,idx(k))
          enddo
          call apply_full(nz, W, AW)
          do k = 1, nz
             X0(:,idx(k)) = W(:,k)
             call daxpy(n, -1.0_dp, AW(1,k), 1, R0(1,idx(k)), 1)
          enddo
          deallocate(W, AW)
       endif
    endif
    X = X0
    R = R0
    hist = 0.0_dp
    do j = 1, ncol
       resid(j) = dnrm2(n, R(1,j), 1)
       hist(1,j) = resid(j)
       conv(j) = resid(j) < tol
       niter(j) = 0
    enddo

    m = 0
    nrestart = 0
    do it = 1, maxit
       na = 0
       do j = 1, ncol
          if (.not. conv(j)) then
             na = na + 1
             act(na) = j
          endif
       enddo
       if (na == 0) exit
       ! preconditioned residuals of the active columns
       allocate(W(n,na))
       do k = 1, na
          W(:,k) = R(:,act(k))
       enddo
       call precondition(na, W, d, inner, ninner, cx)
       ! orthogonalisation against the basis (twice) and within the block
       if (m > 0) then
          allocate(G(m,na))
          do pass = 1, 2
             call dgemm('T','N', m, na, n, 1.0_dp, V, n, W, n, 0.0_dp, G, m)
             call dgemm('N','N', n, na, m, -1.0_dp, V, n, G, m, 1.0_dp, W, n)
          enddo
          deallocate(G)
       endif
       nnew = 0
       do k = 1, na
          nrm0 = dnrm2(n, W(1,k), 1)
          if (k /= nnew + 1) W(:,nnew+1) = W(:,k)
          do pass = 1, 2
             do kk = 1, nnew
                h = ddot(n, W(1,kk), 1, W(1,nnew+1), 1)
                call daxpy(n, -h, W(1,kk), 1, W(1,nnew+1), 1)
             enddo
          enddo
          nrm = dnrm2(n, W(1,nnew+1), 1)
          if (nrm > 1.0e-8_dp*nrm0 .and. nrm > 1.0e-300_dp) then
             nnew = nnew + 1
             call dscal(n, 1.0_dp/nrm, W(1,nnew), 1)
          endif
       enddo
       if (nnew == 0) then
          deallocate(W)
          exit
       endif
       ! restart from the current solutions when the basis is full
       if (m + nnew > mmax) then
          X0 = X
          R0 = R
          m = 0
          nrestart = nrestart + 1
       endif
       ! operator products and the incremental QR of AV
       allocate(AW(n,nnew))
       call apply_full(nnew, W, AW)
       do k = 1, nnew
          col = m + k
          V(:,col) = W(:,k)
          Rm(1:col,col) = 0.0_dp
          if (col > 1) then
             do pass = 1, 2
                call dgemv('T', n, col-1, 1.0_dp, Q, n, AW(1,k), 1, 0.0_dp, c, 1)
                call dgemv('N', n, col-1, -1.0_dp, Q, n, c, 1, 1.0_dp, AW(1,k), 1)
                Rm(1:col-1,col) = Rm(1:col-1,col) + c(1:col-1)
             enddo
          endif
          h = dnrm2(n, AW(1,k), 1)
          if (h < 1.0e-300_dp) h = 1.0e-300_dp
          Rm(col,col) = h
          Q(:,col) = AW(:,k)/h
       enddo
       m = m + nnew
       deallocate(W, AW)
       ! minimal-residual update of the active columns
       t1 = wall_time()
       do k = 1, na
          j = act(k)
          call dgemv('T', n, m, 1.0_dp, Q, n, R0(1,j), 1, 0.0_dp, c, 1)
          y(1:m) = c(1:m)
          call dtrsv('U', 'N', 'N', m, Rm, mmax, y, 1)
          X(:,j) = X0(:,j)
          call dgemv('N', n, m, 1.0_dp, V, n, y, 1, 1.0_dp, X(1,j), 1)
          R(:,j) = R0(:,j)
          call dgemv('N', n, m, -1.0_dp, Q, n, c, 1, 1.0_dp, R(1,j), 1)
          resid(j) = dnrm2(n, R(1,j), 1)
          niter(j) = it
          hist(it+1,j) = resid(j)
          if (resid(j) < tol) conv(j) = .true.
       enddo
    enddo

    ! Ritz values of D^-1 H in the final subspace: T y = theta M y with
    ! T = V^T A V = (V^T Q) Rm (symmetrised) and M = V^T D V
    info(4) = real(m, dp)
    if (m > 0) call ritz_diag(n, m, mmax, V, Q, Rm, d, info)
    info(1) = real(nzs_full, dp)
    info(2) = real(nzs_inner, dp)
    info(3) = real(nrestart, dp)
    info(8) = time_zprod
    info(9) = time_zinner
    info(10) = wall_time() - t0
    time_zsolve = time_zsolve + info(10)
    deallocate(d, X0, R0, R, V, Q, Rm, act, idx, conv, c, y)

  contains

    subroutine apply_full(nv, Z, AZ)
      integer(is), intent(in) :: nv
      real(dp), intent(in)    :: Z(n,nv)
      real(dp), intent(out)   :: AZ(n,nv)
      real(dp) :: ta
      ta = wall_time()
      call zvec_apply(nv, cx, have_xc, Z, AZ)
      nzs_full = nzs_full + nv
      time_zprod = time_zprod + wall_time() - ta
    end subroutine apply_full

  end subroutine zvec_solve

!######################################################################
! precondition: W (lz, nv) residuals -> preconditioned directions: the
! diagonal, or ninner steps of Jacobi-PCG on the kernel-free operator
! (from zero; a column stops at non-positive curvature, and falls back
! to the diagonal step if that happens at the first step)
!######################################################################
  subroutine precondition(nv, W, d, inner, ninner, cx)

    integer(is), intent(in) :: nv, ninner
    real(dp), intent(inout) :: W(lz,nv)
    real(dp), intent(in)    :: d(lz), cx
    logical, intent(in)     :: inner

    real(dp), allocatable :: Y(:,:), Rr(:,:), Zz(:,:), P(:,:), AP(:,:), rz(:)
    logical, allocatable  :: done(:), moved(:)
    integer(is) :: k, step
    real(dp)    :: pap, alpha, rzn, ta

    if (.not. inner) then
       do k = 1, nv
          W(:,k) = W(:,k)/d
       enddo
       return
    endif
    ta = wall_time()
    allocate(Y(lz,nv), Rr(lz,nv), Zz(lz,nv), P(lz,nv), AP(lz,nv), rz(nv), done(nv), moved(nv))
    Y = 0.0_dp
    Rr = W
    do k = 1, nv
       Zz(:,k) = Rr(:,k)/d
       rz(k) = ddot(lz, Rr(1,k), 1, Zz(1,k), 1)
    enddo
    P = Zz
    done = .false.
    moved = .false.
    do step = 1, ninner
       call zvec_apply(nv, cx, .false., P, AP)
       nzs_inner = nzs_inner + nv
       do k = 1, nv
          if (done(k)) cycle
          pap = ddot(lz, P(1,k), 1, AP(1,k), 1)
          if (.not. (pap > 0.0_dp) .or. abs(rz(k)) < 1.0e-300_dp) then
             done(k) = .true.
             cycle
          endif
          alpha = rz(k)/pap
          call daxpy(lz, alpha, P(1,k), 1, Y(1,k), 1)
          call daxpy(lz, -alpha, AP(1,k), 1, Rr(1,k), 1)
          moved(k) = .true.
          if (step == ninner) cycle
          Zz(:,k) = Rr(:,k)/d
          rzn = ddot(lz, Rr(1,k), 1, Zz(1,k), 1)
          P(:,k) = Zz(:,k) + (rzn/rz(k))*P(:,k)
          rz(k) = rzn
       enddo
       if (all(done)) exit
    enddo
    do k = 1, nv
       if (moved(k)) then
          W(:,k) = Y(:,k)
       else
          W(:,k) = W(:,k)/d
       endif
    enddo
    deallocate(Y, Rr, Zz, P, AP, rz, done, moved)
    time_zinner = time_zinner + wall_time() - ta

  end subroutine precondition

!######################################################################
! ritz_diag: lowest Ritz values of D^-1 H in the subspace V
!######################################################################
  subroutine ritz_diag(n, m, mmax, V, Q, Rm, d, info)

    integer(is), intent(in) :: n, m, mmax
    real(dp), intent(in)    :: V(n,mmax), Q(n,mmax), Rm(mmax,mmax), d(n)
    real(dp), intent(inout) :: info(16)

    real(dp), allocatable :: VQ(:,:), T(:,:), Mm(:,:), DV(:,:), w(:), work(:), Ru(:,:)
    integer(is) :: k, lwork, ierr

    allocate(VQ(m,m), T(m,m), Mm(m,m), DV(n,m), w(m), Ru(m,m))
    Ru = 0.0_dp
    do k = 1, m
       Ru(1:k,k) = Rm(1:k,k)
    enddo
    call dgemm('T','N', m, m, n, 1.0_dp, V, n, Q, n, 0.0_dp, VQ, m)
    call dgemm('N','N', m, m, m, 1.0_dp, VQ, m, Ru, m, 0.0_dp, T, m)
    T = 0.5_dp*(T + transpose(T))
    do k = 1, m
       DV(:,k) = d*V(:,k)
    enddo
    call dgemm('T','N', m, m, n, 1.0_dp, V, n, DV, n, 0.0_dp, Mm, m)
    lwork = max(1, 3*m)
    allocate(work(lwork))
    call dsygv(1, 'N', 'U', m, T, m, Mm, m, w, work, lwork, ierr)
    if (ierr == 0) then
       info(5) = w(1)
       if (m > 1) info(6) = w(2)
       info(7) = real(count(w < 0.0_dp), dp)
    endif
    deallocate(VQ, T, Mm, DV, w, work, Ru)

  end subroutine ritz_diag

end module mrsf_zsolve

!**********************************************************************
! C-bound entry point
!**********************************************************************
module mrsf_zsolve_interface

  use iso_c_binding
  use mrsf_constants
  use mrsf_global
  use mrsf_zvec, only: lz
  use mrsf_zsolve

  implicit none

contains

  subroutine mrsf_zvec_solve_c(ncol, cx, have_xc, B, X, nx0, tol, maxit, precond, ninner, maxsub, &
       niter, resid, hist, info, ierr) bind(c, name='mrsf_zvec_solve')
    integer(is), intent(in)     :: ncol, nx0, maxit, precond, ninner, maxsub
    real(dp), intent(in)        :: cx, tol, B(lz,ncol)
    logical(c_bool), intent(in) :: have_xc
    real(dp), intent(inout)     :: X(lz,ncol)
    integer(is), intent(out)    :: niter(ncol), ierr
    real(dp), intent(out)       :: resid(ncol), hist(maxit+1,ncol), info(16)
    call zvec_solve(ncol, cx, logical(have_xc), B, X, nx0, tol, maxit, precond, ninner, maxsub, &
         niter, resid, hist, info, ierr)
  end subroutine mrsf_zvec_solve_c

end module mrsf_zsolve_interface
