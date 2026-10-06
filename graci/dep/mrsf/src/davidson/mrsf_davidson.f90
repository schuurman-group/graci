!**********************************************************************
! mrsf_davidson: in-core block Davidson solver for the lowest roots of
! the (symmetric) MRSF-TDDFT A matrix of one (multiplicity, irrep)
! block. Exact-diagonal preconditioner, locking of converged roots,
! re-orthogonalised Gram-Schmidt expansion, Ritz-vector collapse.
!**********************************************************************
module mrsf_davidson

  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_space
  use mrsf_sigma

  implicit none

contains

!######################################################################
! davidson_solve
!   nroots : requested roots (reduced if the active space is smaller)
!   nextra : additional tracked roots
!   maxvec : maximum subspace dimension (raised to >= 2*nsolve)
!   tol    : convergence threshold on the residual norm
!   ener   : eigenvalues (omega, relative to the reference energy)
!   vecs   : compressed eigenvectors
!   iconv  : 1 if converged, 0 otherwise
!######################################################################
  subroutine davidson_solve(mult, irrep, nroots, nextra, maxvec_in, maxiter, tol, &
       ener, vecs, niter, iconv)

    integer(is), intent(in)    :: mult, irrep, nextra, maxvec_in, maxiter
    integer(is), intent(inout) :: nroots
    real(dp), intent(in)       :: tol
    real(dp), intent(out)      :: ener(*), vecs(xdim,*)
    integer(is), intent(out)   :: niter, iconv

    integer(is), allocatable :: act(:), idx(:)
    real(dp), allocatable    :: d(:), V(:,:), AV(:,:), G(:,:), Gc(:,:), theta(:), work(:)
    real(dp), allocatable    :: R(:,:), AR(:,:), res(:,:), rnorm(:), dnew(:,:), ovl(:,:)
    real(dp), allocatable    :: dtmp(:)
    integer(is)              :: nact, nsolve, maxvec, currdim, ist, iend, nnew, nkept
    integer(is)              :: ia, a, i, k, l, iter, lwork, info, kmin
    real(dp)                 :: denom, rmax, dnorm, proj, dmin
    logical                  :: collapsed

    iconv = 0
    niter = 0

    ! active slots
    allocate(act(xdim))
    nact = 0
    do i = 1, nocca
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

    allocate(d(xdim), V(xdim,maxvec), AV(xdim,maxvec), G(maxvec,maxvec), Gc(maxvec,maxvec))
    allocate(theta(maxvec), R(xdim,nsolve), AR(xdim,nsolve), res(xdim,nsolve))
    allocate(rnorm(nsolve), dnew(xdim,nsolve), ovl(maxvec,nsolve), dtmp(nact), idx(nsolve))
    V = 0.0_dp; AV = 0.0_dp; G = 0.0_dp

    ! diagonal and preconditioner
    call diagonal(mult, d)
    do ia = 1, xdim
       if (slot_irrep(ia) /= irrep .and. irrep >= 0) d(ia) = 1.0e20_dp
    enddo

    ! guess: unit vectors on the lowest active diagonal elements
    do k = 1, nact
       dtmp(k) = d(act(k))
    enddo
    do k = 1, nsolve
       kmin = 1
       dmin = dtmp(1)
       do l = 2, nact
          if (dtmp(l) < dmin) then
             dmin = dtmp(l)
             kmin = l
          endif
       enddo
       idx(k) = act(kmin)
       dtmp(kmin) = huge(1.0_dp)
       V(idx(k),k) = 1.0_dp
    enddo

    ! LAPACK workspace
    lwork = max(1_is, 3*maxvec + 64)
    allocate(work(lwork))

    if (verbose) then
       write(6,'(/,2x,a,i0,a,i0,a,i0,a,i0)') 'Davidson: mult = ', mult, &
            ', irrep = ', irrep, ', active dim = ', nact, ', roots = ', nroots
       write(6,'(2x,a)') ' iter   dim   max|r|         lowest eigenvalues'
    endif

    currdim = nsolve
    ist = 1
    iend = nsolve

    do iter = 1, maxiter
       niter = iter
       nnew = iend - ist + 1

       ! sigma vectors for the new subspace vectors
       call sigma_batch(nnew, mult, irrep, V(:,ist:iend), AV(:,ist:iend))

       ! subspace matrix update
       call dgemm('T','N', currdim, nnew, xdim, 1.0_dp, V, xdim, AV(1,ist), xdim, &
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
       call dgemm('N','N', xdim, nsolve, currdim, 1.0_dp, V, xdim, Gc, maxvec, 0.0_dp, R, xdim)
       call dgemm('N','N', xdim, nsolve, currdim, 1.0_dp, AV, xdim, Gc, maxvec, 0.0_dp, AR, xdim)
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

       ! preconditioned corrections for the unconverged tracked roots
       nnew = 0
       do k = 1, nsolve
          if (rnorm(k) < tol) cycle
          nnew = nnew + 1
          dnew(:,nnew) = 0.0_dp
          do l = 1, nact
             ia = act(l)
             denom = d(ia) - theta(k)
             if (abs(denom) < 1.0e-4_dp) denom = sign(1.0e-4_dp, denom)
             dnew(ia,nnew) = -res(ia,k) / denom
          enddo
       enddo
       if (nnew == 0) exit

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
          call dgemm('T','N', currdim, nnew, xdim, 1.0_dp, V, xdim, dnew, xdim, 0.0_dp, ovl, maxvec)
          call dgemm('N','N', xdim, nnew, currdim, -1.0_dp, V, xdim, ovl, maxvec, 1.0_dp, dnew, xdim)
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

    ener(1:nroots) = theta(1:nroots)
    vecs(1:xdim,1:nroots) = R(:,1:nroots)

    if (verbose) then
       if (iconv == 1) then
          write(6,'(2x,a,i0,a)') 'Davidson converged in ', niter, ' iterations'
       else
          write(6,'(2x,a,i0,a)') 'Davidson NOT converged after ', niter, ' iterations'
       endif
    endif

    deallocate(act, idx, d, V, AV, G, Gc, theta, work, R, AR, res, rnorm, dnew, ovl, dtmp)

  end subroutine davidson_solve

end module mrsf_davidson
