!**********************************************************************
! grid_wderiv: nuclear derivatives of the partition weights contracted
! with per-point integrands,
!   grad(:,G,c) = sum_p q(p,c) vol_p dW_i(r_p)/dR_G,
! the point moving rigidly with its parent atom i. With r fixed,
! dW_i/dR_G for G /= i follows from the chain rule through the
! mu_lm of every (cell l, factor m) pair of the weight expression,
!   dmu_lm/dR_l = -(u_l + mu e_lm)/d_lm,  dmu_lm/dR_m = (u_m + mu e_lm)/d_lm,
! and the parent term from translational invariance,
!   dW_i/dR_i = -sum_{G /= i} dW_i/dR_G.
! Cells below cell_cut contribute nothing (as in grid_partition).
!**********************************************************************
module grid_wderiv

  use grid_constants
  use grid_tables
  use grid_step
  use grid_partition

  implicit none

contains

  subroutine wderiv_contract(natm, R, Z, adjust, step, mu, cells, weights, npts, xyz, iatom, vol, &
       nnb, nb, nq, q, grad, ierr)

    integer(is), intent(in)  :: natm, Z(natm), adjust, step, cells, weights, npts, nq
    integer(is), intent(in)  :: iatom(npts), nnb(natm), nb(max(natm-1,1),natm)
    real(dp), intent(in)     :: R(3,natm), mu, xyz(3,npts), vol(npts), q(npts,nq)
    real(dp), intent(out)    :: grad(3,natm,nq)
    integer(is), intent(out) :: ierr
    real(dp), allocatable    :: a(:,:), dinv(:,:), dist(:), gacc(:,:,:), dW(:,:), s(:), pre(:), suf(:), cl(:)
    integer(is), allocatable :: stamp(:), slot(:), loc(:), mlist(:)
    integer(is)              :: ld, p, i, k, ip, l, ll, nl, m, mm, nm, nloc, g, c, maxm
    real(dp)                 :: px, py, pz, ci, ssum, rl, rm, muv, nu, sv, ds, coef, t, dmu, ul(3), um(3), e(3)
    real(dp)                 :: tsum(3), fac
    logical                  :: local

    ierr = 0
    grad = 0.0_dp
    call step_check(step, mu, ierr)
    if (ierr /= 0) return
    if ((cells /= mode_becke .and. cells /= mode_delaunay) .or. &
        (weights /= mode_becke .and. weights /= mode_delaunay)) then
       ierr = 45
       return
    endif
    if (natm == 1) return
    ld = max(natm-1, 1)
    local = (cells == mode_delaunay .or. weights == mode_delaunay)
    if (local) then
       if (any(nnb < 0) .or. any(nnb > natm-1)) then
          ierr = 46
          return
       endif
    endif
    allocate(a(natm,natm), dinv(natm,natm))
    call partition_tables(natm, R, Z, adjust, a, dinv, ierr)
    if (ierr /= 0) return
    if (cells == mode_delaunay) then
       maxm = max(maxval(nnb), 1)
    else
       maxm = natm
    endif

    !$omp parallel private(dist, stamp, slot, loc, mlist, gacc, dW, s, pre, suf, cl, ip, p, i, k, px, py, pz, &
    !$omp   ci, ssum, l, ll, nl, m, mm, nm, nloc, g, c, rl, rm, muv, nu, sv, ds, coef, t, dmu, ul, um, e, tsum, fac)
    allocate(dist(natm), stamp(natm), slot(natm), loc(natm), mlist(maxm), gacc(3,natm,nq), dW(3,natm))
    allocate(s(maxm), pre(0:maxm), suf(maxm+1), cl(natm))
    stamp = 0
    slot = 0
    gacc = 0.0_dp
    ip = 0
    !$omp do schedule(dynamic, 128)
    do p = 1, npts
       i = iatom(p)
       if (i < 1 .or. i > natm) cycle
       if (all(q(p,:) == 0.0_dp)) cycle
       ip = ip + 1
       px = xyz(1,p); py = xyz(2,p); pz = xyz(3,p)
       k = i
       if (local) k = nearest_walk(natm, R, ld, nnb, nb, px, py, pz, i, dist, stamp, ip)
       if (weights == mode_delaunay) then
          if (k /= i) then
             if (.not. in_list(i, nnb(k), nb(1,k))) cycle
          endif
          nl = nnb(k) + 1
       else
          nl = natm
       endif
       ! cell values, the sum and the weight
       ssum = 0.0_dp
       ci = 0.0_dp
       do ll = 1, nl
          if (weights == mode_delaunay) then
             if (ll == 1) then
                l = k
             else
                l = nb(ll-1,k)
             endif
          else
             l = ll
          endif
          cl(ll) = cell_value(l, natm, R, a, dinv, step, mu, cells, ld, nnb, nb, px, py, pz, dist, stamp, ip)
          ssum = ssum + cl(ll)
          if (l == i) ci = cl(ll)
       enddo
       if (ssum <= 0.0_dp) cycle
       ! derivatives at fixed r, accumulated per touched atom (slots)
       nloc = 0
       do ll = 1, nl
          if (weights == mode_delaunay) then
             if (ll == 1) then
                l = k
             else
                l = nb(ll-1,k)
             endif
          else
             l = ll
          endif
          if (cl(ll) <= 0.0_dp) cycle
          coef = -ci/ssum/ssum
          if (l == i) coef = coef + 1.0_dp/ssum
          rl = point_dist(l, natm, R, px, py, pz, dist, stamp, ip)
          ul = ((/ px, py, pz /) - R(:,l))/rl
          ! factor list of cell l
          if (cells == mode_delaunay) then
             nm = nnb(l)
             mlist(1:nm) = nb(1:nm,l)
          else
             nm = 0
             do m = 1, natm
                if (m == l) cycle
                nm = nm + 1
                mlist(nm) = m
             enddo
          endif
          pre(0) = 1.0_dp
          do mm = 1, nm
             m = mlist(mm)
             rm = point_dist(m, natm, R, px, py, pz, dist, stamp, ip)
             muv = (rl - rm)*dinv(m,l)
             nu = muv + a(l,m)*(1.0_dp - muv*muv)
             s(mm) = step_eval(step, mu, nu)
             pre(mm) = pre(mm-1)*s(mm)
          enddo
          suf(nm+1) = 1.0_dp
          do mm = nm, 1, -1
             suf(mm) = suf(mm+1)*s(mm)
          enddo
          do mm = 1, nm
             m = mlist(mm)
             rm = dist(m)
             muv = (rl - rm)*dinv(m,l)
             nu = muv + a(l,m)*(1.0_dp - muv*muv)
             call step_eval_der(step, mu, nu, sv, ds)
             t = coef*pre(mm-1)*suf(mm+1)*ds*(1.0_dp - 2.0_dp*a(l,m)*muv)
             if (t == 0.0_dp) cycle
             e = (R(:,l) - R(:,m))*dinv(m,l)
             um = ((/ px, py, pz /) - R(:,m))/rm
             if (l /= i) then
                if (slot(l) == 0) then
                   nloc = nloc + 1
                   loc(nloc) = l
                   slot(l) = nloc
                   dW(:,nloc) = 0.0_dp
                endif
                dW(:,slot(l)) = dW(:,slot(l)) - t*(ul + muv*e)*dinv(m,l)
             endif
             if (m /= i) then
                if (slot(m) == 0) then
                   nloc = nloc + 1
                   loc(nloc) = m
                   slot(m) = nloc
                   dW(:,nloc) = 0.0_dp
                endif
                dW(:,slot(m)) = dW(:,slot(m)) + t*(um + muv*e)*dinv(m,l)
             endif
          enddo
       enddo
       ! scatter: G /= i directly, the parent by translational invariance
       tsum = 0.0_dp
       do g = 1, nloc
          tsum = tsum + dW(:,g)
       enddo
       do c = 1, nq
          fac = vol(p)*q(p,c)
          if (fac == 0.0_dp) cycle
          do g = 1, nloc
             gacc(:,loc(g),c) = gacc(:,loc(g),c) + fac*dW(:,g)
          enddo
          gacc(:,i,c) = gacc(:,i,c) - fac*tsum
       enddo
       do g = 1, nloc
          slot(loc(g)) = 0
       enddo
    enddo
    !$omp end do
    !$omp critical
    grad = grad + gacc
    !$omp end critical
    deallocate(dist, stamp, slot, loc, mlist, gacc, dW, s, pre, suf, cl)
    !$omp end parallel
    deallocate(a, dinv)

  end subroutine wderiv_contract

end module grid_wderiv
