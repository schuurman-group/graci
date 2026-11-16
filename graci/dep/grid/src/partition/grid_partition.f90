!**********************************************************************
! grid_partition: partition-of-unity weights of the grid points
! (Becke 1988 and the local Delaunay cells / weights of Gill &
! Tsoukalas, JCP 165, 174104 (2026), Algorithm 1):
!   cells   = mode_becke:    C_l = prod_{m /= l} s_lm
!             mode_delaunay: C_l = prod_{m in D_l} s_lm
!   weights = mode_becke:    W_i = C_i / sum_l C_l (all l)
!             mode_delaunay: W_i = C_i / sum_{l in D_k^+} C_l, k the
!                            nearest nucleus of the point, W_i = 0 when
!                            i is not in D_k^+
! w = vol W_i. The cell products stop once below cut (eq. 10 of the
! paper: the product can only decrease). Thread-private distance
! caches are passed explicitly (no host association inside OpenMP).
!**********************************************************************
module grid_partition

  use grid_constants
  use grid_tables
  use grid_step

  implicit none

  real(dp), parameter :: cell_cut = 1.0e-20_dp

contains

!######################################################################
! point_dist: |r - R_m|, computed once per point (stamp cache)
!######################################################################
  function point_dist(m, natm, R, x, y, z, dist, stamp, ip) result(dm)

    integer(is), intent(in)    :: m, natm, ip
    real(dp), intent(in)       :: R(3,natm), x, y, z
    real(dp), intent(inout)    :: dist(natm)
    integer(is), intent(inout) :: stamp(natm)
    real(dp)                   :: dm

    if (stamp(m) /= ip) then
       dist(m) = sqrt((x - R(1,m))**2 + (y - R(2,m))**2 + (z - R(3,m))**2)
       stamp(m) = ip
    endif
    dm = dist(m)

  end function point_dist

!######################################################################
! nearest_walk: nearest nucleus to (x,y,z) by walking the neighbour
! lists from k0
!######################################################################
  function nearest_walk(natm, R, ld, nnb, nb, x, y, z, k0, dist, stamp, ip) result(k)

    integer(is), intent(in)    :: natm, ld, nnb(natm), nb(ld,natm), k0, ip
    real(dp), intent(in)       :: R(3,natm), x, y, z
    real(dp), intent(inout)    :: dist(natm)
    integer(is), intent(inout) :: stamp(natm)
    integer(is)                :: k, kbest, jj, j, it
    real(dp)                   :: dk, dj

    k = k0
    dk = point_dist(k, natm, R, x, y, z, dist, stamp, ip)
    do it = 1, natm
       kbest = k
       do jj = 1, nnb(k)
          j = nb(jj,k)
          dj = point_dist(j, natm, R, x, y, z, dist, stamp, ip)
          if (dj < dk) then
             dk = dj
             kbest = j
          endif
       enddo
       if (kbest == k) exit
       k = kbest
    enddo

  end function nearest_walk

!######################################################################
! in_list: i in nb(1:n)
!######################################################################
  pure function in_list(i, n, nb) result(found)

    integer(is), intent(in) :: i, n, nb(*)
    logical                 :: found
    integer(is)             :: jj

    found = .false.
    do jj = 1, n
       if (nb(jj) == i) then
          found = .true.
          return
       endif
    enddo

  end function in_list

!######################################################################
! cell_value: C_l at the point
!######################################################################
  function cell_value(l, natm, R, a, dinv, step, mu, cells, ld, nnb, nb, x, y, z, dist, stamp, ip) result(cl)

    integer(is), intent(in)    :: l, natm, step, cells, ld, nnb(natm), nb(ld,natm), ip
    real(dp), intent(in)       :: R(3,natm), a(natm,natm), dinv(natm,natm), mu, x, y, z
    real(dp), intent(inout)    :: dist(natm)
    integer(is), intent(inout) :: stamp(natm)
    real(dp)                   :: cl, rl, rm, muv, nu
    integer(is)                :: m, mm, nm

    cl = 1.0_dp
    rl = point_dist(l, natm, R, x, y, z, dist, stamp, ip)
    if (cells == mode_delaunay) then
       nm = nnb(l)
    else
       nm = natm
    endif
    do mm = 1, nm
       if (cells == mode_delaunay) then
          m = nb(mm,l)
       else
          m = mm
          if (m == l) cycle
       endif
       rm = point_dist(m, natm, R, x, y, z, dist, stamp, ip)
       muv = (rl - rm)*dinv(m,l)
       nu = muv + a(l,m)*(1.0_dp - muv*muv)
       cl = cl*step_eval(step, mu, nu)
       if (cl < cell_cut) then
          cl = 0.0_dp
          return
       endif
    enddo

  end function cell_value

!######################################################################
! partition_tables: size adjustment and inverse distances
!######################################################################
  subroutine partition_tables(natm, R, Z, adjust, a, dinv, ierr)

    integer(is), intent(in)  :: natm, Z(natm), adjust
    real(dp), intent(in)     :: R(3,natm)
    real(dp), intent(out)    :: a(natm,natm), dinv(natm,natm)
    integer(is), intent(out) :: ierr
    integer(is)              :: i, j

    call adjust_table(natm, Z, adjust, a, ierr)
    if (ierr /= 0) return
    dinv = 0.0_dp
    do j = 1, natm
       do i = 1, natm
          if (i == j) cycle
          dinv(i,j) = 1.0_dp/sqrt(sum((R(:,i) - R(:,j))**2))
       enddo
    enddo

  end subroutine partition_tables

!######################################################################
! partition_weights
!######################################################################
  subroutine partition_weights(natm, R, Z, adjust, step, mu, cells, weights, npts, xyz, iatom, vol, &
       nnb, nb, w, knear, ierr)

    integer(is), intent(in)  :: natm, Z(natm), adjust, step, cells, weights, npts
    integer(is), intent(in)  :: iatom(npts), nnb(natm), nb(max(natm-1,1),natm)
    real(dp), intent(in)     :: R(3,natm), mu, xyz(3,npts), vol(npts)
    real(dp), intent(out)    :: w(npts)
    integer(is), intent(out) :: knear(npts), ierr
    real(dp), allocatable    :: a(:,:), dinv(:,:), dist(:)
    integer(is), allocatable :: stamp(:)
    integer(is)              :: ld, p, i, k, ip, l, ll, nl
    real(dp)                 :: px, py, pz, ci, ssum, cl
    logical                  :: local

    ierr = 0
    call step_check(step, mu, ierr)
    if (ierr /= 0) return
    if ((cells /= mode_becke .and. cells /= mode_delaunay) .or. &
        (weights /= mode_becke .and. weights /= mode_delaunay)) then
       ierr = 45
       return
    endif
    ld = max(natm-1, 1)
    local = (cells == mode_delaunay .or. weights == mode_delaunay)
    if (local) then
       if (any(nnb < 0) .or. any(nnb > natm-1)) then
          ierr = 46
          return
       endif
    endif
    if (natm == 1) then
       do p = 1, npts
          if (iatom(p) == 1) then
             w(p) = vol(p)
             knear(p) = 1
          else
             w(p) = 0.0_dp
             knear(p) = 0
          endif
       enddo
       return
    endif
    allocate(a(natm,natm), dinv(natm,natm))
    call partition_tables(natm, R, Z, adjust, a, dinv, ierr)
    if (ierr /= 0) return

    !$omp parallel private(dist, stamp, ip, p, i, k, px, py, pz, ci, ssum, cl, l, ll, nl)
    allocate(dist(natm), stamp(natm))
    stamp = 0
    ip = 0
    !$omp do schedule(dynamic, 256)
    do p = 1, npts
       i = iatom(p)
       if (i < 1 .or. i > natm) then
          w(p) = 0.0_dp
          knear(p) = 0
          cycle
       endif
       ip = ip + 1
       px = xyz(1,p); py = xyz(2,p); pz = xyz(3,p)
       k = 0
       if (local) k = nearest_walk(natm, R, ld, nnb, nb, px, py, pz, i, dist, stamp, ip)
       knear(p) = k
       if (weights == mode_delaunay) then
          if (k /= i) then
             if (.not. in_list(i, nnb(k), nb(1,k))) then
                w(p) = 0.0_dp
                cycle
             endif
          endif
          nl = nnb(k) + 1
       else
          nl = natm
       endif
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
          cl = cell_value(l, natm, R, a, dinv, step, mu, cells, ld, nnb, nb, px, py, pz, dist, stamp, ip)
          ssum = ssum + cl
          if (l == i) ci = cl
       enddo
       if (ssum > 0.0_dp) then
          w(p) = vol(p)*ci/ssum
       else
          w(p) = 0.0_dp
       endif
    enddo
    !$omp end do
    deallocate(dist, stamp)
    !$omp end parallel
    deallocate(a, dinv)

  end subroutine partition_weights

!######################################################################
! partition_sigma: sigma(r) = sum_all C_l - sum_{l in D_k^+} C_l
! (eq. 22 of the paper) with the chosen cell definition; O(natm) per
! point, diagnostics only
!######################################################################
  subroutine partition_sigma(natm, R, Z, adjust, step, mu, cells, npts, xyz, iatom, nnb, nb, sigma, ierr)

    integer(is), intent(in)  :: natm, Z(natm), adjust, step, cells, npts
    integer(is), intent(in)  :: iatom(npts), nnb(natm), nb(max(natm-1,1),natm)
    real(dp), intent(in)     :: R(3,natm), mu, xyz(3,npts)
    real(dp), intent(out)    :: sigma(npts)
    integer(is), intent(out) :: ierr
    real(dp), allocatable    :: a(:,:), dinv(:,:), dist(:)
    integer(is), allocatable :: stamp(:)
    integer(is)              :: ld, p, i, k, ip, l, ll
    real(dp)                 :: px, py, pz, sall, sloc, cl

    ierr = 0
    call step_check(step, mu, ierr)
    if (ierr /= 0) return
    ld = max(natm-1, 1)
    sigma = 0.0_dp
    if (natm == 1) return
    if (any(nnb < 0) .or. any(nnb > natm-1)) then
       ierr = 46
       return
    endif
    allocate(a(natm,natm), dinv(natm,natm))
    call partition_tables(natm, R, Z, adjust, a, dinv, ierr)
    if (ierr /= 0) return

    !$omp parallel private(dist, stamp, ip, p, i, k, px, py, pz, sall, sloc, cl, l, ll)
    allocate(dist(natm), stamp(natm))
    stamp = 0
    ip = 0
    !$omp do schedule(dynamic, 256)
    do p = 1, npts
       i = iatom(p)
       if (i < 1 .or. i > natm) cycle
       ip = ip + 1
       px = xyz(1,p); py = xyz(2,p); pz = xyz(3,p)
       k = nearest_walk(natm, R, ld, nnb, nb, px, py, pz, i, dist, stamp, ip)
       sall = 0.0_dp
       sloc = 0.0_dp
       do l = 1, natm
          cl = cell_value(l, natm, R, a, dinv, step, mu, cells, ld, nnb, nb, px, py, pz, dist, stamp, ip)
          sall = sall + cl
          if (l == k) then
             sloc = sloc + cl
          else if (in_list(l, nnb(k), nb(1,k))) then
             sloc = sloc + cl
          endif
       enddo
       sigma(p) = sall - sloc
       ll = 0
    enddo
    !$omp end do
    deallocate(dist, stamp)
    !$omp end parallel
    deallocate(a, dinv)

  end subroutine partition_sigma

end module grid_partition
