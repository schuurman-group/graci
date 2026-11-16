!**********************************************************************
! grid_sort: PySCF's spatial grouping of the grid points
! (gen_grid.arg_group_grids): boxes of 1.2 bohr over the nuclear
! bounding box widened by 4.2 bohr, points ordered by the lexicographic
! box index (x, y, z), stable within a box. perm(k) = original index of
! the k-th sorted point (1-based).
!**********************************************************************
module grid_sort

  use grid_constants

  implicit none

  real(dp), parameter :: group_box_size = 1.2_dp
  real(dp), parameter :: group_boundary_penalty = 4.2_dp

contains

!######################################################################
! round_half_even: NumPy's round() of a real to an integer
!######################################################################
  function round_half_even(x) result(n)

    real(dp), intent(in) :: x
    integer(is)          :: n
    real(dp)             :: f, d

    f = floor(x)
    d = x - f
    if (d > 0.5_dp) then
       n = int(f, is) + 1
    else if (d < 0.5_dp) then
       n = int(f, is)
    else
       n = int(f, is)
       if (mod(n, 2) /= 0) n = n + 1
    endif

  end function round_half_even

!######################################################################
! box_sort
!######################################################################
  subroutine box_sort(natm, R, npts, xyz, perm, ierr)

    integer(is), intent(in)  :: natm, npts
    real(dp), intent(in)     :: R(3,natm), xyz(3,npts)
    integer(is), intent(out) :: perm(npts), ierr
    real(dp)                 :: lo(3), hi(3), bsize(3), binv(3), frac
    integer(is)              :: nbox(3), c, p, id(3), key, nkey
    integer(is), allocatable :: keys(:), count(:)

    ierr = 0
    if (natm < 1) then
       ierr = 21
       return
    endif
    do c = 1, 3
       lo(c) = minval(R(c,:)) - group_boundary_penalty
       hi(c) = maxval(R(c,:)) + group_boundary_penalty
       nbox(c) = round_half_even((hi(c) - lo(c))*(1.0_dp/group_box_size))
       if (nbox(c) < 1) nbox(c) = 1
       bsize(c) = (hi(c) - lo(c))/real(nbox(c), dp)
       binv(c) = 1.0_dp/bsize(c)
    enddo
    nkey = (nbox(1) + 2)*(nbox(2) + 2)*(nbox(3) + 2)
    allocate(keys(npts), count(0:nkey))
    count = 0
    do p = 1, npts
       do c = 1, 3
          frac = (xyz(c,p) - lo(c))*binv(c)
          id(c) = int(floor(frac), is)
          if (id(c) < -1) id(c) = -1
          if (id(c) > nbox(c)) id(c) = nbox(c)
       enddo
       key = ((id(1) + 1)*(nbox(2) + 2) + (id(2) + 1))*(nbox(3) + 2) + (id(3) + 1)
       keys(p) = key
       count(key+1) = count(key+1) + 1
    enddo
    ! prefix sums: count(k) = first position of key k (0-based)
    do key = 1, nkey
       count(key) = count(key) + count(key-1)
    enddo
    do p = 1, npts
       key = keys(p)
       count(key) = count(key) + 1
       perm(count(key)) = p
    enddo
    deallocate(keys, count)

  end subroutine box_sort

end module grid_sort
