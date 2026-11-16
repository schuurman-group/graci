!**********************************************************************
! grid_lebedev: Lebedev angular grids from the orbit tables of
! grid_lebedev_tab (the generator gen_oh expands one octahedral orbit;
! point order identical to PySCF's MakeAngularGrid)
!**********************************************************************
module grid_lebedev

  use grid_constants
  use grid_tables
  use grid_lebedev_tab

  implicit none

contains

!######################################################################
! lebedev_index: position of the size n in leb_sizes (0 if unsupported)
!######################################################################
  function lebedev_index(n) result(k)

    integer(is), intent(in) :: n
    integer(is)             :: k, i

    k = 0
    do i = 1, nleb
       if (leb_sizes(i) == n) then
          k = i
          return
       endif
    enddo

  end function lebedev_index

!######################################################################
! lebedev_size_of_order: number of points of the grid of a given
! algebraic order (0 if the order is not tabulated)
!######################################################################
  function lebedev_size_of_order(order) result(n)

    integer(is), intent(in) :: order
    integer(is)             :: n, i

    n = 0
    do i = 1, nleb
       if (leb_orders(i) == order) then
          n = leb_sizes(i)
          return
       endif
    enddo

  end function lebedev_size_of_order

!######################################################################
! lebedev_grid: the n-point grid, xyzw(4,n) = (x, y, z, weight), the
! weights summing to one; ierr 1: unsupported size, 2: table error
!######################################################################
  subroutine lebedev_grid(n, xyzw, ierr)

    integer(is), intent(in)  :: n
    real(dp), intent(out)    :: xyzw(4,n)
    integer(is), intent(out) :: ierr
    integer(is)              :: k, io, m, num

    ierr = 0
    if (n == 1) then
       xyzw(1:3,1) = 0.0_dp
       xyzw(4,1) = 1.0_dp
       return
    endif
    k = lebedev_index(n)
    if (k == 0) then
       ierr = 1
       return
    endif
    m = 0
    do io = leb_orb0(k), leb_orb0(k+1) - 1
       call gen_oh(leb_code(io), leb_a(io), leb_b(io), leb_v(io), xyzw(1,m+1), n - m, num, ierr)
       if (ierr /= 0) return
       m = m + num
    enddo
    if (m /= n) ierr = 2

  end subroutine lebedev_grid

!######################################################################
! gen_oh: the points of one orbit of the octahedral group with
! inversion (Lebedev-Laikov), code 0: (1,0,0) 6 points; 1: (0,a,a)
! 12; 2: (a,a,a) 8; 3: (a,a,b) 24; 4: (a,b,0) 24; 5: (a,b,c) 48;
! g(4,*) receives at most room points; ierr 3: no room
!######################################################################
  subroutine gen_oh(code, a0, b0, v, g, room, num, ierr)

    integer(is), intent(in)  :: code, room
    real(dp), intent(in)     :: a0, b0, v
    real(dp), intent(inout)  :: g(4,*)
    integer(is), intent(out) :: num, ierr
    real(dp)                 :: a, b, c
    integer(is)              :: p

    ierr = 0
    a = a0
    b = b0
    select case (code)
    case (0)
       num = 6
    case (1)
       num = 12
    case (2)
       num = 8
    case (3, 4)
       num = 24
    case (5)
       num = 48
    case default
       ierr = 2
       num = 0
       return
    end select
    if (num > room) then
       ierr = 3
       return
    endif
    do p = 1, num
       g(4,p) = v
    enddo

    select case (code)
    case (0)
       a = 1.0_dp
       call put(1,  a, 0.0_dp, 0.0_dp)
       call put(2, -a, 0.0_dp, 0.0_dp)
       call put(3, 0.0_dp,  a, 0.0_dp)
       call put(4, 0.0_dp, -a, 0.0_dp)
       call put(5, 0.0_dp, 0.0_dp,  a)
       call put(6, 0.0_dp, 0.0_dp, -a)
    case (1)
       a = sqrt(0.5_dp)
       call put(1, 0.0_dp,  a,  a)
       call put(2, 0.0_dp, -a,  a)
       call put(3, 0.0_dp,  a, -a)
       call put(4, 0.0_dp, -a, -a)
       call put(5,  a, 0.0_dp,  a)
       call put(6, -a, 0.0_dp,  a)
       call put(7,  a, 0.0_dp, -a)
       call put(8, -a, 0.0_dp, -a)
       call put(9,   a,  a, 0.0_dp)
       call put(10, -a,  a, 0.0_dp)
       call put(11,  a, -a, 0.0_dp)
       call put(12, -a, -a, 0.0_dp)
    case (2)
       a = sqrt(1.0_dp/3.0_dp)
       call put(1,  a,  a,  a)
       call put(2, -a,  a,  a)
       call put(3,  a, -a,  a)
       call put(4, -a, -a,  a)
       call put(5,  a,  a, -a)
       call put(6, -a,  a, -a)
       call put(7,  a, -a, -a)
       call put(8, -a, -a, -a)
    case (3)
       b = sqrt(1.0_dp - 2.0_dp*a*a)
       call put(1,  a,  a,  b)
       call put(2, -a,  a,  b)
       call put(3,  a, -a,  b)
       call put(4, -a, -a,  b)
       call put(5,  a,  a, -b)
       call put(6, -a,  a, -b)
       call put(7,  a, -a, -b)
       call put(8, -a, -a, -b)
       call put(9,   a,  b,  a)
       call put(10, -a,  b,  a)
       call put(11,  a, -b,  a)
       call put(12, -a, -b,  a)
       call put(13,  a,  b, -a)
       call put(14, -a,  b, -a)
       call put(15,  a, -b, -a)
       call put(16, -a, -b, -a)
       call put(17,  b,  a,  a)
       call put(18, -b,  a,  a)
       call put(19,  b, -a,  a)
       call put(20, -b, -a,  a)
       call put(21,  b,  a, -a)
       call put(22, -b,  a, -a)
       call put(23,  b, -a, -a)
       call put(24, -b, -a, -a)
    case (4)
       b = sqrt(1.0_dp - a*a)
       call put(1,  a,  b, 0.0_dp)
       call put(2, -a,  b, 0.0_dp)
       call put(3,  a, -b, 0.0_dp)
       call put(4, -a, -b, 0.0_dp)
       call put(5,  b,  a, 0.0_dp)
       call put(6, -b,  a, 0.0_dp)
       call put(7,  b, -a, 0.0_dp)
       call put(8, -b, -a, 0.0_dp)
       call put(9,   a, 0.0_dp,  b)
       call put(10, -a, 0.0_dp,  b)
       call put(11,  a, 0.0_dp, -b)
       call put(12, -a, 0.0_dp, -b)
       call put(13,  b, 0.0_dp,  a)
       call put(14, -b, 0.0_dp,  a)
       call put(15,  b, 0.0_dp, -a)
       call put(16, -b, 0.0_dp, -a)
       call put(17, 0.0_dp,  a,  b)
       call put(18, 0.0_dp, -a,  b)
       call put(19, 0.0_dp,  a, -b)
       call put(20, 0.0_dp, -a, -b)
       call put(21, 0.0_dp,  b,  a)
       call put(22, 0.0_dp, -b,  a)
       call put(23, 0.0_dp,  b, -a)
       call put(24, 0.0_dp, -b, -a)
    case (5)
       c = sqrt(1.0_dp - a*a - b*b)
       call put(1,  a,  b,  c)
       call put(2, -a,  b,  c)
       call put(3,  a, -b,  c)
       call put(4, -a, -b,  c)
       call put(5,  a,  b, -c)
       call put(6, -a,  b, -c)
       call put(7,  a, -b, -c)
       call put(8, -a, -b, -c)
       call put(9,   a,  c,  b)
       call put(10, -a,  c,  b)
       call put(11,  a, -c,  b)
       call put(12, -a, -c,  b)
       call put(13,  a,  c, -b)
       call put(14, -a,  c, -b)
       call put(15,  a, -c, -b)
       call put(16, -a, -c, -b)
       call put(17,  b,  a,  c)
       call put(18, -b,  a,  c)
       call put(19,  b, -a,  c)
       call put(20, -b, -a,  c)
       call put(21,  b,  a, -c)
       call put(22, -b,  a, -c)
       call put(23,  b, -a, -c)
       call put(24, -b, -a, -c)
       call put(25,  b,  c,  a)
       call put(26, -b,  c,  a)
       call put(27,  b, -c,  a)
       call put(28, -b, -c,  a)
       call put(29,  b,  c, -a)
       call put(30, -b,  c, -a)
       call put(31,  b, -c, -a)
       call put(32, -b, -c, -a)
       call put(33,  c,  a,  b)
       call put(34, -c,  a,  b)
       call put(35,  c, -a,  b)
       call put(36, -c, -a,  b)
       call put(37,  c,  a, -b)
       call put(38, -c,  a, -b)
       call put(39,  c, -a, -b)
       call put(40, -c, -a, -b)
       call put(41,  c,  b,  a)
       call put(42, -c,  b,  a)
       call put(43,  c, -b,  a)
       call put(44, -c, -b,  a)
       call put(45,  c,  b, -a)
       call put(46, -c,  b, -a)
       call put(47,  c, -b, -a)
       call put(48, -c, -b, -a)
    end select

  contains

    subroutine put(p, x, y, z)
      integer(is), intent(in) :: p
      real(dp), intent(in)    :: x, y, z
      g(1,p) = x
      g(2,p) = y
      g(3,p) = z
    end subroutine put

  end subroutine gen_oh

end module grid_lebedev
