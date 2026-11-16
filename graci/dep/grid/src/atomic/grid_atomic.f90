!**********************************************************************
! grid_atomic: atomic grids (radial x Lebedev with pruning) in PySCF's
! conventions and point order (gen_atomic_grids): for each distinct
! angular size in increasing order, the radial shells with that size in
! groups of 12, angular point outer and radial shell inner within a
! group. Offsets are relative to the nucleus; vol = 4 pi r^2 dr w_ang.
!**********************************************************************
module grid_atomic

  use grid_constants
  use grid_tables
  use grid_lebedev
  use grid_radial

  implicit none

contains

!######################################################################
! period_of: PySCF's period index (1..7) of nuclear charge z
!######################################################################
  function period_of(z) result(p)

    integer(is), intent(in) :: z
    integer(is)             :: p, i

    p = 1
    do i = 1, 7
       if (z > period_zmax(i)) p = i + 1
    enddo
    p = min(p, 7_is)

  end function period_of

!######################################################################
! default_sizes: radial and angular (points) sizes of the PySCF level
!######################################################################
  subroutine default_sizes(z, level, nrad, nang, ierr)

    integer(is), intent(in)  :: z, level
    integer(is), intent(out) :: nrad, nang, ierr
    integer(is)              :: p

    ierr = 0
    if (level < 0 .or. level > 9) then
       ierr = 11
       nrad = 0
       nang = 0
       return
    endif
    p = period_of(z)
    nrad = rad_grids(level, p)
    nang = lebedev_size_of_order(ang_order(level, p))

  end subroutine default_sizes

!######################################################################
! resolve_sizes: sizes of one atom from the level and the optional
! user overrides (0 = default; an angular override may be a size or an
! algebraic order); ierr 12: unsupported angular size
!######################################################################
  subroutine resolve_sizes(z, level, nrad_user, nang_user, nrad, nang, ierr)

    integer(is), intent(in)  :: z, level, nrad_user, nang_user
    integer(is), intent(out) :: nrad, nang, ierr
    integer(is)              :: n

    call default_sizes(z, level, nrad, nang, ierr)
    if (ierr /= 0) return
    if (nrad_user > 0) nrad = nrad_user
    if (nang_user > 0) then
       if (lebedev_index(nang_user) > 0) then
          nang = nang_user
       else
          n = lebedev_size_of_order(nang_user)
          if (n == 0) then
             ierr = 12
             return
          endif
          nang = n
       endif
    endif

  end subroutine resolve_sizes

!######################################################################
! prune_angs: angular size of every radial shell
!######################################################################
  subroutine prune_angs(scheme, z, nrad, rad, nang, angs, ierr)

    integer(is), intent(in)  :: scheme, z, nrad, nang
    real(dp), intent(in)     :: rad(nrad)
    integer(is), intent(out) :: angs(nrad), ierr
    real(dp), parameter      :: alphas(4,3) = reshape([0.25_dp, 0.5_dp, 1.0_dp, 4.5_dp, &
                                                       0.1667_dp, 0.5_dp, 0.9_dp, 3.5_dp, &
                                                       0.1_dp, 0.4_dp, 0.8_dp, 2.5_dp], [4, 3])
    integer(is), parameter   :: sg1_sizes(5) = [6, 38, 86, 194, 86]
    integer(is)              :: row, i, k, place, idx, leb_l(5), m
    real(dp)                 :: r_atom

    ierr = 0
    select case (scheme)

    case (prune_none)
       angs = nang

    case (prune_nwchem)
       if (z < 0 .or. z > zmax_bragg) then
          ierr = 13
          return
       endif
       if (nang < 50) then
          angs = nang
          return
       else if (nang == 50) then
          leb_l = [1, 2, 2, 2, 1]
       else
          ! position of nang in LEBEDEV_NGRID(4:) (0-based, PySCF)
          idx = lebedev_index(nang) - 5
          if (idx < 0) then
             ierr = 12
             return
          endif
          leb_l = [1, 3, idx-1, idx, idx-1]
       endif
       row = 3
       if (z <= 2) then
          row = 1
       else if (z <= 10) then
          row = 2
       endif
       r_atom = bragg_radii(z+1) + 1.0e-200_dp
       do i = 1, nrad
          place = 0
          do k = 1, 4
             if (rad(i)/r_atom > alphas(k,row)) place = place + 1
          enddo
          m = leb_l(place+1)
          ! leb_ngrid(m) = LEBEDEV_NGRID[4+m] = leb_sizes(5+m)
          angs(i) = leb_sizes(5+m)
       enddo

    case (prune_sg1)
       if (z < 0 .or. z > zmax_sg1) then
          ierr = 14
          return
       endif
       row = 3
       if (z <= 2) then
          row = 1
       else if (z <= 10) then
          row = 2
       endif
       r_atom = sg1_radii(z+1) + 1.0e-200_dp
       do i = 1, nrad
          place = 0
          do k = 1, 4
             if (rad(i)/r_atom > alphas(k,row)) place = place + 1
          enddo
          angs(i) = sg1_sizes(place+1)
       enddo

    case (prune_treutler)
       do i = 1, nrad
          if (i <= nrad/3) then
             angs(i) = 14
          else if (i <= nrad/2) then
             angs(i) = 50
          else
             angs(i) = nang
          endif
       enddo

    case default
       ierr = 15
    end select

  end subroutine prune_angs

!######################################################################
! atomic_shells: radial grid and the pruned angular sizes of one atom
!######################################################################
  subroutine atomic_shells(z, level, radial, prune, nrad_user, nang_user, nrad, rad, rw, angs, ierr)

    integer(is), intent(in)               :: z, level, radial, prune, nrad_user, nang_user
    integer(is), intent(out)              :: nrad, ierr
    real(dp), allocatable, intent(out)    :: rad(:), rw(:)
    integer(is), allocatable, intent(out) :: angs(:)
    integer(is)                           :: nang, i
    real(dp), allocatable                 :: dr(:)

    call resolve_sizes(z, level, nrad_user, nang_user, nrad, nang, ierr)
    if (ierr /= 0) return
    allocate(rad(nrad), rw(nrad), dr(nrad), angs(nrad))
    call radial_grid(radial, nrad, z, rad, dr, ierr)
    if (ierr /= 0) return
    do i = 1, nrad
       rw(i) = 4.0_dp*pi*rad(i)**2*dr(i)
    enddo
    call prune_angs(prune, z, nrad, rad, nang, angs, ierr)
    deallocate(dr)

  end subroutine atomic_shells

!######################################################################
! atomic_size: number of points of one atom
!######################################################################
  subroutine atomic_size(z, level, radial, prune, nrad_user, nang_user, npts, ierr)

    integer(is), intent(in)  :: z, level, radial, prune, nrad_user, nang_user
    integer(is), intent(out) :: npts, ierr
    integer(is)              :: nrad
    real(dp), allocatable    :: rad(:), rw(:)
    integer(is), allocatable :: angs(:)

    npts = 0
    call atomic_shells(z, level, radial, prune, nrad_user, nang_user, nrad, rad, rw, angs, ierr)
    if (ierr /= 0) return
    npts = sum(angs)

  end subroutine atomic_size

!######################################################################
! atomic_grid: offsets xyz(3,npts) and volumes vol(npts) of one atom
! (npts from atomic_size)
!######################################################################
  subroutine atomic_grid(z, level, radial, prune, nrad_user, nang_user, npts, xyz, vol, ierr)

    integer(is), intent(in)  :: z, level, radial, prune, nrad_user, nang_user, npts
    real(dp), intent(out)    :: xyz(3,npts), vol(npts)
    integer(is), intent(out) :: ierr
    integer(is)              :: nrad, n, i0, i1, i, j, p, ng, nsz, k
    real(dp), allocatable    :: rad(:), rw(:), ang(:,:)
    integer(is), allocatable :: angs(:), idx(:), sizes(:)

    call atomic_shells(z, level, radial, prune, nrad_user, nang_user, nrad, rad, rw, angs, ierr)
    if (ierr /= 0) return
    if (sum(angs) /= npts) then
       ierr = 16
       return
    endif
    ! distinct angular sizes in increasing order
    allocate(sizes(nrad), idx(nrad))
    nsz = 0
    do i = 1, nrad
       if (nsz == 0) then
          nsz = 1
          sizes(1) = angs(i)
       else if (all(sizes(1:nsz) /= angs(i))) then
          nsz = nsz + 1
          sizes(nsz) = angs(i)
       endif
    enddo
    call isort(nsz, sizes)
    p = 0
    do k = 1, nsz
       n = sizes(k)
       allocate(ang(4,n))
       call lebedev_grid(n, ang, ierr)
       if (ierr /= 0) return
       ! radial shells with this angular size
       ng = 0
       do i = 1, nrad
          if (angs(i) == n) then
             ng = ng + 1
             idx(ng) = i
          endif
       enddo
       do i0 = 1, ng, 12
          i1 = min(i0 + 11, ng)
          do j = 1, n
             do i = i0, i1
                p = p + 1
                xyz(1,p) = rad(idx(i))*ang(1,j)
                xyz(2,p) = rad(idx(i))*ang(2,j)
                xyz(3,p) = rad(idx(i))*ang(3,j)
                vol(p) = rw(idx(i))*ang(4,j)
             enddo
          enddo
       enddo
       deallocate(ang)
    enddo
    if (p /= npts) ierr = 16

  end subroutine atomic_grid

!######################################################################
! molecular_size / molecular_grid: all atoms, xyz = R_i + offset,
! iatom(p) = parent atom (1-based)
!######################################################################
  subroutine molecular_size(natm, Z, level, radial, prune, nrad_user, nang_user, npts, ierr)

    integer(is), intent(in)  :: natm, Z(natm), level, radial, prune, nrad_user(natm), nang_user(natm)
    integer(is), intent(out) :: npts, ierr
    integer(is)              :: i, n

    npts = 0
    ierr = 0
    do i = 1, natm
       call atomic_size(Z(i), level, radial, prune, nrad_user(i), nang_user(i), n, ierr)
       if (ierr /= 0) return
       npts = npts + n
    enddo

  end subroutine molecular_size

  subroutine molecular_grid(natm, R, Z, level, radial, prune, nrad_user, nang_user, npts, xyz, iatom, vol, ierr)

    integer(is), intent(in)  :: natm, Z(natm), level, radial, prune, nrad_user(natm), nang_user(natm), npts
    real(dp), intent(in)     :: R(3,natm)
    real(dp), intent(out)    :: xyz(3,npts), vol(npts)
    integer(is), intent(out) :: iatom(npts), ierr
    integer(is)              :: i, n, p0, p

    ierr = 0
    p0 = 0
    do i = 1, natm
       call atomic_size(Z(i), level, radial, prune, nrad_user(i), nang_user(i), n, ierr)
       if (ierr /= 0) return
       if (p0 + n > npts) then
          ierr = 16
          return
       endif
       call atomic_grid(Z(i), level, radial, prune, nrad_user(i), nang_user(i), n, xyz(1,p0+1), vol(p0+1), ierr)
       if (ierr /= 0) return
       do p = p0 + 1, p0 + n
          xyz(1,p) = xyz(1,p) + R(1,i)
          xyz(2,p) = xyz(2,p) + R(2,i)
          xyz(3,p) = xyz(3,p) + R(3,i)
          iatom(p) = i
       enddo
       p0 = p0 + n
    enddo
    if (p0 /= npts) ierr = 16

  end subroutine molecular_grid

!######################################################################
! isort: insertion sort of a short integer list
!######################################################################
  subroutine isort(n, a)

    integer(is), intent(in)    :: n
    integer(is), intent(inout) :: a(n)
    integer(is)                :: i, j, t

    do i = 2, n
       t = a(i)
       j = i - 1
       do while (j >= 1)
          if (a(j) <= t) exit
          a(j+1) = a(j)
          j = j - 1
       enddo
       a(j+1) = t
    enddo

  end subroutine isort

end module grid_atomic
