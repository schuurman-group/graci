!**********************************************************************
! grid_interface: bind(C) entry points of libgrid (stateless)
!   grid_size            number of points of the molecular grid
!   grid_atomic          atomic grids (unsorted, by atom), offsets added
!   grid_neighbours      Delaunay neighbour lists
!   grid_partition       partition weights (Becke / Delaunay cells and
!                        weights)
!   grid_sort            PySCF's box ordering (permutation)
!   grid_sigma           diagnostic sigma(r) of the Delaunay weights
!   grid_wderiv_contract weight derivatives contracted with integrands
!   grid_lebedev, grid_radial   (tests)
! Arrays are passed by address in Fortran (column-major) order: float64
! ('dptr') and int32 ('iptr') from Python with no copy; ierr = 0 on
! success. Scheme codes: see grid_constants.
!**********************************************************************
module grid_interface

  use iso_c_binding
  use grid_constants
  use grid_io
  use grid_lebedev
  use grid_radial
  use grid_atomic
  use grid_sort
  use grid_voronoi
  use grid_partition
  use grid_wderiv

  implicit none

contains

  subroutine grid_size_c(natm, Z, level, radial, prune, nrad_user, nang_user, npts, ierr) &
       bind(c, name='grid_size')

    integer(is), intent(in)  :: natm, Z(natm), level, radial, prune, nrad_user(natm), nang_user(natm)
    integer(is), intent(out) :: npts, ierr

    call molecular_size(natm, Z, level, radial, prune, nrad_user, nang_user, npts, ierr)

  end subroutine grid_size_c

  subroutine grid_atomic_c(natm, R, Z, level, radial, prune, nrad_user, nang_user, npts, xyz, iatom, vol, &
       ierr) bind(c, name='grid_atomic')

    integer(is), intent(in)  :: natm, Z(natm), level, radial, prune, nrad_user(natm), nang_user(natm), npts
    real(dp), intent(in)     :: R(3,natm)
    real(dp), intent(out)    :: xyz(3,npts), vol(npts)
    integer(is), intent(out) :: iatom(npts), ierr

    call molecular_grid(natm, R, Z, level, radial, prune, nrad_user, nang_user, npts, xyz, iatom, vol, ierr)

  end subroutine grid_atomic_c

  subroutine grid_neighbours_c(natm, R, tol, nnb, nb, ierr) bind(c, name='grid_neighbours')

    integer(is), intent(in)  :: natm
    real(dp), intent(in)     :: R(3,natm), tol
    integer(is), intent(out) :: nnb(natm), nb(max(natm-1,1),natm), ierr

    call delaunay_neighbours(natm, R, tol, nnb, nb, ierr)

  end subroutine grid_neighbours_c

  subroutine grid_partition_c(natm, R, Z, adjust, step, mu, cells, weights, npts, xyz, iatom, vol, &
       nnb, nb, w, knear, ierr) bind(c, name='grid_partition')

    integer(is), intent(in)  :: natm, Z(natm), adjust, step, cells, weights, npts
    integer(is), intent(in)  :: iatom(npts), nnb(natm), nb(max(natm-1,1),natm)
    real(dp), intent(in)     :: R(3,natm), mu, xyz(3,npts), vol(npts)
    real(dp), intent(out)    :: w(npts)
    integer(is), intent(out) :: knear(npts), ierr

    call partition_weights(natm, R, Z, adjust, step, mu, cells, weights, npts, xyz, iatom, vol, &
         nnb, nb, w, knear, ierr)

  end subroutine grid_partition_c

  subroutine grid_sort_c(natm, R, npts, xyz, perm, ierr) bind(c, name='grid_sort')

    integer(is), intent(in)  :: natm, npts
    real(dp), intent(in)     :: R(3,natm), xyz(3,npts)
    integer(is), intent(out) :: perm(npts), ierr

    call box_sort(natm, R, npts, xyz, perm, ierr)

  end subroutine grid_sort_c

  subroutine grid_sigma_c(natm, R, Z, adjust, step, mu, cells, npts, xyz, iatom, nnb, nb, sigma, ierr) &
       bind(c, name='grid_sigma')

    integer(is), intent(in)  :: natm, Z(natm), adjust, step, cells, npts
    integer(is), intent(in)  :: iatom(npts), nnb(natm), nb(max(natm-1,1),natm)
    real(dp), intent(in)     :: R(3,natm), mu, xyz(3,npts)
    real(dp), intent(out)    :: sigma(npts)
    integer(is), intent(out) :: ierr

    call partition_sigma(natm, R, Z, adjust, step, mu, cells, npts, xyz, iatom, nnb, nb, sigma, ierr)

  end subroutine grid_sigma_c

  subroutine grid_wderiv_contract_c(natm, R, Z, adjust, step, mu, cells, weights, npts, xyz, iatom, vol, &
       nnb, nb, nq, q, grad, ierr) bind(c, name='grid_wderiv_contract')

    integer(is), intent(in)  :: natm, Z(natm), adjust, step, cells, weights, npts, nq
    integer(is), intent(in)  :: iatom(npts), nnb(natm), nb(max(natm-1,1),natm)
    real(dp), intent(in)     :: R(3,natm), mu, xyz(3,npts), vol(npts), q(npts,nq)
    real(dp), intent(out)    :: grad(3,natm,nq)
    integer(is), intent(out) :: ierr

    call wderiv_contract(natm, R, Z, adjust, step, mu, cells, weights, npts, xyz, iatom, vol, &
         nnb, nb, nq, q, grad, ierr)

  end subroutine grid_wderiv_contract_c

  subroutine grid_lebedev_c(n, xyzw, ierr) bind(c, name='grid_lebedev')

    integer(is), intent(in)  :: n
    real(dp), intent(out)    :: xyzw(4,n)
    integer(is), intent(out) :: ierr

    call lebedev_grid(n, xyzw, ierr)

  end subroutine grid_lebedev_c

  subroutine grid_radial_c(scheme, n, z, r, dr, ierr) bind(c, name='grid_radial')

    integer(is), intent(in)  :: scheme, n, z
    real(dp), intent(out)    :: r(n), dr(n)
    integer(is), intent(out) :: ierr

    call radial_grid(scheme, n, z, r, dr, ierr)

  end subroutine grid_radial_c

end module grid_interface
