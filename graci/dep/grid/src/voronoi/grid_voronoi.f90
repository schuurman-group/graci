!**********************************************************************
! grid_voronoi: Delaunay neighbour lists in the broad sense of Gill &
! Tsoukalas (JCP 165, 174104 (2026), Appendix A): j is a neighbour of
! i when the Voronoi cells of i and j share at least one point (facet,
! edge or vertex). The cell of i is built by clipping a large cube by
! the bisector half-spaces of the other nuclei in order of distance,
! keeping only the vertices (each with the set of planes it lies on; an
! edge is a vertex pair sharing two planes). Clipping stops once the
! next nucleus is farther than twice the cell radius. j touches the
! cell if the minimum over the vertices of the signed distance to its
! bisector (a linear function, so the minimum is at a vertex) is at
! most tol. Contacts beyond the cube (half-width = molecular extent +
! 100 bohr) are lost; no grid point lies there.
!**********************************************************************
module grid_voronoi

  use grid_constants

  implicit none

  integer(is), parameter :: maxp = 32      ! planes per vertex
  integer(is), parameter :: maxv0 = 1024   ! initial vertex capacity
  ! half-width of the bounding cube beyond the molecular extent: contacts
  ! between the unbounded cells of hull atoms can lie very far out (two
  ! nearly coplanar hull atoms meet at a distance ~ extent^2 / height
  ! difference); 1e4 bohr keeps double precision comfortable (tol 1e-6)
  real(dp), parameter    :: box_margin = 1.0e4_dp

contains

!######################################################################
! delaunay_neighbours: nnb(i) neighbours nb(1:nnb(i), i) (ascending)
! ierr 31: too many planes on one vertex; 32: too many vertices
!######################################################################
  subroutine delaunay_neighbours(natm, R, tol, nnb, nb, ierr)

    integer(is), intent(in)  :: natm
    real(dp), intent(in)     :: R(3,natm), tol
    integer(is), intent(out) :: nnb(natm), nb(max(natm-1,1),natm), ierr
    real(dp)                 :: L, lo(3), hi(3)
    integer(is)              :: i, ie, c

    ierr = 0
    nnb = 0
    if (natm <= 1) return
    if (natm == 2) then
       nnb = 1
       nb(1,1) = 2
       nb(1,2) = 1
       return
    endif
    do c = 1, 3
       lo(c) = minval(R(c,:))
       hi(c) = maxval(R(c,:))
    enddo
    L = maxval(hi - lo) + box_margin

    !$omp parallel do private(i, ie) schedule(dynamic)
    do i = 1, natm
       call cell_neighbours(natm, R, i, L, tol, nnb(i), nb(1,i), ie)
       if (ie /= 0) then
          !$omp atomic write
          ierr = ie
       endif
    enddo
    !$omp end parallel do

  end subroutine delaunay_neighbours

!######################################################################
! cell_neighbours: the Voronoi cell of nucleus i and its neighbours
!######################################################################
  subroutine cell_neighbours(natm, R, i, L, tol, nn, list, ierr)

    integer(is), intent(in)  :: natm, i
    real(dp), intent(in)     :: R(3,natm), L, tol
    integer(is), intent(out) :: nn, list(natm-1), ierr
    real(dp), allocatable    :: d(:), u(:,:), g(:), un(:,:)
    integer(is), allocatable :: ord(:), np(:), pl(:,:), cls(:), npn(:), pln(:,:), keep(:)
    real(dp)                 :: Dv(3), dn, Dh(3), rho2, t, w(3), gmin, dmax
    integer(is)              :: maxv, nv, j, jj, k, a, b, nin, non, nout, nnew, m, nsh, kk, nkeep, q
    logical                  :: dup

    ierr = 0
    nn = 0
    allocate(d(natm), ord(natm))
    do j = 1, natm
       d(j) = sqrt(sum((R(:,j) - R(:,i))**2))
    enddo
    ! other nuclei in order of distance
    m = 0
    do j = 1, natm
       if (j == i) cycle
       m = m + 1
       ord(m) = j
    enddo
    call sort_index(natm-1, d, ord)

    ! the cube
    maxv = maxv0
    allocate(u(3,maxv), g(maxv), np(maxv), pl(maxp,maxv), cls(maxv), keep(maxv))
    allocate(un(3,maxv), npn(maxv), pln(maxp,maxv))
    nv = 0
    do a = -1, 1, 2
       do b = -1, 1, 2
          do k = -1, 1, 2
             nv = nv + 1
             u(1,nv) = real(a, dp)*L
             u(2,nv) = real(b, dp)*L
             u(3,nv) = real(k, dp)*L
             np(nv) = 3
             pl(1,nv) = merge(-1_is, -2_is, a > 0)
             pl(2,nv) = merge(-3_is, -4_is, b > 0)
             pl(3,nv) = merge(-5_is, -6_is, k > 0)
          enddo
       enddo
    enddo
    rho2 = 3.0_dp*L*L

    ! successive clipping
    do jj = 1, natm - 1
       j = ord(jj)
       if (d(j)*d(j) > 4.0_dp*rho2) exit
       Dv = R(:,j) - R(:,i)
       dn = d(j)
       Dh = Dv/dn
       nin = 0; non = 0; nout = 0
       do k = 1, nv
          g(k) = 0.5_dp*dn - (u(1,k)*Dh(1) + u(2,k)*Dh(2) + u(3,k)*Dh(3))
          if (g(k) > tol) then
             cls(k) = 1
             nin = nin + 1
          else if (g(k) < -tol) then
             cls(k) = -1
             nout = nout + 1
          else
             cls(k) = 0
             non = non + 1
          endif
       enddo
       if (nout == 0) then
          ! the plane at most touches: record it on the touching vertices
          do k = 1, nv
             if (cls(k) == 0) then
                if (np(k) >= maxp) then
                   ierr = 31
                   return
                endif
                np(k) = np(k) + 1
                pl(np(k),k) = j
             endif
          enddo
          cycle
       endif
       if (nin == 0 .and. non == 0) then
          ! cannot happen (the nucleus itself is inside its cell)
          ierr = 33
          return
       endif
       ! new vertices on the inside-outside edges
       nnew = 0
       do a = 1, nv
          if (cls(a) /= 1) cycle
          do b = 1, nv
             if (cls(b) /= -1) cycle
             call shared_planes(np(a), pl(1,a), np(b), pl(1,b), nsh, pln(1,maxv))
             if (nsh < 2) cycle
             t = g(a)/(g(a) - g(b))
             w = u(:,a) + t*(u(:,b) - u(:,a))
             ! merge with an existing touching vertex or a previous new one
             dup = .false.
             do k = 1, nv
                if (cls(k) == 0 .and. sum((u(:,k) - w)**2) <= (2.0_dp*tol)**2) then
                   call merge_planes(np(k), pl(1,k), nsh, pln(1,maxv), j, ierr)
                   if (ierr /= 0) return
                   dup = .true.
                   exit
                endif
             enddo
             if (.not. dup) then
                do k = 1, nnew
                   if (sum((un(:,k) - w)**2) <= (2.0_dp*tol)**2) then
                      call merge_planes(npn(k), pln(1,k), nsh, pln(1,maxv), j, ierr)
                      if (ierr /= 0) return
                      dup = .true.
                      exit
                   endif
                enddo
             endif
             if (dup) cycle
             if (nnew + 1 >= maxv) then
                ierr = 32
                return
             endif
             nnew = nnew + 1
             un(:,nnew) = w
             npn(nnew) = nsh + 1
             pln(1:nsh,nnew) = pln(1:nsh,maxv)
             pln(nsh+1,nnew) = j
          enddo
       enddo
       ! compact: inside and touching vertices, then the new ones
       nkeep = 0
       do k = 1, nv
          if (cls(k) == -1) cycle
          nkeep = nkeep + 1
          keep(nkeep) = k
       enddo
       if (nkeep + nnew > maxv) then
          ierr = 32
          return
       endif
       do q = 1, nkeep
          k = keep(q)
          if (q /= k) then
             u(:,q) = u(:,k)
             np(q) = np(k)
             pl(1:np(q),q) = pl(1:np(k),k)
             cls(q) = cls(k)
          endif
          if (cls(q) == 0) then
             if (np(q) >= maxp) then
                ierr = 31
                return
             endif
             np(q) = np(q) + 1
             pl(np(q),q) = j
          endif
       enddo
       nv = nkeep
       do k = 1, nnew
          nv = nv + 1
          u(:,nv) = un(:,k)
          np(nv) = npn(k)
          pl(1:npn(k),nv) = pln(1:npn(k),k)
       enddo
       rho2 = 0.0_dp
       do k = 1, nv
          rho2 = max(rho2, u(1,k)**2 + u(2,k)**2 + u(3,k)**2)
       enddo
    enddo

    ! neighbours: bisectors touching the final cell
    dmax = 2.0_dp*sqrt(rho2) + 2.0_dp*tol
    do jj = 1, natm - 1
       j = ord(jj)
       if (d(j) > dmax) exit
       Dh = (R(:,j) - R(:,i))/d(j)
       gmin = huge(1.0_dp)
       do k = 1, nv
          gmin = min(gmin, 0.5_dp*d(j) - (u(1,k)*Dh(1) + u(2,k)*Dh(2) + u(3,k)*Dh(3)))
       enddo
       if (gmin <= tol) then
          nn = nn + 1
          list(nn) = j
       endif
    enddo
    ! ascending atom index
    do a = 2, nn
       k = list(a)
       b = a - 1
       do while (b >= 1)
          if (list(b) <= k) exit
          list(b+1) = list(b)
          b = b - 1
       enddo
       list(b+1) = k
    enddo
    deallocate(d, ord, u, g, np, pl, cls, keep, un, npn, pln)

  end subroutine cell_neighbours

!######################################################################
! shared_planes: the plane ids common to two vertices
!######################################################################
  subroutine shared_planes(na, pa, nbp, pb, nsh, sh)

    integer(is), intent(in)  :: na, nbp, pa(na), pb(nbp)
    integer(is), intent(out) :: nsh, sh(maxp)
    integer(is)              :: a, b

    nsh = 0
    do a = 1, na
       do b = 1, nbp
          if (pa(a) == pb(b)) then
             nsh = nsh + 1
             sh(nsh) = pa(a)
             exit
          endif
       enddo
    enddo

  end subroutine shared_planes

!######################################################################
! merge_planes: add the planes sh(1:nsh) and j to a vertex's set
!######################################################################
  subroutine merge_planes(np, pl, nsh, sh, j, ierr)

    integer(is), intent(inout) :: np, pl(maxp)
    integer(is), intent(in)    :: nsh, sh(maxp), j
    integer(is), intent(out)   :: ierr
    integer(is)                :: a

    ierr = 0
    do a = 1, nsh + 1
       if (a <= nsh) then
          if (any(pl(1:np) == sh(a))) cycle
          if (np >= maxp) then
             ierr = 31
             return
          endif
          np = np + 1
          pl(np) = sh(a)
       else
          if (any(pl(1:np) == j)) cycle
          if (np >= maxp) then
             ierr = 31
             return
          endif
          np = np + 1
          pl(np) = j
       endif
    enddo

  end subroutine merge_planes

!######################################################################
! sort_index: heapsort of the index list ord(1:n) by the key d
!######################################################################
  subroutine sort_index(n, d, ord)

    integer(is), intent(in)    :: n
    real(dp), intent(in)       :: d(*)
    integer(is), intent(inout) :: ord(n)
    integer(is)                :: i, j, k, t

    if (n < 2) return
    do i = n/2, 1, -1
       call sift(i, n)
    enddo
    do i = n, 2, -1
       t = ord(1); ord(1) = ord(i); ord(i) = t
       call sift(1, i-1)
    enddo

  contains

    subroutine sift(start, last)
      integer(is), intent(in) :: start, last
      integer(is) :: root, child, tmp
      root = start
      do
         child = 2*root
         if (child > last) exit
         if (child < last) then
            if (d(ord(child)) < d(ord(child+1))) child = child + 1
         endif
         if (d(ord(root)) < d(ord(child))) then
            tmp = ord(root); ord(root) = ord(child); ord(child) = tmp
            root = child
         else
            exit
         endif
      enddo
    end subroutine sift

  end subroutine sort_index

end module grid_voronoi
