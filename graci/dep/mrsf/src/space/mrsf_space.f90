!**********************************************************************
! mrsf_space: orbital classification, amplitude bookkeeping, the
! dimensional (spin-adaptation) transformation and masks
!**********************************************************************
module mrsf_space

  use mrsf_constants
  use mrsf_global
  use mrsf_io

  implicit none

contains

!######################################################################
! classify_orbitals: C/O/V classes from an occupation vector.
! hmap: hole set (C ascending, then O1, O2); pmap: particle set
! (O1, O2, then V ascending). Does not touch module state.
!######################################################################
  subroutine classify_orbitals(n, occv, nc1, nv1, hmap1, pmap1, io1, io2)

    integer(is), intent(in)               :: n
    real(dp), intent(in)                  :: occv(n)
    integer(is), intent(out)              :: nc1, nv1, io1, io2
    integer(is), allocatable, intent(out) :: hmap1(:), pmap1(:)
    integer(is)                           :: p, no, ic, iv

    nc1 = 0; nv1 = 0; no = 0; io1 = 0; io2 = 0
    do p = 1, n
       if (abs(occv(p)-2.0_dp) < 1.0e-6_dp) then
          nc1 = nc1 + 1
       else if (abs(occv(p)-1.0_dp) < 1.0e-6_dp) then
          no = no + 1
          if (no == 1) io1 = p
          if (no == 2) io2 = p
       else if (abs(occv(p)) < 1.0e-6_dp) then
          nv1 = nv1 + 1
       else
          call mrsf_error('fractional occupation in the reference')
       endif
    enddo
    if (no /= 2) call mrsf_error('the reference must have exactly two '&
         //'singly occupied MOs')

    allocate(hmap1(nc1+2), pmap1(nv1+2))
    ic = 0; iv = 0
    do p = 1, n
       if (abs(occv(p)-2.0_dp) < 1.0e-6_dp) then
          ic = ic + 1
          hmap1(ic) = p
       else if (abs(occv(p)) < 1.0e-6_dp) then
          iv = iv + 1
          pmap1(2+iv) = p
       endif
    enddo
    hmap1(nc1+1) = io1
    hmap1(nc1+2) = io2
    pmap1(1) = io1
    pmap1(2) = io2

  end subroutine classify_orbitals

!######################################################################
! setup_space: fills the module-level orbital maps and slot irreps;
! with the extended method the vector has ncol = nocca + nC columns
!######################################################################
  subroutine setup_space()

    integer(is) :: a, i, p, ref_irrep

    if (allocated(Hmap)) deallocate(Hmap)
    if (allocated(Pmap)) deallocate(Pmap)
    call classify_orbitals(nmo, occ, nC, nV, Hmap, Pmap, iO1, iO2)
    nocca = nC + 2
    nvirb = nV + 2
    xdim  = nocca * nvirb
    if (extended) then
       ncol = nocca + nC
       ncv  = nvirb * nC
    else
       ncol = nocca
       ncv  = 0
    endif
    xdim_tot = nvirb * ncol

    if (allocated(Hinv)) deallocate(Hinv)
    if (allocated(Pinv)) deallocate(Pinv)
    allocate(Hinv(nmo), Pinv(nmo))
    Hinv = 0; Pinv = 0
    do i = 1, nocca
       Hinv(Hmap(i)) = i
    enddo
    do a = 1, nvirb
       Pinv(Pmap(a)) = a
    enddo

    ! irrep of a spin-flipped configuration i -> a: the irrep of the
    ! open-shell triplet reference, Gamma(O1) x Gamma(O2), times
    ! Gamma(i) x Gamma(a) (abelian groups: products are XORs of the
    ! PySCF irrep ids). The reference factor is what distinguishes
    ! this from the closed-shell-reference rule used in bitci. The CV
    ! columns of the extended method are excitations of the closed-shell
    ! configuration G (totally symmetric): irrep = Gamma(i) x Gamma(a).
    ref_irrep = ieor(mosym(iO1), mosym(iO2))
    if (allocated(slot_irrep)) deallocate(slot_irrep)
    allocate(slot_irrep(xdim_tot))
    do i = 1, nocca
       do a = 1, nvirb
          slot_irrep(a + nvirb*(i-1)) = &
               ieor(ref_irrep, ieor(mosym(Pmap(a)), mosym(Hmap(i))))
       enddo
    enddo
    do i = nocca+1, ncol
       do a = 1, nvirb
          slot_irrep(a + nvirb*(i-1)) = ieor(mosym(Pmap(a)), mosym(Hmap(i-nocca)))
       enddo
    enddo
    nirrep = 1
    do p = 1, nmo
       nirrep = max(nirrep, mosym(p)+1)
    enddo

  end subroutine setup_space

!######################################################################
! slot_masked: .true. if the packed slot (a,i) is redundant for the
! given target multiplicity (local indices; i > nocca: CV columns, whose
! SOMO rows are always masked)
!######################################################################
  pure function slot_masked(mult, a, i) result(masked)

    integer(is), intent(in) :: mult, a, i
    logical                 :: masked

    masked = .false.
    if (i > nocca) then
       if (a <= 2) masked = .true.
       return
    endif
    if (a == 2 .and. i == nC+2) masked = .true.
    if (mult == 3) then
       if (a == 1 .and. i == nC+2) masked = .true.
       if (a == 2 .and. i == nC+1) masked = .true.
    endif

  end function slot_masked

!######################################################################
! slot_active: active for (mult, irrep); irrep < 0 means no irrep mask
!######################################################################
  pure function slot_active(mult, irrep, a, i) result(active)

    integer(is), intent(in) :: mult, irrep, a, i
    logical                 :: active

    active = .not. slot_masked(mult, a, i)
    if (irrep >= 0) then
       if (slot_irrep(a + nvirb*(i-1)) /= irrep) active = .false.
    endif

  end function slot_active

!######################################################################
! expand: compressed amplitudes X(nvirb,ncol,nvec) -> expanded Xt = U x
! (the CV columns are copied unchanged)
!######################################################################
  subroutine expand(mult, nvec, X, Xt)

    integer(is), intent(in) :: mult, nvec
    real(dp), intent(in)    :: X(nvirb,ncol,nvec)
    real(dp), intent(out)   :: Xt(nvirb,ncol,nvec)
    integer(is)             :: v, hO1, hO2
    real(dp)                :: xoo

    hO1 = nC + 1
    hO2 = nC + 2

    !$omp parallel do private(v, xoo)
    do v = 1, nvec
       Xt(:,:,v) = X(:,:,v)
       xoo = X(1,hO1,v)
       Xt(1,hO1,v) = xoo * isqrt2
       if (mult == 1) then
          Xt(2,hO2,v) = -xoo * isqrt2
       else
          Xt(2,hO2,v) = xoo * isqrt2
          Xt(1,hO2,v) = 0.0_dp
          Xt(2,hO1,v) = 0.0_dp
       endif
    enddo
    !$omp end parallel do

  end subroutine expand

!######################################################################
! fold: expanded sigma St -> compressed sigma (U^T St), in place;
! the redundant slots are zeroed (MRSF columns only)
!######################################################################
  subroutine fold(mult, nvec, St)

    integer(is), intent(in) :: mult, nvec
    real(dp), intent(inout) :: St(nvirb,ncol,nvec)
    integer(is)             :: v, hO1, hO2
    real(dp)                :: s11, s22

    hO1 = nC + 1
    hO2 = nC + 2

    do v = 1, nvec
       s11 = St(1,hO1,v)
       s22 = St(2,hO2,v)
       if (mult == 1) then
          St(1,hO1,v) = (s11 - s22) * isqrt2
       else
          St(1,hO1,v) = (s11 + s22) * isqrt2
          St(1,hO2,v) = 0.0_dp
          St(2,hO1,v) = 0.0_dp
       endif
       St(2,hO2,v) = 0.0_dp
    enddo

  end subroutine fold

!######################################################################
! apply_mask: zero the redundant slots and, if irrep >= 0, all slots
! not belonging to the irrep
!######################################################################
  subroutine apply_mask(mult, irrep, nvec, S)

    integer(is), intent(in) :: mult, irrep, nvec
    real(dp), intent(inout) :: S(nvirb,ncol,nvec)
    integer(is)             :: i, a, v, hO1, hO2

    hO1 = nC + 1
    hO2 = nC + 2

    S(2,hO2,:) = 0.0_dp
    if (mult == 3) then
       S(1,hO2,:) = 0.0_dp
       S(2,hO1,:) = 0.0_dp
    endif
    if (ncol > nocca) then
       S(1:2,nocca+1:ncol,:) = 0.0_dp
    endif

    if (irrep >= 0) then
       !$omp parallel do private(i, a, v)
       do v = 1, nvec
          do i = 1, ncol
             do a = 1, nvirb
                if (slot_irrep(a + nvirb*(i-1)) /= irrep) S(a,i,v) = 0.0_dp
             enddo
          enddo
       enddo
       !$omp end parallel do
    endif

  end subroutine apply_mask

end module mrsf_space
