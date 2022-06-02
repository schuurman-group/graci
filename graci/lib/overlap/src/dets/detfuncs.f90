!**********************************************************************
! Sundry determinant manipulation routines
!**********************************************************************
module detfuncs

  implicit none

contains

!######################################################################
! truncate_wave_functions: truncation of the input wave functions
!                          based on a minimum norm threshold
!######################################################################
  subroutine truncate_wave_functions(n_int,ndet,nroots,det,vec,&
       normthrsh,ndet_new,det_new,vec_new)

    use constants
    use global
    use utils
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: n_int,ndet,nroots

    ! Determinant bit strings
    integer(ib), intent(in)  :: det(n_int,2,ndet)

    ! Eigenvectors
    real(dp), intent(in)     :: vec(ndet,nroots)

    ! Norm-based truncation threshold
    real(dp), intent(in)     :: normthrsh

    ! Truncated determinant and eigenvector arrays
    integer(is), intent(out) :: ndet_new
    integer(ib), allocatable :: det_new(:,:,:)
    real(dp), allocatable    :: vec_new(:,:)
    
    ! Sorting
    integer(is), allocatable :: indx(:)
    real(dp), allocatable    :: cabs(:)

    ! Surviving determinants
    integer(is), allocatable :: idet(:)
    
    ! Everything else
    integer(is)              :: i,k,n
    real(dp)                 :: normsq,targ
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(ndet))
    indx=0

    allocate(cabs(ndet))
    cabs=0.0d0

    allocate(idet(ndet))
    idet=0
    
!----------------------------------------------------------------------
! Determine the indices of the surviving determinants
!----------------------------------------------------------------------
    ! Target squared norm
    targ=normthrsh**2

    ! Loop over roots
    do i=1,nroots

       ! Sort the coefficients by absolute value
       cabs=abs(vec(:,i))
       call dsortindxa1('D',ndet,cabs,indx)

       ! Fill in the surviving determinants for this state
       normsq=0.0d0
       do k=1,ndet

          ! Update the squared norm
          normsq=normsq+vec(indx(k),i)**2

          ! Flag the determinant for survival
          idet(indx(k))=1

          ! Exit if we have hit the target squared norm
          if (normsq >= targ) exit

       enddo
       
    enddo

!----------------------------------------------------------------------
! Set up the truncated determinant and eigenvector arrays
!----------------------------------------------------------------------
    !
    ! Number of surviving determinants
    !
    ndet_new=sum(idet)

    !
    ! Allocate arrays
    !
    allocate(det_new(n_int,2,ndet_new))
    det_new=0_ib
    allocate(vec_new(ndet_new,nroots))
    vec_new=0.0d0

    !
    ! Fill in the determinant bit strings
    !
    n=0
    ! Loop over the original set of determinants
    do k=1,ndet

       ! Are we at a surviving determinant?
       if (idet(k) == 1) then
          n=n+1
          det_new(:,:,n)=det(:,:,k)
       endif
       
    enddo

    !
    ! Fill in the eigenvectors
    !
    ! Loop over roots
    do i=1,nroots

       ! Surviving determinant counter
       n=0
       
       ! Loop over the original set of determinants
       do k=1,ndet

          ! Are we at a surviving determinant?
          if (idet(k) == 1) then
             n=n+1
             vec_new(n,i)=vec(k,i)
          endif
             
       enddo

    enddo
    
    return
    
  end subroutine truncate_wave_functions

!######################################################################
! get_nel: determines the number of electrons (total, alpha and beta)
!          in a given determinant d
!######################################################################
  subroutine get_nel(n_int,d,nel,nel_alpha,nel_beta)

    use constants
    
    implicit none

    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: d(n_int,2)
    integer(is), intent(out) :: nel,nel_alpha,nel_beta

    integer(is)              :: k

    !
    ! Number of alpha electrons
    !
    nel_alpha=0
    do k=1,n_int
       nel_alpha=nel_alpha+popcnt(d(k,1))
    enddo
    
    !
    ! Number of beta electrons
    !
    nel_beta=0
    do k=1,n_int
       nel_beta=nel_beta+popcnt(d(k,2))
    enddo

    !
    ! Total number of electrons
    !
    nel=nel_alpha+nel_beta
    
    return
    
  end subroutine get_nel
  
!######################################################################
  
end module detfuncs