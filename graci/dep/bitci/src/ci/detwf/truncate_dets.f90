!**********************************************************************
! Routines for the truncation of the determinant representation of
! the MRCI wave functions
!**********************************************************************

!######################################################################
! ndet_truncated: Given a wave function scratch file number, determines
!                 the number of determinants required to reach a
!                 norm-based truncation threshold normthrsh
!######################################################################
#ifdef CBINDING
subroutine ndet_truncated(wfscr,nroots,normthrsh,ndet_trunc) &
     bind(c,name="ndet_truncated")
#else
subroutine ndet_truncated(wfscr,nroots,normthrsh,ndet_trunc)
#endif

  use constants
  use bitglobal
  use utils
  use iomod
  
  implicit none

  ! Determinant wave function scratch file number
  integer(is), intent(in)  :: wfscr

  ! No. roots
  integer(is), intent(in)  :: nroots
  
  ! Norm-based truncation threshold
  real(dp), intent(in)     :: normthrsh

  ! Truncated no. determinants
  integer(is), intent(out) :: ndet_trunc

  ! Untruncated determinants
  integer(is)              :: ndet

  ! Wave functions
  integer(ib), allocatable :: det(:,:,:)
  real(dp), allocatable    :: vec(:,:)
  integer(is), allocatable :: iroots(:)
  
  ! Sorting
  integer(is), allocatable :: indx(:)
  real(dp), allocatable    :: cabs(:)

  ! Surviving determinants
  integer(is), allocatable :: idet(:)
  
  ! Everything else
  integer(is)              :: i,k
  real(dp)                 :: targ,normsq
  
!----------------------------------------------------------------------
! Get the untruncated number of determinants
!----------------------------------------------------------------------
  call read_ndet(wfscr,ndet)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(det(n_int,2,ndet))
  det=0_ib

  allocate(vec(ndet,nroots))
  vec=0.0d0

  allocate(iroots(nroots))
  iroots=0

  allocate(indx(ndet))
  indx=0

  allocate(cabs(ndet))
  cabs=0.0d0

  allocate(idet(ndet))
  idet=0
  
!----------------------------------------------------------------------
! Read in the eigenvectors
!----------------------------------------------------------------------
  ! Read all the eigenvectors saved to disk
  do i=1,nroots
     iroots(i)=i
  enddo

  call read_detwf(wfscr,ndet,nroots,n_int,det,vec,iroots)

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
! Number of surviving determinants
!----------------------------------------------------------------------
  ndet_trunc=sum(idet)
  
  return
  
end subroutine ndet_truncated

!######################################################################
! retrieve_det_truncated: Construction of truncated determinant bit
!                         strings and eigenvectors based on a minimum
!                         norm threshold
!
!                         Note that, in order to support calling this
!                         routine from outside of Fortran, the
!                         truncated no. determinants, ndet_trunc,
!                         must be known ahead of time. To facilitate
!                         this, the ndet_truncated subroutine may be
!                         called beforehand.
!######################################################################
#ifdef CBINDING
subroutine retrieve_det_truncated(wfscr,ndet,ndet_trunc,nroots,&
     det_trunc,vec_trunc,normthrsh) &
     bind(c,name="retrieve_det_truncated")
#else
subroutine retrieve_det_truncated(wfscr,ndet,ndet_trunc,nroots,&
     det_trunc,vec_trunc,normthrsh)
#endif
  
  use constants
  use bitglobal
  use utils
  use iomod
  
  implicit none

  ! Wave function scratch file number
  integer(is), intent(in)  :: wfscr
  
  ! Dimensions
  integer(is), intent(in)  :: ndet,ndet_trunc,nroots

  ! Truncated determinant bit strings and eigenvectors
  integer(ib), intent(out) :: det_trunc(n_int,2,ndet_trunc)
  real(dp), intent(out)    :: vec_trunc(ndet_trunc,nroots)
  
  ! Norm-based truncation threshold
  real(dp), intent(in)     :: normthrsh

  ! Untruncated determinant bit strings and eigenvectors
  integer(ib), allocatable :: det(:,:,:)
  real(dp), allocatable    :: vec(:,:)
  integer(is), allocatable :: iroots(:)

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
  allocate(det(n_int,2,ndet))
  det=0_ib

  allocate(vec(ndet,nroots))
  vec=0.0d0

  allocate(iroots(nroots))
  iroots=0

  allocate(indx(ndet))
  indx=0

  allocate(cabs(ndet))
  cabs=0.0d0

  allocate(idet(ndet))
  idet=0

!----------------------------------------------------------------------
! Read in the untruncated eigenvectors
!----------------------------------------------------------------------
  do i=1,nroots
     iroots(i)=i
  enddo

  call read_detwf(wfscr,ndet,nroots,n_int,det,vec,iroots)
  
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
! Sanity check on the number surviving determinants
!----------------------------------------------------------------------
  if (sum(idet) /= ndet_trunc) then
     errmsg='Error in retrieve_det_truncated: '&
          //'incorrect value of ndet_trunc'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Fill in the truncated determinant and eigenvector arrays
!----------------------------------------------------------------------
  !
  ! Fill in the determinant bit strings
  !
  n=0
  ! Loop over the original set of determinants
  do k=1,ndet

     ! Are we at a surviving determinant?
     if (idet(k) == 1) then
        n=n+1
        det_trunc(:,:,n)=det(:,:,k)
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
           vec_trunc(n,i)=vec(k,i)
        endif
             
     enddo
     
  enddo
    
  return
    
end subroutine retrieve_det_truncated
