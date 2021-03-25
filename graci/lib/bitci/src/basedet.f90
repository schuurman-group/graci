!######################################################################
! generate_base_det: generates a determinant corresponding to the base
!                    spatial occupation pattern
!######################################################################
subroutine generate_base_det
  
  use constants
  use bitglobal
  use bitutils
  
  implicit none
  
  integer(is) :: imo,imo1,indx,iab,ishift
  integer(is) :: ndocc,nopen
  
!----------------------------------------------------------------------
! Allocate and initialise the base determinant to 00...0
!----------------------------------------------------------------------
  allocate(det0(n_int,2))
  det0=0_ib

!----------------------------------------------------------------------
! Number of doubly-occupied and singly-occupied spatial orbitals in
! the base determinant(s)
!----------------------------------------------------------------------
  nopen=imult-1
  ndocc=(nel-nopen)/2
  
!----------------------------------------------------------------------
! Set up one of the base determinants following the aufbau principle
!----------------------------------------------------------------------
  !
  ! Doubly-occupied spatial orbitals
  !
  do imo=1,ndocc

     ! Index of the element of the det0 array in which the current
     ! spatial orbital resides
     indx=(imo-1)/64+1
     
     ! Alpha electron
     det0(indx,1)=det0(indx,1)+2_ib**(imo-1-(indx-1)*64)
     
     ! Beta electron
     det0(indx,2)=det0(indx,2)+2_ib**(imo-1-(indx-1)*64)
     
  enddo
  
  !
  ! Singly-occupied orbitals (alpha spin by default)
  !
  ! Loop over the unpaired electrons
  do imo=ndocc+1,ndocc+nopen
     
     ! Index of the element of the det0 array in which the current
     ! spatial orbital resides
     indx=(imo-1)/64+1
     
     ! Add the unpaired electron
     det0(indx,1)=det0(indx,1)+2_ib**(imo-1-(indx-1)*64)
     
  enddo

  return
  
end subroutine generate_base_det
