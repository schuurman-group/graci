!**********************************************************************
! Construction of the base determinant and configuration
!**********************************************************************
module baseconf

  implicit none

contains

!######################################################################
! generate_base_det: Generates a determinant corresponding to the base
!                    spatial occupation pattern. For a given spin
!                    multiplicity, we take the maximum possible m_s
!                    value
!######################################################################
subroutine generate_base_det
  
  use constants
  use bitglobal
    
  implicit none
  
  integer(is) :: imo,i,k
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
! Set up one the base determinants following the aufbau principle and
! taking the maximum possible ms value
!----------------------------------------------------------------------
  !
  ! Doubly-occupied spatial orbitals
  !
  do imo=1,ndocc

     ! Block index
     k=(imo-1)/64+1

     ! Orbital index
     i=imo-1-(k-1)*64
     
     ! Add the alpha electron
     det0(k,1)=ibset(det0(k,1),i)
     
     ! Add the beta electron
     det0(k,2)=ibset(det0(k,2),i)

  enddo
  
  !
  ! Singly-occupied (alpha) orbitals
  !
  ! Loop over the unpaired electrons
  do imo=ndocc+1,ndocc+nopen
     
     ! Block index
     k=(imo-1)/64+1

     ! Orbital index
     i=imo-1-(k-1)*64
     
     ! Add the unpaired alpha electron
     det0(k,1)=ibset(det0(k,1),i)
     
  enddo
  
  return
  
end subroutine generate_base_det
  
!######################################################################
! generate_base_conf: Generates a configuration bit string pair
!                     corresponding to the base spatial occupation
!######################################################################
  subroutine generate_base_conf

    use constants
    use bitglobal
        
    implicit none

    integer(is) :: imo,i,k
    integer(is) :: ndocc,nopen
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    allocate(conf0(n_int,2))
    conf0=0_ib

    allocate(iopen0(nmo))
    iopen0=0
    
!----------------------------------------------------------------------
! Number of doubly-occupied and singly-occupied spatial orbitals in
! the base determinant(s)
!----------------------------------------------------------------------
    nopen=imult-1
    ndocc=(nel-nopen)/2

!----------------------------------------------------------------------
! Construct the base configuration following the aufbau principle
!----------------------------------------------------------------------
    !
    ! Doubly-occupied orbitals
    !
    do imo=1,ndocc

       ! Block index
       k=(imo-1)/64+1

       ! Orbital index
       i=imo-(k-1)*64-1

       ! Set the bits
       conf0(k,1)=ibset(conf0(k,1),i)       
       conf0(k,2)=ibset(conf0(k,2),i)
       
    enddo

    !
    ! Singly-occupied orbitals
    !
    ! Loop over the unpaired electrons
    do imo=ndocc+1,ndocc+nopen

       ! Save the indices of the open-shell MOs
       iopen0(imo)=1
       
       ! Block index
       k=(imo-1)/64+1
       
       ! Orbital index
       i=imo-1-(k-1)*64
       
       ! Set the bit
       conf0(k,1)=ibset(conf0(k,1),i)       
       
    enddo
    
    return
    
  end subroutine generate_base_conf
    
!######################################################################
  
end module baseconf
