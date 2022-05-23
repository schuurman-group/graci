!**********************************************************************
! Construction of the base determinant and configuration
!**********************************************************************
! *** IMPORTANT ***
!----------------------------------------------------------------------
! The base determinant has the same spin multiplicity as the computed
! eigenstates
!----------------------------------------------------------------------
! However, for technical reasons associated with the application of
! DFT/MRCI corrections, the base configuration will always correspond
! to either:
!
! (1) the singlet occupation string 2...22 if mod(imult,2) = 1
!
! or
!
! (2) the doublet occupation string 2...21 if mod(imult,2) = 0
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
     k=(imo-1)/n_bits+1

     ! Orbital index
     i=imo-1-(k-1)*n_bits
     
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
     k=(imo-1)/n_bits+1

     ! Orbital index
     i=imo-1-(k-1)*n_bits
     
     ! Add the unpaired alpha electron
     det0(k,1)=ibset(det0(k,1),i)
     
  enddo
  
  return
  
end subroutine generate_base_det
  
!######################################################################
! generate_base_conf: Generates a configuration bit string pair
!                     corresponding to the base spatial occupation
!
!                     For mod(imult,2) = 0, this is the closed shell
!                     configuration 2...22
!
!                     For mod(imult,2) = 1, this is the open shell
!                     configuration 2...21
!######################################################################
  subroutine generate_base_conf

    use constants
    use bitglobal
        
    implicit none

    integer(is) :: ihomo,imo,i,k
    integer(is) :: ndocc,nopen
    
!----------------------------------------------------------------------
! Allocate and initialise arrays
!----------------------------------------------------------------------
    allocate(conf0(n_int,2))
    conf0=0_ib

    allocate(iopen0(nmo))
    iopen0=0

    allocate(iocc0(nmo))
    iocc0=0
    
!----------------------------------------------------------------------
! Construct the base configuration
!----------------------------------------------------------------------
    select case(mod(imult,2))

    case(1)
       ! mod(imult,2) = 1 <-> |2...22>
       
       ! Loop over occupied MOs
       do imo=1,nel/2

          ! Block index
          k=(imo-1)/n_bits+1

          ! Orbital index
          i=imo-(k-1)*n_bits-1
          
          ! Set the bits
          conf0(k,1)=ibset(conf0(k,1),i)       
          conf0(k,2)=ibset(conf0(k,2),i)

          ! Save the MO occupation
          iocc0(imo)=2
          
       enddo
       
    case(0)
       ! mod(imult,2) = 0 <-> |2...21>

       ! HOMO index
       ihomo=(nel+1)/2

       ! Loop over occupied MOs
       do imo=1,ihomo

          ! Block index
          k=(imo-1)/n_bits+1
             
          ! Orbital index
          i=imo-(k-1)*n_bits-1
             
          if (imo == ihomo) then
             ! Singly-occupied MO
             conf0(k,1)=ibset(conf0(k,1),i)
             iopen0(imo)=1
             iocc0(imo)=1
          else
             ! Doubly-occupied MO
             conf0(k,1)=ibset(conf0(k,1),i)       
             conf0(k,2)=ibset(conf0(k,2),i)
             iocc0(imo)=2
          endif
          
       enddo
       
    end select
    
    return
    
  end subroutine generate_base_conf
    
!######################################################################
  
end module baseconf
