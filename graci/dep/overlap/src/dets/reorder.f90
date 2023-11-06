!**********************************************************************
! Routines for the re-ordering of the unique alpha and beta strings
! s.t. the fixed-occupation orbitals come first
!**********************************************************************
module reorder

  implicit none

contains

!######################################################################
! det_reorder: re-orders a given set of determinants s.t. the fixed-
!               occupation orbitals come first
!######################################################################
  subroutine det_reorder(n_int,ndet,det)

    use constants
    
    implicit none

    ! Determinant bit strings
    integer(is), intent(in)    :: n_int,ndet
    integer(ib), intent(inout) :: det(n_int,2,ndet)

    ! Bit string encoding of the fixed-occupation spin orbitals
    integer(ib)                :: fixed(n_int,2)

    ! Number of fixed-occupation spin orbitals
    integer(is)                :: nfixed(2)
    
    ! Everything else
    integer(is)                :: n,k
    
!----------------------------------------------------------------------
! Construct the bit string encodings of the fixed-occupation alpha-
! and beta-spin orbitals
!----------------------------------------------------------------------
    ! Initialisation
    fixed(:,:)=det(:,:,1)

    ! Loop over determinants
    do n=2,ndet

       ! Alpha spin
       do k=1,n_int
          fixed(k,1)=iand(fixed(k,1),det(k,1,n))
       enddo
          
       ! Beta spin
       do k=1,n_int
          fixed(k,2)=iand(fixed(k,2),det(k,2,n))
       enddo
          
    enddo

!----------------------------------------------------------------------
! Numbers of alpha and beta fixed-occupation spin orbitals
!----------------------------------------------------------------------
    ! Initialisation
    nfixed=0
    
    ! Alpha spin
    do k=1,n_int
       nfixed(1)=nfixed(1)+popcnt(fixed(k,1))
    enddo

    ! Beta spin
    do k=1,n_int
       nfixed(2)=nfixed(2)+popcnt(fixed(k,2))
    enddo

!----------------------------------------------------------------------
! Determine the indices of the fixed occupation spin-orbitals
!----------------------------------------------------------------------
    
    
    return
    
  end subroutine det_reorder
  
!######################################################################
  
end module reorder
