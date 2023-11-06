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
    use detfuncs, only: list_from_bitstring
    
    implicit none

    ! Determinant bit strings
    integer(is), intent(in)    :: n_int,ndet
    integer(ib), intent(inout) :: det(n_int,2,ndet)

    ! Bit string encoding of the fixed-occupation spin orbitals
    integer(ib)                :: fixed(n_int,2)

    ! Number of fixed-occupation spin orbitals
    integer(is)                :: nfixed(2)

    ! Fixed-occupation spin orbital indices
    integer(is), allocatable   :: ifixed_alpha(:),ifixed_beta(:)
    
    ! Everything else
    integer(is)                :: n,k
    integer(ib)                :: string(n_int)
    
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
    ! Allocate arrays
    allocate(ifixed_alpha(nfixed(1)))
    allocate(ifixed_beta(nfixed(2)))
    ifixed_alpha=0
    ifixed_beta=0
    
    ! Alpha spin
    string=fixed(:,1)
    call list_from_bitstring(n_int,string,ifixed_alpha,nfixed(1))

    ! Beta spin
    string=fixed(:,2)
    call list_from_bitstring(n_int,string,ifixed_beta,nfixed(2))

!----------------------------------------------------------------------
! Re-order the alpha and beta strings, computing the associated phase
! factors as we go
!----------------------------------------------------------------------
! Important: We are going to treat the orbital re-ordering as a
!            series of pairwise pi/2 orbital rotations, *not* orbital
!            permutations
!----------------------------------------------------------------------
    
    
    return
    
  end subroutine det_reorder
  
!######################################################################
  
end module reorder
