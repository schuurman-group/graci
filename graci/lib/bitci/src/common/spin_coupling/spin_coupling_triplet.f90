!**********************************************************************
! Routines for the pre-computation of the unique spin-coupling
! coefficients <w omega| T_pq^(1,k) |w' omega'> for triplet spin tensor
! operators T_pq^(1,k), k=0,+1
!**********************************************************************
! Follows the formalism detailed in
!
! R. W. Wetmore and G. A. Segal, Chem. Phys. Lett., 36, 478 (1975)
!
! and
!
! M. Kleinschmidt and C. M. Marian, Chem. Phys, 311, 71 (2005)
!**********************************************************************
module spin_coupling_triplet

  implicit none

contains
  
!######################################################################
! generate_triplet_coupling_coefficients:
!
! Generates the unique triplet spin-coupling coefficients for a given
! pair of bra and ket spin multiplicities
!
! Considers a single component of the triplet spin tensor operator
! T^(1,k) indexed by the input variable kindx 
!######################################################################
  subroutine generate_triplet_coupling_coefficients(kindx,&
       imultB,imultK,&
       nocase1,nocase2,&
       maxcsfB,maxdetB,ncsfsB,ndetsB,csfcoeB,detvecB,&
       maxcsfK,maxdetK,ncsfsK,ndetsK,csfcoeK,detvecK,&
       npattern1,npattern2,nspincp,&
       N1s,verbose,spincp,patternmap,offspincp)

    use constants
    use timing

    ! Component of the triplet spin tensor operator to use
    integer(is), intent(in) :: kindx

    ! Bra and ket spin multiplicities
    integer(is), intent(in) :: imultB,imultK

    ! Maximum number open shells for the Case 1 and Case 2 CSFs
    integer(is), intent(in)  :: nocase1,nocase2

    ! Numbers of bra and ket CSFs and determinants as a function
    ! of the number of open shells
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ndetsB(0:nocase2)
    integer(is), intent(in)  :: ncsfsK(0:nocase2),ndetsK(0:nocase2)
    
    ! Maximum number of bra and ket CSFs/determinants across all
    ! numbers of open shells
    integer(is), intent(in)  :: maxcsfB,maxdetB
    integer(is), intent(in)  :: maxcsfK,maxdetK
    
    ! Bra and ket CSF expansion coefficients
    real(dp), intent(in)     :: csfcoeB(maxcsfB,maxdetB,nocase2)
    real(dp), intent(in)     :: csfcoeK(maxcsfK,maxdetK,nocase2)
    
    ! Bit string encoding of the determinants contributing to the
    ! bra and ket CSFs
    integer(ib), intent(in)  :: detvecB(maxdetB,nocase2)
    integer(ib), intent(in)  :: detvecK(maxdetK,nocase2)
    
    ! Number of Case 1 and Case 2 patterns
    integer(is), intent(out) :: npattern1,npattern2

    ! Number of unique spin coupling coefficients
    integer(is), intent(out) :: nspincp(2)

    ! Bit strings comprised of N 1's
    integer(ib), allocatable :: N1s(:)

    ! Verbose output
    logical, intent(in)      :: verbose

    ! All spin coupling coefficients
    integer(is)              :: spincpdim(3)
    real(dp), allocatable    :: spincp(:)

    ! All pattern -> array index mappings
    integer(is)              :: mapdim
    integer(is), allocatable :: patternmap(:)

    ! Spin coupling coefficient offsets for the various
    ! different cases

    ! THIS WILL NEED RE-DIMENSIONING DEPENDING ON THE
    ! SPIN-TENSOR OPERATOR CONSIDERED
    ! SO, IN GENERAL, WE SHOULD MAKE THIS ARRAY ALLOCATABLE
    ! IN BITGLOBAL
    integer(is), intent(out) :: offspincp(4)


    
    print*,''
    print*,'here...'
    print*,''
    stop
    
    return
    
  end subroutine generate_triplet_coupling_coefficients
  
!######################################################################
  
end module spin_coupling_triplet
