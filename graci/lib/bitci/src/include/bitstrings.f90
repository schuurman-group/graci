!**********************************************************************
! Hard-wired bitstring encodings of the spin-coupling coefficient
! sub-cases and sub-case pairs.
!
! This is a rather horrible way to do things, but it means that these
! variables can be used in select case statements.
!**********************************************************************
module bitstrings
  
  use constants

  implicit none

!**********************************************************************
! First, a bit of documentation
!**********************************************************************
! Case 1a: singly-occupied to unoccupied excitation
!
! Case 1b: doubly-occupied to singly-occupied excitation
!
! Case 2a: doubly-occupied to unoccupied excitation
!
! Case 2b: singly-occupied to singly-occupied excitation
!
! Each spin-coupling coefficient corresponds to a singlet excitation
! operator E_a^i operating on a CSF |w omega>.
! The type of MO (unoccupied, singly-occupied, or doubly-occupied)
! indexed by annihilator index i and creator index a will be encoded
! as follows:
!
! Unoccupied:      00
! Singly-occupied: 10
! Doubly-occupied: 11
!
! The bitstring encodings of the combined annihilator/creator types
! will be the concatenation of the encodings of each.
!
! e.g., 1000 corresponds to excitation from a singly-occupied MO to
! to an unoccupied one (case 1a), and 1110 to excitation from a
! doubly-occupied MO to a singly-occupied MO (case1b).
!
! The numerical values of these bitstrings will be used in the
! Hamiltonian build to select the spin coupling sub-cases.
!**********************************************************************

!----------------------------------------------------------------------
! Bitstring encodings for single case types
!----------------------------------------------------------------------
! Case 1a: 1000 (1)
! Case 1b: 1110 (7)
! Case 2a: 1100 (3)
! Case 2b: 1010 (5)
!----------------------------------------------------------------------
  integer(ib), parameter :: i1a=1
  integer(ib), parameter :: i1b=7
  integer(ib), parameter :: i2a=3
  integer(ib), parameter :: i2b=5

!----------------------------------------------------------------------
! Spin-coupling sub-case labels
!----------------------------------------------------------------------
  character(len=3), dimension(7) :: caselbl=&
       ['1a','  ','2a','  ','2b','  ','1b']
  
!----------------------------------------------------------------------
! Bitstring encodings for pairs of case types:
!----------------------------------------------------------------------
! 1a,1a: 1000|1000 (17)
! 1a,1b: 1000|1110 (113)
! 1a,2a: 1000|1100 (49)
! 1a,2b: 1000|1010 (81)
! 1b,1a: 1110|1000 (23)
! 1b,1b: 1110|1110 (119)
! 1b,2a: 1110|1100 (55)
! 1b,2b: 1110|1010 (87)
! 2a,1a: 1100|1000 (19)
! 2a,1b: 1100|1110 (115)
! 2a,2a: 1100|1100 (51)
! 2a,2b: 1100|1010 (83)
! 2b,1a: 1010|1000 (21)
! 2b,1b: 1010|1110 (117)
! 2b,2a: 1010|1100 (53)
! 2b,2b: 1010|1010 (85)
!----------------------------------------------------------------------
  integer(ib), parameter, dimension(7,7) :: ipair= &
       reshape( &
       [17,0,19,0,21,0,23, &
       0,0,0,0,0,0,0, &
       49,0,51,0,53,0,55, &
       0,0,0,0,0,0,0, &
       81,0,83,0,85,0,87, &
       0,0,0,0,0,0,0, &
       113,0,115,0,117,0,119], &
       [7,7])
  
end module bitstrings
