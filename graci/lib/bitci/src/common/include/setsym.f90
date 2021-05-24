!**********************************************************************
! Loading of the symmetry information
!**********************************************************************
module setsym

  implicit none

contains
  
!######################################################################
! initialise_symmetry: sets point group and irrep labels, and the the
!                      point group multiplication array
!######################################################################
  subroutine initialise_symmetry

    use constants
    use bitglobal
    use iomod
  
    implicit none
  
    integer(is)      :: i
    character(len=3) :: all_labels(8)
  
!----------------------------------------------------------------------
! Set the point group labels
!----------------------------------------------------------------------
    pglbls=(/'c1 ','ci ','c2 ','cs ','c2h','c2v','d2 ','d2h'/)

!----------------------------------------------------------------------
! Set the irrep labels
!----------------------------------------------------------------------
    ! C1
    pgdim(1)=1
    irreplbl(0,1)='A'
  
    ! Ci
    pgdim(2)=2
    irreplbl(0,2)='Ag'
    irreplbl(1,2)='Au'
    
    ! C2
    pgdim(3)=2
    irreplbl(0,3)='A'
    irreplbl(1,3)='B'
    
    ! Cs
    pgdim(4)=2
    irreplbl(0,4)='A'''
    irreplbl(1,4)='A'''''
  
    ! C2h
    pgdim(5)=4
    irreplbl(0,5)='Ag'
    irreplbl(1,5)='Bg'
    irreplbl(2,5)='Au'
    irreplbl(3,5)='Bu'
  
    ! C2v
    pgdim(6)=4
    irreplbl(0,6)='A1'
    irreplbl(1,6)='A2'
    irreplbl(2,6)='B1'
    irreplbl(3,6)='B2'
  
    ! D2
    pgdim(7)=4
    irreplbl(0,7)='A1'
    irreplbl(1,7)='B1'
    irreplbl(2,7)='B2'
    irreplbl(3,7)='B3'
    
    ! D2h
    pgdim(8)=8
    irreplbl(0,8)='A1g'
    irreplbl(1,8)='B1g'
    irreplbl(2,8)='B2g'
    irreplbl(3,8)='B3g'
    irreplbl(4,8)='A1u'
    irreplbl(5,8)='B1u'
    irreplbl(6,8)='B2u'
    irreplbl(7,8)='B3u'
  
    return
  
  end subroutine initialise_symmetry

!######################################################################
  
end module setsym
