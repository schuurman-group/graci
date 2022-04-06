!######################################################################
! bitwf_initialise: Interface to the bitwf library. Sets all the
!                   globally accessible variables that are needed to
!                   evaluate matrix elements <psi_m|O|psi'_n> in the
!                   Slater deteminant representation.
!######################################################################
#ifdef CBINDING
subroutine bitwf_initialise(imultB1,imultK1,nelB1,nelK1,nmoB1,nmoK1,&
     Smat1,calctype_in) bind(c,name='bitwf_initialise')
#else
subroutine bitwf_initialise(imultB1,imultK1,nelB1,nelK1,nmoB1,nmoK1,&
     Smat1,calctype_in)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal

  implicit none

  ! Calculation type
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: calctype_in(*)
  character(len=255)                 :: calctype
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: calctype_in
  character(len=255)                 :: calctype
#endif
  
  ! Bra and Ket spin multiplicities and numbers of electrons
  integer(is), intent(in) :: imultB1,imultK1,nelB1,nelK1

  ! Number of bra and ket MOs
  integer(is), intent(in) :: nmoB1,nmoK1

  ! Bra-ket MO overlap matrix
  real(dp), intent(in)    :: Smat1(nmoB1,nmoK1)

  ! Everything else
  integer(is)             :: i,j

  do i=1,nmoB1
     do j=1,nmoK1
        print*,i,j,Smat1(i,j)
     enddo
  enddo
  
  stop
  
  return
  
end subroutine bitwf_initialise
