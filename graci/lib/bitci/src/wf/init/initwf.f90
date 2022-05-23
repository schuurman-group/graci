!######################################################################
! bitwf_initialise: Interface to the bitwf library. Sets all the
!                   globally accessible variables that are needed to
!                   evaluate matrix elements <psi_m|O|psi'_n> in the
!                   Slater deteminant representation.
!######################################################################
#ifdef CBINDING
subroutine bitwf_initialise(imultB1,imultK1,nelB1,nelK1,nmoB1,nmoK1,&
     Smat1,ipgB1,ipgK1,calctype_in) bind(c,name='bitwf_initialise')
#else
subroutine bitwf_initialise(imultB1,imultK1,nelB1,nelK1,nmoB1,nmoK1,&
     Smat1,ipgB1,ipgK1,calctype_in)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use iomod

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

  ! Bra and ket point group indices
  integer(is), intent(in) :: ipgB1,ipgK1
  
  ! Everything else
  integer(is)             :: i,j
  real(dp)                :: s,smax
  logical                 :: verbose
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the calculation type character
! string from the C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(calctype_in)
  call c2fstr(calctype_in,calctype,length)
#else
  calctype=adjustl(trim(calctype_in))
#endif

!----------------------------------------------------------------------
! Quick sanity check on the numbers of electrons and spin
! multiplicities
!----------------------------------------------------------------------
  ! Bra
  if (mod(imultB1,2) == 0 .and. mod(nelB1,2) == 0 &
       .or. mod(imultB1,2) /= 0.and. mod(nelB1,2) /= 0) then
     write(errmsg,'(a,1x,i0,a,i0)') &
          'Nonsensical bra Nel/multiplicity:',nelB1,'/',imultB1
     call error_control
  endif

  ! Ket
  if (mod(imultK1,2) == 0 .and. mod(nelK1,2) == 0 &
       .or. mod(imultK1,2) /= 0.and. mod(nelK1,2) /= 0) then
     write(errmsg,'(a,1x,i0,a,i0)') &
          'Nonsensical ket Nel/multiplicity:',nelK1,'/',imultK1
     call error_control
  endif

!----------------------------------------------------------------------
! Exit if the requested spin multiplicity is not supported by the
! maximum number of open shells
!----------------------------------------------------------------------
  ! S_max
  smax=dble(nomax)/2.0d0

  ! Bra
  s=dble(imultB1-1)/2.0
  if (s > smax) then
     write(errmsg,'(a)') &
          'The requested bra multiplicity is incompatible'&
          //' with the current value of nomax.'
     call error_control
  endif

  ! Ket
  s=dble(imultK1-1)/2.0
  if (s > smax) then
     write(errmsg,'(a)') &
          'The requested ket multiplicity is incompatible'&
          //' with the current value of nomax.'
     call error_control
  endif

!----------------------------------------------------------------------
! Exit if the given point groups are not recognised
!----------------------------------------------------------------------
  ! Bra
  if (ipgB1 < 1 .or. ipgB1 > 8) then
     write(errmsg,'(a,1x,i0)') 'Illegal bra point group index:',ipgB1
     call error_control
  endif

  ! Ket
  if (ipgK1 < 1 .or. ipgK1 > 8) then
     write(errmsg,'(a,1x,i0)') 'Illegal ket point group index:',ipgB1
     call error_control
  endif

!----------------------------------------------------------------------
! Set the bra and ket spin multiplicities
!----------------------------------------------------------------------
  imultB=imultB1
  imultK=imultK1

!----------------------------------------------------------------------
! Set the bra and ket numbers of electrons
!----------------------------------------------------------------------
  nelB=nelB1
  nelB_beta=(nelB-imultB+1)/2
  nelB_alpha=(nelB-imultB+1)/2+imultB-1
  nelB_spin(1)=nelB_alpha
  nelB_spin(2)=nelB_beta

  nelK=nelK1
  nelK_beta=(nelK-imultK+1)/2
  nelK_alpha=(nelK-imultK+1)/2+imultK-1
  nelK_spin(1)=nelK_alpha
  nelK_spin(2)=nelK_beta

!----------------------------------------------------------------------
! Set the number of spatial orbitals
!----------------------------------------------------------------------
  ! We now need nmoB and nmoK variables
  ! Also: n_int can be different for the bra and ket WFs...
  
  STOP
  
  return
  
end subroutine bitwf_initialise
