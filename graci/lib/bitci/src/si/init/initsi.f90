!######################################################################
! bitsi_initialise: Interface to the bitsi library. Sets all the
!                   globally accessible variables that are needed to
!                   perform a State Interaction calculation.
!######################################################################
#ifdef CBINDING
subroutine bitsi_intialise(imultB1,imultK1,nelB1,nelK1) &
     bind(c,name='bitsi_initialise')
#else
subroutine bitsi_intialise(imultB1,imultK1,nelB1,nelK1)
#endif

  use constants
  use bitglobal
  use iomod
  
  implicit none

  ! Bra and Ket spin multiplicities and numbers of electrons
  integer(is), intent(in) :: imultB1,imultK1,nelB1,nelK1

  ! Everything else
  real(dp)                :: s,smax
  logical                 :: lsoc
  
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
! maximum number of openshells
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
! Is this a spin-orbit coupling calculation?
!----------------------------------------------------------------------
  if (abs(imultB-imultK) == 2) then
     ! SOC calculation
     lsoc=.true.
  else if (imultB == imultK) then
     ! 1-RDM or 1-TDM calculation
     lsoc=.false.
  else
     ! Non-sensical combination of bra and ket spin multiplicities
     write(errmsg,'(a,1x,i0,a,i0)') &
          'Non-sensical combination of bra/ket spin multiplicities:',&
          imultB,'/',imultK
     call error_control
  endif

!----------------------------------------------------------------------
! Generate the bra and ket CSFs
!----------------------------------------------------------------------
  ! Bra CSFs


  ! Ket CSFs
  
!----------------------------------------------------------------------
! Compute the spin-coupling coefficients
!----------------------------------------------------------------------
  if (lsoc) then
     ! Spin-coupling coefficients over spin operators
     errmsg='SOC calculations are not yet supported'
     call error_control
  else
     ! Spin-coupling coefficients over singlet excitation operators
     print*,'HERE'
     stop
  endif
  
  return
  
end subroutine bitsi_intialise
