!######################################################################
! bitsi_initialise: Interface to the bitsi library. Sets all the
!                   globally accessible variables that are needed to
!                   perform a State Interaction calculation.
!######################################################################
#ifdef CBINDING
subroutine bitsi_intialise() &
     bind(c,name='bitsi_initialise')
#else
subroutine bitsi_intialise()
#endif

  use constants
    
  implicit none

  print*,"HERE"
  stop
  
  return
  
end subroutine bitsi_intialise
