!######################################################################
! bitwf_finalise: Finalisation of a bitwd calculation.
!                 Simply deallocates all global dynamically allocated
!                 arrays and re-initialises all other global variables
!######################################################################
#ifdef CBINDING
subroutine bitwf_finalise() bind(c,name="bitwf_finalise")
#else
subroutine bitwf_finalise()
#endif

  use constants
  use bitglobal
  
  implicit none

  !
  ! Deallocate arrays
  !
  if (allocated(scrunit))        deallocate(scrunit)
  if (allocated(scrname))        deallocate(scrname)
  if (allocated(csfcoe))         deallocate(csfcoe)
  if (allocated(csfcoeB))        deallocate(csfcoeB)
  if (allocated(csfcoeK))        deallocate(csfcoeK)
  if (allocated(detvec))         deallocate(detvec)
  if (allocated(detvecB))        deallocate(detvecB)
  if (allocated(detvecK))        deallocate(detvecK)
  if (allocated(smo))            deallocate(smo)
  
  !
  ! To be on the safe side, scrub all other global variables
  !
  scratchdir=''
  maxunits=0
  nscratch=0
  nmo=0
  nmoB=0
  nmoK=0
  nel=0
  nel_alpha=0
  nel_beta=0
  nel_spin=0
  ncsfs=0
  ncsfsB=0
  ncsfsK=0
  maxcsf=0
  maxcsfB=0
  maxcsfK=0
  maxdet=0
  maxdetB=0
  maxdetK=0
  imult=0
  imultB=0
  imultK=0
  nelB=0
  nelB_alpha=0
  nelB_beta=0
  nelB_spin=0
  nelK=0
  nelK_alpha=0
  nelK_beta=0
  nelK_spin=0

  return
  
end subroutine bitwf_finalise
