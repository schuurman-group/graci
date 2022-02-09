!######################################################################
! bitsi_finalise: Finalisation of a bitsi calculation.
!                 Simply deallocates all global dynamically allocated
!                 arrays and re-initialises all other global variables
!######################################################################
#ifdef CBINDING
subroutine bitsi_finalise() bind(c,name="bitsi_finalise")
#else
subroutine bitsi_finalise()
#endif

  use constants
  use bitglobal
  
  implicit none

  !
  ! Deallocate arrays
  !
  if (allocated(scrunit))       deallocate(scrunit)
  if (allocated(scrname))       deallocate(scrname)
  if (allocated(patmap))        deallocate(patmap)
  if (allocated(patmap1))       deallocate(patmap1)
  if (allocated(patmap2a_plus)) deallocate(patmap2a_plus)
  if (allocated(patmap2b_plus)) deallocate(patmap2b_plus)
  if (allocated(N1s))           deallocate(N1s)
  if (allocated(spincp))        deallocate(spincp)
  if (allocated(spincp_plus))   deallocate(spincp_plus)
  if (allocated(spincp_minus))  deallocate(spincp_minus)
  if (allocated(offplus))       deallocate(offplus)
  if (allocated(offminus))      deallocate(offminus)
  if (allocated(det0))          deallocate(det0)
  if (allocated(conf0))         deallocate(conf0)
  if (allocated(csfcoe))        deallocate(csfcoe)
  if (allocated(detvec))        deallocate(detvec)

  !
  ! To be on the safe side, scrub all other global variables
  !
  scratchdir=''
  maxunits=0
  nscratch=0
  nspincp=0
  nmo=0
  nel=0
  nel_alpha=0
  nel_beta=0
  nel_spin=0
  imult=0
  ncsfs=0
  maxcsf=0
  maxdet=0
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
  
end subroutine bitsi_finalise
