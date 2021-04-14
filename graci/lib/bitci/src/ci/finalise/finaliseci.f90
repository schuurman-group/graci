!######################################################################
! bitci_finalise: Finalisation of a bitci calculation.
!                 Simply deallocates all global dynamically allocated
!                 arrays and re-initialises all other global variables
!######################################################################
#ifdef CBINDING
subroutine bitci_finalise() bind(c,name="bitci_finalise")
#else
subroutine bitci_finalise()
#endif

  use constants
  use bitglobal
  use hparam
  
  implicit none

  !
  ! Deallocate arrays
  !
  if (allocated(scrunit))     deallocate(scrunit)
  if (allocated(scrname))     deallocate(scrname)
  if (allocated(patternmap1)) deallocate(patternmap1)
  if (allocated(patternmap2)) deallocate(patternmap2)
  if (allocated(N1s))         deallocate(N1s)
  if (allocated(spincp1))     deallocate(spincp1)
  if (allocated(spincp2))     deallocate(spincp2)
  if (allocated(det0))        deallocate(det0)
  if (allocated(conf0))       deallocate(conf0)
  if (allocated(mosym))       deallocate(mosym)
  if (allocated(moen))        deallocate(moen)
  if (allocated(csfcoe))      deallocate(csfcoe)
  if (allocated(detvec))      deallocate(detvec)
  if (allocated(fii))         deallocate(fii)
  if (allocated(fock))        deallocate(fock)
  if (allocated(Vc))          deallocate(Vc)
  if (allocated(Vx))          deallocate(Vx)
  if (allocated(hpar))        deallocate(hpar)

  !
  ! To be on the safe side, scrub all other global variables
  !
  scratchdir=''
  maxunits=0
  nscratch=0
  nspincp=0
  maxpattern=0
  npattern1=0
  npattern2=0
  nmo=0
  nel=0
  nel_alpha=0
  nel_beta=0
  nel_spin(2)=0
  enuc=0.0d0
  ipg=0
  nirrep=0
  pgdim=0
  pgroup=''
  irreplbl=''
  pglbls=''
  imult=0
  ncsfs=0
  maxcsf=0
  maxdet=0
  Escf=0.0d0
  ldftmrci=.false.
  nhpar=0
  desel=0.0d0
  ihamiltonian=0
  
  return
  
end subroutine bitci_finalise
