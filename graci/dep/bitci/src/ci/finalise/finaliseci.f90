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
  if (allocated(scrunit))   deallocate(scrunit)
  if (allocated(scrname))   deallocate(scrname)
  if (allocated(patmap))    deallocate(patmap)
  if (allocated(N1s))       deallocate(N1s)
  if (allocated(spincp))    deallocate(spincp)
  if (allocated(offspincp)) deallocate(offspincp)
  if (allocated(det0))      deallocate(det0)
  if (allocated(conf0))     deallocate(conf0)
  if (allocated(iopen0))    deallocate(iopen0)
  if (allocated(iocc0))     deallocate(iocc0)
  if (allocated(mosym))     deallocate(mosym)
  if (allocated(moen))      deallocate(moen)
  if (allocated(icvs))      deallocate(icvs)
  if (allocated(csfcoe))    deallocate(csfcoe)
  if (allocated(detvec))    deallocate(detvec)
  if (allocated(fii))       deallocate(fii)
  if (allocated(fock))      deallocate(fock)
  if (allocated(Vc))        deallocate(Vc)
  if (allocated(Vx))        deallocate(Vx)
  if (allocated(hpar))      deallocate(hpar)

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
  nel_spin(2)=0
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
  lcvs=.false.
  nhpar=0
  desel=0.0d0
  ihamiltonian=0

  return
  
end subroutine bitci_finalise

!
!
!
#ifdef CBINDING
subroutine bitci_int_finalize() bind(c,name="bitci_int_finalize")
#else
subroutine bitci_int_finalize()
#endif

  use bitglobal

  ! deallocate integral arrays
  if(allocated(bitci_ints)) then
    call bitci_ints%finalize()
    deallocate(bitci_ints)
  endif

  return

end subroutine bitci_int_finalize
