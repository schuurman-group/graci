!######################################################################
! truncate_mrci_wf: truncation of the FOIS space projection of the
!                   MRCI wave functions for a single irrep.
!
!                   The smallest set of configurations s.t. the
!                   truncated wave functions satisfy
!
!                   ||psi||^2 > thrsh
!
!                   is determined and retained.
!######################################################################
#ifdef CBINDING
subroutine truncate_mrci_wf(irrep,nroots,confscr,vecscr,aviiscr,thrsh,&
     nconf_new) bind(c,name="truncate_mrci_wf")
#else
subroutine truncate_mrci_wf(irrep,nroots,confscr,vecscr,aviiscr,thrsh,&
     nconf_new)
#endif

  use constants
  use bitglobal
  use conftype
  use pspace
  use utils
  use iomod
    
  implicit none

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots

  ! Array of MRCI configuration scratch file numbers
  integer(is), intent(in)  :: confscr

  ! Eigenpair scratch file number
  integer(is), intent(in)  :: vecscr

  ! Spin-coupling averaged on-diagonal Hamiltonian matrix
  ! elements scratch file number
  integer(is), intent(in)  :: aviiscr
  
  ! Wave function truncation threshold
  real(dp), intent(in)     :: thrsh

  ! Number of surviving MRCI configurations
  integer(is), intent(out) :: nconf_new
  
  ! MRCI configuration derived types
  type(mrcfg)              :: cfg,cfg_new

  ! Wave function vectors
  real(dp), allocatable    :: vec(:,:),vec_new(:,:)

  ! Spin-coupling averaged on-diagonal Hamiltonian matrix
  ! elements
  real(dp), allocatable    :: averageii(:), averageii_new(:)
  
  ! Surviving configuration flags
  integer(is), allocatable :: i1I(:),i2I(:),i1E(:),i2E(:),i1I1E(:)

  ! Array of CSFs to include in the new space regardless
  ! of their contribution to the WFs
  ! Needed by pspace_conf_indices, but won't be utilised here
  integer(is)              :: nexplicit
  integer(is), allocatable :: iexplicit(:)
  
  ! Everything else
  integer(is)              :: i,iscratch
  integer(is), allocatable :: iroots(:)
  real(dp), allocatable    :: ener(:)
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(vec(cfg%csfdim,nroots))
  vec=0.0d0

  allocate(ener(nroots))
  ener=0.0d0

  allocate(iroots(nroots))
  iroots=0

  allocate(averageii(cfg%confdim))
  averageii=0.0d0
  
  allocate(i1I(cfg%n1I), i2I(cfg%n2I), i1E(cfg%n1E), i2E(cfg%n2E),&
       i1I1E(cfg%n1I1E))
  i1I=0; i2I=0; i1E=0; i2E=0; i1I1E=0

  allocate(iexplicit(cfg%csfdim))
  iexplicit=0
  
!----------------------------------------------------------------------
! Read in the wave function vectors
!----------------------------------------------------------------------
  ! Roots of interest
  do i=1,nroots
     iroots(i)=i
  enddo

  ! Read the wave functions from disk
  call read_some_eigenpairs(vecscr,vec,ener,cfg%csfdim,nroots,iroots)

!----------------------------------------------------------------------
! Read in the spin-coupling averaged on-diagonal Hamiltonian matrix
! elements
!----------------------------------------------------------------------
  call read_averageii(aviiscr,cfg%confdim,averageii)

!----------------------------------------------------------------------
! Get the indices of the configurations that generate CSFs with
! coefficients above threshold
!----------------------------------------------------------------------
  ! Set the no. explicitly added CSFs to zero
  nexplicit=0

  call pspace_conf_indices(cfg,thrsh,vec,cfg%csfdim,cfg%confdim,&
       nroots,nroots,iroots,i1I,i2I,i1E,i2E,i1I1E,cfg%n1I,cfg%n2I,&
       cfg%n1E,cfg%n2E,cfg%n1I1E,nexplicit,iexplicit)

!----------------------------------------------------------------------
! Fill in the surviving configuration information
!----------------------------------------------------------------------
  call pspace_set_new_confs(cfg,cfg_new,i1I,i2I,i1E,i2E,i1I1E,cfg%n1I,&
       cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,confscr,nconf_new)

!----------------------------------------------------------------------
! Set up the CSF information for the truncated set of configurations
! and allocate the new wave function and spin-coupling averaged Hii
! arrays
!----------------------------------------------------------------------
  call cfg_new%finalise
  call cfg_new%initialise(irrep,confscr)
  
  allocate(vec_new(cfg_new%csfdim,nroots))
  vec_new=0.0d0

  allocate(averageii_new(cfg_new%confdim))
  averageii_new=0.0d0
  
!----------------------------------------------------------------------
! Truncate the wave functions
!----------------------------------------------------------------------
  ! Loop over states
  do i=1,nroots

     ! Fil in the surviving CSFs for this state 
     call set_new_vecs(cfg,cfg_new,cfg%csfdim,cfg_new%csfdim,&
          cfg%n1I,cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,i1I,i2I,i1E,i2E,&
          i1I1E,vec(:,i),vec_new(:,i))

  enddo

!----------------------------------------------------------------------
! Truncated spin-coupling averaged Hii array
!----------------------------------------------------------------------
  call set_new_averageii(cfg,cfg_new,cfg%confdim,cfg_new%confdim,&
       cfg%n1I,cfg%n2I,cfg%n1E,cfg%n2E,cfg%n1I1E,i1I,i2I,i1E,i2E,&
       i1I1E,averageii,averageii_new)
  
!----------------------------------------------------------------------
! Symmetric orthonormalisation the truncated wave functions
!----------------------------------------------------------------------
  call symm_ortho(cfg_new%csfdim,nroots,vec_new)
  
!----------------------------------------------------------------------
! Write the truncated wave functions to disk
!----------------------------------------------------------------------
  ! Open the scratch file
  iscratch=scrunit(vecscr)
  open(iscratch,file=scrname(vecscr),form='unformatted',&
       status='unknown')

  ! No. CSFs
  write(iscratch) cfg_new%csfdim

  ! No. roots
  write(iscratch) nroots

  ! Energies
  write(iscratch) ener

  ! Wave functions
  do i=1,nroots
     write(iscratch) vec_new(:,i)
  enddo
  
  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Write the truncated spin-coupling averaged Hii values to disk
!----------------------------------------------------------------------
  ! Open the scratch file
  iscratch=scrunit(aviiscr)
  open(iscratch,file=scrname(aviiscr),form='unformatted',&
       status='unknown')

  ! No. confs
  write(iscratch) cfg_new%confdim

  ! Spin-coupling averaged Hii array
  write(iscratch) averageii_new
  
  ! Close the scratch file
  close(iscratch)
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(vec, vec_new, ener, iroots, i1I, i2I, i1E, i2E, i1I1E, &
       iexplicit, averageii, averageii_new)
    
  return
  
end subroutine truncate_mrci_wf

!######################################################################
! set_new_vecs: given lists of surviving configurations, truncates a
!               single MRCI wave function
!######################################################################
subroutine set_new_vecs(cfg_old,cfg_new,csfdim_old,csfdim_new,&
     n1I_old,n2I_old,n1E_old,n2E_old,n1I1E_old,i1I,i2I,i1E,i2E,i1I1E,&
     vec_old,vec_new)

  use constants
  use bitglobal
  use conftype
  
  implicit none

  ! MRCI configuration derived types
  type(mrcfg), intent(in) :: cfg_old,cfg_new

  ! Dimensions
  integer(is), intent(in) :: csfdim_old,csfdim_new
  integer(is), intent(in) :: n1I_old,n2I_old,n1E_old,n2E_old,n1I1E_old
  
  ! Surviving configuration flags
  integer(is), intent(in) :: i1I(n1I_old),i2I(n2I_old),i1E(n1E_old),&
                             i2E(n2E_old),i1I1E(n1I1E_old)

  ! CSF coefficients
  real(dp), intent(in)    :: vec_old(csfdim_old)
  real(dp), intent(out)   :: vec_new(csfdim_new)

  ! Everything else
  integer(is)             :: n,icsf,ioff
  integer(is)             :: iold,inew

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  ! CSF counters
  iold=0
  inew=0
  
!----------------------------------------------------------------------
! Reference configurations: all retained
!----------------------------------------------------------------------
  ! Loop over ref confs
  do n=1,cfg_old%n0h

     ! Loop over the CSFs generated by this ref conf
     do icsf=cfg_old%csfs0h(n),cfg_old%csfs0h(n+1)-1

        ! Increment the CSF counters
        iold=iold+1
        inew=inew+1
        
        ! New coefficients
        vec_new(inew)=vec_old(iold)
        
     enddo
     
  enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
  if (cfg_old%n1I > 0) then
  
     ! Loop over 1-hole configurations
     do n=1,cfg_old%n1h

        ! Loop over the 1I configurations generated by the 1-hole
        ! configuration
        do ioff=cfg_old%off1I(n),cfg_old%off1I(n+1)-1

           if (i1I(ioff) == 1) then
              ! Surviving conf
              do icsf=cfg_old%csfs1I(ioff),cfg_old%csfs1I(ioff+1)-1
                 iold=iold+1
                 inew=inew+1
                 vec_new(inew)=vec_old(iold)
              enddo
           else
              ! Removed conf
              iold=iold+cfg_old%csfs1I(ioff+1)-cfg_old%csfs1I(ioff)
           endif

        enddo
        
     enddo
        
  endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
  if (cfg_old%n2I > 0) then
  
     ! Loop over 2-hole configurations
     do n=1,cfg_old%n2h

        ! Loop over the 2I configurations generated by the 2-hole
        ! configuration
        do ioff=cfg_old%off2I(n),cfg_old%off2I(n+1)-1

           if (i2I(ioff) == 1) then
              ! Surviving conf
              do icsf=cfg_old%csfs2I(ioff),cfg_old%csfs2I(ioff+1)-1
                 iold=iold+1
                 inew=inew+1
                 vec_new(inew)=vec_old(iold)
              enddo
           else
              ! Removed conf
              iold=iold+cfg_old%csfs2I(ioff+1)-cfg_old%csfs2I(ioff)
           endif

        enddo
        
     enddo
        
  endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
  if (cfg_old%n1E > 0) then
  
     ! Loop over 1-hole configurations
     do n=1,cfg_old%n1h

        ! Loop over the 1E configurations generated by the 1-hole
        ! configuration
        do ioff=cfg_old%off1E(n),cfg_old%off1E(n+1)-1

           if (i1E(ioff) == 1) then
              ! Surviving conf
              do icsf=cfg_old%csfs1E(ioff),cfg_old%csfs1E(ioff+1)-1
                 iold=iold+1
                 inew=inew+1
                 vec_new(inew)=vec_old(iold)
              enddo
           else
              ! Removed conf
              iold=iold+cfg_old%csfs1E(ioff+1)-cfg_old%csfs1E(ioff)
           endif

        enddo
        
     enddo
        
  endif
  
!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
  if (cfg_old%n2E > 0) then
  
     ! Loop over 2-hole configurations
     do n=1,cfg_old%n2h

        ! Loop over the 2E configurations generated by the 2-hole
        ! configuration
        do ioff=cfg_old%off2E(n),cfg_old%off2E(n+1)-1

           if (i2E(ioff) == 1) then
              ! Surviving conf
              do icsf=cfg_old%csfs2E(ioff),cfg_old%csfs2E(ioff+1)-1
                 iold=iold+1
                 inew=inew+1
                 vec_new(inew)=vec_old(iold)
              enddo
           else
              ! Removed conf
              iold=iold+cfg_old%csfs2E(ioff+1)-cfg_old%csfs2E(ioff)
           endif

        enddo
        
     enddo
        
  endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
  if (cfg_old%n1I1E > 0) then
  
     ! Loop over 2-hole configurations
     do n=1,cfg_old%n2h

        ! Loop over the 1I1E configurations generated by the 2-hole
        ! configuration
        do ioff=cfg_old%off1I1E(n),cfg_old%off1I1E(n+1)-1

           if (i1I1E(ioff) == 1) then
              ! Surviving conf
              do icsf=cfg_old%csfs1I1E(ioff),cfg_old%csfs1I1E(ioff+1)-1
                 iold=iold+1
                 inew=inew+1
                 vec_new(inew)=vec_old(iold)
              enddo
           else
              ! Removed conf
              iold=iold+cfg_old%csfs1I1E(ioff+1)-cfg_old%csfs1I1E(ioff)
           endif

        enddo
        
     enddo
        
  endif

  return
  
end subroutine set_new_vecs

!######################################################################
! set_new_averageii: given lists of surviving configurations,
!                    truncates the array of spin-coupling averaged
!                    on-diagonal Hamiltonian matrix elements
!######################################################################
subroutine set_new_averageii(cfg_old,cfg_new,confdim_old,confdim_new,&
     n1I_old,n2I_old,n1E_old,n2E_old,n1I1E_old,i1I,i2I,i1E,i2E,i1I1E,&
     averageii_old,averageii_new)

  use constants
  use bitglobal
  use conftype
  
  implicit none

  ! MRCI configuration derived types
  type(mrcfg), intent(in) :: cfg_old,cfg_new

  ! Dimensions
  integer(is), intent(in) :: confdim_old,confdim_new
  integer(is), intent(in) :: n1I_old,n2I_old,n1E_old,n2E_old,n1I1E_old
  
  ! Surviving configuration flags
  integer(is), intent(in) :: i1I(n1I_old),i2I(n2I_old),i1E(n1E_old),&
                             i2E(n2E_old),i1I1E(n1I1E_old)

  ! CSF coefficients
  real(dp), intent(in)    :: averageii_old(confdim_old)
  real(dp), intent(out)   :: averageii_new(confdim_new)

  ! Everything else
  integer(is)             :: n,ioff
  integer(is)             :: iold,inew

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  ! CSF counters
  iold=0
  inew=0

!----------------------------------------------------------------------
! Reference configurations: all retained
!----------------------------------------------------------------------
  ! Loop over ref confs
  do n=1,cfg_old%n0h
  
     ! Increment the conf counters
     iold=iold+1
     inew=inew+1

     ! New spin-coupling averaged Hii value
     averageii_new(inew)=averageii_old(iold)
     
  enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
  if (cfg_old%n1I > 0) then
  
     ! Loop over 1-hole configurations
     do n=1,cfg_old%n1h

        ! Loop over the 1I configurations generated by the 1-hole
        ! configuration
        do ioff=cfg_old%off1I(n),cfg_old%off1I(n+1)-1

           iold=iold+1
                      
           if (i1I(ioff) == 1) then
              ! Surviving conf
              inew=inew+1
              averageii_new(inew)=averageii_old(iold)
           endif
              
        enddo

     enddo

  endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
  if (cfg_old%n2I > 0) then
  
     ! Loop over 2-hole configurations
     do n=1,cfg_old%n2h

        ! Loop over the 2I configurations generated by the 2-hole
        ! configuration
        do ioff=cfg_old%off2I(n),cfg_old%off2I(n+1)-1

           iold=iold+1
           
           if (i2I(ioff) == 1) then
              ! Surviving conf
              inew=inew+1
              averageii_new(inew)=averageii_old(iold)
           endif

        enddo
        
     enddo
        
  endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
  if (cfg_old%n1E > 0) then
  
     ! Loop over 1-hole configurations
     do n=1,cfg_old%n1h

        ! Loop over the 1E configurations generated by the 1-hole
        ! configuration
        do ioff=cfg_old%off1E(n),cfg_old%off1E(n+1)-1

           iold=iold+1
           
           if (i1E(ioff) == 1) then
              ! Surviving conf
              inew=inew+1
              averageii_new(inew)=averageii_old(iold)
           endif

        enddo
        
     enddo
        
  endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
  if (cfg_old%n2E > 0) then
  
     ! Loop over 2-hole configurations
     do n=1,cfg_old%n2h

        ! Loop over the 2E configurations generated by the 2-hole
        ! configuration
        do ioff=cfg_old%off2E(n),cfg_old%off2E(n+1)-1

           iold=iold+1

           if (i2E(ioff) == 1) then
              ! Surviving conf
              inew=inew+1
              averageii_new(inew)=averageii_old(iold)
           endif

        enddo
        
     enddo
        
  endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
  if (cfg_old%n1I1E > 0) then
  
     ! Loop over 2-hole configurations
     do n=1,cfg_old%n2h

        ! Loop over the 1I1E configurations generated by the 2-hole
        ! configuration
        do ioff=cfg_old%off1I1E(n),cfg_old%off1I1E(n+1)-1

           iold=iold+1
           
           if (i1I1E(ioff) == 1) then
              ! Surviving conf
              inew=inew+1
              averageii_new(inew)=averageii_old(iold)
           endif

        enddo
        
     enddo
        
  endif

  return
  
end subroutine set_new_averageii
