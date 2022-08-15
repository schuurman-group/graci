!**********************************************************************
! Routines for the output of information about the MRCI eigenstates
!**********************************************************************
module mrci_output

  implicit none
  
  public print_mrci_states

contains

!######################################################################
! print_mrci_states: master routine for the printing of the energies
!                    and dominant CSFs of the MRCI eigenstates
!######################################################################
#ifdef CBINDING
  subroutine print_mrci_states(confscr,vecscr,nroots) &
       bind(c,name="print_mrci_states")
#else
  subroutine print_mrci_states(confscr,vecscr,nroots)
#endif

    use constants
    use bitglobal
    use conftype
    use output_common
    use iomod
    
    implicit none
    
    ! MRCI configuration scratch file numbers
    integer(is), intent(in)  :: confscr(0:nirrep-1)

    ! MRCI eigenpair scratch file numbers
    integer(is), intent(in)  :: vecscr(0:nirrep-1)
    
    ! Number of roots per irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)

    ! MRCI configuration derived type
    type(mrcfg), allocatable :: cfg(:)

    ! Energies
    real(dp), allocatable    :: ener(:,:)

    ! Indices of the irrep/roots in order of ascending energy
    integer(is), allocatable :: sindx(:,:)
    
    ! Everything else
    integer(is)              :: irrep,nroot_tot,nroot_max
    integer(is)              :: n,k,i
    real(dp)                 :: emin,minrnorm,rnorm
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Numbers of roots
    nroot_tot=sum(nroots)
    nroot_max=maxval(nroots)

    ! Energies
    allocate(ener(nroot_max,0:nirrep-1))
    ener=0.0d0
    
    ! Root mapping array
    allocate(sindx(nroot_tot,2))
    sindx=0
    
    ! Configuration derived types
    allocate(cfg(0:nirrep-1))
    
!----------------------------------------------------------------------
! Set up the MRCI configuration derived types
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       call cfg(irrep)%initialise(irrep,confscr(irrep))
    enddo

!----------------------------------------------------------------------
! Deadwood analysis
!----------------------------------------------------------------------
    !call deadwood(vecscr,nroots,cfg)
    
!----------------------------------------------------------------------
! Read in the state energies
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       call read_energies(vecscr(irrep),nroots(irrep),&
            ener(1:nroots(irrep),irrep))
    enddo

!----------------------------------------------------------------------
! Sort the energies
!----------------------------------------------------------------------
    call sort_energies(ener,nroot_max,nroot_tot,nroots,sindx)
    
!----------------------------------------------------------------------
! Output the base conf
!----------------------------------------------------------------------
    if (verbose) call print_base_conf
    
!----------------------------------------------------------------------
! Output the energies and dominant CSFs for each state
!----------------------------------------------------------------------
    minrnorm=1.0d0
    
    ! Minimum energy
    emin=ener(sindx(1,1),sindx(1,2))
    
    ! Loop over roots
    do n=1,nroot_tot

       ! State number and irrep
       k=sindx(n,1)
       irrep=sindx(n,2)

       ! Print the report
       call print_report_1state(n,k,irrep,ener(k,irrep),emin,&
            cfg(irrep),vecscr(irrep),rnorm)

       ! Minimum reference space norm
       minrnorm=min(rnorm,minrnorm)

    enddo

!----------------------------------------------------------------------
! Minimum reference space norm
!----------------------------------------------------------------------
    if (verbose) then
       write(6,'(/,2x,a,x,F6.4)') &
            'Minimum norm in the reference space:',minrnorm
    endif

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ener)
    deallocate(sindx)
    deallocate(cfg)

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
    flush(6)
    
    return
     
  end subroutine print_mrci_states
  
!######################################################################
! print_report_1state: prints the energy and dominant CSFs for a
!                      single MRCI state
!######################################################################
  subroutine print_report_1state(n,k,irrep,ener,emin,cfg,vecscr,rnorm)

    use constants
    use bitglobal
    use conftype
    use output_common
    use mrciutils
    use utils
    
    implicit none

    ! State number and irrep
    integer(is), intent(in)  :: n,k,irrep

    ! State energy
    real(dp), intent(in)     :: ener,emin

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Eigenpair scratch file number
    integer(is), intent(in)  :: vecscr

    ! Projection of the eigenvector onto the reference space
    integer(is)              :: nrefcsf
    real(dp), intent(out)    :: rnorm

    ! Dominant configurations and spin-couplings
    integer(ib), allocatable :: dconf(:,:,:)
    integer(is), allocatable :: domega(:)
    integer(is), allocatable :: indx(:)
    real(dp), allocatable    :: dcoe(:),abscoe(:)
    integer(is)              :: hlist(nexmax),plist(nexmax)
        
    ! Everything else
    integer(is)              :: i,i1,modus,ndom,nexci
    real(dp), parameter      :: csfthrsh=0.005d0
    character(len=50)        :: string

!----------------------------------------------------------------------
! Number of reference space CSFs
!----------------------------------------------------------------------
    nrefcsf=cfg%csfs0h(cfg%n0h+1)-1

!----------------------------------------------------------------------
! Determine the configuration and spin-coupling information for the
! dominant CSFs
!----------------------------------------------------------------------
    ! First pass: determine the number of dominant CSFs
    modus=0
    call dominant_csfs(modus,csfthrsh,cfg,vecscr,k,ndom,rnorm,&
         dconf,domega,dcoe)

    ! Second pass: determine the spatial configurations and
    ! spin-couplings of the dominant configurations
    modus=1
    call dominant_csfs(modus,csfthrsh,cfg,vecscr,k,ndom,rnorm,&
         dconf,domega,dcoe)

!----------------------------------------------------------------------
! Sort the dominant coefficients by absolute value
!----------------------------------------------------------------------
    allocate(abscoe(ndom), indx(ndom))

    abscoe=abs(dcoe)

    call dsortindxa1('D',ndom,abscoe,indx)

!----------------------------------------------------------------------
! Header
!----------------------------------------------------------------------
    if (verbose) then
       write(6,'(/,2x,50a)') ('-',i=1,50)
       write(6,'(3x,a,x,i0,a,x,i0,a3,2(2x,F12.6),x,a)') &
            'State' ,n,':',k,irreplbl(irrep,ipg),ener,&
            (ener-emin)*eh2ev,'eV'
       write(6,'(2x,50a)') ('-',i=1,50)
       write(6,'(3x,a,x,F6.4)') '||Psi_R|| = ',rnorm
       write(6,'(2x,50a)') ('-',i=1,50)
       write(6,'(4x,a)') 'Coeff        omega     Delta w'
       write(6,'(2x,50a)') ('-',i=1,50)
    endif

!----------------------------------------------------------------------
! Output the dominant CSF information
!----------------------------------------------------------------------
    ! Loop over dominant CSFs
    do i1=1,ndom
       
       i=indx(i1)

       ! Get the difference configuration information
       nexci=exc_degree_conf(conf0,dconf(:,:,i),n_int)
       hlist=0; plist=0
       call get_exci_indices(conf0,dconf(:,:,i),n_int,&
            hlist(1:nexci),plist(1:nexci),nexci)
       call exci_to_string(string,hlist(1:nexci),plist(1:nexci),&
            nexci,cfg%m2c)
       
       ! Output the CSF coefficient, spin-coupling and spatial
       ! configuration
       if (verbose) &
            write(6,'(3x,F10.7,6x,i0,6x,a)') dcoe(i),domega(i),string
       
    enddo

    ! Footer
    if (verbose) write(6,'(2x,50a)') ('-',i=1,50)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dconf, domega, dcoe, abscoe, indx)
    
    return
    
  end subroutine print_report_1state
  
!######################################################################
  
end module mrci_output
