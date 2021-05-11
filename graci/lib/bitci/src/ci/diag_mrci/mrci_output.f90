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
! Section header
!----------------------------------------------------------------------
    write(6,'(/,52a)') ('-',i=1,52)
    if (ldftmrci) then
       write(6,'(x,a)') 'DFT/MRCI eigenstates'
    else
       write(6,'(x,a)') 'MRCI eigenstates'
    endif
    write(6,'(52a)') ('-',i=1,52)
    
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
    write(6,'(/,2x,a,x,F6.4)') &
         'Minimum norm in the reference space:',minrnorm
    
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
! sort_energies: returns irrep/state number pairs in order of
!                increasing energy
!######################################################################
  subroutine sort_energies(ener,nroot_max,nroot_tot,nroots,sindx)

    use constants
    use bitglobal
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: nroot_max,nroot_tot
    integer(is), intent(in)  :: nroots(0:nirrep-1)
    
    ! Energies
    real(dp), intent(in)     :: ener(nroot_max,0:nirrep-1)

    ! Mapping array
    integer(is), intent(out) :: sindx(nroot_tot,2)

    ! Everything else
    integer(is)              :: irrep,n,k
    real(dp), allocatable    :: etmp(:,:)
    real(dp)                 :: minval
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(etmp(nroot_max,0:nirrep-1))
    etmp=0.0d0
    
!----------------------------------------------------------------------
! Fill in the mapping array
! This is a pretty dumb approach, but speed doesn't matter here
!----------------------------------------------------------------------
    ! Initialisation
    etmp=ener

    ! Loop over the total number of states
    do n=1,nroot_tot

       minval=1e+90_dp
       
       ! Loop over irreps
       do irrep=0,nirrep-1
          
          ! Loop over roots for this irrep
          do k=1,nroots(irrep)

             if (etmp(k,irrep) < minval) then
                minval=etmp(k,irrep)
                sindx(n,2)=irrep
                sindx(n,1)=k
             endif
          
          enddo
          
       enddo

       etmp(sindx(n,1),sindx(n,2))=1e+90_dp
       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(etmp)
    
    return
    
  end subroutine sort_energies
  
!######################################################################
! print_report_1state: prints the energy and dominant CSFs for a
!                      single MRCI state
!######################################################################
  subroutine print_report_1state(n,k,irrep,ener,emin,cfg,vecscr,rnorm)

    use constants
    use bitglobal
    use conftype
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
    write(6,'(/,2x,50a)') ('-',i=1,50)
    write(6,'(3x,a,x,i0,a,x,i0,a3,2(2x,F12.6),x,a)') &
         'State' ,n,':',k,irreplbl(irrep,ipg),ener,(ener-emin)*eh2ev,&
         'eV'
    write(6,'(2x,50a)') ('-',i=1,50)
    write(6,'(3x,a,x,F6.4)') '||Psi_R|| = ',rnorm
    write(6,'(2x,50a)') ('-',i=1,50)
    write(6,'(4x,a)') 'Coeff        omega     Delta w'
    write(6,'(2x,50a)') ('-',i=1,50)
    
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
       write(6,'(3x,F10.7,6x,i0,6x,a)') dcoe(i),domega(i),string
       
    enddo

    ! Footer
    write(6,'(2x,50a)') ('-',i=1,50)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dconf, domega, dcoe, abscoe, indx)
    
    return
    
  end subroutine print_report_1state
  
!######################################################################
! dominant_csfs: for a given eigenvector, returns the configuration
!                and spin-coupling information corresponding to the
!                dominant CSFs
!######################################################################
  subroutine dominant_csfs(modus,csfthrsh,cfg,vecscr,k,ndom,refnorm,&
       dconf,domega,dcoe)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of dominant
    !                                CSFs
    !                    modus=1 <-> package the configuration and
    !                                spin-coupling information
    integer(is), intent(in)    :: modus
    
    ! CSF coefficient threshold
    real(dp), intent(in)       :: csfthrsh

    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! Eigenpair scratch file number
    integer(is), intent(in)    :: vecscr

    ! Root number
    integer(is), intent(in)    :: k

    ! No. dominant CSFs
    integer(is), intent(inout) :: ndom
    
    ! Norm of the wavefunction projected onto the ref space
    real(dp), intent(out)      :: refnorm

    ! Dominant configurations and spin-couplings
    integer(ib), allocatable   :: dconf(:,:,:)
    integer(is), allocatable   :: domega(:)
    real(dp), allocatable      :: dcoe(:)
    
    ! Eigenvector
    real(dp), allocatable      :: vec(:)
    
    ! Everything else
    integer(is)                :: csfdim,refdim,n,ioff,csf,counter
    integer(is)                :: omega
    integer(is)                :: n_int_I
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    csfdim=cfg%csfdim
    allocate(vec(csfdim))
        
!----------------------------------------------------------------------
! Read in the eigenvector
!----------------------------------------------------------------------
    call read_single_vector(vecscr,vec,csfdim,k)

!----------------------------------------------------------------------
! Norm of the wavefunction projected onto the reference space
!----------------------------------------------------------------------
    refdim=cfg%csfs0h(cfg%n0h+1)-1
    refnorm=sqrt(dot_product(vec(1:refdim),vec(1:refdim)))

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    counter=0
    n_int_I=cfg%n_int_I
    
    if (modus == 1) then
       allocate(dconf(n_int,2,ndom))
       allocate(domega(ndom))
       allocate(dcoe(ndom))
       dconf=0_ib
       domega=0
       dcoe=0.0d0
    endif
       
!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
    ! Loop over ref space configurations
    do n=1,cfg%n0h

       ! Loop over the CSFs generated by this configuration
       omega=0
       do csf=cfg%csfs0h(n),cfg%csfs0h(n+1)-1

          ! Spin-coupling index
          omega=omega+1
          
          ! Cycle if this is not a dominant CSF
          if (vec(csf)**2 < csfthrsh) cycle
          
          ! Increment the dominant CSF counter
          counter=counter+1

          ! Cycle if we are just determinig the no. dominant CSFs
          if (modus == 0) cycle

          ! Save the CSF information
          dconf(1:n_int_I,:,counter)=cfg%conf0h(:,:,n)
          domega(counter)=omega
          dcoe(counter)=vec(csf)
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1I configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs1I(ioff),cfg%csfs1I(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf1I(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo

          enddo
          
       enddo

    endif
       
!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (cfg%n2I > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2I configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs2I(ioff),cfg%csfs2I(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf2I(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (cfg%n1E > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs1E(ioff),cfg%csfs1E(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf1E(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (cfg%n2E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs2E(ioff),cfg%csfs2E(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf2E(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (cfg%n1I1E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 1I1E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs1I1E(ioff),cfg%csfs1I1E(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf1I1E(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! Set the number of dominant CSFs
!----------------------------------------------------------------------
    ndom=counter
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(vec)
    
    return
    
  end subroutine dominant_csfs

!######################################################################
! exci_to_string: given lists of hole and particle indices, returns
!                 a character string of the form i, j,... -> a, b,...
!######################################################################
  subroutine exci_to_string(string,hlist,plist,nexci,m2c)

    use constants
    use bitglobal
    
    implicit none

    ! Creation and annihilation operator indices
    integer(is), intent(in) :: nexci
    integer(is), intent(in) :: hlist(nexci),plist(nexci)

    ! MO mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Ouput character string
    character(len=*)        :: string

    ! Everything else
    integer(is)             :: i,k

    !
    ! Initialisation
    !
    string=''
    k=1

    !
    ! Base configuration case
    !
    if (nexci == 0) then
       string='base conf'
       return
    endif
    
    !
    ! Annihilation operators
    !
    do i=1,nexci
       write(string(k:),'(x,i0)') m2c(hlist(i))
       k=len_trim(string)+1
    enddo

    !
    ! Creation operators
    !
    write(string(k:),'(x,a)') '->'
    k=len_trim(string)+1
    do i=1,nexci
       write(string(k:),'(x,i0)') m2c(plist(i))
       k=len_trim(string)+1
    enddo

    return
    
  end subroutine exci_to_string
  
!######################################################################
! deadwood: deadwood analysis
!######################################################################
  subroutine deadwood(vecscr,nroots,cfg)

    use constants
    use bitglobal
    use conftype
    use iomod
    use utils
    
    implicit none

    ! MRCI eigenpair scratch file numbers
    integer(is), intent(in)  :: vecscr(0:nirrep-1)
    
    ! Number of roots per irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg(0:nirrep-1)
    
    ! Eigenvectors
    real(dp), allocatable    :: vec(:)

    ! Everything else
    integer(is)              :: irrep,k,i,csfdim,unit
    integer(is), allocatable :: indx(:)
    real(dp)                 :: normsq
    character(len=255)       :: filename
    
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! No. CSFs
       csfdim=cfg(irrep)%csfdim
       
       ! Allocate arrays
       allocate(vec(csfdim),indx(csfdim))
       
       ! Loop over roots
       do k=1,nroots(irrep)

          ! Open the output file
          write(filename,'(a,i0,a)') &
               'normsq_',k,trim(irreplbl(irrep,ipg))//'.dat'
          call freeunit(unit)
          open(unit,file=filename,form='formatted',status='unknown')
          
          ! Read in the eigenvector
          call read_single_vector(vecscr(irrep),vec,csfdim,k)
       
          ! Sort the coefficients
          vec=abs(vec)
          call dsortindxa1('D',csfdim,vec,indx)

          ! Squared norm of the truncated wavefunction
          normsq=0.0d0
          do i=1,csfdim
             normsq=normsq+vec(indx(i))**2
             write(unit,'(i0,x,ES10.4)') i,normsq
             !if (normsq > 0.9999d0) exit
          enddo

          ! Close the output file
          close(unit)
          
       enddo
          
       ! Deallocate arrays
       deallocate(vec,indx)
       
    enddo
    
    return
    
  end subroutine deadwood
    
!######################################################################
  
end module mrci_output
