!**********************************************************************
! Routines for the output of information about the eigenstates of a
! pruned MRCI calculation
!**********************************************************************
module mrci_pruned_output

  implicit none
  
  public print_pmrci_states

contains

!######################################################################
! print_pmrci_states: master routine for the printing of the energies
!                     and dominant CSFs of the eigenstates from a
!                     pruned MRCI calculation
!######################################################################
#ifdef CBINDING
  subroutine print_pmrci_states(confscr,vecscr,vec0scr,eqscr,nroots) &
       bind(c,name="print_pmrci_states")
#else
  subroutine print_pmrci_states(confscr,vecscr,vec0scr,eqscr,nroots)
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

    ! Reference space eigenpair scratch file numbers
    integer(is), intent(in)  :: vec0scr(0:nirrep-1)

    ! Q-space energy correction scratch file numbers
    integer(is), intent(in)  :: eqscr(0:nirrep-1)
    
    ! Number of roots per irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)

    ! MRCI configuration derived type
    type(mrcfg), allocatable :: cfg(:)

    ! Energies
    real(dp), allocatable    :: ener(:,:)

    ! Indices of the irrep/roots in order of ascending energy
    integer(is), allocatable :: sindx(:,:)

    ! Total no. ref space eigenvectors saved to disk
    integer(is), allocatable :: nvec(:)

    ! Reference space eigenpairs
    integer(is), allocatable :: csfdim0(:)
    integer(is)              :: maxcsf0
    real(dp), allocatable    :: vec0(:,:,:),E0(:,:)
        
    ! Q-space energy corrections
    integer(is)              :: maxvec
    real(dp), allocatable    :: E2Q(:,:)

    ! Ref-to-MRCI mapping
    integer(is), allocatable :: r2m(:,:)
    real(dp), allocatable    :: maxoverlap(:,:)
    
    ! Everything else
    integer(is)              :: irrep,nroot_tot,nroot_max
    integer(is)              :: n,k,i,iscratch
    real(dp)                 :: emin,minrnorm,rnorm
    logical                  :: lowoverlap
    character(len=3)         :: an
        
!----------------------------------------------------------------------
! Section header
!----------------------------------------------------------------------
    write(6,'(/,52a)') ('-',i=1,52)
    if (ldftmrci) then
       write(6,'(x,a)') 'Pruned DFT/MRCI eigenstates'
    else
       write(6,'(x,a)') 'Pruned MRCI eigenstates'
    endif
    write(6,'(52a)') ('-',i=1,52)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Total no. ref space eigenvectors saved to disk
    allocate(nvec(0:nirrep-1))
    nvec=0

    ! Number of ref space CSFs per irrep
    allocate(csfdim0(0:nirrep-1))
    csfdim0=0.0d0
    
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

    ! Ref-to-MRCI state mapping
    allocate(r2m(nroot_max,0:nirrep-1))
    allocate(maxoverlap(nroot_max,0:nirrep-1))
    r2m=0.0d0
    maxoverlap=0.0d0
    
!----------------------------------------------------------------------
! Set up the MRCI configuration derived types
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       call cfg(irrep)%initialise(irrep,confscr(irrep))
    enddo    
    
!----------------------------------------------------------------------
! Get the total number of ref space eigenvectors saved to disk per
! irrep
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Open the ref space eigenvector file
       iscratch=scrunit(vec0scr(irrep))
       open(iscratch,file=scrname(vec0scr(irrep)),form='unformatted',&
            status='old')

       ! Read past the no. CSFs
       read(iscratch) n

       ! Number eigenvectors
       read(iscratch) nvec(irrep)

       ! Close the eigenvector file
       close(iscratch)
       
    enddo

    ! Maximum no. vectors across all irreps
    maxvec=maxval(nvec)
    
!----------------------------------------------------------------------
! Read in the Q-space energy corrections
!----------------------------------------------------------------------
    ! Allocate the Q-space energy correction array
    allocate(E2Q(maxvec,0:nirrep-1))
    E2Q=0.0d0

    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Open the E2Q file
       iscratch=scrunit(eqscr(irrep))
       open(iscratch,file=scrname(eqscr(irrep)),form='unformatted',&
            status='old')

       ! Number vectors
       read(iscratch) n

       ! Consistency check
       if (n /= nvec(irrep)) then
          errmsg='Error in print_pmrci_states: inconsistent numbers '&
               //'of ref space vectors'
          call error_control
       endif

       ! Q-space energy corrections
       read(iscratch) E2Q(1:nvec(irrep),irrep)
       
       ! Close the eigenvector file
       close(iscratch)
       
    enddo

!----------------------------------------------------------------------
! No. ref space CSFs per irrep
!----------------------------------------------------------------------
    csfdim0=0

    ! Loop over irreps
    do irrep=0,nirrep-1

       do n=1,cfg(irrep)%n0h
          csfdim0(irrep)=csfdim0(irrep)&
               +cfg(irrep)%csfs0h(n+1)-cfg(irrep)%csfs0h(n)
       enddo

    enddo

    maxcsf0=maxval(csfdim0)

!----------------------------------------------------------------------
! Read in all the ref space eigenvectors saved to disk
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(E0(maxvec,0:nirrep-1))
    allocate(vec0(maxcsf0,maxvec,0:nirrep-1))
    E0=0.0d0
    vec0=0.0d0
    
    ! Read in the eigenpairs
    do irrep=0,nirrep-1
       call read_all_eigenpairs(&
            vec0scr(irrep),&
            vec0(1:csfdim0(irrep),1:nvec(irrep),irrep),&
            E0(1:nvec(irrep),irrep),&
            csfdim0(irrep),&
            nvec(irrep))
    enddo

!----------------------------------------------------------------------
! Determine which Q-space energy corrections to use based on the
! overlaps of the reference and MRCI eigenvectors
!----------------------------------------------------------------------
    call get_r2m(cfg,nroots,nroot_max,csfdim0,maxcsf0,nvec,maxvec,&
         vec0,r2m,maxoverlap,vecscr)

!----------------------------------------------------------------------
! Read in the MRCI energies
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       call read_energies(vecscr(irrep),nroots(irrep),&
            ener(1:nroots(irrep),irrep))
    enddo

!----------------------------------------------------------------------
! Add on the Q-space energy corrections
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       do k=1,nroots(irrep)
          n=r2m(k,irrep)
          ener(k,irrep)=ener(k,irrep)+E2Q(n,irrep)
       enddo
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
       call print_report_1state_pruned(n,k,irrep,ener(k,irrep),emin,&
            cfg(irrep),vecscr(irrep),maxoverlap(k,irrep),rnorm)

       ! Minimum reference space norm
       minrnorm=min(rnorm,minrnorm)

    enddo

!----------------------------------------------------------------------
! Output a warning if there wasn't a good overlap of an MRCI
! eigenfunction with one of the ref space eigenfunctions
! In this case, the Q-space corrected energies are of questionable
! quality   
!----------------------------------------------------------------------
    lowoverlap=.false.

    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Loop over states
       do n=1,nroots(irrep)

          ! Low overlap?
          if (maxoverlap(n,irrep) < 0.8d0) then
             lowoverlap=.true.
             write(an,'(i0)') n
             write(6,'(/,2x,a)') &
                  'WARNING: small ref-MRCI overlap found for '&
                  //'State '//trim(an)//irreplbl(irrep,ipg)
          endif
          
       enddo
       
    enddo

    ! Helpful suggestion
    if (lowoverlap) &
         write(6,'(/,2x,a)') 'Small overlaps found: Consider '&
         //'increasing the value of nextra used in ref_diag_mrci'
    
!----------------------------------------------------------------------
! Minimum reference space norm
!----------------------------------------------------------------------
    write(6,'(/,2x,a,x,F6.4)') &
         'Minimum norm in the reference space:',minrnorm
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(nvec)
    deallocate(csfdim0)
    deallocate(ener)
    deallocate(sindx)
    deallocate(cfg)
    deallocate(r2m)
    deallocate(maxoverlap)
    deallocate(E2Q)
    deallocate(E0)
    deallocate(vec0)
    
    return

  end subroutine print_pmrci_states

!######################################################################
! get_r2m: determines the reference-to-MRCI state mapping
!######################################################################
  subroutine get_r2m(cfg,nroots,nroot_max,csfdim0,maxcsf0,nvec,&
       maxvec,vec0,r2m,maxoverlap,vecscr)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! MRCI configuration derived types
    type(mrcfg), intent(in)  :: cfg(0:nirrep-1)

    ! Number of roots per irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)
    integer(is), intent(in)  :: nroot_max

    ! Ref space CSF dimensions
    integer(is), intent(in)  :: csfdim0(0:nirrep-1)
    integer(is), intent(in)  :: maxcsf0

    ! Total no. ref space eigenvectors
    integer(is), intent(in)  :: nvec(0:nirrep-1)
    integer(is), intent(in)  :: maxvec

    ! Ref space eigenvectors
    real(dp), intent(in)     :: vec0(maxcsf0,maxvec,0:nirrep-1)

    ! Ref-to-MRCI state mapping
    integer(is), intent(out) :: r2m(nroot_max,0:nirrep-1)
    real(dp), intent(out)    :: maxoverlap(nroot_max,0:nirrep-1)
    
    ! MRCI eigenpair scratch file numbers
    integer(is), intent(in)  :: vecscr(0:nirrep-1)
    
    ! MRCI eigenpairs
    real(dp), allocatable    :: vec(:,:),ener(:)
    
    ! Everything else
    integer(is)              :: irrep,csfdim,i,j
    real(dp)                 :: overlap

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    r2m=0
    maxoverlap=0.0d0

!----------------------------------------------------------------------
! Determine the ref-to-MRCI mappings
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! No. MRCI CSFs
       csfdim=cfg(irrep)%csfdim

       ! Allocate arrays
       allocate(vec(csfdim,nroots(irrep)), ener(nroots(irrep)))

       ! Read in the MRCI eigenpairs
       call read_all_eigenpairs(vecscr(irrep),vec,ener,csfdim,&
            nroots(irrep))

       print*,''
       
       ! Loop over MRCI roots
       do i=1,nroots(irrep)

          ! Determine the corresponding ref space eigenvector
          do j=1,nvec(irrep)

             ! Overlap between the i'th MRCI eigenvector and
             ! the j'th ref space eigenvector
             overlap=abs(dot_product(vec(1:csfdim0(irrep),i),&
                  vec0(1:csfdim0(irrep),j,irrep)))
             
             ! Update the maximim overlap and mapping arrays
             if (overlap > maxoverlap(i,irrep)) then
                maxoverlap(i,irrep)=overlap
                r2m(i,irrep)=j
             endif

           enddo
             
       enddo
       
       ! Deallocate arrays
       deallocate(vec, ener)
       
    enddo
    
    return
    
  end subroutine get_r2m

!######################################################################
! print_report_1state_pruned: prints the energy and dominant CSFs for
!                             a single pruned MRCI state
!######################################################################
  subroutine print_report_1state_pruned(n,k,irrep,ener,emin,cfg,&
       vecscr,maxoverlap,rnorm)

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

    ! Max. overlap of the eigenvector with a ref space eigenvector
    real(dp), intent(in)     :: maxoverlap
    
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
    write(6,'(2(3x,a,x,F6.4))') '||Psi_R|| = ',rnorm,&
         'max <R|I> = ',maxoverlap
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
    
  end subroutine print_report_1state_pruned
    
!######################################################################
  
end module mrci_pruned_output
