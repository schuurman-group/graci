!**********************************************************************
! Routines for the construction and diagonalisation of the GVVPT2
! effective Hamiltonian
!**********************************************************************
module gvvpt2_hamiltonian
  
  implicit none

contains

!######################################################################
! gvvpt2_heff: Computes the eigenpairs of the GVVPT2 effective
!              Hamiltonian as well as the first-order perturbed
!              model functions
!######################################################################
  subroutine gvvpt2_heff(irrep,cfg,hdiag,averageii,csfdim,confdim,&
       vec0scr,mcoeff,Avec,E2,nroots,nref,shift,dspscr,EQD,mix)

    use constants
    use bitglobal
    use conftype
    use pt2_common
    use iomod

    implicit none

    ! Irrep
    integer(is)              :: irrep
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! No. CSFs and configurations
    integer(is), intent(in)  :: csfdim,confdim
    
    ! On-diagonal Hamiltonian matrix elements and their
    ! spin-coupling averaged values
    real(dp), intent(in)     :: hdiag(csfdim),averageii(confdim)

    ! Reference space eigenpair scratch file number
    integer(is), intent(in)  :: vec0scr

    ! Representation of the model states in the reference space
    ! eigenstate basis
    real(dp), intent(in)     :: mcoeff(nref,nroots)
    
    ! Number of roots/model states
    integer(is), intent(in)  :: nroots

    ! Number of reference space states
    integer(is), intent(in)  :: nref
    
    ! ENPT2 wave function and energy corrections
    real(dp), intent(out)    :: Avec(csfdim,nroots)
    real(dp), intent(out)    :: E2(nroots)

    ! ISA shift
    real(dp), intent(in)     :: shift
    
    ! Damped strong perturber scratch file number
    integer, intent(out)     :: dspscr
    
    ! Eigenpairs of the effective Hamiltonian
    real(dp), intent(out)    :: EQD(nroots),mix(nroots,nroots)

    ! Reference space eigenpairs
    integer(is)              :: refdim
    integer(is), allocatable :: iroots(:)
    real(dp), allocatable    :: e0(:),vec0(:,:)

    ! Model states
    real(dp), allocatable    :: mvec(:,:)
    
    ! Damped strong perturbers
    integer(is)              :: ndsp
    integer(is), allocatable :: idsp(:)

    ! Temporary Hij array
    integer(is)              :: harr2dim
    real(dp), allocatable    :: harr2(:)
    
    ! I/O
    integer(is)              :: iscratch
    character(len=60)        :: dspfile
    character(len=2)         :: amult,airrep
        
    ! Everything else
    integer(is)              :: i

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Number of reference space CSFs
    refdim=cfg%csfs0h(cfg%n0h+1)-1

    ! Reference space eigenvalues
    allocate(e0(nref))
    e0=0.0d0

    ! Reference space eigenvectors
    allocate(iroots(nref))
    allocate(vec0(refdim,nref))
    iroots=0
    vec0=0.0d0

    ! Model states
    allocate(mvec(refdim,nroots))
    mvec=0.0d0
    
    ! Hij working array
    harr2dim=maxval(ncsfs(0:nomax))**2
    allocate(harr2(harr2dim))

    ! Indices of the damped strong perturbers
    allocate(idsp(csfdim))
    idsp=0
    
!----------------------------------------------------------------------
! Reference space eigenpairs
!----------------------------------------------------------------------
! Note that we have to use read_some_eigenpairs here as different
! numbers of extra ref space eigenvectors may be used for
! different tasks. We will, however, assume for now that the 1st
! nroots ref space eigenvectors are needed.
!----------------------------------------------------------------------
    ! Read in the eigenpairs
    do i=1,nref
       iroots(i)=i
    enddo
    call read_some_eigenpairs(vec0scr,vec0,e0,refdim,nref,iroots)

    ! Subtract off E_SCF from the energies to get the true eigenvalues
    e0=e0-escf

!----------------------------------------------------------------------
! Construct the model states as linear combination of the reference
! space eigenstates
!----------------------------------------------------------------------
    mvec=matmul(vec0,mcoeff)

!----------------------------------------------------------------------
! (1) 1-hole configurations -> 1I and 1E configurations
!----------------------------------------------------------------------
    call avec_1h(cfg,Avec,averageii,mvec,csfdim,confdim,refdim,nroots,&
         harr2,harr2dim)
    
!----------------------------------------------------------------------
! (2) 2-hole configurations -> 2I, 2E and 1I1E configurations
!----------------------------------------------------------------------
    call avec_2h(cfg,Avec,averageii,mvec,csfdim,confdim,refdim,nroots,&
         harr2,harr2dim)

!----------------------------------------------------------------------
! Construct and diagonalise the GVVPT2 effective Hamiltonian
!----------------------------------------------------------------------
    call diag_heff(cfg,Avec,hdiag,e0,mcoeff,EQD,mix,csfdim,&
         nroots,nref,refdim,shift)
    
!----------------------------------------------------------------------
! Compute the first-order perturbed model states
!----------------------------------------------------------------------
    call model_states1(cfg,Avec,E2,hdiag,e0,mcoeff,csfdim,nroots,nref,&
         refdim,shift,idsp)
   
!----------------------------------------------------------------------
! Write the damped strong perturber array to disk
!----------------------------------------------------------------------
    ! Register the scratch file
    write(amult,'(i0)') imult
    write(airrep,'(i0)') irrep
    call scratch_name('dsp'//'.mult'//trim(amult)//&
         '.sym'//trim(airrep),dspfile)
    call register_scratch_file(dspscr,dspfile)

    ! Open scratch file
    iscratch=scrunit(dspscr)
    open(iscratch,file=scrname(dspscr),form='unformatted',&
         status='unknown')

    ! Number of damped strong perturbers
    ndsp=sum(idsp)
    write(iscratch) ndsp

    ! Damped strong perturber flags
    write(iscratch) idsp
    
    ! Close scratch file
    close(iscratch)

!----------------------------------------------------------------------
! Output the number of damped strong perturbers
!----------------------------------------------------------------------
    if (verbose) write(6,'(/,x,a,x,i0)') &
         'Number of damped strong perurbers:',ndsp
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(e0)
    deallocate(iroots)
    deallocate(vec0)
    deallocate(mvec)
    deallocate(harr2)
    deallocate(idsp)

    return
    
  end subroutine gvvpt2_heff
    
!######################################################################
! model_states1: Application of the (H_nn - E^0_I) to get the 1st-order
!                perturbed model states 
!######################################################################
  subroutine model_states1(cfg,Avec,E2,hdiag,e0,mcoeff,csfdim,nroots,&
       nref,refdim,shift,idsp)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Dimensions
    integer(is), intent(in)  :: csfdim,nroots,nref,refdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hdiag(csfdim)

    ! ENPT2 wave function and energy corrections
    real(dp), intent(inout)  :: Avec(csfdim,nroots)
    real(dp), intent(out)    :: E2(nroots)
    
    ! Reference space eigenvalues
    real(dp), intent(in)     :: e0(nref)

    ! Model space expansion coefficients
    real(dp), intent(in)     :: mcoeff(nref,nroots)
    
    ! ISA shift
    real(dp), intent(in)     :: shift

    ! Damped strong perturbers
    integer(is), intent(out) :: idsp(csfdim)    
    real(dp), parameter      :: cthrsh=0.055d0

    ! Zeroth-erder effective Hamiltonian matrix
    real(dp), allocatable   :: href(:,:),heff0(:,:)
    
    ! Everything else
    integer(is)              :: j,icsf
    real(dp)                 :: dj,Aold,ediff

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(href(nref,nref))
    href=0.0d0

    allocate(heff0(nroots,nroots))
    heff0=0.0d0

!----------------------------------------------------------------------
! Zeroth-order effective Hamiltonian
!----------------------------------------------------------------------
    ! H_eff in the ref space eigenstate basis
    do j=1,nref
       href(j,j)=e0(j)
    enddo

    ! Transform to the model space basis
    heff0=matmul(matmul(transpose(mcoeff),href),mcoeff)
    
!----------------------------------------------------------------------
! Compute the 1st-order perturbed model states (projected onto the
! Q-space)
!----------------------------------------------------------------------
    ! Initialisation
    E2=0.0d0
    idsp=0
    
    ! Loop over roots
    do j=1,nroots

       ! Loop over CSFs (excluding the reference space ones)
       do icsf=refdim+1,csfdim

          ! E^(0) - H_ii
          ediff=heff0(j,j)-hdiag(icsf)
          
          ! ISA factor
          dj=shift/ediff

          ! Unshifted A-vector element
          Aold=Avec(icsf,j)/ediff
          
          ! Energy correction
          E2(j)=E2(j)+Avec(icsf,j)**2/(ediff+dj)
          
          ! A-vector element
          Avec(icsf,j)=Avec(icsf,j)/(ediff+dj)

          ! Make sure that all strong perturbers are captured
          if (abs(Aold) >= cthrsh &
               .and. abs(Avec(icsf,j)) < cthrsh) idsp(icsf)=1
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(href)
    deallocate(heff0)
        
    return
    
  end subroutine model_states1

!######################################################################
! diag_heff: constructs and diagonalises the GVVPT2 effective
!            Hamiltonian using the ENPT2 Hamiltonian partitioning  
!######################################################################
  subroutine diag_heff(cfg,Avec,hdiag,e0,mcoeff,EQD,mix,csfdim,&
       nroots,nref,refdim,shift)

    use constants
    use bitglobal
    use conftype
    use utils
    use iomod
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Dimensions
    integer(is), intent(in) :: csfdim,nroots,nref,refdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(csfdim)

    ! A-vector
    real(dp), intent(inout) :: Avec(csfdim,nroots)
    
    ! Reference space eigenvalues
    real(dp), intent(in)    :: e0(nref)

    ! Model space expansion coefficients
    real(dp), intent(in)    :: mcoeff(nref,nroots)

    ! QDPT2 energies and mixing coefficients
    real(dp), intent(out)   :: EQD(nroots)
    real(dp), intent(out)   :: mix(nroots,nroots)

    ! ISA shift
    real(dp), intent(in)    :: shift
    
    ! Effective Hamiltonian matrix
    real(dp), allocatable   :: href(:,:),heff0(:,:),heff(:,:)

    ! Everything else
    integer(is)             :: i,j,icsf
    real(dp)                :: fac,di,dj
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(href(nref,nref))
    href=0.0d0

    allocate(heff0(nroots,nroots))
    heff0=0.0d0

    allocate(heff(nroots,nroots))
    heff=0.0d0
    
!----------------------------------------------------------------------
! Zeroth-order contribution to the GVVPT2 effective Hamiltonian
!----------------------------------------------------------------------
    ! H_eff in the ref space eigenstate basis
    do i=1,nref
       href(i,i)=e0(i)
    enddo

    ! Transform to the model space basis
    heff0=matmul(matmul(transpose(mcoeff),href),mcoeff)
    
!----------------------------------------------------------------------
! Second-order contribution to the GVVPT2 effective Hamiltonian
!----------------------------------------------------------------------
    ! Initialisation
    heff=heff0
    
    ! Loop over pairs of roots
    do i=1,nroots
       do j=i,nroots

          ! Loop over FOIS CSFs
          do icsf=refdim+1,csfdim

             ! ISA factors
             di=shift/(heff0(i,i)-hdiag(icsf))
             dj=shift/(heff0(j,j)-hdiag(icsf))

             ! ISA-shifted denominators
             fac=0.5d0/(heff0(i,i)-hdiag(icsf)+di)
             fac=fac+0.5d0/(heff0(j,j)-hdiag(icsf)+dj)

             ! <I|H|i> <i|H|J> / (E_I^0 - H_ii + Delta_Ii) + (I<->J)
             heff(i,j)=heff(i,j)+Avec(icsf,i)*Avec(icsf,j)*fac
             
          enddo
          
          heff(j,i)=heff(i,j)

       enddo
    enddo

!----------------------------------------------------------------------
! Diagonalise the effective Hamiltonian
!----------------------------------------------------------------------
    ! Eigenpairs of H_eff
    call diag_matrix_real(heff,EQD,mix,nroots)

    ! Add E_SCF to the eigenvalues
    EQD=EQD+escf
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(href)
    deallocate(heff0)
    
    return
    
  end subroutine diag_heff
  
!######################################################################
  
end module gvvpt2_hamiltonian
