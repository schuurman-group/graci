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
       vec0scr,Avec,E2,nroots,ireg,regfac,dspscr,EQD,mix,heff)

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

    ! Number of roots/model states
    integer(is), intent(in)  :: nroots

    ! ENPT2 wave function and energy corrections
    real(dp), intent(out)    :: Avec(csfdim,nroots)
    real(dp), intent(out)    :: E2(nroots)

    ! Regularizer index: 1 <-> ISA
    !                    2 <-> sigma^p
    integer(is), intent(in)  :: ireg
    
    ! Regularization factor
    real(dp), intent(in)     :: regfac
    
    ! Damped strong perturber scratch file number
    integer, intent(out)     :: dspscr
    
    ! Eigenpairs of the effective Hamiltonian
    real(dp), intent(out)    :: EQD(nroots),mix(nroots,nroots)

    ! Effective Hamiltonian matrix
    real(dp), intent(out)    :: heff(nroots,nroots)
    
    ! Reference space eigenpairs
    integer(is)              :: refdim
    integer(is), allocatable :: iroots(:)
    real(dp), allocatable    :: e0(:),vec0(:,:)

    ! Damped strong perturbers
    integer(is)              :: ndsp
    integer(is), allocatable :: idsp(:)

    ! Temporary Hij array
    integer(is)              :: harr2dim
    real(dp), allocatable    :: harr2(:)
    
    ! I/O
    integer(is)              :: iscratch
    character(len=250)       :: dspfile
    character(len=2)         :: amult,airrep
        
    ! Everything else
    integer(is)              :: i

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Number of reference space CSFs
    refdim=cfg%csfs0h(cfg%n0h+1)-1

    ! Reference space eigenvalues
    allocate(e0(nroots))
    e0=0.0d0

    ! Reference space eigenvectors
    allocate(iroots(nroots))
    allocate(vec0(refdim,nroots))
    iroots=0
    vec0=0.0d0

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
    do i=1,nroots
       iroots(i)=i
    enddo
    call read_some_eigenpairs(vec0scr,vec0,e0,refdim,nroots,iroots)

    ! Subtract off E_SCF from the energies to get the true eigenvalues
    e0=e0-escf

!----------------------------------------------------------------------
! (1) 1-hole configurations -> 1I and 1E configurations
!----------------------------------------------------------------------
    call avec_1h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,nroots,&
         harr2,harr2dim)

!----------------------------------------------------------------------
! (2) 2-hole configurations -> 2I, 2E and 1I1E configurations
!----------------------------------------------------------------------
    call avec_2h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,nroots,&
         harr2,harr2dim)

!----------------------------------------------------------------------
! Construct and diagonalise the GVVPT2 effective Hamiltonian
!----------------------------------------------------------------------
    call diag_heff(cfg,Avec,hdiag,e0,EQD,mix,csfdim,nroots,refdim,&
         ireg,regfac,heff)
    
!----------------------------------------------------------------------
! Compute the first-order perturbed model states
!----------------------------------------------------------------------
    call model_states1(cfg,Avec,E2,hdiag,e0,csfdim,nroots,refdim,&
         ireg,regfac,idsp)
   
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
         'Number of damped strong perturbers:',ndsp
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(e0)
    deallocate(iroots)
    deallocate(vec0)
    deallocate(harr2)
    deallocate(idsp)

    return
    
  end subroutine gvvpt2_heff
    
!######################################################################
! model_states1: Application of the (H_nn - E^0_I) to get the 1st-order
!                perturbed model states 
!######################################################################
  subroutine model_states1(cfg,Avec,E2,hdiag,e0,csfdim,nroots,refdim,&
       ireg,regfac,idsp)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Dimensions
    integer(is), intent(in)  :: csfdim,nroots,refdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)     :: hdiag(csfdim)

    ! ENPT2 wave function and energy corrections
    real(dp), intent(inout)  :: Avec(csfdim,nroots)
    real(dp), intent(out)    :: E2(nroots)
    
    ! Reference space eigenvalues
    real(dp), intent(in)     :: e0(nroots)

    ! Regularizer index: 1 <-> ISA
    !                    2 <-> sigma^p
    integer(is), intent(in)  :: ireg
    
    ! Regularization factor
    real(dp), intent(in)     :: regfac

    ! Damped strong perturbers
    integer(is), intent(out) :: idsp(csfdim)    
    real(dp), parameter      :: cthrsh=0.055d0
    
    ! Everything else
    integer(is)              :: j,icsf
    real(dp)                 :: dj,Aold,ediff
    real(dp)                 :: deltai,Vi,fi
    
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
          deltai=e0(j)-hdiag(icsf)

          ! <I|H|i>
          Vi=Avec(icsf,j)

          ! Regularized denominators
          if (ireg == 1) then
             ! ISA
             fi=deltai/(deltai**2+regfac)
          else if (ireg == 2) then
             ! sigma^p
             fi=(1.0d0-exp(-abs(deltai/regfac)**2))/deltai
          endif

          ! Unregularized A-vector element
          Aold=Vi/deltai
          
          ! A-vector element
          Avec(icsf,j)=Vi*fi

          ! Energy correction
          E2(j)=E2(j)+Vi**2*fi
          
          ! Make sure that all strong perturbers are captured
          if (abs(Aold) >= cthrsh &
               .and. abs(Avec(icsf,j)) < cthrsh) idsp(icsf)=1
          
       enddo
       
    enddo

    return
    
  end subroutine model_states1

!######################################################################
! diag_heff: constructs and diagonalises the GVVPT2 effective
!            Hamiltonian using the ENPT2 Hamiltonian partitioning  
!######################################################################
  subroutine diag_heff(cfg,Avec,hdiag,e0,EQD,mix,csfdim,nroots,&
       refdim,ireg,regfac,heff)

    use constants
    use bitglobal
    use conftype
    use utils
    use iomod
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfg

    ! Dimensions
    integer(is), intent(in) :: csfdim,nroots,refdim

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(csfdim)

    ! A-vector
    real(dp), intent(inout) :: Avec(csfdim,nroots)
    
    ! Reference space eigenvalues
    real(dp), intent(in)    :: e0(nroots)

    ! QDPT2 energies and mixing coefficients
    real(dp), intent(out)   :: EQD(nroots)
    real(dp), intent(out)   :: mix(nroots,nroots)

    ! Regularizer index: 1 <-> ISA
    !                    2 <-> sigma^p
    integer(is), intent(in) :: ireg
    
    ! Regularization factor
    real(dp), intent(in)    :: regfac
    
    ! Effective Hamiltonian
    real(dp), intent(out)   :: heff(nroots,nroots)

    ! Everything else
    integer(is)             :: i,j,icsf
    real(dp)                :: fac,di,dj,deltai,deltaj,Vi,Vj
    real(dp)                :: fi,fj,epsilon
    
!----------------------------------------------------------------------
! Second-order contribution to the GVVPT2 effective Hamiltonian
!----------------------------------------------------------------------
    ! Initialisation
    heff=0.0d0
    do i=1,nroots
       heff(i,i)=e0(i)
    enddo

    ! Loop over pairs of roots
    do i=1,nroots
       do j=i,nroots

          ! Loop over FOIS CSFs
          do icsf=refdim+1,csfdim

             ! <I|H|i> and <i|H|J>
             Vi=Avec(icsf,i)
             Vj=Avec(icsf,j)
             
             ! Regularized denominators
             deltai=e0(i)-hdiag(icsf)
             deltaj=e0(j)-hdiag(icsf)
             if (ireg == 1) then
                ! ISA
                fi=deltai/(deltai**2+regfac)
                fj=deltaj/(deltaj**2+regfac)
             else if (ireg == 2) then
                ! sigma^p
                fi=(1.0d0-exp(-abs(deltai/regfac)**2))/deltai
                fj=(1.0d0-exp(-abs(deltaj/regfac)**2))/deltaj
             endif
             
             ! <I|H|i> <i|H|J> / (E_I^0 - H_ii) * f_I + (I<->J)
             heff(i,j)=heff(i,j)+0.5d0*Vi*Vj*(fi+fj)
             
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
    
    return
    
  end subroutine diag_heff
  
!######################################################################
  
end module gvvpt2_hamiltonian
