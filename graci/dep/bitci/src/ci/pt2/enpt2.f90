!**********************************************************************
! Routines for the calculation of the ENPT2 energy and wave function
! corrections:
!
! E_w,omega = -|<w omega|H|Psi^0_I>|^2/(H_{w omega,w omega} - E^0_I),
!
! A_w,omega = <w omega|H|Psi^0_I>/(H_{w omega,w omega} - E^0_I),
!
! where {|Psi^0_I>, E^0_I} is the set of reference space eigenpairs.
!**********************************************************************
module epstein_nesbet

  implicit none
  
contains

!######################################################################
! enpt2: Computes a batch of ENPT2 energy and wave function corrections
!######################################################################
  subroutine enpt2(irrep,cfg,hdiag,averageii,csfdim,confdim,vec0scr,&
       Avec,E2,nroots)

    use constants
    use bitglobal
    use conftype
    use pt2_common
    use iomod

    implicit none

    ! Irrep
    integer(is)                     :: irrep
    
    ! MRCI configuration derived type
    type(mrcfg), intent(in)         :: cfg

    ! No. CSFs and configurations
    integer(is), intent(in)         :: csfdim,confdim
    
    ! On-diagonal Hamiltonian matrix elements and their
    ! spin-coupling averaged values
    real(dp), intent(in)            :: hdiag(csfdim),averageii(confdim)

    ! Reference space eigenpair scratch file number
    integer(is), intent(in)         :: vec0scr

    ! Number of roots
    integer(is), intent(in)         :: nroots

    ! ENPT2 wave function and energy corrections
    real(dp), intent(out)           :: Avec(csfdim,nroots)
    real(dp), intent(out)           :: E2(nroots)

    ! Reference space eigenpairs
    integer(is)                     :: refdim
    integer(is), allocatable        :: iroots(:)
    real(dp), allocatable           :: e0(:),vec0(:,:)

    ! Temporary Hij array
    integer(is)                     :: harr2dim
    real(dp), allocatable           :: harr2(:)
    
    ! Everything else
    integer(is)                     :: i

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
! (2)  2-hole configurations -> 2I, 2E and 1I1E configurations
!----------------------------------------------------------------------
    call avec_2h(cfg,Avec,averageii,vec0,csfdim,confdim,refdim,nroots,&
         harr2,harr2dim)
    
!----------------------------------------------------------------------
! Divide by (H_nn - E^0_I)
!----------------------------------------------------------------------
    call apply_denominator(cfg,Avec,E2,hdiag,e0,csfdim,nroots,refdim)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(e0)
    deallocate(iroots)
    deallocate(vec0)
    deallocate(harr2)

    return
    
  end subroutine enpt2
    
!######################################################################
! apply_denominator: Application of the (H_nn - E^0_I) to get
!                    the A-vector element values
!######################################################################
  subroutine apply_denominator(cfg,Avec,E2,hdiag,e0,csfdim,nroots,&
       refdim)

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

    ! Everything else
    integer(is)              :: j,icsf
    real(dp)                 :: Aold,ediff
    
!----------------------------------------------------------------------
! ENPT2 energy and wave function corrections
!----------------------------------------------------------------------
    ! Initialisation
    E2=0.0d0
    
    ! Loop over roots
    do j=1,nroots

       ! Loop over CSFs (excluding the reference space ones)
       do icsf=refdim+1,csfdim

          ! E^(0) - H_ii
          ediff=e0(j)-hdiag(icsf)
          
          ! Energy correction
          E2(j)=E2(j)+Avec(icsf,j)**2/ediff
          
          ! A-vector element
          Avec(icsf,j)=Avec(icsf,j)/ediff

       enddo
       
    enddo

    return
    
  end subroutine apply_denominator

!######################################################################
  
end module epstein_nesbet
