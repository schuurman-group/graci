!**********************************************************************
! Routines for the calculation of wave function overlaps using
! pre-computed beta factors, etc.
!**********************************************************************
module wfoverlap

  implicit none

contains

!######################################################################
! get_overlaps: semi-direct calculation of wave function overlaps using
!               pre-computed beta factors
!######################################################################
  subroutine get_overlaps(npairs,ipairs,Sij)

    use constants
    use global
    use detfuncs
    use factors
    
    implicit none

    ! Indices of the pairs of states for which overlaps are requested
    integer(is), intent(in)  :: npairs
    integer(is), intent(in)  :: ipairs(npairs,2)
    
    ! Wave function overlaps
    real(dp), intent(out)    :: Sij(npairs)
    
    ! Transposed eigenvector arrays
    real(dp), allocatable    :: vecTB(:,:),vecTK(:,:)

    ! Occupied MOs
    integer(is)              :: noccB,noccK
    integer(is), allocatable :: occB(:),occK(:)

    ! Work arrays
    real(dp), allocatable    :: fwork(:,:)
    integer(is), allocatable :: iwork(:)

    ! Everything else
    integer(is)              :: n,isB,isK,iaB,iaK,ibB,ibK,idB,idK
    real(dp)                 :: afac,bfac

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(vecTB(nrootsB,ndetB))
    vecTB=0.0d0

    allocate(vecTK(nrootsK,ndetK))
    vecTK=0.0d0

    ! Note here that the no. bra and ket alpha electrons are equal
    ! in a wave function overlap calculation
    allocate(occB(nel_alphaB), occK(nel_alphaB))
    occB=0; occK=0

    allocate(fwork(nel_alphaB,nel_alphaB))
    fwork=0.0d0

    allocate(iwork(nel_alphaB))
    iwork=0
    
!----------------------------------------------------------------------
! Transposes of the eigenvector arrays
!----------------------------------------------------------------------
    vecTB=transpose(vecB)
    vecTK=transpose(vecK)
    
!----------------------------------------------------------------------
! Compute the wave function overlaps
!----------------------------------------------------------------------
    ! Initialisation
    Sij=0.0d0
    
    ! Loop over ket alpha strings
    do iaK=1,nalphaK
    
       ! Get the ket occupied MO indices
       call mo_occ_string(n_intK,alphaK(:,iaK),nel_alphaK,noccK,occK)
       
       ! Loop over bra alpha strings
       do iaB=1,nalphaB
    
          ! Get the ket occupied MO indices
          call mo_occ_string(n_intB,alphaB(:,iaB),nel_alphaB,&
               noccB,occB)

          ! Compute the alpha factor for this pair of strings
          call get_one_factor(nel_alphaB,alphaB(:,iaB),alphaK(:,iaK),&
               occB,occK,fwork,iwork,afac)

          ! Cycle if the alpha factor is below threshold
          if (abs(afac) < fthrsh) cycle
    
          ! Loop over determinants in the ket block
          do idK=offsetK(iaK),offsetK(iaK+1)-1
    
             ! Ket beta string index
             ibK=det2betaK(idK)
             
             ! Loop over determinants in the bra block
             do idB=offsetB(iaB),offsetB(iaB+1)-1
    
                ! Bra beta string index
                ibB=det2betaB(idB)
                
                ! Beta factor
                bfac=betafac(ibB,ibK)
    
                ! Cycle if the beta factor is below threshold
                if (abs(bfac) < fthrsh) cycle
                
                ! Loop over bra-ket state pairs
                do n=1,npairs
    
                   ! Bra and ket state indices
                   isB=ipairs(n,1)
                   isK=ipairs(n,2)
    
                   ! Contributions to the overlap
                   Sij(n)=Sij(n) &
                        +afac*bfac*vecTB(isB,idB)*vecTK(isK,idK)
                   
                enddo
    
             enddo
                
          enddo
             
       enddo
    
    enddo

    return
    
  end subroutine get_overlaps

!######################################################################
! get_overlaps_schur: semi-direct calculation of wave function overlaps
!                     using pre-computed beta factors and alpha factors
!                     computed using Schur's determinant identity
!######################################################################
  subroutine get_overlaps_schur(npairs,ipairs,Sij)

    use constants
    use global
    use detfuncs
    use factors
    
    implicit none

    ! Indices of the pairs of states for which overlaps are requested
    integer(is), intent(in)  :: npairs
    integer(is), intent(in)  :: ipairs(npairs,2)
    
    ! Wave function overlaps
    real(dp), intent(out)    :: Sij(npairs)

    ! Transposed eigenvector arrays
    real(dp), allocatable    :: vecTB(:,:),vecTK(:,:)

    ! Occupied MOs
    integer(is)              :: noccB,noccK
    integer(is), allocatable :: occB(:),occK(:)

    ! Work arrays
    real(dp), allocatable    :: S(:,:),Svv(:,:),Svf(:,:),Sfv(:,:),&
                                invSffSfv(:,:)
    real(dp), allocatable    :: fwork(:,:)
    integer(is), allocatable :: iwork(:)

    ! Everything else
    integer(is)              :: n,isB,isK,iaB,iaK,ibB,ibK,idB,idK
    integer(is)              :: fdim
    real(dp)                 :: afac,bfac
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(vecTB(nrootsB,ndetB))
    vecTB=0.0d0

    allocate(vecTK(nrootsK,ndetK))
    vecTK=0.0d0

    ! Note here that the no. bra and ket alpha electrons are equal
    ! in a wave function overlap calculation
    allocate(occB(nel_alphaB), occK(nel_alphaB))
    occB=0; occK=0

    allocate(S(nel_alphaB,nel_alphaB),Svv(nvar_alphaB,nvar_alphaB),&
         Svf(nvar_alphaB,nfixed), &
         Sfv(nfixed,nvar_alphaB), invSffSfv(nfixed,nvar_alphaB))
    Svv=0.0d0; Svf=0.0d0; Sfv=0.0d0; invSffSfv=0.0d0
    
    allocate(fwork(nvar_alphaB,nvar_alphaB))
    fwork=0.0d0

    allocate(iwork(nvar_alphaB))
    iwork=0
    
!----------------------------------------------------------------------
! Transposes of the eigenvector arrays
!----------------------------------------------------------------------
    vecTB=transpose(vecB)
    vecTK=transpose(vecK)

!----------------------------------------------------------------------
! Compute the wave function overlaps
!----------------------------------------------------------------------
    ! Initialisation
    Sij=0.0d0

    ! Loop over ket alpha strings
    do iaK=1,nalphaK
    
       ! Get the ket occupied MO indices
       call mo_occ_string(n_intK,alphaK(:,iaK),nel_alphaK,noccK,occK)
       
       ! Loop over bra alpha strings
       do iaB=1,nalphaB
    
          ! Get the ket occupied MO indices
          call mo_occ_string(n_intB,alphaB(:,iaB),nel_alphaB,&
               noccB,occB)

          ! Compute the alpha factor for this pair of strings
          call get_one_factor_schur(nel_alphaB,nfixed,nvar_alphaB,&
               alphaB(:,iaB),alphaK(:,iaK),&
               occB,occK,&
               S,Svv,Svf,Sfv,invSffSfv,fwork,iwork,&
               afac)

          ! Cycle if the alpha factor is below threshold
          if (abs(afac) < fthrsh) cycle

          ! Loop over determinants in the ket block
          do idK=offsetK(iaK),offsetK(iaK+1)-1
    
             ! Ket beta string index
             ibK=det2betaK(idK)
             
             ! Loop over determinants in the bra block
             do idB=offsetB(iaB),offsetB(iaB+1)-1
    
                ! Bra beta string index
                ibB=det2betaB(idB)
                
                ! Beta factor
                bfac=betafac(ibB,ibK)
    
                ! Cycle if the beta factor is below threshold
                if (abs(bfac) < fthrsh) cycle
                
                ! Loop over bra-ket state pairs
                do n=1,npairs
    
                   ! Bra and ket state indices
                   isB=ipairs(n,1)
                   isK=ipairs(n,2)
    
                   ! Contributions to the overlap
                   Sij(n)=Sij(n) &
                        +afac*bfac*vecTB(isB,idB)*vecTK(isK,idK)
                   
                enddo
    
             enddo

          enddo
          
       enddo

    enddo
    
    return
    
  end subroutine get_overlaps_schur
    
!######################################################################
  
end module wfoverlap
