!**********************************************************************
! Routines for the calculation of wave function overlaps using
! pre-computed beta factors, etc.
!**********************************************************************
module wfoverlap

  implicit none

contains

!######################################################################
! get_overlaps: Top level routine for the calculation of the wave
!               function overlaps
!######################################################################
  subroutine get_overlaps(npairs,ipairs,Sij)

    use constants
    use global
    use detfuncs
    use factors
    use timing
        
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

    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end

    ! Everything else
    integer(is)              :: n,isB,isK,iaB,iaK,ibB,ibK,idB,idK
    real(dp)                 :: afac,bfac

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)

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

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'get_overlaps')
    
    return
    
  end subroutine get_overlaps
  
!######################################################################
  
end module wfoverlap
