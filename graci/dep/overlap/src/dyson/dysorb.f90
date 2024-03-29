module dyson_builder

  implicit none

contains

!######################################################################
! get_dysorbs: semi-direct calculation of Dyson orbitals using
!              pre-computed tau factors
!######################################################################
  subroutine get_dysorbs(n_basis,npairs,ipairs,dysorb)

    use constants
    use global
    use detfuncs
    use factors
    
    implicit none

    ! Indices of the pairs of states for which overlaps are requested
    integer(is), intent(in) :: npairs
    integer(is), intent(in) :: ipairs(npairs,2)
    
    ! Dyson orbitals
    integer(is), intent(in) :: n_basis
    real(dp), intent(inout) :: dysorb(n_basis,npairs)

    ! Transposed eigenvector arrays
    real(dp), allocatable    :: vecTB(:,:),vecTK(:,:)

    ! Occupied MOs
    integer(is)              :: noccB,noccK
    integer(is), allocatable :: occB(:),occK(:)

    ! Work arrays
    real(dp), allocatable    :: fwork(:,:)
    integer(is), allocatable :: iwork(:)
    
    ! Everything else
    integer(is)              :: isHK,isK,isB,itB,itK,idB,idK
    integer(is)              :: iloc,imo,iphase,n,istaB,istaK
    real(dp)                 :: sfac,tfac,prefac
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(vecTB(nrootsB,ndetB))
    vecTB=0.0d0

    allocate(vecTK(nrootsK,ndetK))
    vecTK=0.0d0

    allocate(occB(nel_sigmaB), occK(nel_sigmaK))
    occB=0; occK=0

    allocate(fwork(nel_sigmaB,nel_sigmaK))
    fwork=0.0d0

    allocate(iwork(nel_sigmaB))
    iwork=0

!----------------------------------------------------------------------
! Transposes of the eigenvector arrays
!----------------------------------------------------------------------
    vecTB=transpose(vecB)
    vecTK=transpose(vecK)

!----------------------------------------------------------------------
! Compute the Dyson orbitals
!----------------------------------------------------------------------
    ! Initialisation
    dysorb=0.0d0

    ! Loop over ket sigma-hole strings
    do isHK=1,nsigmaHK

       ! Get the ket occupied MO indices
       call mo_occ_string(n_intK,sigmaHK(:,isHK),nel_sigmaK,noccK,occK)

       ! Loop over bra sigma strings
       do isB=1,nsigmaB

          ! Get the ket occupied MO indices
          call mo_occ_string(n_intB,sigmaB(:,isB),nel_sigmaB,&
               noccB,occB)

          ! Compute the sigma factor for this pair of strings
          call get_one_factor(nel_sigmaB,sigmaB(:,isB),&
               sigmaHK(:,isHK),occB,occK,fwork,iwork,sfac)

          ! Cycle if the sigma factor is below threshold
          if (abs(sfac) < fthrsh) cycle

          ! Loop over the ket sigma strings corresponding to
          ! the ket sigma-hole string
          do iloc=offHinfo(isHK),offHinfo(isHK+1)-1

             ! Annihilation operator MO index
             imo=Hinfo(1,iloc)

             ! Ket sigma string index
             isK=Hinfo(2,iloc)

             ! Phase factor
             iphase=Hinfo(3,iloc)

             ! Loop over the determinants generated by the ket
             ! sigma string
             do idK=offsetK(isK),offsetK(isK+1)-1

                ! Ket tau string index
                itK=det2tauK(idK)

                ! Loop over determinants in the bra block
                do idB=offsetB(isB),offsetB(isB+1)-1

                   ! Bra beta string index
                   itB=det2tauB(idB)
                
                   ! tau factor
                   tfac=taufac(itB,itK)

                   ! Cycle if the tau factor is below threshold
                   if (abs(tfac) < fthrsh) cycle

                   ! sigma-factor * tau-factor * phase-factor
                   prefac=iphase*sfac*tfac
                   
                   ! Loop over bra-ket state pairs
                   do n=1,npairs

                      ! Bra and ket state indices
                      istaB=ipairs(n,1)
                      istaK=ipairs(n,2)
                      
                      ! Contribution to the Dyson orbital
                      dysorb(imo,n)=dysorb(imo,n) &
                           +prefac*vecTB(istaB,idB)*vecTK(istaK,idK)
                      
                   enddo
                   
                enddo
                   
             enddo
             
          enddo
          
       enddo
       
    enddo

    return
    
  end subroutine get_dysorbs

!######################################################################
  
end module dyson_builder
