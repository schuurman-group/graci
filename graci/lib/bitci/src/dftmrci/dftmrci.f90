!***********************************************************************
! Routines for the application of the DFT/MRCI corrections to the
! Hamiltonian matrix
!***********************************************************************
module dftmrci

  implicit none

contains
  
!######################################################################
! hii_dftmrci: applies DFT/MRCI corrections to a batch of on-diagonal
!              Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci(harr,nsp,Dw,ndiff,nopen,m2c)
    
    use constants
    use bitglobal
    use hparam

    implicit none

    ! Array of on-diagonal Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: harr(nsp)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    select case(ihamiltonian)
       
    case(2:3)
       ! Grimme's parameterisation
       call hii_dftmrci_grimme(harr,nsp,Dw,ndiff,nopen,m2c)
       
    case default
       print*,'Your Hamiltonian choice has not been implemented yet'
       stop
       
    end select
    
    return
    
  end subroutine hii_dftmrci

!######################################################################
! hij_dftmrci: applies DFT/MRCI corrections to a batch of off-diagonal
!              Hamiltonian matrix element
!######################################################################
  subroutine hij_dftmrci_batch(hij,bdim,kdim,bav,kav)

    use constants
    use bitglobal
    use hparam
    
    implicit none

    ! Hamiltonian matrix elements
    integer(is), intent(in) :: bdim,kdim
    real(dp), intent(inout) :: hij(:,:)
    real(dp), intent(in)    :: bav,kav
    real(dp)                :: damp
    
!----------------------------------------------------------------------
! Compute the damping factor
!----------------------------------------------------------------------
    select case(ihamiltonian)
       
    case(2:3)
       ! Grimme's parameterisation
       damp=damping_grimme(bav,kav)
       
    case default
       print*,'Your Hamiltonian choice has not been implemented yet'
       stop
       
    end select

!----------------------------------------------------------------------
! Apply the damping factor
!----------------------------------------------------------------------
    hij(1:bdim,1:kdim)=hij(1:bdim,1:kdim)*damp
    
    return
    
  end subroutine hij_dftmrci_batch
    
!######################################################################
! hii_dftmrci_grimme: applies Grimme's DFT/MRCI correction to a batch
!                     of on-diagonal Hamiltonian matrix elements
!######################################################################
  subroutine hii_dftmrci_grimme(harr,nsp,Dw,ndiff,nopen,m2c)
    
    use constants
    use bitglobal
    use hparam

    implicit none

    ! Array of on-diagonal Hamiltonian matrix elements
    integer(is), intent(in) :: nsp
    real(dp), intent(inout) :: harr(nsp)

    ! Difference configuration information
    integer(is), intent(in) :: ndiff
    integer(is), intent(in) :: Dw(nmo,2)

    ! Number of open shells
    integer(is), intent(in) :: nopen
    
    ! MO index mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Everything else
    integer(is)             :: i,j,i1,j1,Dwi,Dwj,ipos
    integer(is)             :: nexci
    real(dp)                :: Viijj,Vijji
    real(dp)                :: pJ,pN0
    real(dp)                :: contrib
    
!----------------------------------------------------------------------
! Return if we are at the base configuration
!----------------------------------------------------------------------
    if (ndiff == 0) return

!----------------------------------------------------------------------
! Parameter values
!----------------------------------------------------------------------
    pJ=hpar(3)
    pN0=hpar(4)+nopen*hpar(5)
    
!----------------------------------------------------------------------
! Sum_i F_ii^KS - F_ii^HF Delta w_i
!----------------------------------------------------------------------
    ! Loop over non-zero Delta w_i values
    do i=1,ndiff
    
       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
    
       ! Sum the contribution
       harr=harr+(moen(i1)-Fii(i1))*Dwi
       
    enddo

!----------------------------------------------------------------------
!  1/nexc Sum_i Sum_j (pJ V_iijj - p[N0] V_ijji) |Delta w_i| Delta w_j
!         w_i<0 w_j<0
!----------------------------------------------------------------------
    !
    ! Find the start of the postive Delta w_i values
    !
    do i=1,ndiff
       if (Dw(i,2) > 0) then
          ipos=i
          exit
       endif
    enddo

    !
    ! Excitation degree
    !
    nexci=sum(Dw(ipos:ndiff,2))

    !
    ! Sum_i Sum_j (pJ V_iijj - p[N0] V_ijji) |Delta w_i| Delta w_j
    ! w_i<0 w_j<0
    !
    contrib=0.0d0
    ! Loop over negative Delta w_i values
    do i=1,ipos-1

       ! MO index
       i1=m2c(Dw(i,1))
    
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Loop over positive Delta w_j values
       do j=ipos,ndiff

          ! MO index
          j1=m2c(Dw(j,1))
          
          ! Delta w_j value
          Dwj=Dw(j,2)

          ! V_iijj
          Viijj=Vc(i1,j1)
          
          ! V_ijji
          Vijji=Vx(i1,j1)

          ! Sum the contibution
          contrib=contrib-(pJ*Viijj-pN0*Vijji)*Dwi*Dwj
          
       enddo
    enddo

    !
    ! Make the correction
    !
    harr=harr+contrib/nexci
    
    return
    
  end subroutine hii_dftmrci_grimme

!######################################################################
! damping_grimme: for two CSF-averaged on-diagonal matrix element
!                 values, returns the value of Grimme's original
!                 DFT/MRCI damping function
!######################################################################
  function damping_grimme(av1,av2) result(func)

    use constants
    use bitglobal
    use hparam
    
    implicit none

    ! Function result
    real(dp)             :: func

    ! CSF-averaged on-diagonal matrix elements
    real(dp), intent(in) :: av1,av2

    ! Everything else
    real(dp)             :: DE4
    
    !
    !  p1 exp(-p2 DeltaE^4)
    !
    DE4=(av1-av2)**4
    func=hpar(1)*exp(-hpar(2)*DE4)
    
    return
    
  end function damping_grimme

!######################################################################
  
end module dftmrci
