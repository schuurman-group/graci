!**********************************************************************
! Routines for the calculation of the diabatic potential matrix
! as the 2nd-order GVVPT2 effective Hamiltonian for a set of
! quasi-diabatic model states
!**********************************************************************
module diabpot

  implicit none

contains

!######################################################################
! heff_diab: construction of the GVVPT2 effective Hamiltonian for
!            a quasi-diabatic model space
!######################################################################
  subroutine heff_diab(refdim,csfdim,nroots,hdiag,Avec,H01,ireg,regfac)

    use constants
    use bitglobal
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: refdim,csfdim,nroots
    
    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(csfdim)

    ! A-vectors
    real(dp), intent(in)    :: Avec(csfdim,nroots)

    ! Zeroth- plus first-order effecive Hamiltonian
    real(dp), intent(in)    :: H01(nroots,nroots)

    ! Regularizer index: 1 <-> ISA
    !                    2 <-> sigma^p
    integer(is), intent(in) :: ireg
  
    ! Regularization factor
    real(dp), intent(in)    :: regfac
    
    ! Zeroth-order effective Hamiltonian (diagonal)
    real(dp), allocatable   :: H0(:)
    
    ! Effective Hamiltonian up to second order
    real(dp), allocatable   :: Heff(:,:)
    
    ! Everything else
    integer(is)             :: i,j,icsf
    real(dp)                :: deltai,deltaj,Vi,Vj
    real(dp)                :: fi,fj,epsilon

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(Heff(nroots,nroots))
    Heff=0.0d0

    allocate(H0(nroots))
    H0=0.0d0

!----------------------------------------------------------------------
! Zeroth-order effective Hamiltonian elements
!----------------------------------------------------------------------
    do i=1,nroots
       H0(i)=H01(i,i)
    enddo

!----------------------------------------------------------------------
! Zeroth- and first-order contributions
!----------------------------------------------------------------------
    Heff=H01

!----------------------------------------------------------------------
! Second-order contributions
!----------------------------------------------------------------------

    print*,''

    ! Loop over pairs of roots
    do i=1,nroots
       do j=i,nroots

          ! Loop over FOIS CSFs
          do icsf=refdim+1,csfdim

             ! <I|H|i> and <i|H|J>
             Vi=Avec(icsf,i)
             Vj=Avec(icsf,j)

             ! Regularized denominators
             deltai=H0(i)-hdiag(icsf)
             deltaj=H0(j)-hdiag(icsf)
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
             Heff(i,j)=Heff(i,j)+0.5d0*Vi*Vj*(fi+fj)
             
          enddo

          Heff(j,i)=Heff(i,j)

          print*,i,j,Heff(i,j)
          
       enddo
    enddo

    STOP
    
    return
    
  end subroutine heff_diab
  
!######################################################################
  
end module diabpot
