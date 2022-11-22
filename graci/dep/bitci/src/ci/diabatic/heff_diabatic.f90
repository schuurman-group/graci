!**********************************************************************
! Routines for the calculation of the diabatic potential matrix
! as the 2nd-order GVVPT2 effective Hamiltonian for a set of
! quasi-diabatic model states
!**********************************************************************
module heff_diabatic

  implicit none

contains

!######################################################################
! heff_diab: construction of the GVVPT2 effective Hamiltonian for
!            a quasi-diabatic model space
!######################################################################
  subroutine heff_diab(refdim,csfdim,nroots,hdiag,Avec,H01,ireg,&
       regfac,Heff)

    use constants
    use bitglobal
    
    implicit none

    ! Dimensions
    integer(is), intent(in) :: refdim,csfdim,nroots
    
    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(in)    :: hdiag(csfdim)

    ! A-vectors
    real(dp), intent(inout) :: Avec(csfdim,nroots)

    ! Zeroth- plus first-order effecive Hamiltonian
    real(dp), intent(in)    :: H01(nroots,nroots)

    ! Regularizer index: 1 <-> ISA
    !                    2 <-> sigma^p
    integer(is), intent(in) :: ireg
  
    ! Regularization factor
    real(dp), intent(in)    :: regfac

    ! Effective Hamiltonian up to second order
    real(dp), intent(out)   :: Heff(nroots,nroots)
    
    ! Zeroth-order effective Hamiltonian (diagonal)
    real(dp), allocatable   :: H0(:)
    
    ! Everything else
    integer(is)             :: i,j,icsf
    real(dp)                :: deltai,deltaj,Vi,Vj
    real(dp)                :: fi,fj,epsilon

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(H0(nroots))
    H0=0.0d0

!----------------------------------------------------------------------
! Effective Hamiltonian
!----------------------------------------------------------------------
    !
    ! Zeroth-order Hamiltonian elements
    !
    do i=1,nroots
       H0(i)=H01(i,i)
    enddo

    !
    ! Zeroth- and first-order contributions
    !
    Heff=H01

    !
    ! Second-order contributions
    !
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

       enddo
    enddo

!----------------------------------------------------------------------
! First-order perturbed model states (quasi-diabatic states)
!----------------------------------------------------------------------
    ! Loop over roots
    do j=1,nroots

       ! Loop over CSFs (excluding the reference space ones)
       do icsf=refdim+1,csfdim

          ! E^(0) - H_ii
          deltai=H0(j)-hdiag(icsf)

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

          ! First-order perturbed model state coefficient
          Avec(icsf,j)=Vi*fi
          
       enddo
       
    enddo
    
    return
    
  end subroutine heff_diab
  
!######################################################################
  
end module heff_diabatic
