!**********************************************************************
! Routines for the calculation of reduced matrix elements
! < S_bra I || T_ij^(1,.) || S_ket J >
!**********************************************************************
module redmat

  implicit none

  private

  public :: reduced_matrix
  
contains

!######################################################################
! reduced_matrix: calculates the reduced matrix elements
! < S_bra I || T_ij^(1) || S_ket K > for a given pair of set of triplet
! TDMs < S_bra M_bra=S_bra I | T_ij^(1,k) | S_ket M_ket=S_ket J >
!######################################################################
  subroutine reduced_matrix(kindx,npairs,uij)

    use constants
    use bitglobal
    use clebsch_gordan
    use iomod
    
    implicit none

    ! Triplet spin tensor component
    integer(is), intent(in) :: kindx

    ! uij: on input contains the triplet TDMs,
    !      on output contains the reduced matrix elements
    integer(is), intent(in) :: npairs
    real(dp), intent(inout) :: uij(nmo,nmo,npairs)
    
    ! Bra and ket total spin angular momentum quantum numbers
    real(dp)                :: SB,SK

    ! Clebsch-Gordan coefficients
    integer(is)             :: dim1,dim2,dim12,dimtot
    real(dp), allocatable   :: cg(:,:)
    real(dp)                :: cgcoe
    
    ! Everything else
    integer(is)             :: i,i1,i2,i12
    real(dp)                :: m1,m2,M
    real(dp)                :: j1,j2,J
    real(dp), parameter     :: tiny=1e-10_dp
    
!----------------------------------------------------------------------
! Compute the Clebsch-Gordan coefficient
! < S_ket S_ket; 1 k | S_bra S_bra >
!----------------------------------------------------------------------
    ! Bra and ket total spins
    SB=0.5d0*(dble(imultB)-1.0d0)
    SK=0.5d0*(dble(imultK)-1.0d0)

    ! Angular momentum to be coupled
    j1=SK
    j2=1.0d0

    ! Total angular momentum
    J=SB

    ! Size of the direct product basis |j1 m1> x |j2 m2>
    dim1=int(2*j1)+1
    dim2=int(2*j2)+1
    dim12=dim1*dim2

    ! Number of |J M> states
    dimtot=int(2*J)+1
    
    ! Allocate the coefficient array
    allocate(cg(dim12,dimtot))
    cg=0.0d0
    
    ! Compute the full set of Clebsch-Gordan coefficients
    ! < j1m1; j2m2 | JM > for all ((m1,m2),M) pairs
    call cgcoeff(j1,j2,J,cg)

    ! Value of < S_ket S_ket; 1 k | S_bra S_bra >
    !                  ^        ^         ^
    !                  m1       m2        M
    m1=SK
    m2=real(kindx)
    M=SB
    i=int(M+J)+1
    i1=int(m1+j1)+1
    i2=int(m2+j2)+1
    i12=(i1-1)*dim2+i2
    cgcoe=cg(i12,i)

    ! Check on the Clebsch-Gordan coefficient value
    if (abs(cgcoe) < tiny) then
       errmsg='Error in reduced_matrix: zero Clebsch-Gordan ' &
            //'coefficient'
       call error_control
    endif

!----------------------------------------------------------------------
! Compute the reduced matrix elements in place
!----------------------------------------------------------------------
    uij=uij/cgcoe
    
    return
    
  end subroutine reduced_matrix
    
!######################################################################
  
end module redmat
