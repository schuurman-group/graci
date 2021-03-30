!**********************************************************************
! Routines for the direct calculation of the sigma vector
!**********************************************************************

module sigma

contains

!######################################################################
! h_ondiag: Computes the on-diagonal elements of the Hamiltonian
!           matrix
!######################################################################
  subroutine h_ondiag(irrep,hdiag,ndet,offdim,nsym,d,nunique,offset)
    
    use constants
    use bitglobal
    use utils
    use slater_condon
    
    implicit none

    ! Determinant arrays
    integer(is), intent(in) :: irrep
    integer(is), intent(in) :: ndet,offdim
    integer(is), intent(in) :: nsym(0:nirrep-1)
    integer(ib), intent(in) :: d(n_int,2,ndet)
    integer(is), intent(in) :: nunique(0:nirrep-1)
    integer(is), intent(in) :: offset(offdim,0:nirrep-1)
    
    ! On-diagonal elements of the Hamiltonian matrix
    real(dp), intent(out)   :: hdiag(nsym(irrep))

    ! Hamiltonian matrix build variables
    integer(is)              :: i,a,b
    integer(is)              :: idet
    integer(is)              :: occlist(nmo,2)
    integer(is)              :: nsym_sum(0:nirrep-1)
    
!----------------------------------------------------------------------
! Cumulative numbers of determinants in the symmetry blocks
!----------------------------------------------------------------------
    nsym_sum=0
    nsym_sum(1:nirrep-1)=nsym(0:nirrep-2)
    do i=1,nirrep-1
       nsym_sum(i)=nsym_sum(i)+nsym_sum(i-1)
    enddo

!----------------------------------------------------------------------
! Calculate the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
    ! Loop over unique alpha strings for the current irrep
    do a=1,nunique(irrep)

       ! Loop over determinants with the current alpha
       ! string
       do b=offset(a,irrep),offset(a+1,irrep)-1
          
          ! Determinant index
          idet=b-nsym_sum(irrep)
          
          ! Get the indices occupied alpha and beta spin orbitals in
          ! the ket determinant
          call occ_orbs(d(:,:,b),occlist)
          
          ! On-diagonal element of H
          hdiag(idet)=hii(occlist)

       enddo
       
    enddo
    
    return
    
  end subroutine h_ondiag

!######################################################################
  
end module sigma
