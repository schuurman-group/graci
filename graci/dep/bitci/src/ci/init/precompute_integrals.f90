!**********************************************************************
! Pre-computation of Fock matrix elements and Coulomb and exchange
! integrals
!**********************************************************************
module precompute

  implicit none

contains

!######################################################################
  
  subroutine precompute_integrals

    use constants
    use bitglobal
    use detutils
    
    implicit none

    integer(is) :: i,j,k,n
    integer(is) :: ablist(nmo,2)
    integer(is) :: iocc(nmo)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! On-diagonal elements of the Fock matrix
    allocate(fii(nmo))
    fii=0.0d0

    ! Fock-matrix
    allocate(fock(nmo,nmo))
    fock=0.0d0
    
    ! Coulomb integrals V_iijj
    allocate(Vc(nmo,nmo))
    Vc=0.0d0
    
    ! Exchange integrals V_ijji
    allocate(Vx(nmo,nmo))
    Vx=0.0d0

!----------------------------------------------------------------------
! Get the list of occupied alpha and beta MOs in the base determinant
!----------------------------------------------------------------------
    ! Lists of occupied alpha and beta orbitals
    call get_occ(det0,ablist)

    ! MO occupations
    iocc=0
    do i=1,nmo
       if (ablist(i,1) /= 0) iocc(i)=iocc(i)+1
       if (ablist(i,2) /= 0) iocc(i)=iocc(i)+1
    enddo

!----------------------------------------------------------------------
! Coulomb and exchange integrals
!----------------------------------------------------------------------
    do i=1,nmo
       do j=i,nmo
          Vc(i,j)=bitci_ints%mo_int(i,i,j,j)
          Vx(i,j)=bitci_ints%mo_int(i,j,j,i)
          Vc(j,i)=Vc(i,j)
          Vx(j,i)=Vx(i,j)
       enddo
    enddo
    
!----------------------------------------------------------------------
! On-diagonal Fock matrix elements in the HF/DFT ordering
!----------------------------------------------------------------------
    ! Loop over MOs
    do i=1,nmo
       
       ! Core Hamiltonian contribution
       fii(i) = bitci_ints%h_1e(i,i)
       
       ! Coulomb and exchange integral contributions
       do j=1,nmo
          if (iocc(j) == 0) cycle
          fii(i)=fii(i)+iocc(j)*(Vc(i,j)-0.5d0*Vx(i,j))
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Complete Fock matrix in the HF/DFT ordering
!----------------------------------------------------------------------
    ! Loop over pairs of MOs
    do i=1,nmo
       do j=i,nmo

          ! If this is a DFT/MRCI run, then the off-diagonal elements
          ! are not needed

          ! TEST
          if (ldftmrci .and. i /= j) cycle
          ! TEST
          
          ! Core Hamiltonian contribution
          fock(i,j) = bitci_ints%h_1e(i,j)

          ! Coulomb and exchange contributions
          do k=1,nmo
             if (iocc(k) == 0) cycle
             fock(i,j)=fock(i,j)+iocc(k)*(bitci_ints%mo_int(i,j,k,k) &
                  -0.5d0*bitci_ints%mo_int(i,k,k,j))
          enddo

          ! F_ji=F_ij
          fock(j,i)=fock(i,j)

       enddo
    enddo
    
    return
    
  end subroutine precompute_integrals

!######################################################################
  
end module precompute
