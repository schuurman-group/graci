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

!----------------------------------------------------------------------
! [INTCHK] TEMPORARY DIAGNOSTIC -- remove once the C1-after-C1 bug is
! closed.
!
! Fingerprints the integrals this calculation was actually built from.
! Two CI calculations on identical input must print identical values
! here. c1c1.inp gives identical configuration counts (8586, 1525) and
! identical saved-element counts (4753686, 300898) but different
! energies, so the Hamiltonian VALUES differ. This says whether that
! starts in the integrals or downstream of them.
!----------------------------------------------------------------------
    write(6,'(/,x,a)')         '[INTCHK] ---------------------------------'
    write(6,'(x,a,i0)')        '[INTCHK] nmo       : ',nmo
    write(6,'(x,a,i0)')        '[INTCHK] sum(iocc) : ',sum(iocc)
    write(6,'(x,a,es25.17)')   '[INTCHK] sum|fii|  : ',sum(abs(fii))
    write(6,'(x,a,es25.17)')   '[INTCHK] sum|fock| : ',sum(abs(fock))
    write(6,'(x,a,es25.17)')   '[INTCHK] sum|Vc|   : ',sum(abs(Vc))
    write(6,'(x,a,es25.17)')   '[INTCHK] sum|Vx|   : ',sum(abs(Vx))
    write(6,'(x,a,es25.17)')   '[INTCHK] h_1e(1,1) : ',bitci_ints%h_1e(1,1)
    write(6,'(x,a,es25.17)')   '[INTCHK] h_1e(N,N) : ',bitci_ints%h_1e(nmo,nmo)
    write(6,'(x,a,es25.17)')   '[INTCHK] (11|11)   : ',&
         bitci_ints%mo_int(1,1,1,1)
    write(6,'(x,a,es25.17)')   '[INTCHK] (NN|NN)   : ',&
         bitci_ints%mo_int(nmo,nmo,nmo,nmo)
    write(6,'(x,a)')           '[INTCHK] ---------------------------------'
    flush(6)

    return
    
  end subroutine precompute_integrals

!######################################################################
  
end module precompute
