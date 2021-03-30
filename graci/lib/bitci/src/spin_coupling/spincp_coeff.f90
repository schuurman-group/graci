!**********************************************************************
! Functions for the retrieval of individual spin-coupling coefficients
! for a given pair of bra and ket CSFs and pairs of
! creation/annihilation operator indices.
!**********************************************************************
! This is intended to be used mainly for debugging purposes
! For Hamiltonian and density matrix builds, it is more efficient
! to precompute batches of pattern indices and spin-coupling sub-case
! bit strings, and to use them to retrieve spin-coupling coefficients
! on-the-fly.
!**********************************************************************
module spincp_coeff

  implicit none

contains

!######################################################################
! spincp_1electron: For a given ket configuration (ksop), bra and ket
!                   spin couplings (indexed by bomega and komega), and
!                   creation and annihilation operator indices (ac and
!                   ia, respectively), returns a single one-electron
!                   spin-coupling coefficient
!                   < w' omega' | E_a^i | w omega >
!######################################################################
  function spincp_1electron(ksop,bomega,komega,ac,ia,knopen) &
       result(spincp)

    use constants
    use bitglobal
    use bitstrings
    use pattern_indices
    use mrciutils
    
    implicit none

    ! Function result
    real(dp)                :: spincp

    ! Ket SOP
    integer(ib), intent(in) :: ksop(n_int,2)

    ! Bra and ket spin-coupling indices
    integer(is), intent(in) :: bomega,komega

    ! Creation and annihilation operator indices
    integer(is), intent(in) :: ac,ia

    ! Number of open shells in the ket configuration
    integer(is), intent(in) :: knopen
    
    ! Spin-coupling sub-case index
    integer(ib)             :: icase

    ! Pattern index
    integer(is)             :: indx

    ! No. open shells preceding the created/annhilated MOs
    integer(is)             :: nc,na
    
!----------------------------------------------------------------------
! Spin-coupling sub-case index
!----------------------------------------------------------------------
    icase=get_icase(ksop,ac,ia)

!----------------------------------------------------------------------
! Number of open-shells before the creation and annihilation operator
! indices
!----------------------------------------------------------------------
    nc=n_bits_set_before(ksop(:,1),n_int,ac)
    na=n_bits_set_before(ksop(:,1),n_int,ia)
    
!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    indx=pattern_index(ksop,ac,ia,nc,na,knopen,icase)

!----------------------------------------------------------------------
! Retrieve the spin-coupling coefficient
!----------------------------------------------------------------------
    select case(icase)
    case(i1a)
       spincp=spincp1(komega,bomega,indx)
    case(i1b)
       spincp=-spincp1(komega,bomega,indx)
    case(i2a)
       spincp=spincp2(komega,bomega,indx)
    case(i2b)
       spincp=spincp2(bomega,komega,indx)
    end select
    
    return
    
  end function spincp_1electron

!######################################################################
! spincp_2electron: For a given pair of bra and ket configuration
!                   (bsop and ksop), bra and ket
!                   spin couplings (indexed by bomega and komega), and
!                   creation and annihilation operator indices (ac and
!                   ia, respectively), returns a single one-electron
!                   spin-coupling coefficient
!                   < w' omega' | E_a^i | w omega >
!######################################################################
  function spincp_2electron(bsop,ksop,bomega,komega,plist,hlist,&
       bnopen,knopen) result(spincp)

    use constants
    use bitglobal
    use bitstrings
    use pattern_indices
    use mrciutils
    use mrci_integrals
    
    implicit none

    ! Function result
    real(dp)                :: spincp

    ! Bra and ket SOPs
    integer(ib), intent(in) :: bsop(n_int,2),ksop(n_int,2)

    ! Bra and ket spin-coupling indices
    integer(is), intent(in) :: bomega,komega

    ! Creation and annihilation operator indices
    integer(is), intent(in) :: plist(2),hlist(2)
    integer(is)             :: ac,bc,ia,ja

    ! Number of open shells in the bra and ket configurations
    integer(is), intent(in) :: bnopen,knopen

    ! Pattern indices
    integer(is)             :: bindx,kindx
    
    ! Spin-coupling sub-case bitstrings
    integer(ib)             :: bicase,kicase

    ! Spin-coupling sub-case pair bitstring
    integer(ib)             :: pairindx

    ! No. open shells preceding the created/annhilated MOs
    integer(is)             :: nc,na
    
!----------------------------------------------------------------------    
! Creation and annihilation operator indices
!----------------------------------------------------------------------
    ! Creation operators
    ac=plist(1)
    bc=plist(2)

    ! Annihilation operators
    ia=hlist(1)
    ja=hlist(2)

!----------------------------------------------------------------------
! Bitstring encodings of the spin-coupling sub-cases
!----------------------------------------------------------------------
    ! Bra
    bicase=get_icase(bsop,ia,ac)

    ! Ket
    kicase=get_icase(ksop,bc,ja)

    print*,caselbl(bicase),caselbl(kicase)
    
!----------------------------------------------------------------------
! Bitstring encodings of the spin-coupling sub-case pair
!----------------------------------------------------------------------
    pairindx=ipair(bicase,kicase)

!----------------------------------------------------------------------
! Pattern indices
!----------------------------------------------------------------------
    ! Bra
    nc=n_bits_set_before(bsop(:,1),n_int,ia)
    na=n_bits_set_before(bsop(:,1),n_int,ac)
    bindx=pattern_index(bsop,ia,ac,nc,na,bnopen,bicase)
    
    ! Ket
    nc=n_bits_set_before(ksop(:,1),n_int,bc)
    na=n_bits_set_before(ksop(:,1),n_int,ja)
    kindx=pattern_index(ksop,bc,ja,nc,na,knopen,kicase)

!----------------------------------------------------------------------
! Compute the spin-coupling coefficient
!----------------------------------------------------------------------
    spincp=contract_spincp(bomega,komega,bindx,kindx,pairindx,&
         bnopen,knopen)

    return
    
  end function spincp_2electron
    
!######################################################################
  
end module spincp_coeff
