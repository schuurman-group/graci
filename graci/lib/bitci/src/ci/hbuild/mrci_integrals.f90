!**********************************************************************
! Routines for the packaging and retrieval of all the two-electron and
! spin-coupling coefficient integrals needed to evaluate the
! Hamiltonian matrix elements
!**********************************************************************
module mrci_integrals

  implicit none

contains

!######################################################################
! package_confinfo: returns all the configuration information needed
!                   to evaluate a diagonal matrix element
!######################################################################
  subroutine package_confinfo(sop,conf,unocc,socc,docc,nunocc,nsocc,&
       ndocc,Dw,ndiff,nbefore)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    ! SOP and conf bitstrings in the full MO space
    integer(ib), intent(in)  :: sop(n_int,2)
    integer(ib), intent(in)  :: conf(n_int,2)
    
    ! Lists of MOs in the different occupation classes
    integer(is), intent(out) :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is), intent(out) :: nsocc,ndocc,nunocc

    ! Difference configuration information
    integer(is), intent(out) :: ndiff
    integer(is), intent(out) :: Dw(nmo,2)

    ! Number of open shells preceding each MO
    integer(is), intent(out) :: nbefore(nmo)

    !
    ! Indices of the singly-occupied MOs
    !
    call sop_socc_list(sop,n_int,socc,nmo,nsocc)

    !
    ! Indices of the doubly-occupied MOs
    !
    call sop_docc_list(sop,n_int,docc,nmo,ndocc)

    !
    ! Indices of the unoccupied MOs
    !
    call sop_unocc_list(sop,n_int,unocc,nmo,nunocc)

    !
    ! Difference configuration relative to the base configuration
    !
    call diffconf(conf,n_int,Dw,nmo,ndiff)

    !
    ! Number of open shells preceding each MO
    !
    call nobefore(sop,nbefore)
    
    return
    
  end subroutine package_confinfo

!######################################################################
! package_confinfo_offdiag: returns all the configuration information
!                           needed to evaluate an off-diagonal matrix
!                            element  
!######################################################################
  subroutine package_confinfo_offdiag(sop,conf,socc,nsocc,Dw,ndiff,&
       nbefore)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    ! SOP and conf bitstrings in the full MO space
    integer(ib), intent(in)  :: sop(n_int,2)
    integer(ib), intent(in)  :: conf(n_int,2)
    
    ! Lists of sigly-occupied MOs
    integer(is), intent(out) :: socc(nmo)
    integer(is), intent(out) :: nsocc

    ! Difference configuration information
    integer(is), intent(out) :: ndiff
    integer(is), intent(out) :: Dw(nmo,2)

    ! Number of open shells preceding each MO
    integer(is), intent(out) :: nbefore(nmo)

    !
    ! Indices of the singly-occupied MOs
    !
    call sop_socc_list(sop,n_int,socc,nmo,nsocc)
    
    !
    ! Difference configuration relative to the base configuration
    !
    call diffconf(conf,n_int,Dw,nmo,ndiff)
    
    !
    ! Number of open shells preceding each MO
    !
    call nobefore(sop,nbefore)
    
    return
    
  end subroutine package_confinfo_offdiag

!######################################################################
! Packaging of the two-electron integrals and spin-coupling
! coefficient information needed to evaluate a Hamiltonian matrix
! element between two CSFs linked by a single excitation
!######################################################################
  subroutine package_integrals_nexci1(bsop,ksop,hindx,pindx,bnopen,&
       knopen,bpattern,kpattern,Vpqrs,m2c,socc,nsocc,knbefore,Dw,&
       ndiff,icase,insp)

    use constants
    use bitglobal
    use pattern_indices
    use bitstrings
    use mrciutils
    use int_pyscf
    use hparam
    
    implicit none

    ! Bra and ket SOPs
    integer(ib), intent(in)  :: bsop(n_int,2),ksop(n_int,2)

    ! Spin-coupling sub-case bitstring encodings
    integer(ib), intent(out) :: icase
    
    ! Indices of the annihilation and creation operators
    integer(is), intent(in)  :: hindx,pindx
    integer(is)              :: ia,ac,ia1,ac1
    integer(is)              :: i1
    
    ! Numbers of open shells
    integer(is), intent(in)  :: bnopen,knopen

    ! Pattern indices
    integer(is), intent(out) :: bpattern(nmo+1),kpattern(nmo+1)
    
    ! Spin-coupling sub-case bitstrings
    integer(ib)              :: icase_b,icase_k

    ! Integrals and functions of integrals
    real(dp), intent(out)    :: Vpqrs(nmo)

    ! MO mapping array
    integer(is), intent(in)  :: m2c(nmo)

    ! Indices of the singly-occupied MOs in the ket configuration
    integer(is), intent(in)  :: nsocc
    integer(is), intent(in)  :: socc(nmo)

    ! Number of open shells preceding each MO
    integer(is), intent(in)  :: knbefore(nmo)
    integer(is)              :: nc_k,na_k,nc_b,na_b

    ! No. CSFs for the intermediate configuration in the case
    ! of spin-coupling coefficients
    ! <w' omega' | E_i^j E_k^l | w omega>
    integer(is), intent(out) :: insp(nmo)
    
    ! Difference configuration information
    integer(is)             :: ndiff
    integer(is)             :: Dw(nmo,2)
    
    ! Everything else
    integer(is)             :: i,k,k1,count
    integer(is)             :: wac,wia
    integer(is)             :: bnsp,knsp
    
!----------------------------------------------------------------------
! Indices of the annihilation and creation operators
!----------------------------------------------------------------------
    ! Annihilation operator
    ia=hindx

    ! Creation operator
    ac=pindx

!----------------------------------------------------------------------
! Occupations of the MOs corresponding to the annihilation and
! creation operators
!----------------------------------------------------------------------
    ! Creation operator
    k=(ac-1)/64+1
    i=ac-(k-1)*64-1
    if (btest(ksop(k,1),i)) then
       wac=1
    else
       wac=0
    endif

    ! Annihilation operator
    k=(ia-1)/64+1
    i=ia-(k-1)*64-1
    if (btest(ksop(k,1),i)) then
       wia=1
    else
       wia=2
    endif 
    
!----------------------------------------------------------------------
! Bitstring encoding of the sub-case for the spin-coupling coefficients
! < w' omega'| E_a^k E_k^i | w omega >
!----------------------------------------------------------------------
    if (nsocc > 0) then
       do k=1,nsocc
          k1=socc(k)
          if (k1 == ia) cycle
          if (k1 == ac) cycle
          icase_b=get_icase(bsop,k1,ac)
          icase_k=get_icase(ksop,k1,ia)
          exit
       enddo
    endif

!----------------------------------------------------------------------
! Number of CSFs for the intermediate configuration obtained by acting
! on the ket CSF with with first singlet excitation operator for the
! spin-coupling coefficients < w' omega'| E_a^k E_k^i | w omega >
!----------------------------------------------------------------------
    ! No. ket CSFs
    knsp=ncsfs(knopen)

    ! Set the number of intermediate CSFs
    if (nsocc > 0) then
       do k=1,nsocc
          k1=socc(k)
          if (k1 == ia) cycle
          if (k1 == ac) cycle
          select case(icase_k)
          case(i1a)
             ! Case 1a
             insp(k)=knsp
          case(i1b)
             ! Case 1b
             insp(k)=knsp
          case(i2a)
             ! Case 2a
             insp(k)=ncsfs(knopen+2)
          case(i2b)
             ! Case 2b
             insp(k)=ncsfs(knopen-2)
          end select
       enddo
    endif

!----------------------------------------------------------------------
! Pattern indices for the spin-coupling coefficients
! < w' omega'| E_a^k E_k^i | w omega > (bpattern and kpattern)
!----------------------------------------------------------------------    
    ! No. open shells before the annihilation operator index in the
    ! ket configuration
    na_k=knbefore(ia)

    ! No. open shells before the creation operator index in the
    ! ket configuration
    na_b=n_bits_set_before(bsop(:,1),n_int,ac)
    
    ! Loop over singly-occupied MOs in the ket configuration
    do k=1,nsocc
       
       ! MO index
       k1=socc(k)

       ! Cycle if k=a or k=i
       if (k1 == ia) cycle
       if (k1 == ac) cycle
       
       ! bra
       nc_b=n_bits_set_before(bsop(:,1),n_int,k1)
       bpattern(k)=pattern_index(bsop,k1,ac,nc_b,na_b,bnopen,icase_b)
       
       ! ket
       nc_k=knbefore(k1)
       kpattern(k)=pattern_index(ksop,k1,ia,nc_k,na_k,knopen,icase_k)
       
    enddo

!----------------------------------------------------------------------
! Pattern index for the spin-coupling coefficients
! < w' omega'| E_a^i | w omega > (icase)
!----------------------------------------------------------------------
    nc_k=knbefore(ac)
    na_k=knbefore(ia)
    icase=get_icase(ksop,ac,ia)
    kpattern(nsocc+1)=pattern_index(ksop,ac,ia,nc_k,na_k,knopen,icase)
    
!----------------------------------------------------------------------
! Two-electron integrals
!----------------------------------------------------------------------
    ! Initialise the integral counter
    count=0

    ! HF/DFT MO indices
    ia1=m2c(ia)
    ac1=m2c(ac)
    
    !
    ! V_ikka, k indexing a singly-occupied MO in the ket configuration
    !
    do k=1,nsocc

       ! Increment the integral counter
       count=count+1

       ! MO index
       k1=m2c(socc(k))

       ! V_ikka
       Vpqrs(count)=mo_integral(ia1,k1,k1,ac1)
       
    enddo

    !
    ! F_ia + Sum_k (V_iakk - 1/2 V_ikka) Delta w_k
    !
    count=count+1
    if (ihamiltonian == 1) then
       ! Canonical Hamiltonian
       Vpqrs(count)=fock(ia1,ac1)
    else
       ! DFT/MRCI Hamiltonians - F_ia contribution is neglected
       Vpqrs(count)=0.0d0
    endif
    do k=1,ndiff

       ! MO index
       k1=m2c(Dw(k,1))

       ! (V_iakk - 1/2 V_ikka) Delta w_k
       Vpqrs(count)=Vpqrs(count) &
            +(mo_integral(ia1,ac1,k1,k1) &
            -0.5d0*mo_integral(ia1,k1,k1,ac1))*Dw(k,2)
       
    enddo

    !
    ! 1/2 [ V_aaai w_a + V_aiii (w_i -2) ]
    !
    count=count+1
    Vpqrs(count)=0.5d0*(mo_integral(ac1,ac1,ac1,ia1)*wac &
         +mo_integral(ac1,ia1,ia1,ia1)*(wia-2))
    
    return
    
  end subroutine package_integrals_nexci1  
  
!######################################################################
! Packaging of the two-electron integrals and spin-coupling
! coefficient information needed to evaluate a Hamiltonian matrix
! element between two CSFs linked by a double excitation
!**********************************************************************
! IMPORTANT: We scale the two-electron integrals by the
!            1 / [(1+delta_ab) * (1+delta_ij)] prefactor to save
!            effort down the road
!######################################################################
  subroutine package_integrals_nexci2(bsop,ksop,hlist,plist,bnopen,&
       knopen,bpattern,kpattern,Vpqrs,m2c,knbefore,insp)

    use constants
    use bitglobal
    use pattern_indices
    use bitstrings
    use mrciutils
    use int_pyscf
    
    implicit none

    ! Bra and ket SOPs
    integer(ib), intent(in)  :: bsop(n_int,2),ksop(n_int,2)

    ! Indices of the annihilation and creation operators
    integer(is), intent(in)  :: hlist(2),plist(2)
    integer(is)              :: ia,ja,ac,bc
    integer(is)              :: i1,j1,a1,b1
    
    ! Numbers of open shells
    integer(is), intent(in)  :: bnopen,knopen

    ! Pattern indices
    integer(is), intent(out) :: bpattern(2),kpattern(2)
    
    ! Spin-coupling sub-case bitstrings
    integer(ib)              :: icase_b(2),icase_k(2)

    ! Integrals
    real(dp), intent(out)    :: Vpqrs(2)

    ! MO mapping array
    integer(is), intent(in)  :: m2c(nmo)

    ! Numbers of open shells preceding the annihilation
    ! and creation operator indices
    integer(is), intent(in)  :: knbefore(nmo)
    integer(is)              :: nc,na

    ! Number of CSFs for the intermediate configuration obtained
    ! by acting on the ket CSF with the first singlet excitation
    ! operator
    integer(is), intent(out) :: insp(2)
    
    ! Prefactor
    real(dp)                 :: fac1,fac2

    ! Number of bra and ket CSFs
    integer(is)              :: bnsp,knsp

    ! Everything else
    integer(is)              :: m
    
!----------------------------------------------------------------------
! Indices of the annihilation and creation operators
!----------------------------------------------------------------------
    ! Annihilation operators
    ia=hlist(1)
    ja=hlist(2)

    ! Creation operators
    ac=plist(1)
    bc=plist(2)

!----------------------------------------------------------------------
! Bitstring encodings of the spin-coupling sub-cases
!----------------------------------------------------------------------
    ! bra, V_aibj
    icase_b(1)=get_icase(bsop,ia,ac)

    ! bra, V_ajbi
    icase_b(2)=get_icase(bsop,ja,ac)
    
    ! ket, V_aibj
    icase_k(1)=get_icase(ksop,bc,ja)

    ! ket, V_ajbi
    icase_k(2)=get_icase(ksop,bc,ia)

!----------------------------------------------------------------------
! Number of CSFs for the intermediate configuration obtained by acting
! on the ket CSF with with first singlet excitation operator
!----------------------------------------------------------------------
    ! No. ket CSFs
    knsp=ncsfs(knopen)

    ! Set the no. intermediate CSFs
    do m=1,2
       select case(icase_k(m))
       case(i1a)
          ! Case 1a
          insp(m)=knsp
       case(i1b)
          ! Case 1b
          insp(m)=knsp
       case(i2a)
          ! Case 2a
          insp(m)=ncsfs(knopen+2)
       case(i2b)
          ! Case 2b
          insp(m)=ncsfs(knopen-2)
       end select
    enddo
          
!----------------------------------------------------------------------
! Pattern indices
!----------------------------------------------------------------------
    ! bra, V_aibj
    nc=n_bits_set_before(bsop(:,1),n_int,ia)
    na=n_bits_set_before(bsop(:,1),n_int,ac)
    bpattern(1)=pattern_index(bsop,ia,ac,nc,na,bnopen,icase_b(1))
    
    ! bra, V_ajbi
    nc=n_bits_set_before(bsop(:,1),n_int,ja)
    bpattern(2)=pattern_index(bsop,ja,ac,nc,na,bnopen,icase_b(2))
    
    ! ket, V_aibj
    nc=knbefore(bc)
    na=knbefore(ja)
    kpattern(1)=pattern_index(ksop,bc,ja,nc,na,knopen,icase_k(1))
    
    ! ket, V_ajbi
    na=knbefore(ia)
    kpattern(2)=pattern_index(ksop,bc,ia,nc,na,knopen,icase_k(2))
    
!----------------------------------------------------------------------
! Two-electron integrals scaled by the 1/[(1+delta_ab)*(1+delta_ij)]
! prefactor
!----------------------------------------------------------------------
    ! MO indices
    a1=m2c(ac)
    b1=m2c(bc)
    i1=m2c(ia)
    j1=m2c(ja)
    
    ! V_aibj
    Vpqrs(1)=mo_integral(a1,i1,b1,j1)
    
    ! V_ajbi
    Vpqrs(2)=mo_integral(a1,j1,b1,i1)

    ! Scaling by the prefactor
    if (ac == bc) then
       fac1=2.0d0
    else
       fac1=1.0d0
    endif
    if (ia == ja) then
       fac2=2.0d0
    else
       fac2=1.0d0
    endif
    Vpqrs=Vpqrs/(fac1*fac2)
    
    return
    
  end subroutine package_integrals_nexci2
  
!######################################################################
  
end module mrci_integrals
