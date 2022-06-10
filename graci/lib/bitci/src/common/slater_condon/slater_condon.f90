!**********************************************************************
! Routines for the evaluation of the Slater-Condon rules using the bit
! string representation of Slater determinants
!**********************************************************************
! Makes heavy use of the functions given in Chapter 3 of Yann
! Garniron's PhD Thesis:
!
! `Development and parallel implementation of selected configuration
! interaction methods', PhD thesis, Université de Toulouse,
! February 2019. URL: https://doi.org/10.5281/zenodo.2558127,
! doi:10.5281/zenodo.2558127
!
!**********************************************************************

module slater_condon

contains

!######################################################################
! exc_degree_det: Computes the excitation degree between two
!                 determinants d1 and d2
!######################################################################
  function exc_degree_det(d1,d2) result(nexci)

    use constants
    use bitglobal
    
    implicit none

    integer(is)             :: nexci
    
    integer(ib), intent(in) :: d1(n_int,2),d2(n_int,2)
    integer(is)             :: ispin,k

    ! Compute the number of creation and annihilation operators
    ! connecting d1 and d2
    nexci=popcnt(ieor(d1(1,1),d2(1,1))) &
         +popcnt(ieor(d1(1,2),d2(1,2)))
    do k=2,n_int
       nexci=nexci+popcnt(ieor(d1(k,1),d2(k,1))) &
            +popcnt(ieor(d1(k,2),d2(k,2)))
    enddo

    ! Divide by 2 to get the no. excitations connecting d1 and d2
    ! Note that this is equivalent to performing a right bit shift
    ! by 1 place
    nexci=shiftr(nexci,1)

    return
    
  end function exc_degree_det
    
!######################################################################
! exc_degree_string: Computes the excitation degree between two alpha
!                    or beta strings s1 and s2
!######################################################################  
  function exc_degree_string(s1,s2) result(nexci)

    use constants
    use bitglobal
    
    implicit none

    integer(is)             :: nexci
    
    integer(ib), intent(in) :: s1(n_int),s2(n_int)
    integer(is)             :: ispin,k

    ! Compute the number of creation and annihilation operators
    ! connecting s1 and s2
    nexci=popcnt(ieor(s1(1),s2(1)))
    do k=2,n_int
       nexci=nexci+popcnt(ieor(s1(k),s2(k)))
    enddo

    ! Divide by 2 to get the no. excitations connecting s1 and s2
    ! Note that this is equivalent to performing a right bit shift
    ! by 1 place
    nexci=shiftr(nexci,1)

    return
    
  end function exc_degree_string
    
!######################################################################
! exc: Given two determinants d1 and d2, returns bitstrings p and h
!      encoding, respectively, the indices of the particles and holes
!      created by the excitation linking d1 and d2
!######################################################################
  subroutine exc(d1,d2,p,h)

    use constants
    use bitglobal
    
    implicit none

    integer(ib), intent(in)  :: d1(n_int,2),d2(n_int,2)
    integer(ib), intent(out) :: p(n_int,2),h(n_int,2)
    integer(ib)              :: c
    integer(is)              :: ispin,k
    
    ! Loop over spins
    do ispin=1,2
       ! Loop over blocks
       do k=1,n_int
          c=ieor(d1(k,ispin),d2(k,ispin))
          p(k,ispin)=iand(c,d2(k,ispin))
          h(k,ispin)=iand(c,d1(k,ispin))
       enddo
    enddo
    
    return
    
  end subroutine exc
  
!######################################################################
! phasemask: Computes the phase mask associated with the determinant d
!######################################################################
  
  subroutine phasemask(d,mask)

    use constants
    use bitglobal
    
    implicit none

    integer(ib), intent(in)  :: d(n_int,2)
    integer(ib), intent(out) :: mask(n_int,2)
    integer(is)              :: ispin,k
    integer(ib)              :: r,n,cnt
    
    ! Loop over spins
    do ispin=1,2

       r=0_ib

       ! Loop over blocks
       do k=1,n_int

          mask(k,ispin)=ieor(d(k,ispin),1_ib)

          do n=0,5
             mask(k,ispin)=ieor(mask(k,ispin),&
                  ishft(mask(k,ispin),2_ib**n))
          enddo

          mask(k,ispin)=ieor(mask(k,ispin),r)

          cnt=popcnt(d(k,ispin))
          
          if (iand(cnt,1_ib)==1) then
             r=not(r)
          endif
          
       enddo
       
    enddo
    
    return
    
  end subroutine phasemask

!######################################################################
! phase_phasemask: For a given phase mask for some determinant |I>
!                  and the single excitation i_sigma -> j_sigma,
!                  returns the corresponding phase factor.
!
!                  Note that here mask is the part of the phase mask
!                  for the single spin, sigma, of the excitation
!                  operator.
!######################################################################
  function phase_phasemask(mask,i,j) result(phase)

    use constants
    use bitglobal
    
    implicit none

    integer(is)             :: phase

    integer(is), intent(in) :: i,j
    integer(ib), intent(in) :: mask(n_int)
    integer(is)             :: c,in,jn
    integer(ib)             :: ibb,jbb,B

    if (j<i) then
       c=0
    else
       c=1
    endif

    in=(i-1)/n_bits+1
    jn=(j-1)/n_bits+1

    ibb=modulo(i-1,n_bits)
    jbb=modulo(j-1,n_bits)

    B=ieor(ishft(mask(in),-ibb),ishft(mask(jn),-jbb))

    if (iand(B,1_ib)==c) then
       phase=-1
    else
       phase=1
    endif
    
    return
    
  end function phase_phasemask

!######################################################################
! phase_pure_exc: Computes the phase factor for a pure alpha or beta
!                 single or double excitation.
!######################################################################
  function phase_pure_exc(phase_mask,hlist,plist,maxex,nexci) &
       result(phase)

    use constants
    use bitglobal
    
    implicit none

    integer(is)             :: phase

    integer(is), intent(in) :: maxex,nexci
    integer(ib), intent(in) :: phase_mask(n_int)
    integer(is), intent(in) :: hlist(maxex),plist(maxex)
    integer(is)             :: p,q,r,s
    
    select case(nexci)

    case(1) ! Single alpha or beta excitation

       phase=phase_phasemask(phase_mask,hlist(1),plist(1))

    case(2) ! Double alpha,alpha or beta,beta excitation

       p=hlist(1)
       q=hlist(2)
       r=plist(1)
       s=plist(2)
       
       phase=phase_phasemask(phase_mask,p,r) &
            *phase_phasemask(phase_mask,q,s)
       
       if (max(p,r) > min(q,s)) phase=-phase

    end select
    
    return
    
  end function phase_pure_exc

!######################################################################
! phase_slow: Simple, naive calculation of phase factors based on the
!             sequential application of creation and annihilation
!             operators on the occupation number vector corresponding
!             to the determinant dI and the accumulation of the
!             associated phase factors.
!             This is only intended to be used as a check on the more
!             complicated phase mask-based routines.
!######################################################################
  function phase_slow(dI,hlist,plist,maxex,nexci) result(phase)

    use constants
    use bitglobal
    use detutils
    
    implicit none

    integer(is)             :: phase

    integer(is), intent(in) :: maxex,nexci
    integer(ib), intent(in) :: dI(n_int,2)
    integer(is), intent(in) :: hlist(maxex,2),plist(maxex,2)
    integer(ib)             :: onv(n_int,2)
    integer(is)             :: nalpha,ia,ic
    integer(is)             :: k,i
    integer(is)             :: k1,i1,nbset
    
!----------------------------------------------------------------------
! Make a copy of the input determinant to operate on.
! This will be our representation of the ONV.
!----------------------------------------------------------------------
    onv=dI

!----------------------------------------------------------------------
! Operate on the ONV with the creation and annihilation operators in
! the order: (1) beta then alpha spin operators;
!            (2) within each spin block, annihilator with the highest
!                index, then the creator with the highest index, then
!                the annihilator with the lowest index, then the
!                creator with the lowest index.
!----------------------------------------------------------------------
    !
    ! Initialisation
    !
    phase=1

    !
    ! Number of alpha electrons
    !
    nalpha=n_alpha_electrons(dI)
    
    !
    ! Alpha spin operators
    !
    if (hlist(2,1) /= 0) then
       
       !
       ! Highest-index alpha annihilator
       !
       ia=hlist(2,1)
       ! Block and orbital indices corresponding to the annihilator
       k=(ia-1)/n_bits+1
       i=ia-1-(k-1)*n_bits
       ! Number of bits set before the orbital being annihilated
       nbset=0
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,1),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,1),i1)) nbset=nbset+1
       enddo
       ! Clear the bit corresponding to the annihilator
       onv(k,1)=ibclr(onv(k,1),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset

       !
       ! Highest-index alpha creator
       !
       ic=plist(2,1)
       ! Block and orbital indices corresponding to the creator
       k=(ic-1)/n_bits+1
       i=ic-1-(k-1)*n_bits
       ! Number of bits set before the orbital being created
       nbset=0
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,1),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,1),i1)) nbset=nbset+1
       enddo
       ! Set the bit corresponding to the creator          
       onv(k,1)=ibset(onv(k,1),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset
       
    endif

    if (hlist(1,1) /= 0) then

       !
       ! Lowet-index alpha annihilator
       !
       ia=hlist(1,1)
       ! Block and orbital indices corresponding to the annihilator
       k=(ia-1)/n_bits+1
       i=ia-1-(k-1)*n_bits
       ! Number of bits set before the orbital being annihilated
       nbset=0
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,1),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,1),i1)) nbset=nbset+1
       enddo
       ! Clear the bit corresponding to the annihilator
       onv(k,1)=ibclr(onv(k,1),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset
       
       !
       ! Lowest-index alpha creator
       !
       ic=plist(1,1)
       ! Block and orbital indices corresponding to the creator
       k=(ic-1)/n_bits+1
       i=ic-1-(k-1)*n_bits
       ! Number of bits set before the orbital being created
       nbset=0
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,1),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,1),i1)) nbset=nbset+1
       enddo
       ! Set the bit corresponding to the creator
       onv(k,1)=ibset(onv(k,1),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset
       
    endif

    !
    ! Beta spin operators
    !
    if (hlist(2,2) /= 0) then
       
       !
       ! Highest-index beta annihilator
       !
       ia=hlist(2,2)
       ! Block and orbital indices corresponding to the annihilator
       k=(ia-1)/n_bits+1
       i=ia-1-(k-1)*n_bits
       ! Number of bits set before the orbital being annihilated
       nbset=nalpha
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,2),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,2),i1)) nbset=nbset+1
       enddo
       ! Clear the bit corresponding to the annihilator
       onv(k,2)=ibclr(onv(k,2),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset

       !
       ! Highest-index beta creator
       !
       ic=plist(2,2)
       ! Block and orbital indices corresponding to the creator
       k=(ic-1)/n_bits+1
       i=ic-1-(k-1)*n_bits
       ! Number of bits set before the orbital being created
       nbset=nalpha
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,2),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,2),i1)) nbset=nbset+1
       enddo
       ! Set the bit corresponding to the creator          
       onv(k,2)=ibset(onv(k,2),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset
       
    endif

    if (hlist(1,2) /= 0) then

       !
       ! Lowet-index beta annihilator
       !
       ia=hlist(1,2)
       ! Block and orbital indices corresponding to the annihilator
       k=(ia-1)/n_bits+1
       i=ia-1-(k-1)*n_bits
       ! Number of bits set before the orbital being annihilated
       nbset=nalpha
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,2),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,2),i1)) nbset=nbset+1
       enddo
       ! Clear the bit corresponding to the annihilator
       onv(k,2)=ibclr(onv(k,2),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset
       
       !
       ! Lowest-index beta creator
       !
       ic=plist(1,2)
       ! Block and orbital indices corresponding to the creator
       k=(ic-1)/n_bits+1
       i=ic-1-(k-1)*n_bits
       ! Number of bits set before the orbital being created
       nbset=nalpha
       do k1=1,k-1
          do i1=0,n_bits-1
             if (btest(onv(k1,2),i1)) nbset=nbset+1
          enddo
       enddo
       do i1=0,i-1
          if (btest(onv(k,2),i1)) nbset=nbset+1
       enddo
       ! Set the bit corresponding to the creator
       onv(k,2)=ibset(onv(k,2),i)
       ! Accumulate the phase factor
       phase=phase*(-1)**nbset
          
    endif
    
    return
    
  end function phase_slow
  
!######################################################################
! occ_orbs: determines the indices of the occupied alpha and beta
!           spin orbitals in the determinant d
!######################################################################
  subroutine occ_orbs(d,occlist)

    use constants
    use bitglobal
    use bitutils
    
    implicit none

    integer(ib), intent(in)  :: d(n_int,2)
    integer(is), intent(out) :: occlist(nmo,2)
    integer(ib)              :: I(n_int)
    integer(is)              :: k
    
    !
    ! Get the list of occupied alpha and beta spin orbital indices
    !
    I=d(:,1)
    call list_from_bitstring(I,occlist(:,1),nmo)
    I=d(:,2)
    call list_from_bitstring(I,occlist(:,2),nmo)
    
    return
    
  end subroutine occ_orbs
  
!######################################################################
! hii: Computes a single on-diagonal element of the Hamiltonian
!      matrix
!######################################################################
  function hii(occlist)

    use constants
    use bitglobal

    implicit none

    real(dp)                :: hii

    integer(is), intent(in) :: occlist(nmo,2)
    integer(is)             :: i,j,i1,j1

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hii=0.0d0
    
!----------------------------------------------------------------------
! Core Hamiltonian contributions
!----------------------------------------------------------------------
    do i=1,nel_alpha
       i1=occlist(i,1)
       hii=hii + bitci_ints%h_1e(i1,i1)
    enddo

    do i=1,nel_beta
       i1=occlist(i,2)
       hii=hii + bitci_ints%h_1e(i1,i1)
    enddo
    
!----------------------------------------------------------------------
! Coulomb integral contributions
!----------------------------------------------------------------------
    do i=1,nel_alpha
       i1=occlist(i,1)       
       do j=1,nel_alpha
          j1=occlist(j,1)
          hii=hii+0.5d0*bitci_ints%mo_int(i1,i1,j1,j1)
       enddo
    enddo

    do i=1,nel_beta
       i1=occlist(i,2)       
       do j=1,nel_beta
          j1=occlist(j,2)
          hii=hii+0.5d0*bitci_ints%mo_int(i1,i1,j1,j1)
       enddo
    enddo
    
    do i=1,nel_alpha
       i1=occlist(i,1)       
       do j=1,nel_beta
          j1=occlist(j,2)
          hii=hii+0.5d0*bitci_ints%mo_int(i1,i1,j1,j1)
       enddo
    enddo

    do i=1,nel_beta
       i1=occlist(i,2)       
       do j=1,nel_alpha
          j1=occlist(j,1)
          hii=hii+0.5d0*bitci_ints%mo_int(i1,i1,j1,j1)
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange integral contributions
!----------------------------------------------------------------------
    do i=1,nel_alpha
       i1=occlist(i,1)       
       do j=1,nel_alpha
          j1=occlist(j,1)
          hii=hii-0.5d0*bitci_ints%mo_int(i1,j1,i1,j1)
       enddo
    enddo

    do i=1,nel_beta
       i1=occlist(i,2)       
       do j=1,nel_beta
          j1=occlist(j,2)
          hii=hii-0.5d0*bitci_ints%mo_int(i1,j1,i1,j1)
       enddo
    enddo

!----------------------------------------------------------------------
! Nuclear repulsion energy
!----------------------------------------------------------------------
    hii=hii+enuc
    
    return
    
  end function hii

!######################################################################
! hij_single: Computes a single off-diagonal Hamiltonian matrix
!             element between two determinants linked by a single
!             excitation of spin ispin
!######################################################################
  function hij_single(ispin,phase,occlist,hlist,plist,maxex) &
       result(hij)

    use constants
    use bitglobal
    
    implicit none

    real(dp)                :: hij

    integer(is), intent(in) :: ispin,phase,maxex
    integer(is), intent(in) :: occlist(nmo,2)
    integer(is), intent(in) :: plist(maxex,2),hlist(maxex,2)
    integer(is)             :: ia,ic
    integer(is)             :: jspin
    integer(is)             :: i,i1
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hij=0.0d0

!----------------------------------------------------------------------
! Indices of the annihilation and creation operators
!----------------------------------------------------------------------
    ia=hlist(1,ispin)
    ic=plist(1,ispin)
    
!----------------------------------------------------------------------
! Core Hamiltonian contribution
!----------------------------------------------------------------------
    hij=hij + bitci_ints%h_1e(ia,ic)

!----------------------------------------------------------------------
! Coulomb integral contributions
!----------------------------------------------------------------------
    do jspin=1,2
       do i=1,nel_spin(jspin)
          i1=occlist(i,jspin)
          hij=hij+bitci_ints%mo_int(ia,ic,i1,i1)
       enddo
    enddo

!----------------------------------------------------------------------
! Exchange integral contributions
!----------------------------------------------------------------------
    do i=1,nel_spin(ispin)
       i1=occlist(i,ispin)
       hij=hij-bitci_ints%mo_int(ia,i1,ic,i1)
    enddo

!----------------------------------------------------------------------
! Phase factor
!----------------------------------------------------------------------
    hij=hij*phase
    
    return
    
  end function hij_single

!######################################################################
! hij_double_same_spin: Computes a single off-diagonal Hamiltonian
!                       matrix element between two determinants linked
!                       double alpha,alpha or beta,beta excitation
!######################################################################
  function hij_double_same_spin(ispin,phase,occlist,hlist,plist,&
       maxex) result(hij)

    use constants
    use bitglobal
    
    implicit none

    real(dp)                :: hij

    integer(is), intent(in) :: ispin,phase,maxex
    integer(is), intent(in) :: occlist(nmo,2)
    integer(is), intent(in) :: plist(maxex,2),hlist(maxex,2)
    integer(is)             :: ia1,ic1,ia2,ic2
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hij=0.0d0

!----------------------------------------------------------------------
! Indices of the annihilation and creation operators
!----------------------------------------------------------------------
    ia1=hlist(1,ispin)
    ia2=hlist(2,ispin)
    ic1=plist(1,ispin)
    ic2=plist(2,ispin)
    
!----------------------------------------------------------------------
! Coulomb integral contribution
!----------------------------------------------------------------------
    hij=hij+bitci_ints%mo_int(ia1,ic1,ia2,ic2)

!----------------------------------------------------------------------
! Exchange integral contribution
!----------------------------------------------------------------------
    hij=hij-bitci_ints%mo_int(ia1,ic2,ia2,ic1)
    
!----------------------------------------------------------------------
! Phase factor
!----------------------------------------------------------------------
    hij=hij*phase

    return
    
  end function hij_double_same_spin

!######################################################################
! hij_double_diff_spin: Computes a single off-diagonal Hamiltonian
!                       matrix element between two determinants linked
!                       double alpha,beta excitation
!######################################################################
  function hij_double_diff_spin(phase,occlist,hlist,plist,maxex) &
       result(hij)

    use constants
    use bitglobal
    
    implicit none

    real(dp)                :: hij

    integer(is), intent(in) :: phase,maxex
    integer(is), intent(in) :: occlist(nmo,2)
    integer(is), intent(in) :: plist(maxex,2),hlist(maxex,2)
    integer(is)             :: ia1,ic1,ia2,ic2

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hij=0.0d0

!----------------------------------------------------------------------
! Indices of the annihilation and creation operators
!----------------------------------------------------------------------
    ia1=hlist(1,1)
    ia2=hlist(1,2)
    ic1=plist(1,1)
    ic2=plist(1,2)

!----------------------------------------------------------------------
! Coulomb integral contribution
!----------------------------------------------------------------------
    hij=hij+bitci_ints%mo_int(ia1,ic1,ia2,ic2)
    
!----------------------------------------------------------------------
! Phase factor
!----------------------------------------------------------------------
    hij=hij*phase
    
    return
    
  end function hij_double_diff_spin
  
!######################################################################
! s2ii: Computes a single on-diagonal element of the S2 matrix
!######################################################################
  function s2ii(d)

    use constants
    use bitglobal
    
    implicit none

    real(dp)                :: s2ii
    
    integer(ib), intent(in) :: d(n_int,2)
    integer(ib)             :: sop1(n_int)
    integer(is)             :: na,nb,k

    !
    ! Get the indices of the singly-occupied spatial orbitals
    !
    do k=1,n_int
       sop1(k)=ieor(d(k,1),d(k,2))
    enddo

    !
    ! Determine the number of singly-occupied alpha and bets
    ! spin orbitals
    !
    na=0
    nb=0
    do k=1,n_int
       na=na+popcnt(iand(sop1(k),d(k,1)))
       nb=nb+popcnt(iand(sop1(k),d(k,2)))
    enddo

    !
    ! Value of the on-diagonal element of the S2 matrix
    !
    s2ii=0.25d0*(na-nb)**2 + 0.5d0*(na+nb)

    return
    
  end function s2ii

!######################################################################
! s2ij: Computes a single off-diagonal element of the S2 matrix.
!######################################################################
  function s2ij(phase,hlist,plist)

    use constants
    use bitglobal
    
    implicit none

    real(dp)                :: s2ij
    
    integer(is), intent(in) :: phase
    integer(is), intent(in) :: plist(2,2),hlist(2,2)
    
    if (hlist(1,1) == plist(1,2) &
         .and. plist(1,1) == hlist(1,2)) then
       s2ij=-phase
    else
       s2ij=0.0d0
    endif
    
    return
    
  end function s2ij

!######################################################################
  
end module slater_condon
