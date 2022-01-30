!**********************************************************************
! Routines for the pre-computation of the unique spin-coupling
! coefficients <w omega| T_pq^(1,k) |w' omega'> for triplet spin tensor
! operators T_pq^(1,k), k=+1
!**********************************************************************
! Follows the formalism detailed in
!
! R. W. Wetmore and G. A. Segal, Chem. Phys. Lett., 36, 478 (1975)
!
! and
!
! M. Kleinschmidt and C. M. Marian, Chem. Phys, 311, 71 (2005)
!**********************************************************************
! Imporant: We are assuming that the spin-coupling coefficients will
!           only be used in the evaluation of matrix elements over
!           1-electron operators. i.e., the Case 2 spin-coupling
!           coefficients between nomax and nomax+2 open shells are not
!           currently computed.
!**********************************************************************
module spin_coupling_k1

  implicit none

  private
  
  public :: scc_k1
  
contains
  
!######################################################################
! scc_k1: Generates the unique k=+1 component triplet spin-coupling
!         coefficients for a given pair of bra and ket spin
!         multiplicities
!######################################################################
  subroutine scc_k1(imultB,imultK,nocase1,nocase2,maxcsfB,maxdetB,&
       ncsfsB,ndetsB,csfcoeB,detvecB,maxcsfK,maxdetK,ncsfsK,ndetsK,&
       csfcoeK,detvecK,nspincp,N1s,verbose,spincp,patternmap,offspincp)

    use constants
    use iomod
    use timing

    ! Bra and ket spin multiplicities
    integer(is), intent(in) :: imultB,imultK

    ! Maximum number open shells for the Case 1 and Case 2 CSFs
    integer(is), intent(in)  :: nocase1,nocase2

    ! Numbers of bra and ket CSFs and determinants as a function
    ! of the number of open shells
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ndetsB(0:nocase2)
    integer(is), intent(in)  :: ncsfsK(0:nocase2),ndetsK(0:nocase2)
    
    ! Maximum number of bra and ket CSFs/determinants across all
    ! numbers of open shells
    integer(is), intent(in)  :: maxcsfB,maxdetB
    integer(is), intent(in)  :: maxcsfK,maxdetK
    
    ! Bra and ket CSF expansion coefficients
    real(dp), intent(in)     :: csfcoeB(maxcsfB,maxdetB,nocase2)
    real(dp), intent(in)     :: csfcoeK(maxcsfK,maxdetK,nocase2)
    
    ! Bit string encoding of the determinants contributing to the
    ! bra and ket CSFs
    integer(ib), intent(in)  :: detvecB(maxdetB,nocase2)
    integer(ib), intent(in)  :: detvecK(maxdetK,nocase2)
    
    ! Number of unique spin coupling coefficients
    integer(is), intent(out) :: nspincp(2)

    ! Bit strings comprised of N 1's
    integer(ib), allocatable :: N1s(:)

    ! Verbose output
    logical, intent(in)      :: verbose

    ! All spin coupling coefficients
    integer(is)              :: spincpdim(3)
    real(dp), allocatable    :: spincp(:)

    ! All pattern -> array index mappings
    integer(is)              :: mapdim
    integer(is), allocatable :: patternmap(:)

    ! Spin coupling coefficient offsets for the various
    ! different cases
    integer(is), allocatable :: offspincp(:)

    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,twall_end
    
    ! Everything else
    integer(is)              :: i
    real(dp)                 :: SB,SK
    character(len=10)        :: ak
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    if (verbose) then
       write(6,'(/,52a)') ('-',i=1,52)
       write(6,'(x,a)') 'Triplet spin-coupling coefficient generation'
       write(6,'(52a)') ('-',i=1,52)
    endif

!----------------------------------------------------------------------
! Check on the bra and ket spin multiplicities
!----------------------------------------------------------------------
    ! k=+1 only makes sense if the bra and ket total spins differ
    ! by one
    SB=0.5d0*(dble(imultB)-1.0d0)
    SK=0.5d0*(dble(imultK)-1.0d0)
    if (abs(SB-SK) /= 1.0d0) then
       errmsg='|S_bra - S_ket| != 0: k=+1 makes no sense'
       call error_control
    endif

!----------------------------------------------------------------------
! Compute the number of unique spin coupling coefficients for the
! spin multiplicity under consideration
!----------------------------------------------------------------------
    call get_nunique_k1(imultB,imultK,nocase1,nocase2,ncsfsB,ncsfsK,&
         nspincp)
    
!----------------------------------------------------------------------
! Output some information about what we are doing
!----------------------------------------------------------------------
    if (verbose) call print_spincp_info_k1(imultB,imultK,nspincp)

!----------------------------------------------------------------------
! Allocate the spin coupling coefficient arrays
!----------------------------------------------------------------------
    call init_spincp_k1(imultB,imultK,nocase1,nocase2,ncsfsB,ncsfsK,&
         nspincp,verbose,spincp,spincpdim,offspincp)

!----------------------------------------------------------------------
! Allocate the pattern value -> array index mapping array
!----------------------------------------------------------------------
    call init_patternmap_k1(imultB,imultB,nocase1,nocase2,ncsfsB,&
         ncsfsK,patternmap,mapdim)

!----------------------------------------------------------------------
! Compute the Case 1 spin coupling coefficients
!----------------------------------------------------------------------
    call case1_coeffs_k1(nocase1,nocase2,maxcsfB,maxcsfK,maxdetB,&
         maxdetK,ncsfsB,ncsfsK,ndetsB,ndetsK,csfcoeB,csfcoeK,detvecB,&
         detvecK,spincpdim,spincp,nspincp,mapdim,patternmap)
    
    STOP

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) call report_times(twall_end-twall_start,&
         tcpu_end-tcpu_start,'generate_coupling_coefficients_triplet')

!----------------------------------------------------------------------    
! Flush stdout
!----------------------------------------------------------------------
    flush(6)
    
    return
    
  end subroutine scc_k1

!######################################################################
! get_nunique_k1: Computes the total number of *unique* spin
!                 coupling coefficients as a function of the bra
!                 and ket spin multiplicities and the maximum number
!                 of open shells
!######################################################################
  subroutine get_nunique_k1(imultB,imultK,nocase1,nocase2,ncsfsB,&
       ncsfsK,nspincp)

    use constants
        
    implicit none

    integer(is), intent(in)  :: imultB,imultK,nocase1,nocase2
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ncsfsK(0:nocase2)
    integer(is), intent(out) :: nspincp(2)
    
    integer(is)              :: nopen,n,npat

!----------------------------------------------------------------------
! Case 1: N_Bra = N_Ket
!----------------------------------------------------------------------
! This covers excitations from a singly-occupied orbital into an
! unoccupied orbital and from a doubly-occupied orbital into a
! singly-occupied orbital
!----------------------------------------------------------------------
    nspincp(1)=0

    ! Loop over numbers of open shells
    do nopen=1,nocase1
       
       ! Number of patterns
       npat=(nopen+1)*nopen/2

       ! Number of unique spin-coupling coefficients corresponding to
       ! this pattern
       nspincp(1)=nspincp(1)+2*npat*ncsfsB(nopen)*ncsfsK(nopen)
       
    enddo

!----------------------------------------------------------------------
! Case 2: N_Bra = N_Ket +/- 2
!----------------------------------------------------------------------
! This covers excitations from a doubly-occupied orbital into an
! unoccupied orbital or from a singly-occupied orbital into another
! singly-occupied orbital
!---------------------------------------------------------------------- 
    nspincp(2)=0
    
    ! Loop over numbers of open shells
    do nopen=2,nocase1

       ! Number of patterns
       npat=(nopen)*(nopen-1)/2
       
       ! Number of unique spin-coupling coefficients corresponding to
       ! this pattern
       ! Case 2a
       nspincp(2)=nspincp(2)+npat*ncsfsB(nopen)*ncsfsK(nopen-2)
       ! Case 2b
       nspincp(2)=nspincp(2)+npat*ncsfsB(nopen-2)*ncsfsK(nopen)
              
    enddo
    
    return
    
  end subroutine get_nunique_k1

!######################################################################
! print_spincp_info_k1: Printing of some information about the
!                       k=+1 component triplet spin coupling
!                       coefficients being calculated
!######################################################################
  subroutine print_spincp_info_k1(imultB,imultK,nspincp)
    
    use constants
        
    implicit none

    integer(is), intent(in) :: imultB,imultK
    integer(is), intent(in) :: nspincp(2)
    
    integer(is)             :: n,i
    
!----------------------------------------------------------------------
! Triplet spin tensor component
!----------------------------------------------------------------------
    write(6,'(/,x,a)') 'Spin tensor component: +1'
    
!----------------------------------------------------------------------
! Spin multiplicities
!----------------------------------------------------------------------
    write(6,'(/,x,a,x,i0)') 'Bra spin multiplicity:',imultB
    write(6,'(x,a,x,i0)') 'Ket spin multiplicity:',imultK
    
!----------------------------------------------------------------------
! Number of unique spin-coupling coefficients
!----------------------------------------------------------------------
    ! Total number of spin-coupling coefficients
    write(6,'(/,x,a,2x,i0)') &
         'Total number of spin coupling coefficients:',&
         sum(nspincp)
    
    ! Number of Case 1 spin-coupling coefficients
    write(6,'(x,a,x,i0)') &
         'Number of Case 1 spin coupling coefficients:',nspincp(1)
    
    ! Number of Case 2 spin-coupling coefficients
    write(6,'(x,a,x,i0)') &
         'Number of Case 2 spin coupling coefficients:',nspincp(2)
        
    return
    
  end subroutine print_spincp_info_k1
    
!######################################################################
! init_spincp_k1: Allocation and initialisation of the various spin
!                 coupling coefficient arrays for the k=+1 component
!                 of the triplet spin tensor operator
!######################################################################
  subroutine init_spincp_k1(imultB,imultK,nocase1,nocase2,ncsfsB,&
       ncsfsK,nspincp,verbose,spincp,spincpdim,offspincp)

    use constants
    
    implicit none

    integer(is), intent(in)  :: imultB,imultK,nocase1,nocase2
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ncsfsK(0:nocase2)
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(out) :: spincpdim(3)
    integer(is), allocatable :: offspincp(:)
    real(dp), allocatable    :: spincp(:)
    logical, intent(in)      :: verbose
        
    integer(is)              :: nopen,npat,n1,n2a,n2b
    real(dp)                 :: mem

!----------------------------------------------------------------------
! Numbers of spin coupling coefficients to be stored in the spincp
! array
!----------------------------------------------------------------------
    ! Case 1 coefficients
    ! We will explicitly store: (i) 1a i>j, (ii) 1a i<j,
    !                         (iii) 1b i>j, (iv) 1b i<j
    spincpdim(1)=2*nspincp(1)
       
    ! Case 1 coefficients
    ! We will explicitly store: (i) 2a i>j, (ii) 2a i<j,
    !                         (iii) 2b i>j, (iv) 2b i<j
    spincpdim(2)=2*nspincp(2)
    
    ! Total number of coefficients
    spincpdim(3)=sum(spincpdim(1:2))

!----------------------------------------------------------------------
! Allocate the spin coupling coefficient and offset arrays
!----------------------------------------------------------------------    
    ! Spin coupling coefficients
    allocate(spincp(spincpdim(3)))
    spincp=0.0d0

    ! Offsets
    allocate(offspincp(6))
    offspincp=0
    
!----------------------------------------------------------------------
! Number of unique Case 1a, 1b, 2a and 2b coefficients
!----------------------------------------------------------------------
    ! Case 1a and 1b
    n1=nspincp(1)/2

    ! Case 2a and 2b
    n2a=0
    n2b=0
    ! Loop over numbers of open shells
    do nopen=2,nocase1
       ! Number of patterns
       npat=(nopen)*(nopen-1)/2
       ! Case 2a: N_bra = N_ket + 2
       n2a=n2a+npat*ncsfsB(nopen)*ncsfsK(nopen-2)
       ! Case 2b: N_bra = N_ket - 2
       n2b=n2b+npat*ncsfsB(nopen-2)*ncsfsK(nopen)
    enddo
       
!----------------------------------------------------------------------
! Offsets for the various different spin coupling coefficients within
! the spincp array
!----------------------------------------------------------------------
! The spin coupling coefficients will be stored in the order
! ([1a i>j], [1a i<j], [1b i>j], [1b i<j], [2a i>j], [2a i<j],
!  [2b i>j], [2b i<j])
!----------------------------------------------------------------------
! The offsets correspond to the shifts that will be needed to be added
! to the pattern indices to reach each block of coefficients in the
! spincp array
!----------------------------------------------------------------------
    ! Case 1a, i<j
    offspincp(1)=n1
    
    ! Case 1b, i>j
    offspincp(2)=2*n1
    
    ! Case 1b, i<j
    offspincp(3)=3*n1
    
    ! Case 2a, i<j
    offspincp(4)=n2a
    
    ! Case 2b, i>j
    offspincp(5)=2*n2a
    
    ! Case 2b, i<j
    offspincp(6)=2*n2a+n2b

!----------------------------------------------------------------------
! Output the memory used
!----------------------------------------------------------------------
    mem=spincpdim(3)*8/1024.0d0**2
    
    if (verbose) write(6,'(/,x,a,F8.2,x,a)') 'Memory used:',&
         mem,'MB'
    
    return
    
  end subroutine init_spincp_k1

!######################################################################
! init_patternmap_k1: Allocation and initialisation of the array
!                     holding the pattern value -> array index mapping
!######################################################################
  subroutine init_patternmap_k1(imultB,imultK,nocase1,nocase2,ncsfsB,&
       ncsfsK,patternmap,mapdim)

    use constants
    
    implicit none

    integer(is), intent(in)  :: imultB,imultK
    integer(is), intent(in)  :: nocase1,nocase2
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ncsfsK(0:nocase2)
    integer(is), intent(out) :: mapdim
    integer(is), allocatable :: patternmap(:)

!----------------------------------------------------------------------
! Maximum possible pattern value
!----------------------------------------------------------------------
    ! Case 1 pattern values
    if (ncsfsK(nocase1) /= 0) then
       mapdim=2**(nocase1+1)-4
    else
       mapdim=2**nocase1-4
    endif

    ! Case 2 pattern values
    if (modulo(imultK,2) == 0) then
       mapdim=max(mapdim,2**(nocase1-1)-4)
    else
       mapdim=max(mapdim,2**nocase1-4)
    endif

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(patternmap(0:mapdim))
    patternmap=0
    
    return
    
  end subroutine init_patternmap_k1
    
!######################################################################
! case1_coeffs: Calculation of the Case 1 spin coupling coefficients.
!               These are the coefficients corresponding to pairs of
!               initial and final spatial occupations with equal
!               numbers of open shells
!######################################################################
  subroutine case1_coeffs_k1(nocase1,nocase2,maxcsfB,maxcsfK,maxdetB,&
       maxdetK,ncsfsB,ncsfsK,ndetsB,ndetsK,csfcoeB,csfcoeK,detvecB,&
       detvecK,spincpdim,spincp,nspincp,mapdim,patternmap)

    use constants
    use bitutils
        
    implicit none

    integer(is), intent(in)  :: nocase1,nocase2
    integer(is), intent(in)  :: maxcsfB,maxcsfK,maxdetB,maxdetK
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ncsfsK(0:nocase2)
    integer(is), intent(in)  :: ndetsB(0:nocase2),ndetsK(0:nocase2)
    real(dp), intent(in)     :: csfcoeB(maxcsfB,maxdetB,nocase2)
    real(dp), intent(in)     :: csfcoeK(maxcsfK,maxdetK,nocase2)
    integer(ib), intent(in)  :: detvecB(maxdetB,nocase2)
    integer(ib), intent(in)  :: detvecK(maxdetB,nocase2)
    integer(is), intent(in)  :: spincpdim(3)
    real(dp), intent(out)    :: spincp(spincpdim(3))
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(in)  :: mapdim
    integer(is), intent(out) :: patternmap(0:mapdim)
    
    integer(is)              :: nopen,is1,is2,icsf1,icsf2
    integer(is)              :: nspB,nspK,ndetB,ndetK
    integer(is)              :: ioff,istart,iend
    integer(ib), allocatable :: ws(:,:)
    integer(ib)              :: pattern
    real(dp)                 :: coeff
    real(dp), allocatable    :: work(:),workT(:),work2(:,:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Simplified spatial occupation vectors
    allocate(ws(nocase1+1,nocase1))
    ws=0_ib

!----------------------------------------------------------------------
! Generate the unique Case 1 simplified spatial occupation vectors
! i.e, all possible permutations of N 1's and one 0 for
! N=1,...,nocase1
!----------------------------------------------------------------------
    ! Loop over number of open shells
    do nopen=1,nocase1

       ! Generate the permutations
       call get_permutations(nopen,1,ws(1:nopen+1,nopen),nopen+1)
       
    enddo

!----------------------------------------------------------------------
! Compute the Case 1 a triplet spin coupling coefficients
! < w' omega' | T_ij^(1,k=+1) | w omega > for i>j
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=1

    ! Loop over numbers of open shells
    do nopen=1,nocase1

       ! No. bra and ket CSFs and determinants
       nspB=ncsfsB(nopen)
       ndetB=ndetsB(nopen)
       nspK=ncsfsK(nopen)
       ndetK=ndetsK(nopen)

       ! Cycle if there are no CSFs for this number of open shells
       if (nspB == 0 .or. nspK==0) cycle

       ! Allocate the spin-coupling coefficient work array
       allocate(work(nspB*nspK))

       ! Loop over the unique pairs of simplified spatial occupation
       ! vectors
       do is1=1,nopen
          do is2=is1+1,nopen+1

             ! Pattern number and the corresponding array index
             pattern=iand(ws(is1,nopen),ws(is2,nopen))
             patternmap(pattern)=ioff

             ! Evaluate the Case 1a i>j spin coupling coefficients
             ! for all bra and ket CSF pairs
             call spincp_1a_igtj(ws(is1,nopen),ws(is2,nopen),&
                  nopen,nocase2,maxcsfB,maxcsfK,maxdetB,maxdetK,nspB,&
                  nspK,ndetB,ndetK,detvecB,detvecK,csfcoeB,csfcoeK,&
                  work)
             
          enddo
       enddo
          
       ! Deallocate the spin-coupling coefficient work array
       deallocate(work)
       
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ws)
    
    return

  end subroutine case1_coeffs_k1

!######################################################################
! spincp_1a_igtj: Computes a single batch of Case 1a i>j spin coupling
!                 coefficients for the k=+1 component of the triplet
!                 spin tensor operator
!######################################################################
  subroutine spincp_1a_igtj(ws1,ws2,nopen,nocase2,maxcsfB,maxcsfK,&
       maxdetB,maxdetK,nspB,nspK,ndetB,ndetK,detvecB,detvecK,csfcoeB,&
       csfcoeK,spincoe)

    use constants
    use bitutils
    use slater_condon
    
    implicit none

    integer(ib), intent(in)  :: ws1,ws2
    integer(is), intent(in)  :: nopen
    integer(is), intent(in)  :: nocase2
    integer(is), intent(in)  :: maxcsfB,maxcsfK,maxdetB,maxdetK
    integer(is), intent(in)  :: nspB,nspK,ndetB,ndetK
    integer(ib), intent(in)  :: detvecB(maxdetB,nocase2)
    integer(ib), intent(in)  :: detvecK(maxdetK,nocase2)
    real(dp), intent(in)     :: csfcoeB(maxcsfB,maxdetB,nocase2)
    real(dp), intent(in)     :: csfcoeK(maxcsfK,maxdetK,nocase2)
    real(dp), intent(out)    :: spincoe(nspB,nspK)
    integer(is)              :: icsf1,icsf2
    integer(is)              :: idet,n,vecindx,iket,ibra,ispin
    integer(is)              :: ic,ia,ih,ip
    integer(ib), allocatable :: d1(:,:),d2(:,:)
    integer(is)              :: n_int_save
    integer(ib)              :: phase_mask(2)
    integer(is)              :: nexci
    integer(ib)              :: p(2),h(2)
    integer(is), parameter   :: maxex=2
    integer(is)              :: plist(maxex,2),hlist(maxex,2)
    integer(is)              :: phase
    real(dp), allocatable    :: phasemat(:,:),pcT(:,:)
    real(dp), allocatable    :: coeB(:,:),coeK(:,:)

!----------------------------------------------------------------------
! Save the actual value of n_int and then set this to 1 for use in the
! following. This allows us to use the slater_condon module to
! calculate the spin coupling coefficients.
!----------------------------------------------------------------------
    n_int_save=n_int
    n_int=1

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Determinant bit strings
    allocate(d1(2,ndetK))
    allocate(d2(2,ndetB))
    d1=0_ib
    d2=0_ib

    ! Phase factor matrix
    allocate(phasemat(ndetB,ndetK))
    phasemat=0.0d0

    ! Phase factor matrix contracted with the transpose of the ket CSF
    ! coefficient matrix
    allocate(pcT(ndetk,nspk))
    pcT=0.0d0
    
    ! Working arrays
    allocate(coeB(nspB,ndetB))
    allocate(coeK(nspK,ndetK))
    coeB=0.0d0
    coeK=0.0d0
    
!----------------------------------------------------------------------
! Generate the determinants for CSF 1 (the ket CSF)
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    d1=0_ib

    ! Initialise the detvec open shell orbital counter
    vecindx=0
    
    ! Loop over orbitals in the spatial occupation
    do n=1,nopen+1

       ! Cycle if we are at an unset bit
       if (.not. btest(ws1,n-1)) cycle

       ! Increment the detvec open shell orbital counter
       vecindx=vecindx+1

       ! Loop over determinants
       do idet=1,ndetK
          
          ! Add the next spin orbital
          if (btest(detvecK(idet,nopen),vecindx-1)) then
             ! Occupied alpha spin-orbital
             d1(1,idet)=ibset(d1(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             d1(2,idet)=ibset(d1(2,idet),n-1)
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Generate the determinants for CSF 2 (the bra CSF)
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    d2=0_ib

    ! Initialise the detvec open shell orbital counter
    vecindx=0
    
    ! Loop over orbitals in the spatial occupation
    do n=1,nopen+1

       ! Cycle if we are at an unset bit
       if (.not. btest(ws2,n-1)) cycle

       ! Increment the detvec open shell orbital counter
       vecindx=vecindx+1
       
       ! Loop over determinants
       do idet=1,ndetB
          
          ! Add the next spin orbital
          if (btest(detvecB(idet,nopen),vecindx-1)) then
             ! Occupied alpha spin-orbital
             d2(1,idet)=ibset(d2(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             d2(2,idet)=ibset(d2(2,idet),n-1)
          endif
          
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Determine the indices of the triple excitation operator T_ij^(1,k=+1)
! Remember that these are given by the positions of the 0's in
! ws1 and ws2, i.e., the indices of the first unset bits
!----------------------------------------------------------------------
    ! Creation operator index
    ic=trailz(not(ws1))+1

    ! Annihilation operator index
    ia=trailz(not(ws2))+1

!----------------------------------------------------------------------
! Compute the matrix of phase factors for all pairs of determinants
! for the spin-flip excitation operator a_ic,alpha^\dagger a_ia,beta
!----------------------------------------------------------------------
! phasemat(I,J) = < det_I | a_ic,alpha^\dagger a_ia,beta | det_J >
!----------------------------------------------------------------------
    ! Loop over ket determinants
    do iket=1,ndetK

       ! Loop over bra determinants
       do ibra=1,ndetB

          ! Get the excitation degree
          nexci=exc_degree_det(d2(:,ibra),d1(:,iket))

          ! Cycle if the excitation degree is not equal to 1
          if (nexci /= 1) cycle

          ! Get the indices of the spin-orbitals involved in the
          ! excitations linking the bra and ket determinants
          call exc(d1(:,iket),d2(:,ibra),p,h)
          call list_from_bitstring(p(ialpha),plist(:,ialpha),maxex)
          call list_from_bitstring(h(ialpha),hlist(:,ialpha),maxex)
          call list_from_bitstring(p(ibeta),plist(:,ibeta),maxex)
          call list_from_bitstring(h(ibeta),hlist(:,ibeta),maxex)

          ! Make sure that the hole is in a beta spin orbital
          if (hlist(1,2) == 0) cycle

          ! Cycle if the indices do not match those of the creation
          ! and annihilation operators
          ip=plist(1,1)
          ih=hlist(1,2)
          if (ip /= ic .or. ih /= ia) cycle

          ! Phase factor
          phase=phase_slow(d1(:,iket),hlist,plist,maxex,nexci)
          phasemat(ibra,iket)=dble(phase)
          
       enddo

    enddo

!----------------------------------------------------------------------
! Compute the matrix of spin coupling coefficients for the current
! pattern index
!----------------------------------------------------------------------
    ! Working coefficient arrays
    coeB=csfcoeB(1:nspB,1:ndetB,nopen)
    coeK=csfcoeK(1:nspK,1:ndetK,nopen)
    
    ! Contract the phase matrix with the ket CSF expansion
    ! coefficients: pcT = phasemat * transpose(csfcoeK)
    call dgemm('N','T',ndetB,nspK,ndetK,1.0d0,phasemat,ndetB,coeK,&
         nspK,0.0d0,pcT,ndetB)

    ! Contract the bra CSF expansion coefficients with the intermediate
    ! matrix pcT to yield the spin coupling coefficients
    call dgemm('N','N',nspB,nspK,ndetB,1.0d0,coeB,nspB,pcT,ndetB,&
         0.0d0,spincoe,nspB)

    ! Up to now, we have computed the matrix elements
    ! < CSF_B | a_ialpha^dagger a_jbeta | CSF_K >
    !
    ! Multiply by -1 to get the matrix elements
    ! < CSF_B | T_ij^(1,k=+1) | CSF_K >
    spincoe=-spincoe
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(d1)
    deallocate(d2)
    deallocate(phasemat)
    deallocate(pcT)
    deallocate(coeB)
    deallocate(coeK)
    
!----------------------------------------------------------------------
! Reset n_int
!----------------------------------------------------------------------
    n_int=n_int_save
    
    return
    
  end subroutine spincp_1a_igtj

!######################################################################
  
end module spin_coupling_k1
