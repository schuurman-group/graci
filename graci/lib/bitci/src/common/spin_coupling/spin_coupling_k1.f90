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
       csfcoeK,detvecK,nspincp,N1s,verbose,spincp,patternmap1,&
       patternmap2a,patternmap2b,offspincp)

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
    integer(is)              :: mapdim1,mapdim2
    integer(is), allocatable :: patternmap1(:)
    integer(is), allocatable :: patternmap2a(:)
    integer(is), allocatable :: patternmap2b(:)
    
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
    ! k=+1 only makes sense if the bra minus ket total spins equals 1
    SB=0.5d0*(dble(imultB)-1.0d0)
    SK=0.5d0*(dble(imultK)-1.0d0)
    if (SB-SK /= 1.0d0) then
       errmsg='Error in scc_k1: S_bra - S_ket != 1.0'
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
         ncsfsK,patternmap1,patternmap2a,patternmap2b,mapdim1,&
         mapdim2)

!----------------------------------------------------------------------
! Compute the Case 1 spin coupling coefficients
!----------------------------------------------------------------------
    call case1_coeffs_k1(nocase1,nocase2,maxcsfB,maxcsfK,maxdetB,&
         maxdetK,ncsfsB,ncsfsK,ndetsB,ndetsK,csfcoeB,csfcoeK,detvecB,&
         detvecK,spincpdim,spincp,nspincp,mapdim1,patternmap1,offspincp)

!----------------------------------------------------------------------
! Compute the Case 2 spin coupling coefficients
!----------------------------------------------------------------------
    call case2_coeffs_k1(imultB,imultK,nocase1,nocase2,maxcsfB,&
         maxcsfK,maxdetB,maxdetK,ncsfsB,ncsfsK,ndetsB,ndetsK,csfcoeB,&
         csfcoeK,detvecB,detvecK,spincpdim,spincp,nspincp,mapdim2,&
         patternmap2a,patternmap2b,offspincp)
    
!----------------------------------------------------------------------
! Fill in the array of bit strings with N set bits, from which the
! pattern numbers will be derived by clearing bits
!----------------------------------------------------------------------
    call fill_N1s(nocase2,N1s)

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) call report_times(twall_end-twall_start,&
         tcpu_end-tcpu_start,'scc_k1')
    
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
    
    ! Case 2b, i<j
    offspincp(5)=n2b
    
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
       ncsfsK,patternmap1,patternmap2a,patternmap2b,mapdim1,mapdim2)

    use constants
    
    implicit none

    integer(is), intent(in)  :: imultB,imultK
    integer(is), intent(in)  :: nocase1,nocase2
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ncsfsK(0:nocase2)
    integer(is), intent(out) :: mapdim1,mapdim2
    integer(is), allocatable :: patternmap1(:),patternmap2a(:),&
                                patternmap2b(:)

!----------------------------------------------------------------------
! Maximum possible pattern value
!----------------------------------------------------------------------
    ! Case 1 pattern values
    if (ncsfsK(nocase1) /= 0) then
       mapdim1=2**(nocase1+1)-4
    else
       mapdim1=2**nocase1-4
    endif

    ! Case 2 pattern values
    if (modulo(imultK,2) == 0) then
       mapdim2=2**(nocase1-1)-4
    else
       mapdim2=2**nocase1-4
    endif

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(patternmap1(0:mapdim1))
    patternmap1=0

    allocate(patternmap2a(0:mapdim2))
    patternmap2a=0

    allocate(patternmap2b(0:mapdim2))
    patternmap2b=0
    
    return
    
  end subroutine init_patternmap_k1
    
!######################################################################
! case1_coeffs_k1: Calculation of the Case 1 spin coupling
!                  coefficients. These are the coefficients
!                  corresponding to pairs of initial and final spatial
!                  occupations with equal numbers of open shells
!######################################################################
  subroutine case1_coeffs_k1(nocase1,nocase2,maxcsfB,maxcsfK,maxdetB,&
       maxdetK,ncsfsB,ncsfsK,ndetsB,ndetsK,csfcoeB,csfcoeK,detvecB,&
       detvecK,spincpdim,spincp,nspincp,mapdim1,patternmap1,offspincp)

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
    integer(is), intent(in)  :: mapdim1
    integer(is), intent(out) :: patternmap1(0:mapdim1)
    integer(is), intent(in)  :: offspincp(6)
    
    integer(is)              :: nopen,is1,is2,icsf1,icsf2
    integer(is)              :: nspB,nspK,ndetB,ndetK
    integer(is)              :: ioff,istart,iend
    integer(ib), allocatable :: ws(:,:)
    integer(ib)              :: pattern
    real(dp)                 :: coeff
    real(dp), allocatable    :: work(:),workgt(:),worklt(:)

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
! < w' omega' | T_ij^(1,k=+1) | w omega > for i>j and i<j
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

       ! Loop over the *upper triangle* of unique pairs of simplified
       ! spatial occupation vectors
       do is1=1,nopen
          do is2=is1+1,nopen+1

             ! Pattern number and the corresponding array index
             pattern=iand(ws(is1,nopen),ws(is2,nopen))
             patternmap1(pattern)=ioff

             ! Compute the 1a i>j spin coupling coefficients
             call spincp_1a_k1(ws(is1,nopen),ws(is2,nopen),nopen,&
                  nocase2,maxcsfB,maxcsfK,maxdetB,maxdetK,nspB,nspK,&
                  ndetB,ndetK,detvecB,detvecK,csfcoeB,csfcoeK,work)

             ! Save the 1a i>j spin coupling coefficients
             istart=ioff
             iend=ioff+nspB*nspK-1
             spincp(istart:iend)=work

             ! Compute the 1a i<j spin coupling coefficients
             call spincp_1a_k1(ws(is2,nopen),ws(is1,nopen),nopen,&
                  nocase2,maxcsfB,maxcsfK,maxdetB,maxdetK,nspB,nspK,&
                  ndetB,ndetK,detvecB,detvecK,csfcoeB,csfcoeK,work)

             ! Save the 1a i<j spin coupling coefficients
             istart=istart+offspincp(1)
             iend=iend+offspincp(1)
             spincp(istart:iend)=work
             
             ! Update the spincp offset
             ioff=ioff+nspB*nspK
             
          enddo
       enddo
          
       ! Deallocate the spin-coupling coefficient work array
       deallocate(work)
       
    enddo

!----------------------------------------------------------------------
! Fill in the Case 1b spin-coupling coefficients
!----------------------------------------------------------------------
! Note that for the k=+1 component of the triplet spin tensor operator,
! the following relations hold:
!
! (i)  1b i>j = 1a i<j
!
! (ii) 1b i<j = 1a i>j
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=1

    ! Loop over numbers of open shells
    do nopen=1,nocase1

       ! No. bra and ket CSFs
       nspB=ncsfsB(nopen)
       nspK=ncsfsK(nopen)

       ! Cycle if there are no bra or ket CSFs for this number of
       ! open shells
       if (nspB == 0 .or. nspK == 0) cycle

       ! Allocate work arrays
       allocate(workgt(nspB*nspK),worklt(nspB*nspK))

       ! Loop over the *upper triangle* of unique pairs of simplified
       ! spatial occupation vectors
       do is1=1,nopen
          do is2=is1+1,nopen+1

             ! 1a i>j
             istart=ioff
             iend=ioff+nspB*nspK-1
             workgt=spincp(istart:iend)

             ! 1a i<j
             worklt=spincp(istart+offspincp(1):iend+offspincp(1))

             ! 1b i>j
             spincp(istart+offspincp(2):iend+offspincp(2))=worklt

             ! 1b i<j
             spincp(istart+offspincp(3):iend+offspincp(3))=workgt
             
             ! Update the spincp offset
             ioff=ioff+nspB*nspK
             
          enddo
       enddo
          
       ! Deallocate work arrays
       deallocate(workgt,worklt)
       
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ws)
    
    return

  end subroutine case1_coeffs_k1

!######################################################################
! spincp_1a_k1: Computes a single batch of Case 1a spin coupling
!               coefficients for the k=+1 component of the triplet spin
!               tensor operator
!######################################################################
  subroutine spincp_1a_k1(ws1,ws2,nopen,nocase2,maxcsfB,maxcsfK,&
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
    
  end subroutine spincp_1a_k1

!######################################################################
! case2_coeffs_k1: Calculation of the Case 2 spin coupling
!                  coefficients. These are the coefficients
!                  corresponding to pairs of initial and final spatial
!                  occupations with numbers of open shells differing
!                  by two.
!######################################################################
  subroutine case2_coeffs_k1(imultB,imultK,nocase1,nocase2,maxcsfB,&
       maxcsfK,maxdetB,maxdetK,ncsfsB,ncsfsK,ndetsB,ndetsK,csfcoeB,&
       csfcoeK,detvecB,detvecK,spincpdim,spincp,nspincp,mapdim2,&
       patternmap2a,patternmap2b,offspincp)

    use constants
    use bitutils
        
    implicit none

    integer(is), intent(in)  :: imultB,imultK
    integer(is), intent(in)  :: nocase1,nocase2
    integer(is), intent(in)  :: maxcsfB,maxcsfK,maxdetB,maxdetK
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ncsfsK(0:nocase2)
    integer(is), intent(in)  :: ndetsB(0:nocase2),ndetsK(0:nocase2)
    real(dp), intent(in)     :: csfcoeB(maxcsfB,maxdetB,nocase2)
    real(dp), intent(in)     :: csfcoeK(maxcsfK,maxdetK,nocase2)
    integer(ib), intent(in)  :: detvecB(maxdetB,nocase2)
    integer(ib), intent(in)  :: detvecK(maxdetK,nocase2)
    integer(is), intent(in)  :: spincpdim(3)
    real(dp), intent(out)    :: spincp(spincpdim(3))
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(in)  :: mapdim2
    integer(is), intent(out) :: patternmap2a(0:mapdim2)
    integer(is), intent(out) :: patternmap2b(0:mapdim2)
    integer(is), intent(in)  :: offspincp(6)
    
    integer(is)              :: nopen,is1,is2,icsfB,icsfK,lim,i
    integer(is)              :: ioff,istart,iend,nspB,nspK
    integer(ib), allocatable :: ws(:,:)
    integer(ib)              :: wsp(nocase1)
    integer(ib)              :: pattern
    real(dp)                 :: coeff
    real(dp), allocatable    :: work(:),workT(:),work2(:,:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Simplified spatial occupation vectors with N-2 open shells
    allocate(ws((nocase1)*(nocase1-1)/2,0:nocase1))
    ws=0_ib

!----------------------------------------------------------------------
! Generate the simplified final spatial occupation vectors with N+2
! open shells
!----------------------------------------------------------------------
    ! Initialise to all unset bits
    wsp=0_ib

    ! Loop over numbers of open shells
    do nopen=1,nocase1

       ! Construct the vector of nopen 1's
       do i=1,nopen
          wsp(nopen)=ibset(wsp(nopen),i-1)
       enddo

    enddo

!----------------------------------------------------------------------
! Generate the unique Case 2 simplified initial spatial occupation
! vectors with N-2 open shells. i.e, all possible permutations of
! N 1's and two 0's for N=1,...,nocase1-2
!----------------------------------------------------------------------
    !
    ! N = 0 case: one simplified vector corresponding to 00
    !
    ws(1,0)=0_ib

    !
    ! N > 0 cases
    !
    ! Loop over numbers of open shells
    do nopen=1,nocase1-2

       ! Generate the permutations
       lim=(nopen+2)*(nopen+1)/2
       call get_permutations(nopen,2,ws(1:lim,nopen),lim)
       
    enddo

!----------------------------------------------------------------------
! Compute the Case 2a spin coupling coefficient for N = 0
! This only exists for singlets kets as we are enforcing S_B > S_K
! Note that this will then correspond to a 2a i>j coefficient
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=spincpdim(1)+1

    ! N = 0 spin coupling coefficient
    if (imultK == 1) then

       ! Pattern number and corresponding array index
       pattern=iand(ws(1,0),wsp(2))
       patternmap2a(pattern)=ioff

       ! 2a i>j spin coupling coefficient value
       spincp(ioff)=1.0d0

       ! Update the spincp offset
       ioff=ioff+1
       
    endif

!----------------------------------------------------------------------
! Compute the Case 2a i>j spin coupling coefficients for N > 0
! i.e., Ket: 2 0 -> Bra: 1 1
!----------------------------------------------------------------------
    ! Loop over numbers of ket open shells
    do nopen=1,nocase1-2

       ! Cycle if there are no bra or ket CSFs
       if (ncsfsB(nopen+2) == 0 .or. ncsfsK(nopen) == 0) cycle

       ! Allocate the spin-coupling coefficient work array
       nspK=ncsfsK(nopen)
       nspB=ncsfsB(nopen+2)
       allocate(work(nspB*nspK))

       ! Loop over the ket simplified spatial occupation vectors
       do is1=1,(nopen+2)*(nopen+1)/2

          ! Pattern number and the corresponding array index
          pattern=iand(ws(is1,nopen),wsp(nopen+2))
          patternmap2a(pattern)=ioff

          ! Compute the 2a i>j spin coupling coefficients
          call spincp_2a_k1(ws(is1,nopen),wsp(nopen+2),nopen,&
               nopen+2,nocase2,maxcsfB,maxcsfK,maxdetB,maxdetK,&
               ndetsB,ndetsK,detvecB,detvecK,csfcoeB,csfcoeK,nspK,&
               nspB,work)

          ! Save the 2a i>j spin coupling coefficients
          istart=ioff
          iend=ioff+nspB*nspK-1
          spincp(istart:iend)=work
          
          ! Update the spincp offset
          ioff=ioff+nspB*nspK

       enddo
          
       ! Deallocate the spin-coupling coefficient work array
       deallocate(work)

    enddo

!----------------------------------------------------------------------
! Compute the Case 2b i>j spin coupling coefficients
! i.e., Ket: 1 1 -> Bra: 0 2
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=spincpdim(1)+2*offspincp(4)+1

    ! Loop over numbers of ket open shells
    do nopen=2,nocase1

       ! Cycle if there are no bra or ket CSFs
       if (ncsfsB(nopen-2) == 0 .or. ncsfsK(nopen) == 0) cycle

       ! Allocate the spin-coupling coefficient work array
       nspK=ncsfsK(nopen)
       nspB=ncsfsB(nopen-2)
       allocate(work(nspB*nspK))

       ! Loop over the bra simplified spatial occupation vectors
       do is2=1,(nopen)*(nopen-1)/2

          ! Pattern number and the corresponding array index
          pattern=iand(wsp(nopen),ws(is2,nopen-2))
          patternmap2b(pattern)=ioff

          ! Compute the 2b i>j spin coupling coefficients
          call spincp_2b_k1(wsp(nopen),ws(is2,nopen-2),nopen,&
               nopen-2,nocase2,maxcsfB,maxcsfK,maxdetB,maxdetK,&
               ndetsB,ndetsK,detvecB,detvecK,csfcoeB,csfcoeK,nspK,&
               nspB,work)
          
          ! Save the 2b i>j spin coupling coefficients
          istart=ioff
          iend=ioff+nspB*nspK-1
          spincp(istart:iend)=work
          
          ! Update the spincp offset
          ioff=ioff+nspB*nspK

       enddo
          
       ! Deallocate the spin-coupling coefficient work array
       deallocate(work)
       
    enddo

!----------------------------------------------------------------------
! Fill in the Case 2a i<j and 2b i<j spin-coupling coefficients
!----------------------------------------------------------------------
! Note that for the k=+1 component of the triplet spin tensor operator,
! the following relations hold:
!
! (i)  2a i<j = - 2a i>j
!
! (ii) 2b i<j = - 2b i>j
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Fill in the Case 2a i<j spin-coupling coefficients using the relation
! 2a i<j = - 2a i>j
!----------------------------------------------------------------------
    ! Allocate the work array
    istart=spincpdim(1)+1
    iend=istart+offspincp(4)-1
    allocate(work(iend-istart+1))

    ! Copy of the 2a i>j terms
    work=spincp(istart:iend)

    ! Fill in the 2a i<j terms
    istart=iend+1
    iend=istart+offspincp(4)-1
    spincp(istart:iend)=-work

    ! Deallocate the work array
    deallocate(work)

!----------------------------------------------------------------------
! Fill in the Case 2b i<j spin-coupling coefficients using the relation
! 2b i<j = - 2b i>j
!---------------------------------------------------------------------
    ! Allocate the work array
    istart=spincpdim(1)+2*offspincp(4)+1
    iend=istart+offspincp(5)-1
    allocate(work(iend-istart+1))

    ! Copy of the 2b i>j terms
    work=spincp(istart:iend)

    ! Fill in the 2b i<j terms
    istart=iend+1
    iend=istart+offspincp(5)-1
    spincp(istart:iend)=-work

    ! Deallocate the work array
    deallocate(work)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ws)
    
    return
    
  end subroutine case2_coeffs_k1

!######################################################################
! spincp_2a_k1: Computes a batch of 2a i>j k=+1 triplet spin-coupling
!               coefficients
!######################################################################
  subroutine spincp_2a_k1(wsK,wsB,nopenK,nopenB,nocase2,maxcsfB,&
       maxcsfK,maxdetB,maxdetK,ndetsB,ndetsK,detvecB,detvecK,csfcoeB,&
       csfcoeK,ncsfK,ncsfB,spincoe)

    use constants
    use bitutils
    use slater_condon

    integer(ib), intent(in)  :: wsK,wsB
    integer(is), intent(in)  :: nopenK,nopenB,nocase2
    integer(is), intent(in)  :: maxcsfB,maxcsfK,maxdetB,maxdetK
    integer(is), intent(in)  :: ndetsB(0:nocase2),ndetsK(0:nocase2)
    integer(ib), intent(in)  :: detvecB(maxdetB,nocase2)
    integer(ib), intent(in)  :: detvecK(maxdetK,nocase2)
    real(dp), intent(in)     :: csfcoeB(maxcsfB,maxdetB,nocase2)
    real(dp), intent(in)     :: csfcoeK(maxcsfK,maxdetK,nocase2)
    integer(is), intent(in)  :: ncsfK,ncsfB
    real(dp), intent(out)    :: spincoe(ncsfB,ncsfK)
    
    integer(is)              :: icsfB,icsfK
    integer(is)              :: idet,n,vecindx,nunset
    integer(is)              :: iket,ibra,ispin
    integer(is)              :: ic,ia,ih,ip
    integer(ib)              :: b
    integer(ib), allocatable :: dK(:,:),dB(:,:)
    integer(is)              :: n_int_save
    integer(ib)              :: phase_mask(2)
    integer(is)              :: nexci
    integer(ib)              :: p(2),h(2)
    integer(is), parameter   :: maxex=2
    integer(is)              :: plist(maxex,2),hlist(maxex,2)
    integer(is)              :: phase
    integer(is)              :: ndetK,ndetB
    real(dp), allocatable    :: phasemat(:,:),pcT(:,:)
    real(dp), allocatable    :: coeK(:,:),coeB(:,:)

    ! Phase factor accociated with creating the doubly-occupied
    ! MO in the ket determinants from the determinants with
    ! only unoccupied MOs
    integer(is)              :: docc_phase
    integer(is)              :: n2,n2a,n2b,i
    integer(is)              :: noa

!**********************************************************************
! Here we are considering excitations of the form
! Ket: 2 0 -> Bra: 1 1
!**********************************************************************
    
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
    ! Dimensions
    ndetK=ndetsK(nopenK)
    ndetB=ndetsB(nopenB)

    ! Determinant bit strings
    allocate(dK(2,ndetK))
    allocate(dB(2,ndetB))
    dK=0_ib
    dB=0_ib

    ! Phase factor matrix
    allocate(phasemat(ndetB,ndetK))
    phasemat=0.0d0

    ! Phase factor matrix contracted with the transpose of the CSF
    ! coefficient matrix
    allocate(pcT(ndetB,ncsfK))
    pcT=0.0d0
    
    ! Working arrays
    allocate(coeK(ncsfK,ndetK))
    coeK=0.0d0
    allocate(coeB(ncsfB,ndetB))
    coeB=0.0d0
    
!----------------------------------------------------------------------
! Generate the determinants the ket CSF, with N singly-occupied MOs,
! one doubly occupied MO, and one unoccupied MO
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    dK=0_ib

    ! Initialise the detvec open shell counter
    vecindx=0

    ! Initialise the unset bit counter
    nunset=0

    ! Loop over orbitals in the spatial occupation
    do n=1,nopenK+2

       ! We are considering 2a i>j terms here
       ! So, if we are at the first unset bit, then treat it as a
       ! doubly-occupied MO, else treat it as an unoccupied MO
       if (.not. btest(wsK,n-1)) then
          nunset=nunset+1
          if (nunset == 1) then
             ! Position of the doubly-occupied MO
             n2=n
             ! Doubly-occupied MO: set both alpha and beta string bits
             ! in all determinants
             do idet=1,ndetK
                dK(1,idet)=ibset(dK(1,idet),n-1)
                dK(2,idet)=ibset(dK(2,idet),n-1)
             enddo
             cycle
          else
             ! Unoccupied MO: cycle
             cycle
          endif
       endif

       ! Increment the detvec open shell orbital counter
       vecindx=vecindx+1

       ! Loop over ket determinants
       do idet=1,ndetK
          
          ! Add the next spin orbital
          if (btest(detvecK(idet,nopenK),vecindx-1)) then
             ! Occupied alpha spin-orbital
             dK(1,idet)=ibset(dK(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             dK(2,idet)=ibset(dK(2,idet),n-1)
          endif

       enddo
       
    enddo

!----------------------------------------------------------------------
! Generate the determinants for the Bra CSFs, with N+2 singly-occupied
! MOs
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    dB=0_ib

    ! Loop over bra determinants
    do idet=1,ndetB
    
       ! Loop over orbitals in the simplified spatial occupation
       do n=1,nopenB

          ! Add the next spin orbital
          if (btest(detvecB(idet,nopenB),n-1)) then
             ! Occupied alpha spin-orbital
             dB(1,idet)=ibset(dB(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             dB(2,idet)=ibset(dB(2,idet),n-1)
          endif
          
       enddo

    enddo
    
!----------------------------------------------------------------------
! Determine the phase factors associated with operating on the
! open-shell-only determinants to yield the ket determinants (which
! contain a single doubly-occupied MO).
! Note that the phase factor will be the same for all determinants
! with the same spatial configuration.
!----------------------------------------------------------------------
    ! Number of unpaired alpha electrons
    noa=popcnt(dK(1,1))

    ! Number of electrons before the index of the alpha-spin
    ! creation operator
    n2a=0
    if (n2 > 1) then
       do i=1,n2-1
          if (btest(dK(1,1),i-1)) n2a=n2a+1
       enddo
    endif

    ! Number of electrons before the index of the beta-spin
    ! creation operator
    n2b=noa+1
    if (n2 > 1) then
       do i=1,n2-1
          if (btest(dK(2,1),i-1)) n2b=n2b+1
       enddo
    endif
    
    ! Phase factor
    docc_phase=(-1)**(n2a+n2b)

!----------------------------------------------------------------------
! Determine the indices of the triplet excitation operator
! T_ij^(1,k=+1)
!----------------------------------------------------------------------
! Remember that these are given by the positions of the two 0's in
! wsK, i.e., the indices of the first and second unset bits.
!----------------------------------------------------------------------
! As we are considering 2a i>j coefficients here, we are taking the
! first zero in wsK to correspond to the doubly-occupied MO, i.e., to
! the annihilator index. If we ever change this, then the creation
! operator index will have to be determined before the annihilation
! operator index due to the use of the bit clearing operation.
!----------------------------------------------------------------------
    ! Temporary bit string array
    b=not(wsK)
    
    ! Annihilation operator index
    ia=trailz(b)+1
    
    ! Creation operator index
    b=ibclr(b,ia-1)
    ic=trailz(b)+1

!----------------------------------------------------------------------
! Compute the matrix of phase factors for all pairs of determinants
! for the spin-flip excitation operator a_ic,alpha^\dagger a_ia,beta
!----------------------------------------------------------------------
! phasemat(I,J) = < det_I | a_ic,alpha^\dagger a_ia,beta | det_J >
!----------------------------------------------------------------------
    ! Initialisation
    phasemat=0.0d0
    
    ! Loop over ket determinants
    do iket=1,ndetK

       ! Loop over bra determinants
       do ibra=1,ndetB

          ! Get the excitation degree
          nexci=exc_degree_det(dB(:,ibra),dK(:,iket))
          
          ! Cycle if the excitation degree is not equal to 1
          if (nexci /= 1) cycle

          ! Get the indices of the spin-orbitals involved in the
          ! excitations linking the bra and ket determinants
          call exc(dK(:,iket),dB(:,ibra),p,h)
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
          phase=phase_slow(dK(:,iket),hlist,plist,maxex,nexci)
          phasemat(ibra,iket)=dble(phase)
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Compute the matrix of spin coupling coefficients for the current
! pattern index
!----------------------------------------------------------------------
    ! Working coefficient arrays
    coeB=csfcoeB(1:ncsfB,1:ndetB,nopenB)
    coeK=csfcoeK(1:ncsfK,1:ndetK,nopenK)

    ! Multiply the ket CSF expansion coefficients by the phase factors
    ! associated with creating the doubly-occupied MO
    coeK=coeK*docc_phase

    ! Contract the phase matrix with the ket CSF expansion coefficients,
    ! pcT = phasemat * transpose(csfcoeK)
    call dgemm('N','T',ndetB,ncsfK,ndetK,1.0d0,phasemat,ndetB,coeK,&
         ncsfK,0.0d0,pcT,ndetB)

    ! Contract the bra CSF expansion coefficients with the intermediate
    ! matrix pcT to yield the spin coupling coefficients
    call dgemm('N','N',ncsfB,ncsfK,ndetB,1.0d0,coeB,ncsfB,pcT,ndetB,&
         0.0d0,spincoe,ncsfB)

    ! Up to now, we have computed the matrix elements
    ! < CSF_B | a_ialpha^dagger a_jbeta | CSF_K >
    !
    ! Multiply by -1 to get the matrix elements
    ! < CSF_B | T_ij^(1,k=+1) | CSF_K >
    spincoe=-spincoe
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dK)
    deallocate(dB)
    deallocate(phasemat)
    deallocate(pcT)
    deallocate(coeK)
    deallocate(coeB)
    
!----------------------------------------------------------------------
! Reset n_int
!----------------------------------------------------------------------
    n_int=n_int_save
    
    return
    
  end subroutine spincp_2a_k1

!######################################################################
! spincp_2b_k1: Computes a batch of 2b i>j k=+1 triplet spin-coupling
!               coefficients
!######################################################################
  subroutine spincp_2b_k1(wsK,wsB,nopenK,nopenB,nocase2,maxcsfB,&
       maxcsfK,maxdetB,maxdetK,ndetsB,ndetsK,detvecB,detvecK,csfcoeB,&
       csfcoeK,ncsfK,ncsfB,spincoe)

    use constants
    use bitutils
    use slater_condon

    integer(ib), intent(in)  :: wsK,wsB
    integer(is), intent(in)  :: nopenK,nopenB,nocase2
    integer(is), intent(in)  :: maxcsfB,maxcsfK,maxdetB,maxdetK
    integer(is), intent(in)  :: ndetsB(0:nocase2),ndetsK(0:nocase2)
    integer(ib), intent(in)  :: detvecB(maxdetB,nocase2)
    integer(ib), intent(in)  :: detvecK(maxdetK,nocase2)
    real(dp), intent(in)     :: csfcoeB(maxcsfB,maxdetB,nocase2)
    real(dp), intent(in)     :: csfcoeK(maxcsfK,maxdetK,nocase2)
    integer(is), intent(in)  :: ncsfK,ncsfB
    real(dp), intent(out)    :: spincoe(ncsfB,ncsfK)
    
    integer(is)              :: icsfB,icsfK
    integer(is)              :: idet,n,vecindx,nunset
    integer(is)              :: iket,ibra,ispin
    integer(is)              :: ic,ia,ih,ip
    integer(ib)              :: b
    integer(ib), allocatable :: dK(:,:),dB(:,:)
    integer(is)              :: n_int_save
    integer(ib)              :: phase_mask(2)
    integer(is)              :: nexci
    integer(ib)              :: p(2),h(2)
    integer(is), parameter   :: maxex=2
    integer(is)              :: plist(maxex,2),hlist(maxex,2)
    integer(is)              :: phase
    integer(is)              :: ndetK,ndetB
    real(dp), allocatable    :: phasemat(:,:),pcT(:,:)
    real(dp), allocatable    :: coeK(:,:),coeB(:,:)

    ! Phase factor accociated with creating the doubly-occupied
    ! MO in the ket determinants from the determinants with
    ! only unoccupied MOs
    integer(is)              :: docc_phase
    integer(is)              :: n2,n2a,n2b,i
    integer(is)              :: noa

!**********************************************************************
! Here we are considering excitations of the form
! Ket: 1 1 -> Bra: 0 2
!**********************************************************************
    
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
    ! Dimensions
    ndetK=ndetsK(nopenK)
    ndetB=ndetsB(nopenB)

    ! Determinant bit strings
    allocate(dK(2,ndetK))
    allocate(dB(2,ndetB))
    dK=0_ib
    dB=0_ib

    ! Phase factor matrix
    allocate(phasemat(ndetB,ndetK))
    phasemat=0.0d0

    ! Phase factor matrix contracted with the transpose of the CSF
    ! coefficient matrix
    allocate(pcT(ndetB,ncsfK))
    pcT=0.0d0
    
    ! Working arrays
    allocate(coeK(ncsfK,ndetK))
    coeK=0.0d0
    allocate(coeB(ncsfB,ndetB))
    coeB=0.0d0

!----------------------------------------------------------------------
! Generate the determinants the bra CSF, with N singly-occupied MOs,
! one doubly occupied MO, and one unoccupied MO
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    dB=0_ib

    ! Initialise the detvec open shell counter
    vecindx=0

    ! Initialise the unset bit counter
    nunset=0

    ! Loop over orbitals in the spatial occupation
    do n=1,nopenB+2

       ! We are considering 2b i>j terms here
       ! So, if we are at the first unset bit, then treat it as an
       ! unoccupied MO, else treat it as a doubly-occupied MO
       if (.not. btest(wsB,n-1)) then
          nunset=nunset+1
          if (nunset == 1) then
             ! Unoccupied MO: cycle
             cycle
          else if (nunset == 2) then
             ! Position of the doubly-occupied MO
             n2=n
             ! Doubly-occupied MO: set both alpha and beta string bits
             ! in all determinants
             do idet=1,ndetB
                dB(1,idet)=ibset(dB(1,idet),n-1)
                dB(2,idet)=ibset(dB(2,idet),n-1)
             enddo
             cycle
          endif
       endif

       ! Increment the detvec open shell orbital counter
       vecindx=vecindx+1

       ! Loop over bra determinants
       do idet=1,ndetB
          
          ! Add the next spin orbital
          if (btest(detvecB(idet,nopenB),vecindx-1)) then
             ! Occupied alpha spin-orbital
             dB(1,idet)=ibset(dB(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             dB(2,idet)=ibset(dB(2,idet),n-1)
          endif

       enddo
       
    enddo

!----------------------------------------------------------------------
! Generate the determinants for the ket CSFs, with N+2 singly-occupied
! MOs
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    dK=0_ib

    ! Loop over ket determinants
    do idet=1,ndetK
    
       ! Loop over orbitals in the simplified spatial occupation
       do n=1,nopenK

          ! Add the next spin orbital
          if (btest(detvecK(idet,nopenK),n-1)) then
             ! Occupied alpha spin-orbital
             dK(1,idet)=ibset(dK(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             dK(2,idet)=ibset(dK(2,idet),n-1)
          endif
          
       enddo

    enddo

!----------------------------------------------------------------------
! Determine the phase factors associated with operating on the
! open-shell-only determinants to yield the bra determinants (which
! contain a single doubly-occupied MO).
! Note that the phase factor will be the same for all determinants
! with the same spatial configuration.
!----------------------------------------------------------------------
    ! Number of unpaired alpha electrons
    noa=popcnt(dB(1,1))

    ! Number of electrons before the index of the alpha-spin
    ! creation operator
    n2a=0
    if (n2 > 1) then
       do i=1,n2-1
          if (btest(dB(1,1),i-1)) n2a=n2a+1
       enddo
    endif

    ! Number of electrons before the index of the beta-spin
    ! creation operator
    n2b=noa+1
    if (n2 > 1) then
       do i=1,n2-1
          if (btest(dB(2,1),i-1)) n2b=n2b+1
       enddo
    endif
    
    ! Phase factor
    docc_phase=(-1)**(n2a+n2b)

!----------------------------------------------------------------------
! Determine the indices of the triplet excitation operator
! T_ij^(1,k=+1)
!----------------------------------------------------------------------
! Remember that these are given by the positions of the two 0's in
! wsB, i.e., the indices of the first and second unset bits.
!----------------------------------------------------------------------
! As we are considering 2b i>j coefficients here, we are taking the
! first zero in wsB to correspond to the unoccupied MO, i.e., to the
! annihilator index. If we ever change this, then the creation operator
! index will have to be determined before the annihilation operator
! index due to the use of the bit clearing operation.
!----------------------------------------------------------------------
    ! Temporary bit string array
    b=not(wsB)

    ! Annihilation operator index
    ia=trailz(b)+1
    
    ! Creation operator index
    b=ibclr(b,ia-1)
    ic=trailz(b)+1

!----------------------------------------------------------------------
! Compute the matrix of phase factors for all pairs of determinants
! for the spin-flip excitation operator a_ic,alpha^\dagger a_ia,beta
!----------------------------------------------------------------------
! phasemat(I,J) = < det_I | a_ic,alpha^\dagger a_ia,beta | det_J >
!----------------------------------------------------------------------
    ! Initialisation
    phasemat=0.0d0

    ! Loop over ket determinants
    do iket=1,ndetK

       ! Loop over bra determinants
       do ibra=1,ndetB

          ! Get the excitation degree
          nexci=exc_degree_det(dB(:,ibra),dK(:,iket))
          
          ! Cycle if the excitation degree is not equal to 1
          if (nexci /= 1) cycle

          ! Get the indices of the spin-orbitals involved in the
          ! excitations linking the bra and ket determinants
          call exc(dK(:,iket),dB(:,ibra),p,h)
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
          phase=phase_slow(dK(:,iket),hlist,plist,maxex,nexci)
          phasemat(ibra,iket)=dble(phase)
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Compute the matrix of spin coupling coefficients for the current
! pattern index
!----------------------------------------------------------------------
    ! Working coefficient arrays
    coeB=csfcoeB(1:ncsfB,1:ndetB,nopenB)
    coeK=csfcoeK(1:ncsfK,1:ndetK,nopenK)

    ! Multiply the bra CSF expansion coefficients by the phase factors
    ! associated with creating the doubly-occupied MO
    coeB=coeB*docc_phase

    ! Contract the phase matrix with the ket CSF expansion coefficients,
    ! pcT = phasemat * transpose(csfcoeK)
    call dgemm('N','T',ndetB,ncsfK,ndetK,1.0d0,phasemat,ndetB,coeK,&
         ncsfK,0.0d0,pcT,ndetB)

    ! Contract the bra CSF expansion coefficients with the intermediate
    ! matrix pcT to yield the spin coupling coefficients
    call dgemm('N','N',ncsfB,ncsfK,ndetB,1.0d0,coeB,ncsfB,pcT,ndetB,&
         0.0d0,spincoe,ncsfB)

    ! Up to now, we have computed the matrix elements
    ! < CSF_B | a_ialpha^dagger a_jbeta | CSF_K >
    !
    ! Multiply by -1 to get the matrix elements
    ! < CSF_B | T_ij^(1,k=+1) | CSF_K >
    spincoe=-spincoe

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(dK)
    deallocate(dB)
    deallocate(phasemat)
    deallocate(pcT)
    deallocate(coeK)
    deallocate(coeB)
    
!----------------------------------------------------------------------
! Reset n_int
!----------------------------------------------------------------------
    n_int=n_int_save
    
    return
    
  end subroutine spincp_2b_k1

!######################################################################
! fill_N1s: Fills in the array of bit strings with the first N bits
!           set. The spin coupling coefficient pattern numbers will
!           be obtained by unsetting bits in these strings.
!######################################################################
  subroutine fill_N1s(nocase2,N1s)

    use constants
    
    implicit none

    integer(is), intent(in)  :: nocase2
    integer(ib), allocatable :: N1s(:)
    
    integer(is)              :: i
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(N1s(nocase2+2))
    
!----------------------------------------------------------------------
! Fill in the arrays of N set bits
!----------------------------------------------------------------------
    ! Initialisation
    N1s=0_ib
    N1s(1)=ibset(N1s(1),0)
    
    ! Loop over number of bits
    do i=2,nocase2+2
       ! Bit string with only the first i bits set
       N1s(i)=ibset(N1s(i-1),i-1)
    enddo
    
    return
    
  end subroutine fill_N1s
  
!######################################################################
  
end module spin_coupling_k1
