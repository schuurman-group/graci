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

    ! Where did this come from?
    !if (ncsfs(nocase1) /= 0) then
    !   mapdim=2**(nocase1+1)-4
    !else
    !   mapdim=2**nocase1-4
    !endif

    ! Case 2 pattern values

    ! HOW DOES THIS CHANGE NOW THAT WE ARE ONLY CONSIDERING
    ! 1-ELECTRON MATRIX ELEMENTS?
    ! I.E., NOW THAT WE ARE ONLY CONSIDERING NUMBERS OF OPEN SHELLS
    ! UP TO NOCASE1=NOMAX
    
    return
    
  end subroutine init_patternmap_k1
    
!######################################################################
  
end module spin_coupling_k1
