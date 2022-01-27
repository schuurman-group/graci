!**********************************************************************
! Routines for the pre-computation of the unique spin-coupling
! coefficients <w omega| T_pq^(1,k) |w' omega'> for triplet spin tensor
! operators T_pq^(1,k), k=0,+1
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
module spin_coupling_triplet

  implicit none

  private
  
  public :: generate_triplet_coupling_coefficients
  
contains
  
!######################################################################
! generate_triplet_coupling_coefficients:
!
! Generates the unique triplet spin-coupling coefficients for a given
! pair of bra and ket spin multiplicities
!
! Considers a single component of the triplet spin tensor operator
! T^(1,k) indexed by the input variable kindx 
!######################################################################
  subroutine generate_triplet_coupling_coefficients(kindx,imultB,&
       imultK,nocase1,nocase2,maxcsfB,maxdetB,ncsfsB,ndetsB,csfcoeB,&
       detvecB,maxcsfK,maxdetK,ncsfsK,ndetsK,csfcoeK,detvecK,&
       npattern1,npattern2,nspincp,N1s,verbose,spincp,patternmap,&
       offspincp)

    use constants
    use iomod
    use timing

    ! Component of the triplet spin tensor operator to use
    integer(is), intent(in) :: kindx

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
    
    ! Number of Case 1 and Case 2 patterns
    integer(is), intent(out) :: npattern1,npattern2

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

    ! THIS WILL NEED RE-DIMENSIONING DEPENDING ON THE
    ! SPIN-TENSOR OPERATOR CONSIDERED
    ! SO, IN GENERAL, WE SHOULD MAKE THIS ARRAY ALLOCATABLE
    ! IN BITGLOBAL
    integer(is), intent(out) :: offspincp(4)

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
! Check on the component of the triplet spin tensor operator under
! consideration
!----------------------------------------------------------------------
! Note that we only support the calculation of coupling coefficients
! for the k=0 and k=+1 components, as the k=-1 and +1 components are
! linked via T^(1,k=+1) = -[T^(1,k=-1)]^\dagger
!----------------------------------------------------------------------
    if (kindx /= 0 .and. kindx /= 1) then
       write(ak,'(i0)') kindx
       errmsg='Illegal input spin tensor operator component:' &
            //trim(ak)
       call error_control
    endif

!----------------------------------------------------------------------
! Check on the bra and ket spin multiplicities
!----------------------------------------------------------------------
    ! k=0 only makes sense for equal bra and ket spin multiplicities
    if (kindx == 0 .and. imultB /= imultK) then
       errmsg='Bra and ket multiplicities not the same: ' &
            //'k=0 makes no sense'
       call error_control
    endif

    ! k=+1 only makes sense if the bra and ket total spins differ
    ! by one
    if (kindx == 1) then
       SB=0.5d0*(dble(imultB)-1.0d0)
       SK=0.5d0*(dble(imultK)-1.0d0)
       if (abs(SB-SK) /= 1.0d0) then
          errmsg='|S_bra - S_ket| != 0: k=+1 makes no sense'
          call error_control
       endif
    endif
    
!----------------------------------------------------------------------
! Compute the number of unique spin coupling coefficients for the
! spin multiplicity under consideration
!----------------------------------------------------------------------
    select case(kindx)
    case(0)
       ! k=0       
       errmsg='k=0 not yet implemented'
       call error_control
    case(1)
       ! k=1
       call get_nunique_kplus(kindx,imultB,imultK,nocase1,nocase2,&
            ncsfsB,ncsfsK,nspincp)
    end select

!----------------------------------------------------------------------
! Output some information about what we are doing
!----------------------------------------------------------------------
    if (verbose) call print_triplet_spincp_info(kindx,imultB,imultK,&
         nspincp)

!----------------------------------------------------------------------
! Allocate the spin coupling coefficient arrays
!----------------------------------------------------------------------
    call init_triplet_spincp_arrays(kindx,imultB,imultK,nocase1,&
         nocase2,ncsfsB,ncsfsK,nspincp,npattern1,npattern2,verbose,&
         spincp,spincpdim,offspincp)
    
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
    
  end subroutine generate_triplet_coupling_coefficients

!######################################################################
! get_nunique_kplus: Computes the total number of unique spin coupling
!                    coefficients as a function of the bra and ket spin
!                    multiplicities and the maximum number of open
!                    shells
!######################################################################
  subroutine get_nunique_kplus(kindx,imultB,imultK,nocase1,nocase2,&
       ncsfsB,ncsfsK,nspincp)

    use constants
    use math
    
    implicit none

    integer(is), intent(in)  :: kindx,imultB,imultK,nocase1,nocase2
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
    ! *** This is implicitly assuming that k=+1 ***
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
    ! *** This is implicitly assuming that k=+1 ***
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
    
  end subroutine get_nunique_kplus

!######################################################################
! print_triplt_spincp_info: Printing of some information about the
!                           triplet spin coupling coefficients being
!                           calculated
!######################################################################
  subroutine print_triplet_spincp_info(kindx,imultB,imultK,nspincp)
    
    use constants
        
    implicit none

    integer(is), intent(in) :: kindx,imultB,imultK
    integer(is), intent(in) :: nspincp(2)
    
    integer(is)             :: n,i
    character(len=10)       :: ak
    
!----------------------------------------------------------------------
! Triplet spin tensor component
!----------------------------------------------------------------------
    if (kindx == 0) then
       ak='0'
    else if (kindx == 1) then
       ak='+1'
    endif
    
    write(6,'(/,x,a)') 'Spin tensor component: '//trim(ak)
    
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
    
  end subroutine print_triplet_spincp_info
    
!######################################################################
! init_triplet_spincp_arrays: Allocation and initialisation of the
!                             various spin coupling coefficient arrays
!######################################################################
  subroutine init_triplet_spincp_arrays(kindx,imultB,imultK,nocase1,&
       nocase2,ncsfsB,ncsfsK,nspincp,npattern1,npattern2,verbose,&
       spincp,spincpdim,offspincp)

    use constants

    implicit none

    integer(is), intent(in)  :: kindx,imultB,imultK,nocase1,nocase2
    integer(is), intent(in)  :: ncsfsB(0:nocase2),ncsfsK(0:nocase2)
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(out) :: npattern1,npattern2
    integer(is), intent(out) :: spincpdim(3)
    integer(is), intent(out) :: offspincp(4)
    real(dp), allocatable    :: spincp(:)
    logical, intent(in)      :: verbose
        
    integer(is)              :: nopen
    real(dp)                 :: mem

    
    return
    
  end subroutine init_triplet_spincp_arrays

!######################################################################
  
end module spin_coupling_triplet
