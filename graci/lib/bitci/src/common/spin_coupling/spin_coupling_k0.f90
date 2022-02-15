!**********************************************************************
! Routines for the pre-computation of the spin-coupling coefficients
! <w omega| T_pq^(1,k) |w' omega'> for triplet spin tensor operators
! T_pq^(1,k), k=0
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
module spin_coupling_k0

  implicit none

  private

  public :: scc_k0

contains

!######################################################################
! scc_k0: Generates the unique k=0 component triplet spin-coupling
!         coefficients for a given spin multiplicity
!######################################################################
  subroutine scc_k0(imult,nocase1,nocase2,maxcsf,maxdet,ncsfs,ndets,&
       csfcoe,detvec,nspincp,N1s,verbose,spincp,patmap,offspincp)

    use constants
    use iomod
    use timing

    implicit none

    ! Spin multiplicity
    integer(is), intent(in)  :: imult

    ! Maximum number open shells for the Case 1 and Case 2 CSFs
    integer(is), intent(in)  :: nocase1,nocase2

    ! Numbers of CSFs and determinants as a function of the number
    ! of open shells
    integer(is), intent(in)  :: ncsfs(0:nocase2),ndets(0:nocase2)

    ! Maximum number of CSFs/determinants across all numbers of
    ! open shells
    integer(is), intent(in)  :: maxcsf,maxdet

    ! CSF expansion coefficients
    real(dp), intent(in)     :: csfcoe(maxcsf,maxdet,nocase2)

    ! Bit string encoding of the determinants contributing to the
    ! CSFs
    integer(ib), intent(in)  :: detvec(maxdet,nocase2)

    ! Number of unique spin coupling coefficients
    integer(is), intent(out) :: nspincp(2)

    ! All spin coupling coefficients
    integer(is)              :: spincpdim(3)
    real(dp), allocatable    :: spincp(:)

    ! All pattern -> array index mappings
    integer(is)              :: mapdim
    integer(is), allocatable :: patmap(:)

    ! Spin coupling coefficient offsets for the various
    ! different cases
    integer(is), allocatable :: offspincp(:)
    
    ! Bit strings comprised of N 1's
    integer(ib), allocatable :: N1s(:)

    ! Verbose output
    logical, intent(in)      :: verbose

    ! Everything else
    integer(is) :: i
    real(dp)    :: tcpu_start,tcpu_end,twall_start,twall_end

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
! Compute the number of unique spin coupling coefficients for the
! spin multiplicity under consideration
!----------------------------------------------------------------------
    call get_nunique_k0(imult,nocase1,nocase2,ncsfs,nspincp)

!----------------------------------------------------------------------
! Output some information about what we are doing
!----------------------------------------------------------------------
    if (verbose) call print_spincp_info_k0(imult,nspincp)

!----------------------------------------------------------------------
! Allocate the spin coupling coefficient arrays
!----------------------------------------------------------------------
    call init_spincp_k0(imult,nocase1,nocase2,ncsfs,nspincp,verbose,&
         spincp,spincpdim,offspincp)

!----------------------------------------------------------------------
! Allocate the pattern value -> array index mapping array
!----------------------------------------------------------------------
    call init_patternmap_k0(imult,nocase1,nocase2,ncsfs,patmap,mapdim)

!----------------------------------------------------------------------
! Compute the Case 1 spin coupling coefficients
!----------------------------------------------------------------------
    call case1_coeffs_k0(nocase1,nocase2,maxcsf,maxdet,ncsfs,ndets,&
         csfcoe,detvec,spincpdim,spincp,nspincp,mapdim,patmap)

    STOP
    
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) call report_times(twall_end-twall_start,&
         tcpu_end-tcpu_start,'scc_k0')

!----------------------------------------------------------------------    
! Flush stdout
!----------------------------------------------------------------------
    flush(6)
    
    return
    
  end subroutine scc_k0

!######################################################################
! get_nunique_k0: Computes the total number of *unique* spin coupling
!                 coefficients as a function of the spin multiplicity
!                 and the maximum number of open shells
!######################################################################
  subroutine get_nunique_k0(imult,nocase1,nocase2,ncsfs,nspincp)

    use constants
    use math
    
    implicit none

    integer(is), intent(in)  :: imult,nocase1,nocase2
    integer(is), intent(in)  :: ncsfs(0:nocase2)
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
       nspincp(1)=nspincp(1)+npat*ncsfs(nopen)**2
       
    enddo

!----------------------------------------------------------------------
! Case 2: N_Bra = N_Ket +/- 2
!----------------------------------------------------------------------
! This covers excitations from a doubly-occupied orbital into an
! unoccupied orbital or from a singly-occupied orbital into another
! singly-occupied orbital
!---------------------------------------------------------------------- 
    ! Loop over numbers of open shells
    do nopen=2,nocase1

       ! Number of patterns
       npat=(nopen)*(nopen-1)/2

       ! Number of unique spin-coupling coefficients corresponding to
       ! this pattern
       nspincp(2)=nspincp(2)+npat*ncsfs(nopen)*ncsfs(nopen-2)
       
    enddo

    return
    
  end subroutine get_nunique_k0

!######################################################################
! print_spincp_info_k0: Printing of some information about the
!                       k=0 component triplet spin coupling
!                       coefficients being calculated
!######################################################################
  subroutine print_spincp_info_k0(imult,nspincp)
    
    use constants
        
    implicit none

    integer(is), intent(in) :: imult
    integer(is), intent(in) :: nspincp(2)
    
    integer(is)             :: n,i

!----------------------------------------------------------------------
! Spin tensor operator component    
!----------------------------------------------------------------------
    write(6,'(/,x,a)') 'Spin tensor component: 0'

!----------------------------------------------------------------------
! Spin multiplicity
!----------------------------------------------------------------------
    write(6,'(/,x,a,x,i0)') 'Bra/ket spin multiplicity:',imult

!----------------------------------------------------------------------
! Number of unique spin-coupling coefficients
!----------------------------------------------------------------------
    ! Total number of spin-coupling coefficients
    write(6,'(/,x,a,2x,i0)') &
         'Total number of spin coupling coefficients:',sum(nspincp)
    
    ! Number of Case 1 spin-coupling coefficients
    write(6,'(x,a,x,i0)') &
         'Number of Case 1 spin coupling coefficients:',nspincp(1)
    
    ! Number of Case 2 spin-coupling coefficients
    write(6,'(x,a,x,i0)') &
         'Number of Case 2 spin coupling coefficients:',nspincp(2)
    
    return
    
  end subroutine print_spincp_info_k0

!######################################################################
! init_spincp_k0: Allocation and initialisation of the various spin
!                 coupling coefficient arrays for the k=0 component
!                 of the triplet spin tensor operator
!######################################################################
  subroutine init_spincp_k0(imult,nocase1,nocase2,ncsfs,nspincp,&
       verbose,spincp,spincpdim,offspincp)

    use constants
        
    implicit none

    integer(is), intent(in)  :: imult,nocase1,nocase2
    integer(is), intent(in)  :: ncsfs(0:nocase2)
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(out) :: spincpdim(3)
    integer(is), allocatable :: offspincp(:)
    real(dp), allocatable    :: spincp(:)
    logical, intent(in)      :: verbose
        
    integer(is)              :: nopen
    real(dp)                 :: mem

!----------------------------------------------------------------------
! Numbers of spin coupling coefficients
!----------------------------------------------------------------------
    ! Case 1 coefficients: (i) 1a i>j, (ii) 1a i<j, (iii) 1b i>j,
    !                     (iv) 1b i<j
    spincpdim(1)=4*nspincp(1)

    ! Case 2 coefficients: (i) 2a i>j, (ii) 2a i<j, (iii) 2b i>j,
    !                     (iv) 2b i<j
    spincpdim(2)=4*nspincp(2)

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
! The k=0 spin coupling coefficients will be stored in the order
! ([1a i>j], [1a i<j], [1b i>j], [1b i<j], [2a i>j], [2a i<j],
!  [2b i>j], [2b i<j])
!----------------------------------------------------------------------
! The offsets correspond to the shifts that will be needed to be added
! to the pattern indices to reach each block of coefficients in the
! spincp array
!----------------------------------------------------------------------
    ! Case 1a, i<j
    offspincp(1)=nspincp(1)

    ! Case 1b, i>j
    offspincp(2)=2*nspincp(1)

    ! Case 1b, i<j
    offspincp(3)=3*nspincp(1)

    ! Case 2a, i<j
    offspincp(4)=nspincp(2)

    ! Case 2b, i>j
    offspincp(5)=2*nspincp(2)

    ! Case 2b, i<j
    offspincp(6)=3*nspincp(2)

!----------------------------------------------------------------------
! Output the memory used
!----------------------------------------------------------------------
    if (verbose) then
       mem=spincpdim(3)*8/1024.0d0**2
       write(6,'(/,x,a,F8.2,x,a)') 'Memory used:',mem,'MB'
    endif
    
    return
    
  end subroutine init_spincp_k0

!######################################################################
! init_patternmap_k0: Allocation and initialisation of the array
!                     holding the pattern value -> array index mapping
!######################################################################
  subroutine init_patternmap_k0(imult,nocase1,nocase2,ncsfs,patmap,&
       mapdim)

    use constants
    
    implicit none

    integer(is), intent(in)  :: imult
    integer(is), intent(in)  :: nocase1,nocase2
    integer(is), intent(in)  :: ncsfs(0:nocase2)
    integer(is), intent(out) :: mapdim
    integer(is), allocatable :: patmap(:)

!----------------------------------------------------------------------
! Maximum possible pattern value
!----------------------------------------------------------------------
    ! Case 1 pattern values
    if (ncsfs(nocase1) /= 0) then
       mapdim=2**(nocase1+1)-4
    else
       mapdim=2**nocase1-4
    endif
    
    ! Case 2 pattern values
    if (modulo(imult,2) == 0) then
       mapdim=max(mapdim,2**(nocase1-1)-4)
    else
       mapdim=max(mapdim,2**nocase1-4)
    endif
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(patmap(0:mapdim))
    patmap=0
    
    return
    
  end subroutine init_patternmap_k0

!######################################################################
! case1_coeffs_k0: Calculation of the k=0 Case 1 triplet spin coupling
!                  coefficients. These are the coefficients
!                  corresponding to pairs of initial and final spatial
!                  occupations with equal numbers of open shells
!######################################################################
  subroutine case1_coeffs_k0(nocase1,nocase2,maxcsf,maxdet,ncsfs,&
       ndets,csfcoe,detvec,spincpdim,spincp,nspincp,mapdim,patmap)

    use constants
    use bitutils
        
    implicit none

    integer(is), intent(in)  :: nocase1,nocase2,maxcsf,maxdet
    integer(is), intent(in)  :: ncsfs(0:nocase2),ndets(0:nocase2)
    real(dp), intent(in)     :: csfcoe(maxcsf,maxdet,nocase2)
    integer(ib), intent(in)  :: detvec(maxdet,nocase2)
    integer(is), intent(in)  :: spincpdim(3)
    real(dp), intent(out)    :: spincp(spincpdim(3))
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(in)  :: mapdim
    integer(is), intent(out) :: patmap(0:mapdim)
    
    integer(is)              :: nopen,is1,is2,icsf1,icsf2,nsp,ndet
    integer(is)              :: ioff,istart,iend
    integer(ib), allocatable :: ws(:,:)
    integer(ib)              :: pattern
    real(dp)                 :: coeff
    real(dp), allocatable    :: work(:),workT(:),work2(:,:)
    
    print*,''
    print*,'here'
    stop
    
    return
    
  end subroutine case1_coeffs_k0
    
!######################################################################
  
end module spin_coupling_k0
