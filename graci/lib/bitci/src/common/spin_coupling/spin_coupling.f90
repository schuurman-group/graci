!**********************************************************************
! Routines for the pre-computation of the unique spin coupling
! coefficients for the singlet spin tensor operators E_p^q
!**********************************************************************
! Follows the formalism detailed in R. W. Wetmore and G. A. Segal,
! Chem. Phys. Lett., 36, 478 (1975).
!**********************************************************************
module spin_coupling

  implicit none

contains

!######################################################################
! generate_coupling_coefficients: Does what it says on the tin - 
!                                 generates the unique spin coupling
!                                 coefficients
!######################################################################
  subroutine generate_coupling_coefficients(imult,nocase1,nocase2,&
       maxcsf,maxdet,ncsfs,ndets,csfcoe,detvec,npattern1,&
       npattern2,nspincp,N1s,verbose,spincp,patternmap,offspincp)

    use constants
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

    ! Number of Case 1 and Case 2 patterns
    integer(is), intent(out) :: npattern1,npattern2

    ! Number of unique spin coupling coefficients
    integer(is), intent(out) :: nspincp(2)

    ! All spin coupling coefficients
    integer(is)              :: spincpdim(3)
    real(dp), allocatable    :: spincp(:)

    ! All pattern -> array index mappings
    integer(is)              :: mapdim
    integer(is), allocatable :: patternmap(:)

    ! Spin coupling coefficient offsets for the various
    ! different cases
    integer(is), intent(out) :: offspincp(4)
    
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
       write(6,'(x,a)') 'Spin-coupling coefficient generation'
       write(6,'(52a)') ('-',i=1,52)
    endif
       
!----------------------------------------------------------------------
! Compute the number of unique spin coupling coefficients for the
! spin multiplicity under consideration
!----------------------------------------------------------------------
    call get_nunique(imult,nocase1,nocase2,ncsfs,nspincp)

!----------------------------------------------------------------------
! Output some information about what we are doing
!----------------------------------------------------------------------
    call print_spincp_info(imult,nspincp,verbose)

!----------------------------------------------------------------------
! Allocate the spin coupling coefficient arrays
!----------------------------------------------------------------------
    call init_spincp_arrays(imult,nocase1,nocase2,ncsfs,nspincp,&
         npattern1,npattern2,verbose,spincp,spincpdim,offspincp)

!----------------------------------------------------------------------
! Allocate the pattern value -> array index mapping array
!----------------------------------------------------------------------
    call init_patternmap(imult,nocase1,nocase2,ncsfs,patternmap,mapdim)
    
!----------------------------------------------------------------------
! Compute the Case 1 spin coupling coefficients
!----------------------------------------------------------------------
    call case1_coeffs(nocase1,nocase2,maxcsf,maxdet,ncsfs,ndets,&
         csfcoe,detvec,npattern1,spincpdim,spincp,nspincp,mapdim,&
         patternmap)
    
!----------------------------------------------------------------------
! Compute the Case 2 spin coupling coefficients
!----------------------------------------------------------------------
    call case2_coeffs(imult,nocase2,maxcsf,maxdet,ncsfs,ndets,csfcoe,&
         detvec,spincpdim,spincp,nspincp,mapdim,patternmap)

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
         tcpu_end-tcpu_start,'generate_coupling_coefficients')

    return
  
  end subroutine generate_coupling_coefficients

!######################################################################
! get_nunique: Computes the total number of unique spin coupling
!              coefficients as a function of the spin multiplicity and
!              the maximum number of open shells
!######################################################################
  subroutine get_nunique(imult,nocase1,nocase2,ncsfs,nspincp)

    use constants
    use math
    
    implicit none

    integer(is), intent(in)  :: imult,nocase1,nocase2
    integer(is), intent(in)  :: ncsfs(0:nocase2)
    integer(is), intent(out) :: nspincp(2)
    
    integer(is)              :: nopen,n
    
!----------------------------------------------------------------------
! Case 1: N'=N
!----------------------------------------------------------------------
! This covers excitations from a singly-occupied orbital into an
! unoccupied orbital and from a doubly-occupied orbital into a
! singly-occupied orbital
!----------------------------------------------------------------------
    nspincp(1)=0

    ! Loop over numbers of open shells
    do nopen=1,nocase1
       
       ! Cycle if there are no CSFs for nopen open shells
       if (ncsfs(nopen) == 0) cycle
       
       ! Number of unique pairs of permutations of N 1's and one 0
       ! multiplied by the number of pairs of CSFs
       nspincp(1)=nspincp(1)+(nopen+1)*nopen/2 &
            *ncsfs(nopen)*ncsfs(nopen)
       
    enddo

!----------------------------------------------------------------------
! Case 2: N'=N+-2
!----------------------------------------------------------------------
! This covers excitations from a doubly-occupied orbital into an
! unoccupied orbital or from a singly-occupied orbital into another
! singly-occupied orbital
!---------------------------------------------------------------------- 
    ! N=0 contribution
    if (imult == 1) then
       nspincp(2)=1
    else
       nspincp(2)=0
    endif

    ! Loop over numbers of open shells
    do nopen=1,nocase2
    
       ! Cycle if there are no CSFs for nopen open shells
       if (ncsfs(nopen) == 0) cycle
       
       ! Cycle if N+2 > nocase2
       if (nopen+2 > nocase2) cycle
    
       ! Number of unique pairs of permutations of N 1's and two 0's
       ! multiplied by the number of pairs of CSFs
       nspincp(2)=nspincp(2)+(nopen+2)*(nopen+1)/2 &
            *ncsfs(nopen)*ncsfs(nopen+2)

    enddo
    
    return
    
  end subroutine get_nunique

!######################################################################
! print_spincp_info: Printing of some information about the spin
!                    coupling coefficients being calculated
!######################################################################
  subroutine print_spincp_info(imult,nspincp,verbose)

    use constants
        
    implicit none

    integer(is), intent(in) :: imult
    integer(is), intent(in) :: nspincp(2)
    logical, intent(in)     :: verbose
    
    integer(is)             :: n,i

!----------------------------------------------------------------------
! Spin multiplicity
!----------------------------------------------------------------------
    if (verbose) write(6,'(/,x,a,x,i0)') 'Spin multiplicity:',imult
    
!----------------------------------------------------------------------
! Number of unique spin-coupling coefficients
!----------------------------------------------------------------------
    if (verbose) then
       
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

    endif
       
    return
    
  end subroutine print_spincp_info

!######################################################################
! init_spincp_arrays: Allocation and initialisation of the various
!                     spin coupling coefficient arrays
!######################################################################
  subroutine init_spincp_arrays(imult,nocase1,nocase2,ncsfs,nspincp,&
       npattern1,npattern2,verbose,spincp,spincpdim,offspincp)

    use constants
        
    implicit none

    integer(is), intent(in)  :: imult,nocase1,nocase2
    integer(is), intent(in)  :: ncsfs(0:nocase2)
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(out) :: npattern1,npattern2
    integer(is), intent(out) :: spincpdim(3)
    integer(is), intent(out) :: offspincp(4)
    real(dp), allocatable    :: spincp(:)
    logical, intent(in)      :: verbose
        
    integer(is)              :: nopen
    real(dp)                 :: mem
    
!----------------------------------------------------------------------
! Compute the number of unique Case 1 and 2 patterns
!----------------------------------------------------------------------
    !
    ! Case 1 patterns
    !
    npattern1=0
    ! Loop over numbers of open shells
    do nopen=1,nocase1
       ! Cycle if there are no CSFs for nopen open shells
       if (ncsfs(nopen) == 0) cycle
       ! Number of unique pairs of permutations of N 1's and one 0
       npattern1=npattern1+(nopen+1)*nopen/2
    enddo
    
    !
    ! Case 2 patterns
    !
    ! N=0 contribution
    if (imult == 1) then
       npattern2=1
    else
       npattern2=0
    endif
    ! Loop over numbers of open shells
    do nopen=1,nocase2
       ! Cycle if there are no CSFs for nopen open shells
       if (ncsfs(nopen) == 0) cycle
       ! Cycle if N+2 > nocase2
       if (nopen+2 > nocase2) cycle
       ! Number of unique pairs of permutations of N 1's and two 0's
       npattern2=npattern2+(nopen+2)*(nopen+1)/2
    enddo

!----------------------------------------------------------------------
! Numbers of spin coupling coefficients
!----------------------------------------------------------------------
    ! Case 1 coefficients: (i) 1a i>j, (ii) 1a i<j, (iii) 1b i>j,
    !                     (iv) 1b i<j
    spincpdim(1)=4*nspincp(1)

    ! Case 2 coefficients: (i) 2a, (ii) 2b
    spincpdim(2)=2*nspincp(2)

    ! Total number of coefficients
    spincpdim(3)=sum(spincpdim(1:2))

!----------------------------------------------------------------------
! Offsets for the various different spin coupling coefficients within
! the spincp array
!----------------------------------------------------------------------
    ! Case 1a, i<j
    offspincp(1)=nspincp(1)

    ! Case 1b, i>j
    offspincp(2)=2*nspincp(1)

    ! Case 1b, i<j
    offspincp(3)=3*nspincp(1)

    ! Case 2b
    offspincp(4)=nspincp(2)
    
!----------------------------------------------------------------------
! Spin coupling coefficient arrays
!----------------------------------------------------------------------    
    ! All spin coupling coefficients
    allocate(spincp(spincpdim(3)))
    spincp=0.0d0
    
!----------------------------------------------------------------------
! Output the memory used
!----------------------------------------------------------------------
    mem=spincpdim(3)*8/1024.0d0**2

    if (verbose) write(6,'(/,x,a,F8.2,x,a)') 'Memory used:',&
         mem,'MB'

    return
    
  end subroutine init_spincp_arrays

!######################################################################
! init_patternmap: Allocation and initialisation of the array holding
!                  the pattern value -> array index mapping
!######################################################################
  subroutine init_patternmap(imult,nocase1,nocase2,ncsfs,patternmap,&
       mapdim)

    use constants
    
    implicit none

    integer(is), intent(in)  :: imult
    integer(is), intent(in)  :: nocase1,nocase2
    integer(is), intent(in)  :: ncsfs(0:nocase2)
    integer(is), intent(out) :: mapdim
    integer(is), allocatable :: patternmap(:)

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
       mapdim=max(mapdim,2**(nocase2-1)-4)
    else
       mapdim=max(mapdim,2**nocase2-4)
    endif
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(patternmap(0:mapdim))
    patternmap=0

    return
    
  end subroutine init_patternmap
  
!######################################################################
! case1_coeffs: Calculation of the Case 1 spin coupling coefficients.
!               These are the coefficients corresponding to pairs of
!               initial and final spatial occupations with equal
!               numbers of open shells
!######################################################################
  subroutine case1_coeffs(nocase1,nocase2,maxcsf,maxdet,ncsfs,ndets,&
       csfcoe,detvec,npattern1,spincpdim,spincp,nspincp,mapdim,&
       patternmap)

    use constants
    use bitutils
        
    implicit none

    integer(is), intent(in)  :: nocase1,nocase2,maxcsf,maxdet
    integer(is), intent(in)  :: ncsfs(0:nocase2),ndets(0:nocase2)
    real(dp), intent(in)     :: csfcoe(maxcsf,maxdet,nocase2)
    integer(ib), intent(in)  :: detvec(maxdet,nocase2)
    integer(is), intent(in)  :: npattern1
    integer(is), intent(in)  :: spincpdim(3)
    real(dp), intent(out)    :: spincp(spincpdim(3))
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(in)  :: mapdim
    integer(is), intent(out) :: patternmap(0:mapdim)
    
    integer(is)              :: nopen,is1,is2,icsf1,icsf2,nsp,ndet
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
! Compute the Case 1 a spin coupling coefficients
! < w' omega' | E_i^j | w omega > for i>j
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=1
    
    ! Loop over numbers of open shells
    do nopen=1,nocase1

       ! No. CSFs and determinants
       nsp=ncsfs(nopen)
       ndet=ndets(nopen)
       
       ! Cycle if there are no CSFs for this number of open shells
       if (nsp == 0) cycle

       ! Allocate the spin-coupling coefficient work array
       allocate(work(nsp**2))
       
       ! Loop over the unique pairs of simplified spatial occupation
       ! vectors
       do is1=1,nopen
          do is2=is1+1,nopen+1

             ! Pattern number and the corresponding array index
             pattern=iand(ws(is1,nopen),ws(is2,nopen))
             patternmap(pattern)=ioff
             
             ! Evaluate the spin coupling coefficients for all
             ! CSF pairs
             call spincp_coeff_case1(ws(is1,nopen),ws(is2,nopen),&
                  nopen,nocase2,maxcsf,maxdet,nsp,ndet,detvec,csfcoe,&
                  work)
             
             ! Fill in the spincp array
             spincp(ioff:ioff+nsp**2-1)=work

             ! Update the spincp offset
             ioff=ioff+nsp**2
             
          enddo
       enddo

       ! Deallocate the spin-coupling coefficient work array
       deallocate(work)
       
    enddo

!----------------------------------------------------------------------
! Fill in the Case 1a i<j and Case 1b spin-coupling coefficients
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=1
    
    ! Loop over numbers of open shells
    do nopen=1,nocase1

       ! No. CSFs
       nsp=ncsfs(nopen)
       
       ! Cycle if there are no CSFs for this number of open shells
       if (nsp == 0) cycle

       ! Allocate work arrays
       allocate(work(nsp**2), work2(nsp,nsp), workT(nsp**2))
       
       ! Loop over the unique pairs of simplified spatial occupation
       ! vectors
       do is1=1,nopen
          do is2=is1+1,nopen+1

             ! Case 1a, i>j spin-coupling coefficient matrix and its
             ! transpose
             work=spincp(ioff:ioff+nsp**2-1)
             work2=transpose(reshape(work,(/nsp,nsp/)))
             workT=reshape(work2,(/nsp**2/))

             ! Case 1a, i<j coefficients
             istart=ioff+nspincp(1)
             iend=istart+nsp**2-1
             spincp(istart:iend)=workT
             
             ! Case 1b, i>j coefficients
             istart=istart+nspincp(1)
             iend=iend+nspincp(1)
             spincp(istart:iend)=-workT
             
             ! Case 1b, i<j coefficients
             istart=istart+nspincp(1)
             iend=iend+nspincp(1)
             spincp(istart:iend)=-work

             ! Update the spincp offset
             ioff=ioff+nsp**2
             
          enddo
       enddo

       ! Deallocate work arrays
       deallocate(work, workT, work2)
       
    enddo
       
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ws)
    
    return
    
  end subroutine case1_coeffs

!######################################################################
! spincp_coeff_case1: Computes a single Case 1 spin coupling
!                     coefficient
!######################################################################
  subroutine spincp_coeff_case1(ws1,ws2,nopen,nocase2,maxcsf,maxdet,&
       nsp,ndet,detvec,csfcoe,spincoe)
    
    use constants
    use bitutils
    use slater_condon
    
    implicit none

    integer(ib), intent(in)  :: ws1,ws2
    integer(is), intent(in)  :: nopen
    integer(is), intent(in)  :: nocase2,maxcsf,maxdet
    integer(is), intent(in)  :: nsp,ndet
    integer(ib), intent(in)  :: detvec(maxdet,nocase2)
    real(dp), intent(in)     :: csfcoe(maxcsf,maxdet,nocase2)
    real(dp), intent(out)    :: spincoe(nsp,nsp)
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
    real(dp), allocatable    :: coe(:,:)
    
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
    allocate(d1(2,ndet))
    allocate(d2(2,ndet))
    d1=0_ib
    d2=0_ib

    ! Phase factor matrix
    allocate(phasemat(ndet,ndet))
    phasemat=0.0d0

    ! Phase factor matrix contracted with the transpose of the CSF
    ! coefficient matrix
    allocate(pcT(ndet,nsp))
    pcT=0.0d0

    ! Working arrays
    allocate(coe(nsp,ndet))
    coe=0.0d0

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
       if (.not.btest(ws1,n-1)) cycle

       ! Increment the detvec open shell orbital counter
       vecindx=vecindx+1

       ! Loop over determinants
       do idet=1,ndet
          
          ! Add the next spin orbital
          if (btest(detvec(idet,nopen),vecindx-1)) then
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
       if (.not.btest(ws2,n-1)) cycle

       ! Increment the detvec open shell orbital counter
       vecindx=vecindx+1
       
       ! Loop over determinants
       do idet=1,ndet
          
          ! Add the next spin orbital
          if (btest(detvec(idet,nopen),vecindx-1)) then
             ! Occupied alpha spin-orbital
             d2(1,idet)=ibset(d2(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             d2(2,idet)=ibset(d2(2,idet),n-1)
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Determine the indices of the singlet excitation operator E_i^j
! Remember that these are given by the positions of the 0's in
! ws1 and ws2, i.e., the indices of the first unset bits
!----------------------------------------------------------------------
    ! Creation operator index
    ic=trailz(not(ws1))+1

    ! Annihilation operator index
    ia=trailz(not(ws2))+1

!----------------------------------------------------------------------
! Compute the matrix of phase factors for all pairs of determinants
!----------------------------------------------------------------------
! Note that phasemat(I,J) = < det_I | E_i^j | det_J >, where E_i^j is
! the singlet excitation operator for the spin coupling coefficients
! under consideration. We just use the name 'phasemat' because these
! elements are either zero or else equal to the phase factors between
! the determinants.
!----------------------------------------------------------------------
    ! Loop over ket determinants
    do iket=1,ndet

       ! Compute the phase mask
       call phasemask(d1(:,iket),phase_mask)
       
       ! Loop over bra determinants
       do ibra=1,ndet
          
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

          ! Cycle if the indices do not match those of the creation
          ! and annihilation operators
          if (plist(1,ialpha) /= 0) then
             ispin=1
          else
             ispin=2
          endif
          ip=plist(1,ispin)
          ih=hlist(1,ispin)
          if (ip /= ic .or. ih /= ia) cycle
          
          ! Phase factor
          !phase=phase_slow(d1(:,iket),hlist,plist,maxex,nexci)
          phase=phase_pure_exc(phase_mask(ispin),hlist(:,ispin),&
                plist(:,ispin),maxex,nexci)
          phasemat(ibra,iket)=dble(phase)
          
       enddo

    enddo
    
!----------------------------------------------------------------------
! Compute the matrix of spin coupling coefficients for the current
! pattern index
!----------------------------------------------------------------------
    ! Working array to avoid run-time warnings about temporary arrays
    ! being created
    coe=csfcoe(1:nsp,1:ndet,nopen)
    
    ! Contract the phase matrix with the CSF expansion coefficients,
    ! pcT = phasemat * transpose(csfcoe)
    call dgemm('N','T',ndet,nsp,ndet,1.0d0,phasemat,ndet,coe,nsp,&
         0.0d0,pcT,ndet)

    ! Contract the CSF expansion coefficients with the intermediate
    ! matrix pcT to yield the spin coupling coefficients
    call dgemm('N','N',nsp,nsp,ndet,1.0d0,coe,nsp,pcT,ndet,0.0d0,&
         spincoe,nsp)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(d1)
    deallocate(d2)
    deallocate(phasemat)
    deallocate(pcT)
    deallocate(coe)

!----------------------------------------------------------------------
! Reset n_int
!----------------------------------------------------------------------
    n_int=n_int_save
    
    return
    
  end subroutine spincp_coeff_case1
  
!######################################################################
! case2_coeffs: Calculation of the Case 2 spin coupling coefficients.
!               These are the coefficients corresponding to pairs of
!               initial and final spatial occupations with numbers of
!               open shells differing by two.
!######################################################################
  subroutine case2_coeffs(imult,nocase2,maxcsf,maxdet,ncsfs,ndets,&
       csfcoe,detvec,spincpdim,spincp,nspincp,mapdim,patternmap)

    use constants
    use bitutils
        
    implicit none

    integer(is), intent(in)  :: imult,nocase2,maxcsf,maxdet
    integer(is), intent(in)  :: ncsfs(0:nocase2),ndets(0:nocase2)
    real(dp), intent(in)     :: csfcoe(maxcsf,maxdet,nocase2)
    integer(ib), intent(in)  :: detvec(maxdet,nocase2)
    integer(is), intent(in)  :: spincpdim(3)
    real(dp), intent(out)    :: spincp(spincpdim(3))
    integer(is), intent(in)  :: nspincp(2)
    integer(is), intent(in)  :: mapdim
    integer(is), intent(out) :: patternmap(0:mapdim)
    
    integer(is)              :: nopen,is1,is2,icsf1,icsf2,lim,i
    integer(is)              :: ioff,istart,iend,nsp1,nsp2
    integer(ib), allocatable :: ws(:,:)
    integer(ib)              :: wsp(nocase2)
    integer(ib)              :: pattern
    real(dp)                 :: coeff
    real(dp), allocatable    :: work(:),workT(:),work2(:,:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Simplified spatial occupation vectors
    allocate(ws((nocase2+2)*(nocase2+1)/2,0:nocase2))
    ws=0_ib

!----------------------------------------------------------------------
! Generate the Case 2 simplified final spatial occupation vectors,
! i.e., the vectors of all 1's
!----------------------------------------------------------------------
    ! Initialise to all unset bits
    wsp=0_ib

    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Construct the vector of nopen 1's
       do i=1,nopen
          wsp(nopen)=ibset(wsp(nopen),i-1)
       enddo

    enddo
    
!----------------------------------------------------------------------
! Generate the unique Case 2 simplified initial spatial occupation
! vectors, i.e, all possible permutations of N 1's and two 0's
! for N=1,...,nocase2
!----------------------------------------------------------------------
    !
    ! N = 0 case: one simplified vector corresponding to 00
    !
    ws(1,0)=0_ib
    
    !
    ! N > 0 cases
    !
    ! Loop over numbers of open shells
    do nopen=1,nocase2
       
       ! Cycle if there are no CSFs for nopen open shells
       if (ncsfs(nopen) == 0) cycle
       
       ! Cycle if N+2 > nocase2
       if (nopen+2 > nocase2) cycle

       ! Generate the permutations
       lim=(nopen+2)*(nopen+1)/2
       call get_permutations(nopen,2,ws(1:lim,nopen),lim)
              
    enddo
    
!----------------------------------------------------------------------
! Compute the Case 2 spin coupling coefficient for N = 0 (only exists
! for singlets)
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=spincpdim(1)+1
    
    ! N = 0 spin coupling coefficient
    if (imult == 1) then

       ! Pattern number and corresponding array index
       pattern=iand(ws(1,0),wsp(2))
       patternmap(pattern)=ioff
       
       ! Evaluate the spin coupling coefficient
       call spincp_coeff_n0(coeff)
       spincp(ioff)=coeff

       ! Update the spincp offset
       ioff=ioff+1
       
    endif

!----------------------------------------------------------------------
! Compute the Case 2a spin coupling coefficients for N > 0
!----------------------------------------------------------------------
    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Cycle if there are no CSFs for nopen open shells
       if (ncsfs(nopen) == 0) cycle
       
       ! Cycle if N+2 > nocase2
       if (nopen+2 > nocase2) cycle

       ! Allocate the spin-coupling coefficient work array
       nsp1=ncsfs(nopen)
       nsp2=ncsfs(nopen+2)
       allocate(work(nsp1*nsp2))
       
       ! Loop over the simplified spatial occupation vectors
       do is1=1,(nopen+2)*(nopen+1)/2

          ! Pattern number and the corresponding array indx
          pattern=iand(ws(is1,nopen),wsp(nopen+2))
          patternmap(pattern)=ioff

          ! Evaluate the spin coupling coefficient for all CSF
          ! pairs
          call spincp_coeff_case2(ws(is1,nopen),wsp(nopen+2),nopen,&
               nopen+2,nocase2,maxcsf,maxdet,ndets,detvec,&
               csfcoe,nsp1,nsp2,work)

          ! Fill in the spincp array
          spincp(ioff:ioff+nsp1*nsp2-1)=work
          
          ! Update the spincp offset
          ioff=ioff+nsp1*nsp2
          
       enddo

       ! Deallocate the spin-coupling coefficient work array
       deallocate(work)
       
    enddo

!----------------------------------------------------------------------
! Fill in the Case 2b spin coupling coefficients
!----------------------------------------------------------------------
    ! Initialise the spincp offset
    ioff=spincpdim(1)+1

    ! N = 0 case
    if (imult == 1) then
       spincp(ioff+nspincp(2))=spincp(ioff)
       ioff=ioff+1
    endif

    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Cycle if there are no CSFs for nopen open shells
       if (ncsfs(nopen) == 0) cycle
       
       ! Cycle if N+2 > nocase2
       if (nopen+2 > nocase2) cycle

       ! Allocate work arrays
       nsp1=ncsfs(nopen)
       nsp2=ncsfs(nopen+2)
       allocate(work(nsp1*nsp2), workT(nsp1*nsp2), work2(nsp1,nsp2))
       
       ! Loop over the simplified spatial occupation vectors
       do is1=1,(nopen+2)*(nopen+1)/2

          ! Transpose of the Case 2a spin-coupling coefficient matrix
          work=spincp(ioff:ioff+nsp1*nsp2-1)
          work2=transpose(reshape(work, (/nsp2,nsp1/)))
          workT=reshape(work2, (/nsp1*nsp2/))
          
          ! Case 2b coefficients
          istart=ioff+nspincp(2)
          iend=istart+nsp1*nsp2-1
          spincp(istart:iend)=workT
          
          ! Update the spincp offset
          ioff=ioff+nsp1*nsp2
          
       enddo

       ! Deallocate work arrays
       deallocate(work, workT, work2)
       
    enddo

    return
    
  end subroutine case2_coeffs

!######################################################################
! spincp_coeff_case2: Computes a single Case 2 spin coupling
!                     coefficient
!######################################################################
! Note that here ws1 is the simplified occupation vector with
! N=nopen1 open shells, and ws2 is the simplified occupation vector
! with N+2=nopen2 open shells
!######################################################################
  subroutine spincp_coeff_case2(ws1,ws2,nopen1,nopen2,nocase2,maxcsf,&
       maxdet,ndets,detvec,csfcoe,ncsf1,ncsf2,spincoe)

    use constants
    use bitutils
    use slater_condon
    
    implicit none

    integer(ib), intent(in)  :: ws1,ws2
    integer(is), intent(in)  :: nopen1,nopen2
    integer(is), intent(in)  :: nocase2,maxcsf,maxdet
    integer(is), intent(in)  :: ndets(0:nocase2)
    integer(ib), intent(in)  :: detvec(maxdet,nocase2)
    real(dp), intent(in)     :: csfcoe(maxcsf,maxdet,nocase2)
    integer(is), intent(in)  :: ncsf1,ncsf2
    real(dp), intent(out)    :: spincoe(ncsf1,ncsf2)
    
    integer(is)              :: icsf1,icsf2
    integer(is)              :: idet,n,vecindx,nunset
    integer(is)              :: iket,ibra,ispin
    integer(is)              :: ic,ia,ih,ip
    integer(ib)              :: b
    integer(ib), allocatable :: d1(:,:),d2(:,:)
    integer(is)              :: n_int_save
    integer(ib)              :: phase_mask(2)
    integer(is)              :: nexci
    integer(ib)              :: p(2),h(2)
    integer(is), parameter   :: maxex=2
    integer(is)              :: plist(maxex,2),hlist(maxex,2)
    integer(is)              :: phase
    integer(is)              :: ndet1,ndet2
    real(dp), allocatable    :: phasemat(:,:),pcT(:,:)
    real(dp), allocatable    :: coe1(:,:),coe2(:,:)


    ! Phase factor accociated with creating the doubly-occupied
    ! MO in the ket determinants from the determinants with
    ! only unoccupied MOs
    integer(is)              :: docc_phase
    integer(is)              :: n2,n2a,n2b,i
    integer(is)              :: noa
    
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
    ndet1=ndets(nopen1)
    ndet2=ndets(nopen2)

    ! Determinant bit strings
    allocate(d1(2,ndet1))
    allocate(d2(2,ndet2))
    d1=0_ib
    d2=0_ib

    ! Phase factor matrix
    allocate(phasemat(ndet2,ndet1))
    phasemat=0.0d0

    ! Phase factor matrix contracted with the transpose of the CSF
    ! coefficient matrix
    allocate(pcT(ndet2,ncsf1))
    pcT=0.0d0
    
    ! Working arrays
    allocate(coe1(ncsf1,ndet1))
    coe1=0.0d0
    allocate(coe2(ncsf2,ndet2))
    coe2=0.0d0
    
!----------------------------------------------------------------------
! Generate the determinants for CSF 1 (the ket CSF, with N singly-
! occupied MOs, one doubly occupied MO, and one unoccupied MO)
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    d1=0_ib

    ! Initialise the detvec open shell counter
    vecindx=0

    ! Initialise the unset bit counter
    nunset=0

    ! Loop over orbitals in the spatial occupation
    do n=1,nopen1+2

       ! If we are at the first unset bit, then treat it as a
       ! doubly-occupied MO, else treat it as an unoccupied MO
       if (.not.btest(ws1,n-1)) then
          nunset=nunset+1
          if (nunset == 1) then
             ! Position of the doubly-occupied MO
             n2=n
             ! Doubly-occupied MO: set both alpha and beta string bits
             ! in all determinants
             do idet=1,ndet1
                d1(1,idet)=ibset(d1(1,idet),n-1)
                d1(2,idet)=ibset(d1(2,idet),n-1)
             enddo             
             cycle
          else
             ! Unoccupied MO: cycle
             cycle
          endif
       endif
          
       ! Increment the detvec open shell orbital counter
       vecindx=vecindx+1

       ! Loop over determinants
       do idet=1,ndet1
          
          ! Add the next spin orbital
          if (btest(detvec(idet,nopen1),vecindx-1)) then
             ! Occupied alpha spin-orbital
             d1(1,idet)=ibset(d1(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             d1(2,idet)=ibset(d1(2,idet),n-1)
          endif

       enddo

    enddo
    
!----------------------------------------------------------------------
! Generate the determinants for CSF 2 (the bra CSF, with N+2 singly-
! occupied MOs)
!----------------------------------------------------------------------
    ! Initialise the bit string to all unset bits
    d2=0_ib
    
    ! Loop over determinants
    do idet=1,ndet2
    
       ! Loop over orbitals in the simplified spatial occupation
       do n=1,nopen2

          ! Add the next spin orbital
          if (btest(detvec(idet,nopen2),n-1)) then
             ! Occupied alpha spin-orbital
             d2(1,idet)=ibset(d2(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             d2(2,idet)=ibset(d2(2,idet),n-1)
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
    noa=popcnt(d1(1,1))
    
    ! Number of electrons before the index of the alpha-spin
    ! creation operator
    n2a=0
    if (n2 > 1) then
       do i=1,n2-1
          if (btest(d1(1,1),i-1)) n2a=n2a+1
       enddo
    endif
    
    ! Number of electrons before the index of the beta-spin
    ! creation operator
    n2b=noa+1
    if (n2 > 1) then
       do i=1,n2-1
          if (btest(d1(2,1),i-1)) n2b=n2b+1
       enddo
    endif
    
    ! Phase factor
    docc_phase=(-1)**(n2a+n2b)
    
!----------------------------------------------------------------------
! Determine the indices of the singlet excitation operator E_i^j
! Remember that these are given by the positions of the two 0's in
! ws1, i.e., the indices of the first and second unset bits.
!----------------------------------------------------------------------
! Important: We are taking the first zero in ws1 to correspond to the
!            doubly-occupied MO. If we ever change this, then the
!            creation operator index will have to be determined before
!            the annihilation operator index due to the use of the
!            bit clearing operation.
!----------------------------------------------------------------------
    ! Temporary bit string array
    b=not(ws1)
    
    ! Annihilation operator index
    ia=trailz(b)+1
    
    ! Creation operator index
    b=ibclr(b,ia-1)
    ic=trailz(b)+1
    
!----------------------------------------------------------------------
! Compute the matrix of phase factors for all pairs of determinants
!----------------------------------------------------------------------
! Note that phasemat(I,J) = < det_I | E_i^j | det_J >, where E_i^j is
! the singlet excitation operator for the spin coupling coefficients
! under consideration. We just use the name 'phasemat' because these
! elements are either zero or else equal to the phase factors between
! the determinants.
!----------------------------------------------------------------------
    ! Initialisation
    phasemat=0.0d0
    
    ! Loop over ket determinants
    do iket=1,ndet1

       ! Compute the phase mask
       call phasemask(d1(:,iket),phase_mask)

       ! Loop over bra determinants
       do ibra=1,ndet2

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

          ! Cycle if the indices do not match those of the creation
          ! and annihilation operators
          if (plist(1,ialpha) /= 0) then
             ispin=1
          else
             ispin=2
          endif
          ip=plist(1,ispin)
          ih=hlist(1,ispin)
          if (ip /= ic .or. ih /= ia) cycle
          
          ! Phase factor
          !phase=phase_slow(d1(:,iket),hlist,plist,maxex,nexci)
          phase=phase_pure_exc(phase_mask(ispin),hlist(:,ispin),&
                plist(:,ispin),maxex,nexci)
          phasemat(ibra,iket)=dble(phase)
          
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Compute the matrix of spin coupling coefficients for the current
! pattern index
!----------------------------------------------------------------------
    ! Working arrays to avoid run-time warnings about temporary arrays
    ! being created
    coe1=csfcoe(1:ncsf1,1:ndet1,nopen1)
    coe2=csfcoe(1:ncsf2,1:ndet2,nopen2)

    ! Multiply the ket CSF expansion coefficients by the phase factors
    ! associated with creating the doubly-occupied MO
    coe1=coe1*docc_phase
    
    ! Contract the phase matrix with the CSF 1 expansion coefficients,
    ! pcT = phasemat * transpose(csfcoe2)
    call dgemm('N','T',ndet2,ncsf1,ndet1,1.0d0,phasemat,ndet2,coe1,&
         ncsf1,0.0d0,pcT,ndet2)
    
    ! Contract the CSF 2 expansion coefficients with the intermediate
    ! matrix pcT to yield the spin coupling coefficients
    call dgemm('N','N',ncsf2,ncsf1,ndet2,1.0d0,coe2,ncsf2,pcT,ndet2,&
         0.0d0,spincoe,ncsf2)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(d1)
    deallocate(d2)
    deallocate(phasemat)
    deallocate(pcT)
    deallocate(coe1)
    deallocate(coe2)

!----------------------------------------------------------------------
! Reset n_int
!----------------------------------------------------------------------
    n_int=n_int_save
    
    return
    
  end subroutine spincp_coeff_case2

!######################################################################
! spincp_coeff_n0: Computes the N=0 Case 2 spin coupling coefficient.
!                  Note that this could, of course, be hard-coded, but
!                  that this served as a sanity check when writing
!                  the code and was just kept
!######################################################################
  subroutine spincp_coeff_n0(coeff)

    use constants
    use bitglobal
    use bitutils
    use slater_condon
    
    implicit none

    integer(ib)            :: d1(2)
    integer(ib)            :: d2(2,2)
    integer(is)            :: n_int_save
    integer(is)            :: idet,ibra,ispin,n
    integer(is)            :: ic,ia,ip,ih
    integer(ib)            :: phase_mask(2)
    integer(is)            :: nexci
    integer(ib)            :: p(2),h(2)
    integer(is), parameter :: maxex=2
    integer(is)            :: plist(maxex,2),hlist(maxex,2)
    integer(is)            :: phase
    real(dp), intent(out)  :: coeff
    
!----------------------------------------------------------------------
! Save the actual value of n_int and then set this to 1 for use in the
! following. This allows us to use the slater_condon module to
! calculate the spin coupling coefficients.
!----------------------------------------------------------------------
    n_int_save=n_int
    n_int=1
    
!----------------------------------------------------------------------
! Construct the determinants
!----------------------------------------------------------------------
    ! Determinant for the ket CSF
    d1=0_ib
    d1(1)=ibset(d1(1),0)
    d1(2)=ibset(d1(2),0)

    ! Determinants for the bra CSF
    d2=0_ib
    do idet=1,2

       ! Loop over orbitals
       do n=1,2

          ! Add the next spin-orbital
          if (btest(detvec(idet,2),n-1)) then
             ! Occupied alpha spin-orbital
             d2(1,idet)=ibset(d2(1,idet),n-1)
          else
             ! Occupied beta spin-orbital
             d2(2,idet)=ibset(d2(2,idet),n-1)
          endif

       enddo

    enddo

!----------------------------------------------------------------------
! Compute the spin coupling coefficient
!----------------------------------------------------------------------
    ! Initialisation
    coeff=0.0d0

    ! Creation and annihilation operator indices
    ia=1
    ic=2

    ! Compute the phase mask
    call phasemask(d1(:),phase_mask)

    ! Excitation degree
    nexci=1
    
    ! Loop over bra determinants
    do ibra=1,2

       ! Get the indices of the spin-orbitals involved in the
       ! excitations linking the bra and ket determinants
       call exc(d1(:),d2(:,ibra),p,h)
       call list_from_bitstring(p(ialpha),plist(:,ialpha),maxex)
       call list_from_bitstring(h(ialpha),hlist(:,ialpha),maxex)
       call list_from_bitstring(p(ibeta),plist(:,ibeta),maxex)
       call list_from_bitstring(h(ibeta),hlist(:,ibeta),maxex)

       ! Phase factor
       if (plist(1,ialpha) /= 0) then
          ispin=1
       else
          ispin=2
       endif
       phase=phase_pure_exc(phase_mask(ispin),hlist(:,ispin),&
                plist(:,ispin),maxex,nexci)

       ! Contribution to the spin coupling coefficient
       coeff=coeff+phase*csfcoe(1,ibra,2)
       
    enddo

!----------------------------------------------------------------------
! Reset n_int
!----------------------------------------------------------------------
    n_int=n_int_save
    
    return
    
  end subroutine spincp_coeff_n0
  
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
  
end module spin_coupling
