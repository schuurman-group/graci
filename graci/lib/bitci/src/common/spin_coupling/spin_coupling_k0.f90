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

!----------------------------------------------------------------------
! Compute the Case 2 spin coupling coefficients
!----------------------------------------------------------------------
    !call case2_coeffs_k0
    
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
! < w' omega' | T_ij^(1,k=0) | w omega > for i>j
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
             patmap(pattern)=ioff
             
             ! Compute the 1a i>j spin coupling coefficients
             call spincp_1a_k0(ws(is1,nopen),ws(is2,nopen),nopen,&
                  nocase2,maxcsf,maxdet,nsp,ndet,detvec,csfcoe,work)
                  
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
! Fill in the remaining Case 1a k=0 triplet spin coupling coefficients
! using the following relations:
!
! [1a i<j] = [1a i>j]^T
! [1b i>j] = [1a i>j]^T
! [1b i<j] = [1a i>j]
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
             spincp(istart:iend)=workT
             
             ! Case 1b, i<j coefficients
             istart=istart+nspincp(1)
             iend=iend+nspincp(1)
             spincp(istart:iend)=work

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
    
  end subroutine case1_coeffs_k0

!######################################################################
! spincp_1a_k0: Computes a single batch of Case 1a spin coupling
!               coefficients for the k=0 component of the triplet spin
!               tensor operator
!######################################################################
  subroutine spincp_1a_k0(ws1,ws2,nopen,nocase2,maxcsf,maxdet,nsp,&
       ndet,detvec,csfcoe,spincoe)

    use constants
    use bitutils
    use slater_condon
    
    implicit none

    integer(ib), intent(in)           :: ws1,ws2
    integer(is), intent(in)           :: nopen
    integer(is), intent(in)           :: nocase2,maxcsf,maxdet
    integer(is), intent(in)           :: nsp,ndet
    integer(ib), intent(in)           :: detvec(maxdet,nocase2)
    real(dp), intent(in)              :: csfcoe(maxcsf,maxdet,nocase2)
    real(dp), intent(out)             :: spincoe(nsp,nsp)
    integer(is)                       :: icsf1,icsf2
    integer(is)                       :: idet,n,vecindx,iket,ibra,ispin
    integer(is)                       :: ic,ia,ih,ip
    integer(ib), allocatable          :: d1(:,:),d2(:,:)
    integer(is)                       :: n_int_save
    integer(ib)                       :: phase_mask(2)
    integer(is)                       :: nexci
    integer(ib)                       :: p(2),h(2)
    integer(is), parameter            :: maxex=2
    integer(is)                       :: plist(maxex,2),hlist(maxex,2)
    integer(is)                       :: phase
    real(dp), allocatable             :: phasemat(:,:),pcT(:,:)
    real(dp), allocatable             :: coe(:,:)

    ! Spin factors                             alpha   beta
    real(dp), parameter, dimension(2) :: sfac=[1.0d0, -1.0d0]

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
       if (.not. btest(ws1,n-1)) cycle

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
       if (.not. btest(ws2,n-1)) cycle

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
! Determine the indices of the triplet excitation operator T_ij^(1,0)
! Remember that these are given by the positions of the 0's in
! ws1 and ws2, i.e., the indices of the first unset bits
!----------------------------------------------------------------------
    ! Creation operator index
    ic=trailz(not(ws1))+1

    ! Annihilation operator index
    ia=trailz(not(ws2))+1

!----------------------------------------------------------------------
! Compute the matrix of phase factors (multiplied by the appropriate
! spin factors) for all pairs of determinants for the triplet
! excitation operator
!
! T_ic,ia^(1,0) = 1/sqrt(2) (a_ic,alpha^\dagger a_ia,alpha
!                           -a_ic,beta^\dagger a_ia,beta)
!----------------------------------------------------------------------
! phasemat(I,J) = < det_I | T_ic,ia^(1,0) | det_J >
!----------------------------------------------------------------------
    ! Loop over ket determinants
    do iket=1,ndet

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

          ! Spin of the created/annihilated electrons
          if (plist(1,ialpha) /= 0) then
             ispin=1
          else
             ispin=2
          endif
          
          ! Cycle if the indices do not match those of the creation
          ! and annihilation operators
          ip=plist(1,ispin)
          ih=hlist(1,ispin)
          if (ip /= ic .or. ih /= ia) cycle

          ! Product of the phase and spin factors
          phase=phase_slow(d1(:,iket),hlist,plist,maxex,nexci)
          phasemat(ibra,iket)=dble(phase)*sfac(ispin)

       enddo

    enddo

    ! 1/sqrt(2) prefactor
    phasemat=phasemat/sqrt(2.0d0)

!----------------------------------------------------------------------
! Compute the matrix of spin coupling coefficients for the current
! pattern index
!----------------------------------------------------------------------
    ! Working array
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
    
  end subroutine spincp_1a_k0
    
!######################################################################
  
end module spin_coupling_k0
