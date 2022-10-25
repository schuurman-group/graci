!######################################################################
! dyson: Dyson orbital interface with Fortran bindings
!######################################################################
subroutine dyson(nmoL,nmoR,n_intL,n_intR,ndetL,ndetR,&
     nrootsL,nrootsR,detL,detR,vecL,vecR,smo1,normthrsh,&
     npairs,n_basis,dysorb,ipairs,verbose1)

  use constants
  use global
  use detfuncs
  use timing
  
  ! No. MOs
  integer(is), intent(in)     :: nmoL,nmoR

  ! No. (n_bits)-bit integers needed to represent each alpha/beta
  ! string
  integer(is), intent(in)     :: n_intL,n_intR
  
  ! No. determinants
  integer(is), intent(in)     :: ndetL,ndetR

  ! No. roots
  integer(is), intent(in)     :: nrootsL,nrootsR

  ! Untruncated determinant bit strings
  integer(ib), intent(in)     :: detL(n_intL,2,ndetL)
  integer(ib), intent(in)     :: detR(n_intR,2,ndetR)

  ! Untruncated eigenvectors
  real(dp), intent(in)        :: vecL(ndetL,nrootsL)
  real(dp), intent(in)        :: vecR(ndetR,nrootsR)

  ! Matrix of bra-ket MO overlaps
  real(dp), intent(in)        :: smo1(nmoL,nmoR)
  
  ! Norm-based truncation threshold
  real(dp), intent(in)        :: normthrsh
  
  ! Indices of the pairs of states for which overlaps are requested
  integer(is), intent(in)     :: npairs
  integer(is), intent(in)     :: ipairs(npairs,2)

  ! Dyson orbitals
  integer(is), intent(in)     :: n_basis
  real(dp), intent(inout)     :: dysorb(n_basis,npairs)

  ! Verbosity flag
  logical, intent(in)         :: verbose1

  ! Numbers of electrons
  integer(is)                 :: nelL,nel_alphaL,nel_betaL
  integer(is)                 :: nelR,nel_alphaR,nel_betaR
  
  ! Transposition of the input bra and ket WFs
  logical                     :: transposed

  ! Are we annihilating alpha- or beta-spin electrons?
  !
  ! ap_spin = 1 <-> alpha-spin annihilation
  !         = 2 <-> beta-spin annihilation
  !
  integer(is)                 :: ap_spin
  character(len=5), parameter :: ap_label(2)=['alpha', 'beta ']
  
  ! Timing variables
  real(dp)                    :: tcpu_start,tcpu_end,twall_start,&
                                 twall_end

  ! Everything else
  integer(is)                 :: i
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! Make sure that all globally accessible allocatable arrays are
! not allocated
!----------------------------------------------------------------------
  if (allocated(smo))       deallocate(smo)
  if (allocated(detB))      deallocate(detB)
  if (allocated(detK))      deallocate(detK)
  if (allocated(vecB))      deallocate(vecB)
  if (allocated(vecK))      deallocate(vecK)
!  if (allocated(alphaB))    deallocate(alphaB)
!  if (allocated(betaB))     deallocate(betaB)
!  if (allocated(alphaK))    deallocate(alphaK)
!  if (allocated(betaK))     deallocate(betaK)
!  if (allocated(offsetB))   deallocate(offsetB)
!  if (allocated(offsetK))   deallocate(offsetK)
!  if (allocated(det2betaB)) deallocate(det2betaB)
!  if (allocated(det2betaK)) deallocate(det2betaK)
!  if (allocated(betafac))   deallocate(betafac)

!----------------------------------------------------------------------  
! Determine the numbers of electrons
!----------------------------------------------------------------------  
  call get_nel(n_intL,detL(:,:,1),nelL,nel_alphaL,nel_betaL)
  call get_nel(n_intR,detR(:,:,1),nelR,nel_alphaR,nel_betaR)

!----------------------------------------------------------------------  
! Set the verbosity flag
!----------------------------------------------------------------------  
  verbose=verbose1
  
!----------------------------------------------------------------------  
! First, some book kepping
!
! The input variables ending in R/L correspond to the WFs that the
! user wants in the bra/ket positions
!
! However, we require the working bra and ket arrays to correspond to
! the (N-1)- and N-electron WFs, respectively
!
! Hence, there may be a bit of shuffling here
!----------------------------------------------------------------------
  if (nelL == nelR - 1) then
     !
     ! L<->bra and R<->ket
     !

     ! Set the bra-ket transposition flag
     transposed=.false.
     
     ! Numbers of electrons
     nelB=nelL
     nel_alphaB=nel_alphaL
     nel_betaB=nel_betaL
     nelK=nelR
     nel_alphaK=nel_alphaR
     nel_betaK=nel_betaR

     ! Numbers of MOs
     nmoB=nmoL
     nmoK=nmoR

     ! Numbers of (n_bits)-bit integers needed to represent each
     ! determinant
     n_intB=n_intL
     n_intK=n_intR

     ! Numbers of roots
     nrootsB=nrootsL
     nrootsK=nrootsR

     ! Bra-ket MO overlap matrix
     allocate(smo(nmoB,nmoK))
     smo=smo1

  else if (nelL == nelR + 1) then
     !
     ! R<->bra and L<->ket
     !

     ! Set the bra-ket transposition flag
     transposed=.true.

     ! Numbers of electrons
     nelB=nelR
     nel_alphaB=nel_alphaR
     nel_betaB=nel_betaR
     nelK=nelL
     nel_alphaK=nel_alphaL
     nel_betaK=nel_betaL

     ! Numbers of MOs
     nmoB=nmoR
     nmoK=nmoL

     ! Numbers of (n_bits)-bit integers needed to represent each
     ! determinant
     n_intB=n_intR
     n_intK=n_intL

     ! Numbers of roots
     nrootsB=nrootsR
     nrootsK=nrootsL

     ! Bra-ket MO overlap matrix
     allocate(smo(nmoB,nmoK))
     smo=transpose(smo1)

  else
     !
     ! L and R WFs differ by 0 or >1 electrons: die here
     !
     write(6,'(/,x,a)') 'Error in dyson: |nel_bra - nel_ket| != 1'
     stop
     
  endif

!----------------------------------------------------------------------
! Truncate the wave functions
!----------------------------------------------------------------------
  if (.not. transposed) then
     !
     ! L<->bra and R<->ket
     !
     
     ! Bra
     call truncate_wave_functions(n_intL,ndetL,nrootsL,detL,vecL,&
          normthrsh,ndetB,detB,vecB)
     
     ! Ket
     call truncate_wave_functions(n_intR,ndetR,nrootsR,detR,vecR,&
          normthrsh,ndetK,detK,vecK)
     
     ! Ouput the numbers of determinants
     if (verbose) then
        write(6,'(/,x,a,x,i0)') &
             'No. bra determinants before truncation:',ndetL
        write(6,'(x,a,x,i0)') &
             'No. ket determinants before truncation:',ndetR
        write(6,'(/,x,a,x,F10.7)') 'Truncation threshold:',normthrsh
        write(6,'(/,x,a,x,i0)') &
             'No. bra determinants after truncation:',ndetB
        write(6,'(x,a,x,i0)') &
             'No. ket determinants after truncation:',ndetK
     endif
     
  else
     !
     ! R<->bra and L<->ket
     !

     ! Bra
     call truncate_wave_functions(n_intR,ndetR,nrootsR,detR,vecR,&
          normthrsh,ndetB,detB,vecB)
     
     ! Ket
     call truncate_wave_functions(n_intL,ndetL,nrootsL,detL,vecL,&
          normthrsh,ndetK,detK,vecK)

     ! Ouput the numbers of determinants
     if (verbose) then
        write(6,'(/,x,a,x,i0)') &
             'No. bra determinants before truncation:',ndetR
        write(6,'(x,a,x,i0)') &
             'No. ket determinants before truncation:',ndetL
        write(6,'(/,x,a,x,F10.7)') 'Truncation threshold:',normthrsh
        write(6,'(/,x,a,x,i0)') &
             'No. bra determinants after truncation:',ndetB
        write(6,'(x,a,x,i0)') &
             'No. ket determinants after truncation:',ndetK
     endif
     
  endif

!----------------------------------------------------------------------
! Symmetric orthogonalisation the truncated wave functions
!----------------------------------------------------------------------
  ! Bra
  call symm_ortho(n_intB,ndetB,nrootsB,vecB)
  
  ! Ket
  call symm_ortho(n_intK,ndetK,nrootsK,vecK)

!----------------------------------------------------------------------
! Output the numbers of electrons
!----------------------------------------------------------------------
  if (verbose) then
     write(6,'(/,x,28a)') ('-', i=1,28)
     write(6,'(7x,a)') 'N_el  | alpha | beta'
     write(6,'(x,28a)') ('-', i=1,28)
     write(6,'(x,a,x,i5,x,a,x,i5,x,a,x,i5)') &
          'bra |',nelB,'|',nel_alphaB,'|',nel_betaB
     write(6,'(x,a,x,i5,x,a,x,i5,x,a,x,i5)') &
          'ket |',nelK,'|',nel_alphaK,'|',nel_betaK
     write(6,'(x,28a)') ('-', i=1,28)
  endif
  
!----------------------------------------------------------------------
! Determine whether we are annihilating alpha- or beta-spin ket
! electrons
!----------------------------------------------------------------------
  if (nel_alphaK > nel_alphaB) then
     ! Annihilation of alpha-spin electrons
     ap_spin=1
  else
     ! Annihilation of beta-spin electrons
     ap_spin=2
  endif

  write(6,'(/,x,a)') 'Annihilation of '//trim(ap_label(ap_spin)) &
       //' electrons'
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  if (verbose) &
       call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'dyson')
  
  STOP
  
  return
  
end subroutine dyson
