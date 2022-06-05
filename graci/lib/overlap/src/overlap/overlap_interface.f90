!**********************************************************************
! Public-facing interface to the overlap library
!**********************************************************************
subroutine overlap(nmoB1,nmoK1,n_intB1,n_intK1,ndetB1,ndetK1,nrootsB1,&
     nrootsK1,detB1,detK1,vecB1,vecK1,smo1,normthrsh,&
     npairs,Sij,ipairs)

  use constants
  use global
  use detfuncs
  use detsort
  use factors
  use wfoverlap
  use timing
  
  implicit none

  ! No. MOs
  integer(is), intent(in) :: nmoB1,nmoK1

  ! No. (n_bits)-bit integers needed to represent each alpha/beta
  ! string
  integer(is), intent(in) :: n_intB1,n_intK1
  
  ! No. determinants
  integer(is), intent(in) :: ndetB1,ndetK1

  ! No. roots
  integer(is), intent(in) :: nrootsB1,nrootsK1
  
  ! Untruncated determinant bit strings
  integer(ib), intent(in) :: detB1(n_intB1,2,ndetB1)
  integer(ib), intent(in) :: detK1(n_intK1,2,ndetK1)

  ! Untruncated eigenvectors
  real(dp), intent(in)    :: vecB1(ndetB1,nrootsB1)
  real(dp), intent(in)    :: vecK1(ndetK1,nrootsK1)

  ! Matrix of bra-ket MO overlaps
  real(dp), intent(in)    :: smo1(nmoB1,nmoK1)
  
  ! Norm-based truncation threshold
  real(dp), intent(in)    :: normthrsh

  ! Indices of the pairs of states for which overlaps are requested
  integer(is), intent(in) :: npairs
  integer(is), intent(in) :: ipairs(npairs,2)

  ! Wave function overlaps
  real(dp), intent(out)   :: Sij(npairs)

  ! Timing variables
  real(dp)                :: tcpu_start,tcpu_end,twall_start,twall_end

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
  if (allocated(alphaB))    deallocate(alphaB)
  if (allocated(betaB))     deallocate(betaB)
  if (allocated(alphaK))    deallocate(alphaK)
  if (allocated(betaK))     deallocate(betaK)
  if (allocated(offsetB))   deallocate(offsetB)
  if (allocated(offsetK))   deallocate(offsetK)
  if (allocated(det2betaB)) deallocate(det2betaB)
  if (allocated(det2betaK)) deallocate(det2betaK)
  if (allocated(betafac))   deallocate(betafac)
  
!----------------------------------------------------------------------
! Set some globally accessible variables
!----------------------------------------------------------------------
  ! No. MOs
  nmoB=nmoB1
  nmoK=nmoK1

  ! No. integers needed to represent each alpha- and beta-string
  n_intB=n_intB1
  n_intK=n_intK1

  ! No. states
  nrootsB=nrootsB1
  nrootsK=nrootsK1

  ! Bra-ket overlaps
  allocate(smo(nmoB,nmoK))
  smo=smo1
  
  ! No. electrons
  call get_nel(n_intB,detB1(:,:,1),nelB,nel_alphaB,nel_betaB)
  call get_nel(n_intK,detK1(:,:,1),nelK,nel_alphaK,nel_betaK)

!----------------------------------------------------------------------
! Truncate the wave functions
!----------------------------------------------------------------------
  ! Bra
  call truncate_wave_functions(n_intB,ndetB1,nrootsB,detB1,vecB1,&
       normthrsh,ndetB,detB,vecB)
  
  ! Ket
  call truncate_wave_functions(n_intK,ndetK1,nrootsK,detK1,vecK1,&
       normthrsh,ndetK,detK,vecK)

  ! Ouput the numbers of determinants
  write(6,'(/,x,a,x,i0)') &
       'No. bra determinants before truncation:',ndetB1
  write(6,'(x,a,x,i0)') &
       'No. ket determinants before truncation:',ndetK1
  write(6,'(/,x,a,x,F10.7)') 'Truncation threshold:',normthrsh
  write(6,'(/,x,a,x,i0)') &
       'No. bra determinants after truncation:',ndetB
  write(6,'(x,a,x,i0)') &
       'No. ket determinants after truncation:',ndetK

!----------------------------------------------------------------------
! Symmetric orthogonalisation the truncated wave functions
!----------------------------------------------------------------------
  ! Bra
  call symm_ortho(n_intB,ndetB,nrootsB,vecB)

  ! Ket
  call symm_ortho(n_intK,ndetK,nrootsK,vecK)
  
!----------------------------------------------------------------------
! Sorting of the bra and ket determinants, as well as the
! determination of the unique alpha and beta strings
!----------------------------------------------------------------------
  ! Bra
  call det_sorting(n_intB,ndetB,nrootsB,detB,vecB,nalphaB,nbetaB,&
       alphaB,betaB,offsetB,det2betaB)

  ! Ket
  call det_sorting(n_intK,ndetK,nrootsK,detK,vecK,nalphaK,nbetaK,&
       alphaK,betaK,offsetK,det2betaK)

  ! Ouput the number of unique alpha and beta strings
  write(6,'(/,x,a,x,i0,a,i0)') 'No. bra alpha/beta strings:',&
       nalphaB,'/',nbetaB
  write(6,'(x,a,x,i0,a,i0)')   'No. ket alpha/beta strings:',&
       nalphaK,'/',nbetaK
    
!----------------------------------------------------------------------
! Pre-computation of the unique beta factors
!----------------------------------------------------------------------
  allocate(betafac(nbetaB,nbetaK))
  betafac=0.0d0

  call get_all_factors(nel_betaB,nbetaB,nbetaK,betaB,betaK,betafac)

!----------------------------------------------------------------------
! Calculate the wave function overlaps
!----------------------------------------------------------------------
  call get_overlaps(npairs,ipairs,Sij)

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'overlap')
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine overlap