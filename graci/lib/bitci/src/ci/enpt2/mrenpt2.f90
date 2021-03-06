!######################################################################
! Top-level routine to compute the MR-ENPT2 energy and wave function
! corrections
!######################################################################
#ifdef CBINDING
subroutine mrenpt2(irrep,nroots,nextra,shift,multistate,confscr,&
     vecscr,vec0scr,Qscr,dspscr) bind(c,name="mrenpt2")
#else
subroutine mrenpt2(irrep,nroots,nextra,shift,multistate,confscr,&
       vecscr,vec0scr,Qscr,dspscr)
#endif

  use constants
  use bitglobal
  use conftype
  use hii
  use epstein_nesbet
  use utils
  use iomod
  use timing

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots

  ! Number of extra roots to include in the ENPT2 calculation
  integer(is), intent(in)  :: nextra

  ! ISA shift
  real(dp), intent(in)     :: shift
  
  ! Multistate calculation?
  logical, intent(in)      :: multistate
  
  ! Array of MRCI configuration scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1)

  ! Array of reference space eigenvector scratch file numbers
  integer(is), intent(in)  :: vec0scr(0:nirrep-1)

  ! Eigenpair scratch file number
  integer(is), intent(out) :: vecscr

  ! Q-space info scratch file number
  integer(is), intent(out) :: Qscr

  ! Damped strong perturber scratch file number
  integer(is), intent(out) :: dspscr
  
  ! MRCI configuration derived type
  type(mrcfg)              :: cfg

  ! On-diagonal Hamiltonian matrix elements
  real(dp), allocatable    :: hdiag(:)

  ! On-diagonal Hamiltonian matrix elements
  ! averaged over spin couplings
  real(dp), allocatable    :: averageii(:)
  
  ! ENPT2 energy and wavefunction corrections
  real(dp), allocatable    :: Avec(:,:),Avec_ortho(:,:)
  real(dp), allocatable    :: E2(:)

  ! Zeroth-order eigenpairs
  integer(is)              :: refdim
  real(dp), allocatable    :: E0(:)
  real(dp), allocatable    :: vec0(:,:)

  ! 2nd-order corrected energies
  real(dp), allocatable    :: EPT2(:)

  ! QDPT2 energies and mixing coefficients
  real(dp), allocatable    :: EQD(:),mix(:,:),work(:,:)
  
  ! I/O variables
  integer(is)              :: iscratch
  character(len=60)        :: vecfile,Qfile
  character(len=2)         :: amult,airrep
  
  ! Everything else
  integer(is)              :: i,j
  integer(is)              :: nvec
  integer(is), allocatable :: indx(:)
  real(dp)                 :: norm
  real(dp), allocatable    :: Smat(:,:),Sinvsq(:,:)
  real(dp), allocatable    :: Elow(:)
  real(dp), allocatable    :: Qnorm(:),Qener(:)
  
  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  ! Section header
  write(6,'(/,52a)') ('-',i=1,52)
  if (multistate) then
     write(6,'(3(x,a))') 'MS-MR-ENPT2 corrections for the',&
          trim(irreplbl(irrep,ipg)),&
          'subspace'
  else
     write(6,'(3(x,a))') 'MR-ENPT2 corrections for the',&
          trim(irreplbl(irrep,ipg)),&
          'subspace'
  endif
  write(6,'(52a)') ('-',i=1,52)

!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

  write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! Total number of roots for which the ENPT2 corrections will be
  ! calculated
  nvec=nroots+nextra

  ! Number of reference space CSFs
  refdim=cfg%csfs0h(cfg%n0h+1)-1
  
  allocate(hdiag(cfg%csfdim))
  hdiag=0.0d0

  allocate(averageii(cfg%confdim))
  averageii=0.0d0

  allocate(Avec(cfg%csfdim,nvec))
  Avec=0.0d0

  allocate(Avec_ortho(cfg%csfdim,nroots))
  Avec_ortho=0.0d0
  
  allocate(E0(nvec))
  E0=0.0d0

  allocate(vec0(refdim,nvec))
  vec0=0.0d0
  
  allocate(E2(nvec))
  E2=0.0d0

  allocate(EPT2(nvec))
  EPT2=0.0d0

  allocate(Elow(nroots))
  Elow=0.0d0

  allocate(Qnorm(nroots))
  Qnorm=0.0d0

  allocate(Qener(nroots))
  Qener=0.0d0
  
  allocate(indx(nvec))
  indx=0

  allocate(Smat(nroots,nroots))
  Smat=0.0d0

  allocate(Sinvsq(nroots,nroots))
  Sinvsq=0.0d0

  if (multistate) then
     allocate(EQD(nvec))
     allocate(mix(nvec,nvec))
     allocate(work(cfg%csfdim,nvec))
  endif
  
!----------------------------------------------------------------------
! Read in the zeroth-order eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vec0scr(irrep),vec0,E0,refdim,nvec)

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call hmat_diagonal(hdiag,cfg%csfdim,averageii,cfg%confdim,cfg)

!----------------------------------------------------------------------
! Compute the ENPT2 energy and wave function corrections
!----------------------------------------------------------------------
  if (multistate) then
     call enpt2(irrep,cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
          vec0scr(irrep),Avec,E2,nvec,shift,dspscr,multistate,EQD,mix)
  else
     call enpt2(irrep,cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
          vec0scr(irrep),Avec,E2,nvec,shift,dspscr,multistate)
  endif

!----------------------------------------------------------------------
! Add in the zeroth-order wave functions
!----------------------------------------------------------------------
  do i=1,nvec
     Avec(1:refdim,i)=vec0(:,i)
  enddo
  
!----------------------------------------------------------------------
! If this a multi-state calculation, then mix the 1st-order corrected
! wave functions using the QDPT2 coefficients
!----------------------------------------------------------------------
  if (multistate) then
     work=Avec
     Avec=matmul(work,mix)
  endif
  
!----------------------------------------------------------------------
! Sort the 2nd-order corrected energies
!----------------------------------------------------------------------
  ! 2nd-order corrected energies
  if (multistate) then
     ! Take the multi-state energies
     EPT2=EQD
  else
     ! Take the single-state energies
     EPT2=E0+E2
  endif
  
  ! Sorting
  call dsortindxa1('A',nvec,EPT2,indx)

  ! Lowest-lying energies
  do i=1,nroots
     Elow(i)=EPT2(indx(i))
  enddo

!----------------------------------------------------------------------
! Write the Q-space information to disk
!----------------------------------------------------------------------
  ! Norms of the 1st-order wave corrections
  do i=1,nroots
     Qnorm(i)=sqrt(dot_product(Avec(refdim:cfg%csfdim,indx(i)),&
          Avec(refdim:cfg%csfdim,indx(i))))
  enddo

  ! Energy corrections
  do i=1,nroots
     Qener(i)=E2(indx(i))
  enddo
  
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('Qinfo'//'.mult'//trim(amult)//&
       '.sym'//trim(airrep),Qfile)
  call register_scratch_file(Qscr,Qfile)

  ! Open the scratch file
  iscratch=scrunit(Qscr)
  open(iscratch,file=scrname(Qscr),form='unformatted',&
       status='unknown')

  ! Number of roots
  write(iscratch) nroots
  
  ! Norms of the wave function corrections
  write(iscratch) Qnorm

  ! Energy corrections
  write(iscratch) Qener
  
  ! Close the scratch file
  close(iscratch)
  
!----------------------------------------------------------------------
! 1st-order corrected wave functions for the nroots lowest energy
! roots only
!----------------------------------------------------------------------
  do i=1,nroots
     Avec_ortho(:,i)=Avec(:,indx(i))
  enddo
  
!----------------------------------------------------------------------
! Normalisation
!----------------------------------------------------------------------
  do i=1,nroots
     norm=dot_product(Avec_ortho(:,i),Avec_ortho(:,i))
     norm=sqrt(norm)
     Avec_ortho(:,i)=Avec_ortho(:,i)/norm
  enddo
  
!----------------------------------------------------------------------
! Orthonogonalise using Lowdin's symmetric orthogonalisation scheme
!----------------------------------------------------------------------
! Note that this ensures that the resulting vectors are the closest
! possible in the least squares sense to the original ones
!----------------------------------------------------------------------
  ! Overlap matrix: this should be replaced with a BLAS level operation
  ! in the future
  do i=1,nroots
     Smat(i,i)=1.0d0
  enddo
  do i=1,nroots-1
     do j=i+1,nroots
        Smat(i,j)=dot_product(Avec_ortho(:,i),Avec_ortho(:,j))
        Smat(j,i)=Smat(i,j)
     enddo
  enddo

  ! Inverse square root of the overlap matrix
  call invsqrt_matrix(Smat,Sinvsq,nroots)

  ! Orthogonalise
  Avec(:,1:nroots)=Avec_ortho
  Avec_ortho=matmul(Avec(:,1:nroots),Sinvsq)

!----------------------------------------------------------------------
! Write the 2nd-order corrected energies and 1st-order corrected
! wave functions to disk
!----------------------------------------------------------------------
  ! Register the scratch file
   write(amult,'(i0)') imult
   write(airrep,'(i0)') irrep
   call scratch_name('mrenpt2vec'//'.mult'//trim(amult)//&
        '.sym'//trim(airrep),vecfile)
   call register_scratch_file(vecscr,vecfile)

   ! Open the scratch file
   iscratch=scrunit(vecscr)
   open(iscratch,file=scrname(vecscr),form='unformatted',&
        status='unknown')
   
   ! No. CSFs
   write(iscratch) cfg%csfdim
    
   ! No. roots
   write(iscratch) nroots
    
   ! 2nd-order corrected energies
   write(iscratch) Elow
    
   ! 1st-order corrected wave functions
   do i=1,nroots
      write(iscratch) Avec_ortho(:,i)
   enddo

   ! Close the scratch file
   close(iscratch)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(hdiag)
  deallocate(averageii)
  deallocate(Avec)
  deallocate(Avec_ortho)
  deallocate(E0)
  deallocate(vec0)
  deallocate(E2)
  deallocate(EPT2)
  deallocate(Elow)
  deallocate(Qnorm)
  deallocate(Qener)
  deallocate(indx)
  deallocate(Smat)
  deallocate(Sinvsq)
  if (allocated(EQD)) deallocate(EQD)
  if (allocated(mix)) deallocate(mix)
  if (allocated(work)) deallocate(work)
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'mrenpt2')
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
    
  return
  
end subroutine mrenpt2

