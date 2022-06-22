!######################################################################
! Top-level routine to compute the 2nd-order generalised van Vleck
! perturbation theory energies and wave functions using the Epstein-
! Nesbet Hamiltonian partitioning ***with root following***
!######################################################################
#ifdef CBINDING
subroutine gvvpt2_follow(irrep,nroots,nextra,shift,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,confscr,&
     vecscr,vec0scr,Qscr,dspscr) bind(c,name="gvvpt2_follow")
#else
subroutine gvvpt2_follow(irrep,nroots,nextra,shift,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,confscr,&
     vecscr,vec0scr,Qscr,dspscr)
#endif

  use constants
  use bitglobal
  use conftype
  use hii
  use heff
  use csf2det
  use mrciutils
  use utils
  use iomod
  use timing

  implicit none
  
  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots

  ! Number of extra roots to include in the GVVPT2 calculation
  integer(is), intent(in)  :: nextra

  ! ISA shift
  real(dp), intent(in)     :: shift

  ! Eigenvectors to follow, expressed in a Slater determinant basis
  integer(is), intent(in)  :: n_intR0,ndetR0,nrootsR0
  integer(ib), intent(in)  :: detR0(n_intR0,2,ndetR0)
  real(dp), intent(in)     :: vecR0(ndetR0,nrootsR0)

   ! MO overlaps
  integer(is), intent(in)  :: nmoR0
  real(dp), intent(in)     :: smoR0(nmoR0,nmo)

  ! Frozen core orbital info
  integer(is), intent(in)  :: ncore
  integer(is), intent(in)  :: icore(ncore)
  logical(is), intent(in)  :: lfrzcore
  
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
  real(dp), allocatable    :: Avec(:,:)
  real(dp), allocatable    :: E2(:)

  ! Zeroth-order eigenpairs
  integer(is)              :: refdim
  real(dp), allocatable    :: E0(:)
  real(dp), allocatable    :: vec0(:,:)

  ! QDPT2 energies and mixing coefficients
  real(dp), allocatable    :: EQD(:),mix(:,:),work(:,:)

  ! Determinant representation of the 1st-order corrected
  ! wave functions
  integer(is)              :: ndet
  integer(ib), allocatable :: det(:,:,:)
  real(dp), allocatable    :: Avec_det(:,:)

  ! Wave function overlaps
  integer(is)              :: npairs
  integer(is), allocatable :: ipairs(:,:)
  real(dp), allocatable    :: Sij(:)
  real(dp)                 :: normthrsh

  ! Wave function selection
  integer(is)              :: imin(1),imax(1)
  integer(is), allocatable :: isel(:)
  real(dp), allocatable    :: Sij1(:,:),maxSij(:)
  
  ! I/O variables
  integer(is)              :: iscratch
  character(len=60)        :: vecfile,Qfile
  character(len=2)         :: amult,airrep
  
  ! Everything else
  integer(is)              :: i,j,n
  integer(is)              :: nvec
  real(dp), allocatable    :: selarr(:)
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
  write(6,'(x,a,/,x,a)') 'Root-following GVVPT2 calculation for the ' &
       //trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)

!----------------------------------------------------------------------
! For now, we will only support nroots = nrootsR0
! i.e., same numbers of roots and roots to follow
!----------------------------------------------------------------------
  if (nroots /= nrootsR0) then
     errmsg='Error in gvvpt2_follow: nroots != nrootsR0'
     call error_control
  endif
  
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

  allocate(E0(nvec))
  E0=0.0d0

  allocate(vec0(refdim,nvec))
  vec0=0.0d0
  
  allocate(E2(nvec))
  E2=0.0d0

  allocate(Qnorm(nvec))
  Qnorm=0.0d0

  allocate(Qener(nvec))
  Qener=0.0d0
  
  allocate(EQD(nvec))
  EQD=0.0d0
  
  allocate(mix(nvec,nvec))
  mix=0.0d0
  
  allocate(work(cfg%csfdim,nvec))
  work=0.0d0

  allocate(selarr(nroots))
  selarr=0.0d0
  
!----------------------------------------------------------------------
! Read in the zeroth-order eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vec0scr(irrep),vec0,E0,refdim,nvec)

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call hmat_diagonal(hdiag,cfg%csfdim,averageii,cfg%confdim,cfg)

!----------------------------------------------------------------------
! Compute the eigenpairs of the GVVPT2 effective Hamiltonian
! Also returns the 1st-order perturbed model functions in the Avec
! array
!----------------------------------------------------------------------
  call gvvpt2_heff(irrep,cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
       vec0scr(irrep),Avec,E2,nvec,shift,dspscr,EQD,mix)

!----------------------------------------------------------------------
! Add in the zeroth-order wave functions
!----------------------------------------------------------------------
  do i=1,nvec
     Avec(1:refdim,i)=vec0(:,i)
  enddo

!----------------------------------------------------------------------
! Compute the 1st-order wave functions
!----------------------------------------------------------------------
  work=Avec
  Avec=matmul(work,mix)

!----------------------------------------------------------------------
! Q-space information for ***all*** states
!----------------------------------------------------------------------
  ! Norms of the 1st-order wave functions projected onto the Q-space
  do i=1,nvec
     Qnorm(i)=sqrt(dot_product(Avec(refdim:cfg%csfdim,i),&
          Avec(refdim:cfg%csfdim,i)))
  enddo

  ! ENPT2 energy corrections
  do i=1,nvec
     Qener(i)=E2(i)
  enddo
  
!----------------------------------------------------------------------
! Lowdin's symmetric orthonormalisation of the 1st-order corrected
! wave functions
!----------------------------------------------------------------------
  call symm_ortho(cfg%csfdim,nvec,Avec)

!----------------------------------------------------------------------
! Get the determinant representation of the 1st-order corrected
! wave functions
!----------------------------------------------------------------------
  ! Determine the dimension of the determinant basis
  call get_detdim(cfg,ndet)

  ! Allocate arrays
  allocate(det(n_int,2,ndet))
  allocate(Avec_det(ndet,nvec))
  det=0_ib
  Avec_det=0.0d0

  ! Compute the determinant representation of the wave functions
  call det_trans(cfg,cfg%m2c,nvec,cfg%csfdim,ndet,Avec,Avec_det,det)
  
!----------------------------------------------------------------------
! Compute the overlaps with the input/target wave functions
!----------------------------------------------------------------------
  ! Allocate arrays
  npairs=nrootsR0*nvec
  allocate(Sij(npairs))
  allocate(ipairs(npairs,2))
  Sij=0.0d0
  ipairs=0

  ! Truncation threshold
  normthrsh=0.99d0

  ! Fill in the array of bra-ket overlaps required
  n=0
  do i=1,nrootsR0
     do j=1,nvec
        n=n+1
        ipairs(n,1)=i
        ipairs(n,2)=j
     enddo
  enddo

  ! Compute the overlaps
  call overlap(nmoR0,nmo,n_intR0,n_int,ndetR0,ndet,nrootsR0,nvec,&
       detR0,det,vecR0,Avec_det,smoR0,normthrsh,ncore,icore,lfrzcore,&
       npairs,Sij,ipairs)

!----------------------------------------------------------------------
! Determine the 1st-order corrected wave functions of interest
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(isel(nroots))
  allocate(maxSij(nroots))
  allocate(Sij1(nrootsR0,nvec))
  isel=0
  maxSij=0.0d0
  Sij1=0.0d0

  ! Reshape the wave function overlap matrix into a more useful
  ! form
  do n=1,npairs
     i=ipairs(n,1)
     j=ipairs(n,2)
     Sij1(i,j)=Sij(n)
  enddo

  ! Select the reference space eigenfunctions with the greatest
  ! overlaps with the input wave functions
  do i=1,nrootsR0
     imax=maxloc(abs(Sij1(i,:)))
     isel(i)=imax(1)
     maxSij(i)=Sij1(i,imax(1))
     Sij1(:,imax(1))=0.0d0
  enddo

!----------------------------------------------------------------------
! Re-phasing of the states of interest
!----------------------------------------------------------------------
  do i=1,nroots
     if (maxSij(i) < 0.0d0) then
        Avec(:,isel(i))=-Avec(:,isel(i))
     endif
  enddo

!----------------------------------------------------------------------
! Write the Q-space information to disk
!----------------------------------------------------------------------
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
  selarr=0.0d0
  do i=1,nroots
     selarr(i)=Qnorm(isel(i))
  enddo
  write(iscratch) selarr

  ! Energy corrections
  selarr=0.0d0
  do i=1,nroots
     selarr(i)=Qener(isel(i))
  enddo
  write(iscratch) selarr

  ! Close the scratch file
  close(iscratch)
  
!----------------------------------------------------------------------
! Write the 2nd-order energies and 1st-order wave functions to disk
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
    
  ! 2nd-order energies
  do i=1,nroots
     selarr(i)=EQD(isel(i))
  enddo
  write(iscratch) selarr
    
  ! 1st-order wave functions
  do i=1,nroots
     write(iscratch) Avec(:,isel(i))
  enddo

  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(hdiag)
  deallocate(averageii)
  deallocate(Avec)
  deallocate(E0)
  deallocate(vec0)
  deallocate(E2)
  deallocate(Qnorm)
  deallocate(Qener)
  deallocate(EQD)
  deallocate(mix)
  deallocate(work)
  deallocate(selarr)
  deallocate(det)
  deallocate(Avec_det)
  deallocate(Sij)
  deallocate(ipairs)
  deallocate(isel)
  deallocate(maxSij)
  deallocate(Sij1)
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'gvvpt2')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
    
  return
  
end subroutine gvvpt2_follow

