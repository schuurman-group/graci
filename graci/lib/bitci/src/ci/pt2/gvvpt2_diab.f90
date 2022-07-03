!######################################################################
! Top-level routine to compute quasi-diabatic states using 2nd-order
! generalised van Vleck perturbation theory
!######################################################################
#ifdef CBINDING
subroutine gvvpt2_diab(irrep,nroots,nextra,shift,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smoR0,adtR0,ncore,icore,lfrzcore,&
     confscr,vecscr,vec0scr,Qscr,dspscr) bind(c,name="gvvpt2_diab")
#else
subroutine gvvpt2_diab(irrep,nroots,nextra,shift,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smoR0,adtR0,ncore,icore,lfrzcore,&
     confscr,vecscr,vec0scr,Qscr,dspscr)
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

  ! ADT matrix at the previous geometry
  real(dp), intent(in)     :: adtR0(nrootsR0,nrootsR0)
  
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
  type(mrcfg)              :: cfg,cfg_ref

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
  real(dp), allocatable    :: E_ref(:)
  real(dp), allocatable    :: vec_ref(:,:)
  integer(ib), allocatable :: det_ref(:,:,:)
  real(dp), allocatable    :: vec_ref_det(:,:)
  
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

  ! Representation of the prototype diabatic states in the basis
  ! of the ref space eigenfunctions
  real(dp), allocatable    :: Cmat(:,:)

  ! Representation of the prototype diabatic states in the basis
  ! of the ref space eigenfunctions
  real(dp), allocatable    :: vec_proto(:,:)
  
  ! Wave function selection
  integer(is)              :: imin(1),imax(1)
  integer(is), allocatable :: isel(:)
  real(dp), allocatable    :: Sij1(:,:),maxSij(:)
  real(dp), allocatable    :: Avec_sel(:,:)
  real(dp), allocatable    :: EQD_sel(:)
  
  ! I/O variables
  integer(is)              :: iscratch
  character(len=60)        :: vecfile,Qfile
  character(len=2)         :: amult,airrep
  
  ! Everything else
  integer(is)              :: i,j,n
  integer(is)              :: n_int_I,nconf_ref,ndet_ref
  integer(is)              :: nvec
  integer(is), allocatable :: indx(:)
  real(dp), allocatable    :: Qnorm(:),Qener(:),Qnorm1(:),Qener1(:)
  logical                  :: lprint
  
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
  if (verbose) then
     write(6,'(/,52a)') ('-',i=1,52)
     write(6,'(x,a,/,x,a)') &
          'GVVPT2 diabatisation calculation for the ' &
          //trim(irreplbl(irrep,ipg)),'subspace'
     write(6,'(52a)') ('-',i=1,52)
  endif

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

  if (verbose) &
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

  allocate(E_ref(nvec))
  E_ref=0.0d0

  allocate(vec_ref(refdim,nvec))
  vec_ref=0.0d0
  
  allocate(E2(nvec))
  E2=0.0d0

  allocate(Qnorm(nvec))
  Qnorm=0.0d0

  allocate(Qener(nvec))
  Qener=0.0d0

  allocate(Qnorm1(nvec))
  Qnorm=0.0d0

  allocate(Qener1(nvec))
  Qener=0.0d0
  
  allocate(EQD(nvec))
  EQD=0.0d0
  
  allocate(mix(nvec,nvec))
  mix=0.0d0
  
  allocate(work(cfg%csfdim,nvec))
  work=0.0d0

  allocate(indx(nroots))
  indx=0.0d0

!----------------------------------------------------------------------
! Read in the zeroth-order eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vec0scr(irrep),vec_ref,E_ref,refdim,nvec)

!----------------------------------------------------------------------
! Get the determinant representation of the zeroth-order eigenstates
!----------------------------------------------------------------------
  ! Set up a dummy MRCI configuration derived data type to hold
  ! the ref conf info
  n_int_I=cfg%n_int_I
  nconf_ref=cfg%n0h
  allocate(cfg_ref%conf0h(n_int_I,2,nconf_ref))
  allocate(cfg_ref%sop0h(n_int_I,2,nconf_ref))
  cfg_ref%n0h=nconf_ref
  cfg_ref%n_int_I=n_int_I
  cfg_ref%n1I=0
  cfg_ref%n2I=0
  cfg_ref%n1E=0
  cfg_ref%n2E=0
  cfg_ref%n1I1E=0  
  cfg_ref%conf0h=cfg%conf0h
  cfg_ref%sop0h=cfg%sop0h

  ! Determine the dimension of the determinant basis
  call get_detdim(cfg_ref,ndet_ref)

  ! Allocate arrays
  allocate(det_ref(n_int,2,ndet_ref))
  allocate(vec_ref_det(ndet_ref,nvec))
  det_ref=0_ib
  vec_ref_det=0.0d0

  ! Compute the determinant representation of the wave functions
  call det_trans(cfg_ref,cfg%m2c,nvec,refdim,ndet_ref,vec_ref,&
       vec_ref_det,det_ref)
  
!----------------------------------------------------------------------
! Compute the overlaps of the quasi-diabatic states from the
! previous geometry with the zeroth-order eigenstates at the current
! geometry
!----------------------------------------------------------------------
  ! Allocate arrays
  npairs=nrootsR0*nvec
  allocate(Sij(npairs))
  allocate(ipairs(npairs,2))
  allocate(Cmat(nvec,nrootsR0))
  Sij=0.0d0
  ipairs=0
  Cmat=0.0d0

  ! Truncation threshold
  normthrsh=0.999d0

  ! Fill in the array of required bra-ket overlaps
  n=0
  do i=1,nrootsR0
     do j=1,nvec
        n=n+1
        ipairs(n,1)=i
        ipairs(n,2)=j
     enddo
  enddo

  ! Compute the overlaps
  lprint=.false.
  call overlap(nmoR0,nmo,n_intR0,n_int,ndetR0,ndet_ref,nrootsR0,nvec,&
       detR0,det_ref,vecR0,vec_ref_det,smoR0,normthrsh,ncore,icore,&
       lfrzcore,npairs,Sij,ipairs,lprint)

  ! Put the overlaps into a more computationally useful form
  n=0
  do i=1,nrootsR0
     do j=1,nvec
        n=n+1
        Cmat(j,i)=Sij(n)
     enddo
  enddo

!----------------------------------------------------------------------
! Construct the prototype diabatic states
!----------------------------------------------------------------------
  ! Transform using the previous geometry ADT matrix to yield
  ! the prototype diabatic states in the basis of the current geometry
  ! ref space eigenstates
  Cmat=matmul(Cmat,adtR0)

  ! Lowdin's symmetric orthonormalisation of the prototype diabatic
  ! states
  call symm_ortho(nvec,nrootsR0,Cmat)

  ! Prototype diabatic states in the basis of the ref space CSFs
  allocate(vec_proto(refdim,nrootsR0))
  vec_proto=matmul(vec_ref,Cmat)

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call hmat_diagonal(hdiag,cfg%csfdim,averageii,cfg%confdim,cfg)
  
!----------------------------------------------------------------------
! Construct the GVVPT2 effective Hamiltonian and determine its
! eigenpairs
! Also returns the 1st-order perturbed model functions in the Avec
! array
!----------------------------------------------------------------------

  !
  ! Currently, gvvpt2_heff assumes that the model states are identical
  ! to the ref space eigenstates
  !
  ! We need to generalise this by passing a ref-state-to-model-state
  ! transformation matrix
  !
  ! For a normal GVVPT2 calculation, this will simply be the identity
  ! matrix
  !
  ! For a diabatisation calculation, this will be the Cmat array,
  ! i.e., the representation of the prototype diabatic states in the
  ! basis of the ref states
  !
  ! We are also going to have to pass H0, i.e., the model state
  ! representation of the zeroth-order Hamiltonian
  !
  ! In a normal run, this will just be the diagonal matrix of
  ! ref state energies
  !
  ! In a diabatisation run, this will be the transformation of this
  ! matrix using Cmat
  !
  ! We might as well also modify gvvpt2 to take as input the
  ! model space vectors directly rather than the ref space eigenvector
  ! scratch file number
  !
  
  
  !call gvvpt2_heff(irrep,cfg,hdiag,averageii,cfg%csfdim,cfg%confdim,&
  !     vec0scr(irrep),Avec,E2,nvec,shift,dspscr,EQD,mix)

  STOP
  
!!----------------------------------------------------------------------
!! Add in the zeroth-order wave functions
!!----------------------------------------------------------------------
!  do i=1,nvec
!     Avec(1:refdim,i)=vec0(:,i)
!  enddo
!
!!----------------------------------------------------------------------
!! Compute the 1st-order wave functions
!!----------------------------------------------------------------------
!  work=Avec
!  Avec=matmul(work,mix)
!
!!----------------------------------------------------------------------
!! Q-space information for ***all*** states
!!----------------------------------------------------------------------
!  ! Norms of the 1st-order wave functions projected onto the Q-space
!  do i=1,nvec
!     Qnorm(i)=sqrt(dot_product(Avec(refdim+1:cfg%csfdim,i),&
!          Avec(refdim+1:cfg%csfdim,i)))
!  enddo
!
!  ! ENPT2 energy corrections
!  do i=1,nvec
!     Qener(i)=E2(i)
!  enddo
!  
!!----------------------------------------------------------------------
!! Lowdin's symmetric orthonormalisation of the 1st-order corrected
!! wave functions
!!----------------------------------------------------------------------
!  call symm_ortho(cfg%csfdim,nvec,Avec)
!
!!----------------------------------------------------------------------
!! Get the determinant representation of the 1st-order corrected
!! wave functions
!!----------------------------------------------------------------------
!  ! Determine the dimension of the determinant basis
!  call get_detdim(cfg,ndet)
!
!  ! Allocate arrays
!  allocate(det(n_int,2,ndet))
!  allocate(Avec_det(ndet,nvec))
!  det=0_ib
!  Avec_det=0.0d0
!
!  ! Compute the determinant representation of the wave functions
!  call det_trans(cfg,cfg%m2c,nvec,cfg%csfdim,ndet,Avec,Avec_det,det)
!  
!!----------------------------------------------------------------------
!! Compute the overlaps with the input/target wave functions
!!----------------------------------------------------------------------
!  ! Allocate arrays
!  npairs=nrootsR0*nvec
!  allocate(Sij(npairs))
!  allocate(ipairs(npairs,2))
!  Sij=0.0d0
!  ipairs=0
!
!  ! Truncation threshold
!  normthrsh=0.99d0
!
!  ! Fill in the array of bra-ket overlaps required
!  n=0
!  do i=1,nrootsR0
!     do j=1,nvec
!        n=n+1
!        ipairs(n,1)=i
!        ipairs(n,2)=j
!     enddo
!  enddo
!
!  ! Compute the overlaps
!  lprint=.false.
!  call overlap(nmoR0,nmo,n_intR0,n_int,ndetR0,ndet,nrootsR0,nvec,&
!       detR0,det,vecR0,Avec_det,smoR0,normthrsh,ncore,icore,lfrzcore,&
!       npairs,Sij,ipairs,lprint)
!
!!----------------------------------------------------------------------
!! Deallocate the Avec_det array now that it is no longer needed
!!----------------------------------------------------------------------
!  deallocate(Avec_det)
!  
!!----------------------------------------------------------------------
!! Determine the 1st-order corrected wave functions of interest
!!----------------------------------------------------------------------
!  ! Allocate arrays
!  allocate(isel(nroots))
!  allocate(maxSij(nroots))
!  allocate(Sij1(nrootsR0,nvec))
!  isel=0
!  maxSij=0.0d0
!  Sij1=0.0d0
!
!  ! Reshape the wave function overlap matrix into a more useful
!  ! form
!  do n=1,npairs
!     i=ipairs(n,1)
!     j=ipairs(n,2)
!     Sij1(i,j)=Sij(n)
!  enddo
!
!  ! Select the reference space eigenfunctions with the greatest
!  ! overlaps with the input wave functions
!  do i=1,nrootsR0
!     imax=maxloc(abs(Sij1(i,:)))
!     isel(i)=imax(1)
!     maxSij(i)=Sij1(i,imax(1))
!     Sij1(:,imax(1))=0.0d0
!  enddo
!
!!----------------------------------------------------------------------
!! Fill in the selected states
!! We will also fix the phases of the states in this step
!!----------------------------------------------------------------------
!  ! Allocate arrays
!  allocate(Avec_sel(cfg%csfdim,nroots), EQD_sel(nroots))
!  Avec_sel=0.0d0
!  EQD_sel=0.0d0
!
!  ! Fill in the selected states, re-phasing as we go
!  do i=1,nroots
!
!     ! First-order corrected wave function
!     Avec_sel(:,i)=Avec(:,isel(i))
!
!     ! Fix the phase
!     if (maxSij(i) < 0.0d0) Avec_sel(:,i)=-Avec_sel(:,i)
!
!     ! Energy
!     EQD_sel(i)=EQD(isel(i))
!     
!  enddo
!
!  ! Rearrange the Q-space info arrays
!  Qnorm1=Qnorm
!  Qener1=Qener
!  Qnorm=0.0d0
!  Qener=0.0d0
!  do i=1,nroots
!     Qnorm(i)=Qnorm1(isel(i))
!     Qener(i)=Qener1(isel(i))
!  enddo
!  
!!----------------------------------------------------------------------
!! Sort the selected states by energy
!!----------------------------------------------------------------------
!  ! Sort by energy
!  call dsortindxa1('A',nroots,EQD_sel,indx)
!
!  ! Re-order the states, energies and Q-space information arrays
!  ! We will use Avec and EQD arrays as work arrays here
!  do i=1,nroots
!     Avec(:,i)=Avec_sel(:,indx(i))
!     EQD(i)=EQD_sel(indx(i))
!     Qnorm1(i)=Qnorm(indx(i))
!     Qener1(i)=Qener(indx(i))
!  enddo
!
!  Avec_sel(:,1:nroots)=Avec(:,1:nroots)
!  EQD_sel(1:nroots)=EQD(1:nroots)
!  Qnorm(1:nroots)=Qnorm1(1:nroots)
!  Qener(1:nroots)=Qener1(1:nroots)
!  
!  ! Deallocate the Avec and EQD now that they have been 'violated',
!  ! just to be on the safe side
!  deallocate(Avec)
!  deallocate(EQD)
!
!!----------------------------------------------------------------------
!! Write the Q-space information to disk
!!----------------------------------------------------------------------
!  ! Register the scratch file
!  write(amult,'(i0)') imult
!  write(airrep,'(i0)') irrep
!  call scratch_name('Qinfo'//'.mult'//trim(amult)//&
!       '.sym'//trim(airrep),Qfile)
!  call register_scratch_file(Qscr,Qfile)
!
!  ! Open the scratch file
!  iscratch=scrunit(Qscr)
!  open(iscratch,file=scrname(Qscr),form='unformatted',&
!       status='unknown')
!
!  ! Number of roots
!  write(iscratch) nroots
!  
!  ! Norms of the wave function corrections
!  write(iscratch) Qnorm(1:nroots)
!
!  ! Energy corrections
!  write(iscratch) Qener(1:nroots)
!
!  ! Close the scratch file
!  close(iscratch)
!  
!!----------------------------------------------------------------------
!! Write the 2nd-order energies and 1st-order wave functions to disk
!!----------------------------------------------------------------------
!  ! Register the scratch file
!  write(amult,'(i0)') imult
!  write(airrep,'(i0)') irrep
!  call scratch_name('mrenpt2vec'//'.mult'//trim(amult)//&
!       '.sym'//trim(airrep),vecfile)
!  call register_scratch_file(vecscr,vecfile)
!
!  ! Open the scratch file
!  iscratch=scrunit(vecscr)
!  open(iscratch,file=scrname(vecscr),form='unformatted',&
!       status='unknown')
!   
!  ! No. CSFs
!  write(iscratch) cfg%csfdim
!    
!  ! No. roots
!  write(iscratch) nroots
!    
!  ! 2nd-order energies
!  write(iscratch) EQD_sel(1:nroots)
!    
!  ! 1st-order wave functions
!  do i=1,nroots
!     write(iscratch) Avec_sel(:,i)
!  enddo
!
!  ! Close the scratch file
!  close(iscratch)
!
!!----------------------------------------------------------------------
!! Deallocate arrays
!!----------------------------------------------------------------------
!  deallocate(hdiag)
!  deallocate(averageii)
!  deallocate(Avec_sel)
!  deallocate(E0)
!  deallocate(vec0)
!  deallocate(E2)
!  deallocate(Qnorm)
!  deallocate(Qener)
!  deallocate(Qnorm1)
!  deallocate(Qener1)
!  deallocate(EQD_sel)
!  deallocate(mix)
!  deallocate(work)
!  deallocate(det)
!  deallocate(Sij)
!  deallocate(ipairs)
!  deallocate(isel)
!  deallocate(maxSij)
!  deallocate(Sij1)
!  deallocate(indx)
!  
!!----------------------------------------------------------------------
!! Stop timing and print report
!!----------------------------------------------------------------------
!  call get_times(twall_end,tcpu_end)
!  if (verbose) &
!       call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
!       'gvvpt2_diab')
!
!!----------------------------------------------------------------------
!! Flush stdout
!!----------------------------------------------------------------------
!  flush(6)
  
  return
  
end subroutine gvvpt2_diab
  
