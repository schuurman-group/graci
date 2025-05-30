!######################################################################
! Top-level routine to compute GVVPT2 QDPT diabatic states and energies
! using pre-computed A-vectors
!######################################################################
#ifdef CBINDING
subroutine gvvpt2_diab(irrep,nroots,nextra,ireg,regfac,n_intR0,&
     ndetR0,nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,&
     normthrsh,detthrsh,confscr,vec0scr,Ascr,diabpot,diab_vecscr,&
     diab_confscr,diab_nconf,diab_aviiscr) bind(c,name="gvvpt2_diab")
#else
subroutine gvvpt2_diab(irrep,nroots,nextra,ireg,regfac,n_intR0,&
     ndetR0,nrootsR0,detR0,vecR0,nmoR0,smoR0,ncore,icore,lfrzcore,&
     normthrsh,detthrsh,confscr,vec0scr,Ascr,diabpot,diab_vecscr,&
     diab_confscr,diab_nconf,diab_aviiscr)
#endif

  use constants
  use bitglobal
  use conftype
  use hii
  use model_states
  use heff_diabatic
  use utils
  use iomod
  use timing
  
  implicit none

  ! Irrep number
  integer(is), intent(in)  :: irrep

  ! Number of roots requested
  integer(is), intent(in)  :: nroots

  ! Number of buffer states
  integer(is), intent(in)  :: nextra

  ! Regularizer index: 1 <-> ISA
  !                    2 <-> sigma^p
  integer(is), intent(in)  :: ireg
  
  ! Regularization factor
  real(dp), intent(in)     :: regfac

  ! Previous geometry diabatic states expressed in a Slater
  ! determinant basis
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

  ! Norm-based wave function truncation threshold
  real(dp), intent(in)     :: normthrsh

  ! Determinant screening threshold
  real(dp), intent(in)     :: detthrsh
  
  ! Array of MRCI configuration scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1)

  ! Array of reference space eigenvector scratch file numbers
  integer(is), intent(in)  :: vec0scr(0:nirrep-1)

  ! A-vector scratch file number
  integer(is), intent(in)  :: Ascr

  ! Diabatic potential matrix
  real(dp), intent(out)    :: diabpot(nrootsR0,nrootsR0)

  ! Diabatic state vector scratch file number
  integer(is), intent(out) :: diab_vecscr

  ! Diabatic configuration scratch file number
  integer(is), intent(out) :: diab_confscr

  ! Number of diabatic configurations
  integer(is), intent(out) :: diab_nconf

  ! Spin-coupling-averaged on-diagonal Hamiltonian matrix elements
  ! scratch file number
  integer(is), intent(out) :: diab_aviiscr
  
  ! MRCI configuration derived types
  type(mrcfg)              :: cfg
  
  ! Zeroth-order eigenpairs
  integer(is)              :: refdim
  real(dp), allocatable    :: E0(:)
  real(dp), allocatable    :: vec0(:,:)

  ! On-diagonal Hamiltonian matrix elements
  real(dp), allocatable    :: hdiag(:)

  ! On-diagonal Hamiltonian matrix elements
  ! averaged over spin couplings
  real(dp), allocatable    :: averageii(:)
  
  ! A-vectors
  real(dp), allocatable    :: Avec(:,:),Avec_mod(:,:)

  ! Prototype diabatic states
  real(dp), allocatable    :: vec_pds(:,:)

  ! Complement states
  integer(is)              :: ncomp
  real(dp), allocatable    :: vec_comp(:,:)

  ! Model space states
  real(dp), allocatable    :: vec_mod(:,:)
  real(dp), allocatable    :: ref2mod(:,:)

  ! Zeroth- plus first-order effective Hamiltonian
  real(dp), allocatable    :: H01(:,:),H01_mod(:,:)

  ! I/O variables
  integer(is)              :: iscratch
  character(len=250)       :: diab_vecfile,diab_conffile,&
                              diab_aviifile
  character(len=2)         :: amult,airrep
  
  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end
  
  ! Everything else
  integer(is)              :: i,j,n
  integer(is)              :: nvec
  integer(is)              :: dim1,dim2
  real(dp), allocatable    :: diabii(:)
  
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
! i.e., same numbers of reference and displaced geometry states
!----------------------------------------------------------------------
  if (nroots /= nrootsR0) then
     errmsg='Error in gvvpt2_diab: nroots != nrootsR0'
     call error_control
  endif

!----------------------------------------------------------------------
! Output the regularizer being used
!----------------------------------------------------------------------
  select case(ireg)
  case(1)
     if (verbose) write(6,'(/,x,a)') 'Regularizer: ISA'
  case(2)
     if (verbose) write(6,'(/,x,a)') 'Regularizer: sigma^p'
  case default
     errmsg='Error in gvvpt2_diab: unrecognized regularizer index'
     call error_control
  end select
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

  if (verbose) &
       write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! No. reference space states used to construct the prototype
  ! diabatic states (PSDs)
  ! nvec = no. roots + no. buffer states
  nvec=nroots+nextra

  ! Number of reference space CSFs
  refdim=cfg%csfs0h(cfg%n0h+1)-1

  allocate(Avec(cfg%csfdim,nvec))
  Avec=0.0d0

  allocate(E0(nvec))
  E0=0.0d0

  allocate(vec0(refdim,nvec))
  vec0=0.0d0

  allocate(hdiag(cfg%csfdim))
  hdiag=0.0d0

  allocate(averageii(cfg%confdim))
  averageii=0.0d0
  
!----------------------------------------------------------------------
! Read in the zeroth-order eigenpairs
!----------------------------------------------------------------------
  call read_all_eigenpairs(vec0scr(irrep),vec0,E0,refdim,nvec)

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call hmat_diagonal(hdiag,cfg%csfdim,averageii,cfg%confdim,cfg)

!----------------------------------------------------------------------
! Save the spin-coupling-averaged on-diagonal Hamiltonian matrix
! elements to disk
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('diabavhii'//'.mult'//trim(amult)//&
       '.sym'//trim(airrep),diab_aviifile)
  call register_scratch_file(diab_aviiscr,diab_aviifile)

  ! Open the scratch file
  iscratch=scrunit(diab_aviiscr)
  open(iscratch,file=scrname(diab_aviiscr),form='unformatted', &
       status='unknown')

  ! Number of configurations
  write(iscratch) cfg%confdim

  ! Averaged matrix elements
  write(iscratch) averageii
  
  ! Close the scratch file
  close(iscratch)
  
!----------------------------------------------------------------------
! Compute the prototype diabatic states
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vec_pds(refdim,nrootsR0))
  vec_pds=0.0d0

  ! Compute the prototype diabatic states in the ref CSF basis
  call get_pds_basis(cfg,refdim,nvec,vec0,nmoR0,n_intR0,ndetR0,&
       nrootsR0,detR0,vecR0,smoR0,ncore,icore,lfrzcore,vec_pds,&
       normthrsh,detthrsh)

!----------------------------------------------------------------------
! Calculation of the complement space states
!----------------------------------------------------------------------
  ! Allocate arrays
  ncomp=nvec-nrootsR0
  allocate(vec_comp(refdim,ncomp))
  vec_comp=0.0d0

  ! Compute the complement states in the ref CSF basis
  call get_complement_basis(refdim,nrootsR0,nvec,ncomp,vec0,vec_pds,&
       vec_comp)

!----------------------------------------------------------------------
! Calculation of the model space states
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(vec_mod(refdim,nrootsR0))
  allocate(ref2mod(nvec,nrootsR0))
  vec_mod=0.0d0
  ref2mod=0.0d0

  ! Compute the model space states in the ref CSF basis
  call get_model_basis(nvec,nrootsR0,ncomp,refdim,E0,vec0,vec_pds,&
       vec_comp,vec_mod,ref2mod)

!----------------------------------------------------------------------
! Construct the zeroth-order effective Hamiltonian
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(H01(nvec,nvec))
  allocate(H01_mod(nrootsR0,nrootsR0))
  
  ! Hamiltonian in the ref state basis
  ! (subtraction of E_SCF gives us the true eigenvalues)
  H01=0.0d0
  do i=1,nvec
     H01(i,i)=E0(i)-escf
  enddo

  ! Transformation to the model state basis
  H01_mod=matmul(transpose(ref2mod),matmul(H01,ref2mod))

!----------------------------------------------------------------------
! Read in the A-vectors
!----------------------------------------------------------------------
  ! Open the A-vector file
  iscratch=scrunit(Ascr)
  open(iscratch,file=scrname(Ascr),form='unformatted',status='old')

  ! Check on the dimensions
  read(iscratch) dim2
  read(iscratch) dim1
  if (dim1 /= cfg%csfdim .or. dim2 /= nvec) then
     errmsg='Error in gvvpt2_diab: incorrect dimensions found when ' &
          //'parsing the A-vector file'
     call error_control
  endif

  ! Read in the A-vectors
  read(iscratch) Avec
  
  ! Close the A-vector file
  close(iscratch)

!----------------------------------------------------------------------
! Transform the A-vectors to the model state representation
!----------------------------------------------------------------------
  allocate(Avec_mod(cfg%csfdim,nrootsR0))

  call dgemm('N','N',cfg%csfdim,nrootsR0,nvec,1.0d0,Avec,cfg%csfdim,&
       ref2mod,nvec,0.0d0,Avec_mod,cfg%csfdim)
  
!----------------------------------------------------------------------
! Calculate the diabatic potential matrix
!----------------------------------------------------------------------
! On exit, the Avec_mod array will hold the FOIS contribution to the
! diabatic states
!----------------------------------------------------------------------
  call heff_diab(refdim,cfg%csfdim,nrootsR0,hdiag,Avec_mod,H01_mod,&
       ireg,regfac,diabpot)

!----------------------------------------------------------------------
! Add in the model states to obtain the diabatic states and
! orthonormalise  
!----------------------------------------------------------------------
  ! Add in the model states
  do i=1,nrootsR0
     Avec_mod(1:refdim,i)=vec_mod(:,i)
  enddo

  ! Lowdin's symmetric orthonormalisation
  call symm_ortho(cfg%csfdim,nrootsR0,Avec_mod)
  
!----------------------------------------------------------------------
! Add E_SCF to the on-diagonal diabatic potential matrix elements
!----------------------------------------------------------------------
  do i=1,nrootsR0
     diabpot(i,i)=diabpot(i,i)+escf
  enddo
  
!----------------------------------------------------------------------
! On-diagonal diabatic potential matrix elements
!----------------------------------------------------------------------
  allocate(diabii(nrootsR0))

  do i=1,nrootsR0
     diabii(i)=diabpot(i,i)
  enddo
  
!----------------------------------------------------------------------
! Save the diabatic states to disk
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('diabvec'//'.mult'//trim(amult)//&
       '.sym'//trim(airrep),diab_vecfile)
  call register_scratch_file(diab_vecscr,diab_vecfile)

  ! Open the scratch file
  iscratch=scrunit(diab_vecscr)
  open(iscratch,file=scrname(diab_vecscr),form='unformatted',&
       status='unknown')

  ! Dimensions
  write(iscratch) cfg%csfdim
  write(iscratch) nrootsR0

  ! On-diagonal diabatic potentials (included so as to be able to use
  ! the pre-existing eigenpair parsing routines to read these files)
  write(iscratch) diabii

  ! Diabatic state vectors
  do i=1,nrootsR0
     write(iscratch) Avec_mod(:,i)
  enddo
     
  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Write the configuration scratch file for the diabatic states
! At this point, this will simply be the same as for the (adiabatic)
! MRCI configurations; however, due to the use of wave function
! truncation or the like, this may change in later steps of the
! calculation
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('diabconf.mult'//trim(amult)//&
       '.sym'//trim(airrep),diab_vecfile)
  call register_scratch_file(diab_confscr,diab_vecfile)

  ! Open the scratch file
  iscratch=scrunit(diab_confscr)
  open(iscratch,file=scrname(diab_confscr),form='unformatted',&
       status='unknown')
  
  ! Subspace dimensions
  write(iscratch) cfg%n_int_I
  write(iscratch) cfg%nmoI
  write(iscratch) cfg%nmoE
  write(iscratch) cfg%confdim
  write(iscratch) cfg%nR
  write(iscratch) cfg%n1h
  write(iscratch) cfg%n2h
  write(iscratch) cfg%n0h
  write(iscratch) cfg%n1I
  write(iscratch) cfg%n2I
  write(iscratch) cfg%n1E
  write(iscratch) cfg%n2E
  write(iscratch) cfg%n1I1E
  
  ! Configuration information
  write(iscratch) cfg%confR
  write(iscratch) cfg%conf1h
  write(iscratch) cfg%a1h
  write(iscratch) cfg%off1h
  write(iscratch) cfg%conf2h
  write(iscratch) cfg%a2h
  write(iscratch) cfg%off2h
  write(iscratch) cfg%conf0h
  if (cfg%n1I > 0) then
     write(iscratch) cfg%a1I
     write(iscratch) cfg%off1I
  endif
  if (cfg%n2I > 0) then
     write(iscratch) cfg%a2I
     write(iscratch) cfg%off2I
  endif
  if (cfg%n1E > 0) then
     write(iscratch) cfg%a1E
     write(iscratch) cfg%off1E
  endif
  if (cfg%n2E > 0) then
     write(iscratch) cfg%a2E
     write(iscratch) cfg%off2E
  endif
  if (cfg%n1I1E > 0 ) then
     write(iscratch) cfg%a1I1E
     write(iscratch) cfg%off1I1E
  endif
  
  ! MO mapping arrays
  write(iscratch) cfg%m2c
  write(iscratch) cfg%c2m
  
  ! Number of CSFs as a function of the the number of open shells
  write(iscratch) cfg%ncsfs
  
  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Number of diabatic configurations
! Again, this may change in later steps of the calculation  
!----------------------------------------------------------------------
  diab_nconf=cfg%confdim
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  if (verbose) &
       call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'gvvpt2_diab')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine gvvpt2_diab
