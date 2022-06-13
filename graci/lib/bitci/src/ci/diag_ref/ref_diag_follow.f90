!**********************************************************************
! Top level routines for the diagonalisation of the MRCI reference
! space Hamiltonian with root following
!**********************************************************************

!######################################################################
! ref_diag_mrci_follow: Diagonalisation of the reference space
!                       Hamiltonian with root following based on
!                       overlaps with an input set of wave functions,
!                       detR0/vecR0, given in the Slater determinant
!                       representation
!######################################################################
#ifdef CBINDING
subroutine ref_diag_mrci_follow(irrep,nroots,confscr,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smat,nconf,vecscr) &
     bind(c,name="ref_diag_mrci_follow")
#else
subroutine ref_diag_mrci_follow(irrep,nroots,confscr,n_intR0,ndetR0,&
     nrootsR0,detR0,vecR0,nmoR0,smat,nconf,vecscr)
#endif

  use constants
  use bitglobal
  use hbuild_double
  use ref_guess
  use full_diag
  use csf2det
  use mrciutils
  use utils
  use iomod
  use conftype
  
  implicit none

  ! Irrep number and the requested number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(inout) :: nroots
  
  ! Array of reference configuration scratch file numbers
  integer(is), intent(in)    :: confscr(0:nirrep-1)

  ! Eigenvectors to follow, expressed in a Slater determinant basis
  integer(is), intent(in)    :: n_intR0,ndetR0,nrootsR0
  integer(ib), intent(in)    :: detR0(n_intR0,2,ndetR0)
  real(dp), intent(in)       :: vecR0(ndetR0,nrootsR0)

  ! MO overlaps
  integer(is), intent(in)    :: nmoR0
  real(dp), intent(in)       :: smat(nmoR0,nmo)
  
  ! Array of numbers of reference configurations
  integer(is), intent(in)    :: nconf(0:nirrep-1)

  ! Eigenvector scratch file index
  integer(is), intent(out)   :: vecscr

  ! Number of configurations found in the scratch file
  integer(is)                :: nconf1

  ! Number of (n_bits)-bit integers required to represent
  ! the configurations
  integer(is)                :: n_int_I

  ! Number of internal and external MOs
  integer(is)                :: nmoI,nmoE

  ! Reference configurations and SOPs
  integer(ib), allocatable   :: conf(:,:,:),sop(:,:,:)
  
  ! MO mapping arrays
  integer(is)                :: m2c(nmo),c2m(nmo)

  ! Numbers of CSFs
  integer(is)                :: hdim
  integer(is), allocatable   :: offset(:)

  ! Dimension of the CSF basis past which full diagonalisation
  ! will not be used
  integer(is), parameter     :: full_lim=3000
  
  ! Guess vector scratch file number
  integer(is)                :: guessscr

  ! Subspace dimension to be used in the generation of the
  ! guess vectors
  integer(is)                :: guessdim
  
  ! Davidson variables
  integer(is)                :: blocksize,maxvec,ipre,niter
  real(dp)                   :: tol
  
  ! On-diagonal Hamiltonian matrix elements
  real(dp), allocatable      :: hii(:)

  ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
  real(dp), allocatable      :: averageii(:)
    
  ! Hamiltonian scratch file number
  integer(is)                :: hscr
  
  ! Number of Hamiltonian scratch file records
  integer(is)                :: nrec
  
  ! Dummy MRCI configuration derived type: temporarily required
  ! to be passed to the Davidson routines
  type(mrcfg)                :: cfg

  ! Reference space eigenpairs
  integer(is)                :: ndet,ncsf
  integer(ib), allocatable   :: det(:,:,:)
  real(dp), allocatable      :: vec_csf(:,:),vec_det(:,:)
  real(dp), allocatable      :: ener(:)

  ! Norm-based truncation threshold
  real(dp)                   :: normthrsh

  ! Wave function overlaps
  integer(is)                :: npairs
  integer(is), allocatable   :: ipairs(:,:)
  real(dp), allocatable      :: Sij(:)
  integer(is), allocatable   :: isel(:)
  real(dp), allocatable      :: maxSij(:)
  
  ! Everything else
  integer(is)                :: i,j,k,nsave
  integer(is)                :: isigma(3)

  ! TEMPORARY
  integer(is) :: ncore=1
  integer(is) :: icore(1)
  logical(is) :: lfrzcore
  ! TEMPORARY
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,72a)') ('-',i=1,52)
  write(6,'(x,a,/,2(x,a))') &
       'Root-following reference space diagonalisation in',&
       'the '//trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(72a)') ('-',i=1,52)

!----------------------------------------------------------------------
! Return if there are no configurations for the current irrep
!----------------------------------------------------------------------
  if (nconf(irrep) == 0) then
     write(6,'(/,x,a)') 'No reference space configurations of '&
          //trim(irreplbl(irrep,ipg))//' symmetry'
     nroots=0
     return
  endif

!----------------------------------------------------------------------
! Read the configurations from file
!----------------------------------------------------------------------
  call read_ref_confs(confscr(irrep),nconf1,n_int_I,nmoI,nmoE,conf,&
       sop,m2c,c2m)

!----------------------------------------------------------------------
! Sanity check on the number of configurations
!----------------------------------------------------------------------
  if (nconf(irrep) /= nconf1) then
     errmsg='Error in ref_diag_mrci: '&
          //'inconsistent configuration numbers'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Determine the total number of CSFs and the offsets for each
! configuration
!----------------------------------------------------------------------
  allocate(offset(nconf1+1))
  offset=0

  call basis_dimensions(hdim,offset,sop,n_int_I,nconf1)

!----------------------------------------------------------------------
! Sanity check on the number of roots
!----------------------------------------------------------------------
  if (hdim < nroots) then
     errmsg='Error in ref_diag_mrci: N_roots > N_CSF'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(hii(hdim))
  allocate(averageii(nconf1))
  hii=0.0d0
  averageii=0.0d0

!----------------------------------------------------------------------
! Compute the on-diagonal Hamiltonian matrix elements and the
! their spin-coupling averaged values
!----------------------------------------------------------------------
  call hii_double(nconf1,hdim,offset,conf,sop,n_int_I,m2c,irrep,hii,&
       averageii)
  
!----------------------------------------------------------------------
! Save to disk the non-zero off-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call save_hij_double(nconf1,hdim,offset,averageii,conf,sop,n_int_I,&
       m2c,irrep,hscr,nrec,'hij_ref')

!----------------------------------------------------------------------
! In the future, we will implement a root-homing Davidson algorithm,
! but for now, we will only support full diagonalisation here
!----------------------------------------------------------------------
  if (hdim <= full_lim) then
     ! We will request that all roots are saved here
     ! Afterwards, we will find the ones that have the greatest
     ! overlap with the input wave functions
     nsave=hdim
     call diag_full(hscr,nrec,hdim,hii,irrep,&
          nsave,vecscr,'refvec')
  else
     errmsg='Error in ref_diag_mrci_follow: hdim > full_lim and ' &
          //'iterative diagonalisation is not yet supported'
     call error_control
  endif

!----------------------------------------------------------------------
! Get the determinant representations of the reference space wave
! functions
!----------------------------------------------------------------------
  ! Set up a dummy MRCI configuration derived data type to hold
  ! the ref conf info
  allocate(cfg%conf0h(n_int_I,2,nconf1))
  allocate(cfg%sop0h(n_int_I,2,nconf1))
  cfg%n0h=nconf1
  cfg%n_int_I=n_int_I
  cfg%n1I=0
  cfg%n2I=0
  cfg%n1E=0
  cfg%n2E=0
  cfg%n1I1E=0  
  cfg%conf0h=conf
  cfg%sop0h=sop
  
  ! Determine the dimension of the determinant basis
  call get_detdim(cfg,ndet)
  
  ! Allocate arrays
  ncsf=hdim
  allocate(det(n_int,2,ndet))
  allocate(vec_det(ndet,nsave))
  allocate(vec_csf(ncsf,nsave))
  allocate(ener(nsave))
  det=0_ib
  vec_det=0.0d0
  vec_csf=0.0d0
  ener=0.0d0
  
  ! Read in the reference space eigenvectors in the CSF basis
  call read_all_eigenpairs(vecscr,vec_csf,ener,ncsf,nsave)
  
  ! Get the determinant bit strings
  call bitstrings_detbas(cfg,ndet,det)

  ! Put the determinant bit strings into the 'canonical' MO ordering
  call reorder_confs(m2c,det,ndet)

  ! Compute the eigenvectors in the determinant basis
  call eigenvectors_detbas(cfg,nsave,ncsf,ndet,vec_csf,vec_det)

!----------------------------------------------------------------------
! Compute the overlaps with the input wave functions
!----------------------------------------------------------------------
  ! Allocate arrays
  npairs=nrootsR0*nsave
  allocate(Sij(npairs))
  allocate(ipairs(npairs,2))
  Sij=0.0d0
  ipairs=0
  
  ! Temporarily hard-wired truncation threshold
  normthrsh=1.0

  ! Fill in the array of bra-ket overlaps required
  k=0
  do i=1,nrootsR0
     do j=1,nsave
        k=k+1
        ipairs(k,1)=i
        ipairs(k,2)=j
     enddo
  enddo

  ! TEMPORARY HACK
  lfrzcore=.true.
  icore(1)=1
  ! TEMPORARY HACK

  ! Compute the overlaps
  call overlap(nmoR0,nmo,n_intR0,n_int,ndetR0,ndet,nrootsR0,nsave,&
       detR0,det,vecR0,vec_det,smat,normthrsh,ncore,icore,lfrzcore,&
       npairs,Sij,ipairs)

!----------------------------------------------------------------------
! Determine the reference eigenfunctions of interest
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(isel(nroots))
  allocate(maxSij(nroots))
  isel=0
  maxSij=0.0d0
  
  ! Select the reference space eigenfunctions with the greatest
  ! overlaps with the input wave functions
  
  
!----------------------------------------------------------------------
! Re-phasing
!----------------------------------------------------------------------

  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine ref_diag_mrci_follow
