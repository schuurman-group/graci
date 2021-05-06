!**********************************************************************
! Routines for the diagonalisation if the MRCI reference space
! Hamiltonian
!**********************************************************************

!######################################################################
! ref_diag_mrci: Diagonalisation of the reference space Hamiltonian
!######################################################################
#ifdef CBINDING
subroutine ref_diag_mrci(irrep,nroots,confscr,nconf,vecscr) &
     bind(c,name="ref_diag_mrci")
#else
subroutine ref_diag_mrci(irrep,nroots,confscr,nconf,vecscr)
#endif

  use constants
  use bitglobal
  use hbuild_double
  use ref_guess
  use full_diag
  use gendav
  use utils
  use iomod
  use conftype
  
  implicit none

  ! Irrep number and the requested number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(inout) :: nroots
  
  ! Array of reference configuration scratch file numbers
  integer(is), intent(in)    :: confscr(0:nirrep-1)

  ! Array of numbers of referecne configurations
  integer(is), intent(in)    :: nconf(0:nirrep-1)

  ! Eigenvector scratch file index
  integer(is), intent(out)   :: vecscr
  
  ! Number of configurations found in the scratch file
  integer(is)                :: nconf1

  ! Number of 64-bit integers required to represent the configurations
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
  integer(is), parameter     :: full_lim=1200
  
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

  ! Everything else
  integer(is)                :: i
  integer(is)                :: isigma(3)
  
  ! Dummy MRCI configuration derived type: temporarily required
  ! to be passed to the Davidson routines
  type(mrcfg)                :: cfg
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,72a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'Reference space diagonalisation in the',&
       trim(irreplbl(irrep,ipg)),'subspace'
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
  call hii_double(nconf1,hdim,offset,conf,sop,n_int_I,m2c,nmoI,&
       irrep,hii,averageii)
  
!----------------------------------------------------------------------
! Save to disk the non-zero off-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------
  call save_hij_double(nconf1,hdim,offset,averageii,conf,sop,n_int_I,&
       m2c,nmoI,irrep,hscr,nrec,'hij_ref')

!----------------------------------------------------------------------
! Full diagonalisation
!----------------------------------------------------------------------
  if (hdim <= full_lim) call diag_full(hscr,nrec,hdim,hii,irrep,&
       nroots,vecscr,'refvec')

!----------------------------------------------------------------------
! Iterative diagonalisation
!----------------------------------------------------------------------
  if (hdim > full_lim) then

     ! Temporary hard wiring of parameters
     blocksize=nroots*2
     maxvec=4*blocksize
     niter=100
     tol=1e-4_dp
     
     ! Dimension of the subspace used to generate the guess vectors
     guessdim=750
     if (guessdim > hdim) guessdim=int(hdim*0.9d0)
     
     ! Guess vector generation
     call ref_guess_subspace(guessscr,blocksize,hdim,guessdim,hii,&
          nrec,hscr)

     ! Set the sigma-vector algorithm information: disk-based only
     ! for now
     isigma(1)=1
     isigma(2)=hscr
     isigma(3)=nrec

     ! Set the preconditioner: DPR makes sense due to the small
     ! no. CSFs in a reference space diagonalisation
     ipre=1
     
     ! Generalised Davidson diagonalisation
     call generalised_davidson(irrep,isigma,cfg,hdim,hii,guessscr,&
          vecscr,'refvec',nroots,blocksize,maxvec,tol,&
          niter,ipre)

  endif

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(offset)
  deallocate(hii)
  deallocate(averageii)
 
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine ref_diag_mrci

!######################################################################
! basis_dimensions: Determines: (1) the total number of CSFs
!                   generated by the reference space configurations,
!                   (hdim) and; (2) the starting points for the CSFs
!                   generated by each configuration (offset)
!######################################################################
subroutine basis_dimensions(hdim,offset,sop,n_int_I,nconf)

  use constants
  use bitglobal
  
  implicit none

  integer(is), intent(out) :: hdim
  integer(is), intent(in)  :: n_int_I,nconf
  integer(is), intent(out) :: offset(nconf+1)
  integer(ib), intent(in)  :: sop(n_int_I,2,nconf)
  integer(is), allocatable :: nopen(:)
  integer(is)              :: ndet
  integer(is)              :: i,k,sum

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(nopen(nconf))
  
!----------------------------------------------------------------------
! Number of open shells for each configuration
!----------------------------------------------------------------------
  ! Initialisation
  nopen=0

  ! Loop over configurations
  do i=1,nconf

     ! Number of open shells
     do k=1,n_int_I
        nopen(i)=nopen(i)+popcnt(sop(k,1,i))
     enddo
     
  enddo
  
!----------------------------------------------------------------------
! Total number of CSFs and determinants
!----------------------------------------------------------------------
  ! Initialisation
  hdim=0
  ndet=0
  
  ! Loop over configurations
  do i=1,nconf

     ! Number of CSFs generated by the current configuration
     hdim=hdim+ncsfs(nopen(i))

     ! Number of determinants generated by the current configuration
     ndet=ndet+ndets(nopen(i))
     
  enddo

!----------------------------------------------------------------------
! Offsets
!----------------------------------------------------------------------
  sum=1
  offset=0

  do i=1,nconf
     offset(i)=sum
     sum=sum+ncsfs(nopen(i))
  enddo

  offset(nconf+1)=hdim+1

!----------------------------------------------------------------------
! Output the reference space dimensions
!----------------------------------------------------------------------
  write(6,'(/,x,a,x,i0)') &
       'Number of configurations in the reference space:',nconf
  write(6,'(x,a,x,i0)') &
       'Number of CSFs in the reference space:',hdim
  write(6,'(x,a,x,i0)') &
       'Number of determinants in the reference space:',ndet

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(nopen)
  
  return
  
end subroutine basis_dimensions
