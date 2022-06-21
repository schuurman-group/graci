!**********************************************************************
! Extraction of the determinant expansions of the MRCI wave functions
!**********************************************************************

!######################################################################
! wf_mrci: Given configuration and eigenvector scratch file
!          numbers, converts MRCI wave functions from the CSF to
!          determinant basis
!######################################################################
#ifdef CBINDING
subroutine wf_mrci(irrep,nroots,iroots,confscr,vecscr,ndet,wfscr) &
     bind(c,name='wf_mrci')
#else
subroutine wf_mrci(irrep,nroots,iroots,confscr,vecscr,ndet,wfscr)
#endif

  use constants
  use bitglobal
  use conftype
  use csf2det
  use det_expec
  use mrciutils
  use iomod
  use timing
  
  implicit none

  ! Irrep and no. roots
  integer(is), intent(in)  :: irrep,nroots

  ! Indices of the states for which wave functions are requested
  integer(is), intent(in) :: iroots(nroots)
  
  ! MRCI configuration and eigenvector scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1),vecscr(0:nirrep-1)

  ! No. determinants
  integer(is), intent(out) :: ndet

  ! Determinant wave function scratch file numbers
  integer(is), intent(out) :: wfscr
  
  ! MRCI configuration derived type
  type(mrcfg)              :: cfg

  ! Eigenpairs in the CSF basis
  real(dp), allocatable    :: vec_csf(:,:),ener(:)

  ! Determinants and eigenvectors in the determinat basis
  integer(is)              :: detdim
  integer(ib), allocatable :: det(:,:,:)
  real(dp), allocatable    :: vec_det(:,:)

  ! Truncated determinant and eigenvector arrays
  integer(is)              :: detdim_new
  integer(ib), allocatable :: det_new(:,:,:)
  real(dp), allocatable    :: vec_new(:,:)
  
  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end

  ! I/O variables
  integer(is)              :: iscratch
  character(len=60)        :: wffile
  character(len=2)         :: amult,airrep
  
  ! Everything else
  integer(is)              :: i
  real(dp), allocatable    :: S2expec(:)
  real(dp)                 :: S,S2
  real(dp), parameter      :: tiny=5e-12_dp
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'Wave function extraction for the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(52a)') ('-',i=1,52)
  
!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr(irrep))

!----------------------------------------------------------------------
! Determine the dimension of the determinant basis
!----------------------------------------------------------------------
  call get_detdim(cfg,detdim)

  ndet=detdim
  
!----------------------------------------------------------------------  
! Output the basis dimensions
!----------------------------------------------------------------------
  write(6,'(/,x,a,9x,i0)') 'CSF basis dimension:',cfg%csfdim
  write(6,'(x,a,x,i0)') 'Determinant basis dimension:',detdim

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(vec_csf(cfg%csfdim,nroots))
  allocate(vec_det(detdim,nroots))
  allocate(ener(nroots))
  allocate(det(n_int,2,detdim))
  allocate(S2expec(nroots))
  
!----------------------------------------------------------------------
! Read in the eigenvectors in the CSF basis
!----------------------------------------------------------------------
  ! Read in the eigenvectors
  call read_some_eigenpairs(vecscr(irrep),vec_csf,ener,cfg%csfdim,&
       nroots,iroots)

!----------------------------------------------------------------------
! Compute the determinant representation of the wave functions
!----------------------------------------------------------------------
  call det_trans(cfg,nroots,cfg%csfdim,detdim,vec_csf,vec_det,det)
  
!----------------------------------------------------------------------
! Put the determinant bit strings into the 'canonical' MO ordering
! (the reorder_conf subroutine is used for this as the det bit strings
! have the same structure as conf bit strings)
!----------------------------------------------------------------------
  call reorder_confs(cfg%m2c,det,detdim)
  
!----------------------------------------------------------------------
! Debugging: check that the determinant expansions are spin
! eigenfunctions
!----------------------------------------------------------------------
! This should be commented back in if, e.g., nomax is increased
!----------------------------------------------------------------------
!  ! Compute the S^2 expectation values
!  call s2expec_detbas(detdim,nroots,det,vec_det,S2expec)
!
!  ! Correct value
!  S=dble(imult-1)/2.0d0
!  S2=S*(S+1.0d0)
!
!  ! Check the expectation values
!  do i=1,nroots
!     if (abs(S2-S2expec(i)) > 5e-12_dp) then
!        errmsg='Error in wf_mrci: incorrect <S^2> value found'
!        call error_control
!     endif
!  enddo

!----------------------------------------------------------------------
! Write the determinant bit strings and eigenvectors to disk
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('detwf.mult'//trim(amult)//'.sym'//trim(airrep),&
       wffile)
  call register_scratch_file(wfscr,wffile)

  ! Open the scratch file
  iscratch=scrunit(wfscr)
  open(iscratch,file=scrname(wfscr),form='unformatted',status='unknown')

  ! Spin multiplicity
  write(iscratch) imult
  
  ! Dimensions
  write(iscratch) n_int
  write(iscratch) detdim
  write(iscratch) nroots
  write(iscratch) nmo

  ! Eigenvectors in the determinant basis
  do i=1,nroots
     write(iscratch) vec_det(:,i)
  enddo
     
  ! Determinant bit strings
  write(iscratch) det
  
  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  call cfg%finalise
  deallocate(vec_csf)
  deallocate(vec_det)
  deallocate(ener)
  deallocate(det)
  deallocate(S2expec)

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'wf_mrci')
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine wf_mrci
