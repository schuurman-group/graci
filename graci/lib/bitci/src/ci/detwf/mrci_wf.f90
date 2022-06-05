!**********************************************************************
! Extraction of the determinant expansions of the MRCI wave functions
!**********************************************************************

!######################################################################
! wf_mrci: Controls the calculation and writing to an HDF5 file of a
!          a requested subset of the MRCI wave functions for a single
!          irrep
!######################################################################
#ifdef CBINDING
subroutine wf_mrci(irrep,nroots,iroots,confscr,vecscr,h5file_in,&
     grp_name,lbls,ndet) bind(c,name='wf_mrci')
#else
subroutine wf_mrci(irrep,nroots,iroots,confscr,vecscr,h5file_in,&
     grp_name,lbls,ndet)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use conftype
  use csf2det
  use det_expec
  use mrciutils
  use iomod
  use chkpt
  use timing
  
  implicit none

  ! Output HDF5 file name
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: h5file_in
  character(kind=C_CHAR), intent(in) :: grp_name
  character(len=255)                 :: h5file
  character(len=255)                 :: grp
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: h5file_in
  character(len=*), intent(in)       :: grp_name
  character(len=255)                 :: h5file
  character(len=255)                 :: grp
#endif
  
  ! Irrep and no. roots
  integer(is), intent(in)  :: irrep,nroots

  ! Indices of the states for which wave functions are requested
  integer(is), intent(in) :: iroots(nroots)
  
  ! MRCI configuration and eigenvector scratch file numbers
  integer(is), intent(in)  :: confscr(0:nirrep-1),vecscr(0:nirrep-1)

  ! integer labels for the wfs
  integer(is), intent(in)  :: lbls(nroots)

  ! No. determinants
  integer(is), intent(out) :: ndet
  
  ! MRCI configuration derived type
  type(mrcfg)              :: cfg

  ! Eigenpairs in the CSF basis
  real(dp), allocatable    :: vec_csf(:,:),ener(:)

  ! Determinants and eigenvectors in the determinat basis
  integer(is)              :: detdim
  integer(ib), allocatable :: det(:,:,:)
  real(dp), allocatable    :: vec_det(:,:)

  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end
  
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
! If C bindings are on, then convert the output HDF5 filename from the
! C char type to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(h5file_in)
  call c2fstr(h5file_in,h5file,length)
  length=cstrlen(grp_name)
  call c2fstr(grp_name, grp, length)
#else
  h5file = adjustl(trim(h5file_in))
  grp    = adjustl(trim(grp_name))
#endif

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
! Get the determinant bit strings
!----------------------------------------------------------------------
  call bitstrings_detbas(cfg,detdim,det)

!----------------------------------------------------------------------
! Put the determinant bit strings into the 'canonical' MO ordering
! (the reorder_conf subroutine is used for this as the det bit strings
! have the same structure as conf bit strings)
!----------------------------------------------------------------------
  call reorder_confs(cfg%m2c,det,detdim)
  
!----------------------------------------------------------------------
! Compute the eigenvectors in the determinant basis
!----------------------------------------------------------------------
  call eigenvectors_detbas(cfg,nroots,cfg%csfdim,detdim,vec_csf,&
       vec_det)
  
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
! Write the determinant lists to the h5 checkpoint file
!----------------------------------------------------------------------
  do i=1,nroots
     call chkpt_write_wfn(h5file,grp,nmo,lbls(i),detdim,vec_det(:,i),&
          det)
  enddo

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
