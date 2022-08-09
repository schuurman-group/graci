!**********************************************************************
! Routines for the extraction of the Slater determinant representation
! of bitci MRCI-type wave functions
!**********************************************************************

!######################################################################
! detwf: top level routine for the conversion from the CSF to
!        determinant basis of all states of a given irrep, either
!        for the bra or ket manifold (bkstr = 'bra' or 'ket')
!######################################################################
#ifdef CBINDING
subroutine detwf(irrep,conffile_in,vecfile_in,nroots,bkstr_in,wfscr) &
     bind(c,name='detwf')
#else
subroutine detwf(irrep,conffile_in,vecfile_in,nroots,bkstr_in,wfscr)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use conftype
  use csf2det
  use det_expec
  use mrciutils
  use iomod
  use timing

  ! bitCI configuration and eigenvector file names
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: conffile_in(*),vecfile_in(*),&
                                        bkstr_in
  character(len=255)                 :: conffile,vecfile,bkstr
       
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: conffile_in,vecfile_in,bkstr_in
  character(len=255)                 :: conffile,vecfile,bkstr
#endif

  ! Irrep and no. roots
  integer(is), intent(in)  :: irrep,nroots

  ! Determinant wave function scratch file number
  integer(is), intent(out) :: wfscr

  ! MRCI configuration derived type
  type(mrcfg)              :: cfg
  
  ! Scratch file numbers
  integer(is)              :: confscr,vecscr

  ! Eigenpairs in the CSF basis
  real(dp), allocatable    :: vec_csf(:,:),ener(:)

  ! Determinants and eigenvectors in the determinant basis
  integer(is)              :: detdim
  integer(ib), allocatable :: det(:,:,:)
  real(dp), allocatable    :: vec_det(:,:)
  
  ! Timing variables
  real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                              twall_end
  
  ! Everything else
  integer(is)              :: i,iscratch
  real(dp), allocatable    :: S2expec(:)
  real(dp)                 :: S,S2
  real(dp), parameter      :: tiny=5e-12_dp
  character(len=250)       :: wffile
  character(len=2)         :: amult,airrep

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the bitCI configuration and
! eigenvector file names, and bra/ket string from the C char type
! to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(conffile_in)
  call c2fstr(conffile_in,conffile,length)
  length=cstrlen(vecfile_in)
  call c2fstr(vecfile_in,vecfile,length)
  length=cstrlen(bkstr_in)
  call c2fstr(bkstr_in,bkstr,length)
#else
  conffile=adjustl(trim(conffile_in))
  vecfileB=adjustl(trim(vecfileB_in))
  bkstr=adjustl(trim(bkstr_in))
#endif

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  if (verbose) then
     write(6,'(/,52a)') ('-',i=1,52)
     write(6,'(4(x,a))') &
          trim(bkstr),'CSF-to-det conversion for the',&
          trim(irreplbl(irrep,ipg)),'subspace'
     write(6,'(52a)') ('-',i=1,52)
  endif
     
!----------------------------------------------------------------------
! Check on the bra/ket string
!----------------------------------------------------------------------
  if (bkstr /= 'bra' .and. bkstr /= 'ket') then
     errmsg='Error in detwf: bkstr \notin {bra,ket}'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Register the configuration and eigenvector scratch files
!----------------------------------------------------------------------
  call register_scratch_file(confscr,conffile)
  call register_scratch_file(vecscr,vecfile)

!----------------------------------------------------------------------
! Set all globally accessible arrays that are used by the csf2det
! routines
!----------------------------------------------------------------------
! This is a little bit hacky, but these variables are so embedded
! in the bitX routines that we need to set them globally...
!----------------------------------------------------------------------
  ! These arrays _should_ be deallocated, but just in case...
  if (allocated(csfcoe)) deallocate(csfcoe)
  if (allocated(detvec)) deallocate(detvec)

  if (bkstr == 'bra') then
     ! Bra CSFs, etc.
     n_int=n_intB
     nmo=nmoB
     imult=imultB
     maxcsf=maxcsfB
     maxdet=maxdetB
     allocate(csfcoe(maxcsf,maxdet,nocase2))
     allocate(detvec(maxdet,nocase2))
     csfcoe=csfcoeB
     detvec=detvecB
     ncsfs=ncsfsB
     ndets=ndetsB
  else
     ! Ket CSFs, etc.
     n_int=n_intK
     nmo=nmoK
     imult=imultK
     maxcsf=maxcsfK
     maxdet=maxdetK
     allocate(csfcoe(maxcsf,maxdet,nocase2))
     allocate(detvec(maxdet,nocase2))
     csfcoe=csfcoeK
     detvec=detvecK
     ncsfs=ncsfsK
     ndets=ndetsK
  endif

!----------------------------------------------------------------------
! Set up the configuration derived type
!----------------------------------------------------------------------
  call cfg%initialise(irrep,confscr)

!----------------------------------------------------------------------
! Determine the dimension of the determinant basis
!----------------------------------------------------------------------
  call get_detdim(cfg,detdim)
  
!----------------------------------------------------------------------  
! Output the basis dimensions
!----------------------------------------------------------------------
  if (verbose) then
     write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim
     write(6,'(x,a,x,i0)') 'Determinant basis dimension:',detdim
  endif
     
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
  call read_all_eigenpairs(vecscr,vec_csf,ener,cfg%csfdim,nroots)

!----------------------------------------------------------------------
! Compute the determinant representation of the wave functions
!----------------------------------------------------------------------
  call det_trans(cfg,cfg%m2c,nroots,cfg%csfdim,detdim,vec_csf,&
       vec_det,det)
  
!----------------------------------------------------------------------
! Write the determinant representation of the wave functions to disk
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name(trim(bkstr)//'wf.mult'//trim(amult)//&
       '.sym'//trim(airrep),wffile)
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
  
!!----------------------------------------------------------------------
!! Debugging: check that the determinant expansions are spin
!! eigenfunctions
!!----------------------------------------------------------------------
!! This should be commented back in if, e.g., nomax is increased
!!----------------------------------------------------------------------
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
!        errmsg='Error in detwf: incorrect <S^2> value found'
!        call error_control
!     endif
!  enddo
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(csfcoe)
  deallocate(detvec)
  deallocate(vec_csf)
  deallocate(vec_det)
  deallocate(ener)
  deallocate(det)
  deallocate(S2expec)

!----------------------------------------------------------------------
! To be on the safe side, zero all globally accessible arrays that
! were set to the bra/ket values in this subroutine
!----------------------------------------------------------------------
  n_int=0
  nmo=0
  imult=0
  maxcsf=0
  maxdet=0
  ncsfs=0
  ndets=0
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  if (verbose) &
       call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'detwf')
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine detwf
