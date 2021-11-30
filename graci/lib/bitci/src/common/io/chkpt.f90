!######################################################################
! chkpt: a module to write csf and determinant lists to 
!        hdf5 file format
!---------------------------------------------------------------------
!           Reads and writes configurations/unique dets to hdf5 file, 
!           sorted by decreasing coefficient magnitude.
!           This module defines the file format (i.e. group/dataset names)
!           and should be kept consistent with the GRaCI python module
!           chkpt.py. The python module is responsible for writing majority 
!           of information in the checkpoint -- and so 'write' capacity here
!           is limited to det/csf lists. However, It may prove useful to add extensive
!           read capabilities to this module.
!---------------------------------------------------------------------
module chkpt 

  use h5_ops
  use constants
  implicit none

  contains

  !
  ! externally callable routine to write a batch
  ! of csfs or dets to hdf5 file. 
  !  - If the file doesn't exist, it is create
  !  - If the file does exist and newfile=.false., it is appended to
  !  - If the file does exist and newfile=.true., the old file is removed and
  !    a new file is created from the passed det/csf list 
  !
  !    file_name: name of file as character string (in)
  !     wfn_indx: (integer) the wfn label, given as an integer (in)
  !            n: (integer) number of det/csfs to write (in)
  !      cf_list: (double) array of coefficients (in)
  !    conf_list: (integer) array of bitci dets/confs (in)
  subroutine chkpt_write_wfn(file_name, grp_name, nmo, wfn_indx, ndet, cf_list, det_list)

    character(len=255), intent(in) :: file_name
    character(len=255), intent(in) :: grp_name
    integer(is), intent(in)       :: nmo
    integer(is), intent(in)       :: wfn_indx
    integer(is), intent(in)       :: ndet
    real(dp), intent(in)          :: cf_list(ndet)
    integer(ib), intent(in)       :: det_list(n_int, 2, ndet)

    character(len=3)              :: indx_str
    character(len=10)             :: attr
    character(len=255)            :: g_name
    character(len=255)            :: data_name
    integer(is)                   :: dims(2)
    real(dp)                      :: data_set_cf(1, ndet)
    integer(ib)                   :: data_set_det(2*n_int,ndet)

    ! name of data_set
    if (wfn_indx > 999) then
      print *,' trying to write wfn with index > 999'
      call exit(1)
    endif
    write(indx_str, '(i3)')wfn_indx

    ! create the group to put data in if necessary
    g_name = '/'//trim(adjustl(grp_name))
    if(.not.group_exists(file_name, g_name)) call create_group(file_name, g_name)

    ! set buffers for the wfn det_list coefficients
    dims              = (/ 1, ndet /)
    data_name         = trim(adjustl(g_name))//'/wfn_cf'//adjustl(indx_str)
    data_set_cf(1, :) = cf_list(:)    

    ! write wfn coefficients to dedicated dataset
    call write_dataset(file_name, data_name, dims, data_set_cf)

    ! Next write the determinant list -- assumes a single common determinant
    ! basis -- only bother writing if it doesn't exist
    data_name   = trim(adjustl(g_name))//'/wfn_dets'
    if(.not.dataset_exists(file_name, data_name)) then
      dims        = (/ 2*n_int, ndet /)
      data_set_det(:n_int, :ndet)          = det_list(:n_int, 1, :ndet)
      data_set_det(n_int+1:2*n_int, :ndet) = det_list(:n_int, 2, :ndet)
      ! write wfn bit strings to dedicated dataset
      call write_dataset(file_name, data_name, dims, data_set_det)

      attr = 'nmo'
      call write_attribute(file_name, data_name, attr, nmo)
    endif

  end subroutine chkpt_write_wfn

  !
  ! read_wfn: public
  !           Reads the contents of a wfn file into cf_list and conf_list
  !           file
  ! 
  subroutine chkpt_read_wfn(file_name, wfn_indx, ndet, cf_list, det_list)

    character(len=255), intent(in) :: file_name
    integer(is), intent(in)       :: wfn_indx
    integer(is), intent(in)       :: ndet
    real(dp), intent(out)         :: cf_list(ndet)
    integer(ib), intent(out)      :: det_list(n_int, 2, ndet)

    character(len=3)              :: indx_str
    character(len=255)             :: data_name
    integer(is)                   :: dims(2)
    integer(is)                   :: dim_read(2)
    real(dp)                      :: data_set_cf(ndet,1)
    integer(ib)                   :: data_set_det(2*n_int,ndet)

    cf_list  = 0.
    det_list = 0

    ! name of data_set
    if (wfn_indx >= 1000) then
      print *,' trying to read wfn with index > 999'
      call exit(1)
    endif
    write(indx_str, '(i3)')wfn_indx
    data_name        = 'wfn_cf'//adjustl(indx_str)
    dims             = (/ ndet, 1 /)

    ! read the coefficients first
    call read_dataset(file_name, trim(data_name), dims, data_set_cf, dim_read)

    if (any(dims /= dim_read)) then
      print *,'unexpected number of wfn coefficients'
    else
      cf_list(1:ndet) = data_set_cf(1:ndet, 1)
    endif

    ! now read the determinant list
    dims                                 = (/ 2*n_int, ndet /)
    data_name                            = 'wfn_dets'

    call read_dataset(file_name, data_name, dims, data_set_det, dim_read)

    if (any(dims /= dim_read)) then
      print *,'unexpected number of determinants'
    else
      det_list(:n_int, 1, :ndet) = data_set_det(1:n_int,         :ndet)
      det_list(:n_int, 2, :ndet) = data_set_det(n_int+1:2*n_int, :ndet)
    endif

  end subroutine chkpt_read_wfn

end module chkpt
