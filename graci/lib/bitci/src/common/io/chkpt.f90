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

  use hdf5 
  use constants
  implicit none

  interface write_dataset
     module procedure write_dataset_dble
     module procedure write_dataset_int64
  end interface

  interface read_dataset
     module procedure read_dataset_dble
     module procedure read_dataset_int64
  end interface

  interface write_attribute
     module procedure write_attribute_int
  end interface 

  contains

  !
  ! query the dimensions (dims) of that dataset 'data_name'
  ! in the chkpt file 'file_name'
  !
  subroutine chkpt_dset_dims(file_name, data_name, dims)

    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    integer(is), intent(out)       :: dims(2)

    logical                        :: file_exists
    logical                        :: dset_exists
    integer                        :: rank

    integer(hid_t)                 :: file_id, space_id, dset_id
    integer(hsize_t), dimension(2) :: dims_actual, dims_max
    integer                        :: hdferr

    ! initialize to zero
    dims = (/0, 0/)

    ! assume file and data_set don't already exist
    file_exists = .false.
    dset_exists = .false.

    ! initialize the fortran interface
    call h5open_f(hdferr)

    inquire(file=file_name, exist=file_exists)

    if(file_exists) then
      call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, hdferr)

      ! check if dataset exists
      call h5lexists_f(file_id, data_name, dset_exists, hdferr)

      if (dset_exists) then
          call h5dopen_f(file_id, data_name, dset_id, hdferr)
          call h5dget_space_f(dset_id, space_id, hdferr)
          call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
          if (rank /= 2) then
            print *,'error in checking_dataset: existing rank /= 2'
            call exit(1)
          endif
          ! data set is too small, resize    
          call h5sget_simple_extent_dims_f(space_id, dims_actual, dims_max, hdferr)
          dims = dims_actual
      endif

    endif

    ! Close and release resources.
    call h5dclose_f(dset_id , hdferr)
    call h5sclose_f(space_id, hdferr)
    call h5fclose_f(file_id , hdferr)
    call h5close_f(hdferr)

  end subroutine chkpt_dset_dims

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

!#########################################################################3

  function group_exists(file_name, grp_name) result(exists)
    character(len=255), intent(in) :: file_name
    character(len=255), intent(in) :: grp_name
    
    integer(hid_t)                 :: file_id  ! File identifier 
    integer(hid_t)                 :: grp_id ! Group identifier 

    character(len=255)             :: f_name
    character(len=255)             :: g_name
    logical                        :: exists
    integer(is)                    :: error

    f_name    = trim(adjustl(file_name))
    g_name    = trim(adjustl(grp_name))

    ! if file doesn't exist, the group doesn't exist. Not sure if this
    ! default behavior makes most sense
    inquire(file=f_name, exist=exists)

    if(.not.exists) return

    exists    = .false.

    ! else, open the file try to open the group
    call h5open_f(error)

    call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, error)

    ! turn off error reporting -- we're going to explicitly check error code
    call h5eset_auto_f(0, error)
    call h5gopen_f(file_id, g_name, grp_id, error)

    ! if we don't get an error code, group doesn't exist
    if(error == 0) then
      exists = .true.
      call h5gclose_f(grp_id, error)
    else
      exists = .false.
    endif

    ! turn error reporting back on -- even though we're about to shut
    ! everything down
    call h5eset_auto_f(1, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    return
  end function group_exists


  subroutine create_group(file_name, grp_name)
    character(len=255), intent(in) :: file_name
    character(len=255), intent(in) :: grp_name

    integer(hid_t)                 :: file_id  ! File identifier 
    integer(hid_t)                 :: grp_id ! Group identifier 

    character(len=255)             :: f_name
    character(len=255)             :: g_name

    logical                        :: file_exists
    integer(is)                    :: error ! Error flag
    
    
    f_name    = trim(adjustl(file_name))
    g_name    = trim(adjustl(grp_name))

    ! group already exists, nothing to do
    if( group_exists(f_name, grp_name) ) return

    ! else, open the file and create group
    call h5open_f(error)     

    ! if file exists, append, else create a new one
    inquire(file=f_name, exist=file_exists)

    if(file_exists) then
      call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, error)
    else
      call h5fcreate_f(f_name, H5F_ACC_TRUNC_F, file_id, error)
    endif

    !
    ! Create a group named "grp_name" in the file.
    !
    call h5gcreate_f(file_id, grp_name, grp_id, error)

    !
    ! Close the group.
    !
    call h5gclose_f(grp_id, error)

    !
    ! Terminate access to the file.
    !
    call h5fclose_f(file_id, error)

    !
    ! Close FORTRAN interface.
    !
    call h5close_f(error)

    return
  end subroutine create_group

  function dataset_exists(file_name, dset_name) result(exists)
    character(len=255), intent(in) :: file_name
    character(len=255), intent(in) :: dset_name

    integer(hid_t)                 :: file_id  ! File identifier 
    integer(hid_t)                 :: dset_id ! Dataset identifier 

    character(len=255)             :: f_name
    character(len=255)             :: d_name
    logical                        :: exists
    integer(is)                    :: error

    f_name    = trim(adjustl(file_name))
    d_name    = trim(adjustl(dset_name))

    ! if file doesn't exist, the group doesn't exist. Not sure if this
    ! default behavior makes most sense
    inquire(file=f_name, exist=exists)

    if(.not.exists) return

    exists    = .false.

    ! else, open the file try to open the group
    call h5open_f(error)

    ! open file
    call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, error)

    ! check if link exists
    call h5lexists_f(file_id, d_name, exists, error)

    ! close file
    call h5fclose_f(file_id, error)
    call h5close_f(error)

    ! return true if exists and no errors
    exists = exists .and. error == 0

    return
  end function dataset_exists

! # write a 2D array dataset to the specified file_name
  subroutine write_dataset_dble(file_name, data_name, dims, data_set)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    integer(is), intent(in)        :: dims(2)
    real(dp), intent(in)           :: data_set(:,:)

    logical                        :: file_exists
    logical                        :: dset_exists
    integer(is)                    :: rank 
    
    integer(hid_t)                 :: file_id, space_id, dset_id
    integer(hsize_t), dimension(2) :: dset_dims
    integer(hsize_t), dimension(2) :: test_dim, test_max
    character(len=255)             :: f_name
    character(len=255)             :: dset_name
    integer                        :: hdferr

    ! cast dims to the appropriate integer size
    dset_dims = dims
    f_name    = trim(adjustl(file_name))
    dset_name = trim(adjustl(data_name))

    ! assume file and data_set don't already exist
    file_exists = .false.
    dset_exists = .false.

    ! initialize the fortran interface
    call h5open_f(hdferr)

    ! if file exists, append, else create a new one
    inquire(file=f_name, exist=file_exists)

    ! open and change geom dataset if file exists, else create
    if(file_exists) then
      call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, hdferr)
    else
      call h5fcreate_f(f_name, H5F_ACC_TRUNC_F, file_id, hdferr)
    endif

    ! check if dataset exists, if not create it:
    if (file_exists) call h5lexists_f(file_id, dset_name, dset_exists, hdferr)

    ! if data set exists, make sure it's big enough. If not, resize.
    ! It would seem the simpler thing would be to delete the dataset and recreate
    ! if it already exists, but HDF5 doesn't seem to make that straightforward.
    ! Specifically, you can't reclaim the space, so if done repeatedly, this 
    ! would lead to file size bloat.
    ! It may be worth revisiting this in the future.
    if (dset_exists) then
        call h5dopen_f(file_id, dset_name, dset_id, hdferr)
        call h5dget_space_f(dset_id, space_id, hdferr)
        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (rank /= 2) then
          print *,'error in write_dataset: existing rank /= 2'
          call exit(1)  
        endif    
        ! data set is too small, resize    
        call h5sget_simple_extent_dims_f(space_id, test_dim, test_max, hdferr) 
        if (any(dset_dims > test_max)) call h5sset_extent_simple_f(space_id, 2, dset_dims, dset_dims, hdferr)
    else
        ! create the datasapce
        call h5screate_simple_f(2, dset_dims, space_id, hdferr)
        call h5dcreate_f(file_id, dset_name, H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
    endif

    ! write the dataset
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data_set, dset_dims, hdferr)

    ! Close and release resources.
    call h5dclose_f(dset_id , hdferr)
    call h5sclose_f(space_id, hdferr)
    call h5fclose_f(file_id , hdferr)
    call h5close_f(hdferr)

  end subroutine write_dataset_dble

  ! # write a 2D array dataset to the specified file_name
  subroutine write_dataset_int64(file_name, data_name, dims, data_set)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    integer(is), intent(in)        :: dims(2)
    integer(ib), intent(in)        :: data_set(:,:)

    logical                        :: file_exists
    logical                        :: dset_exists
    integer(is)                    :: rank


    integer(hid_t)                 :: file_id, space_id, dset_id
    integer(hsize_t), dimension(2) :: dset_dims
    integer(hsize_t), dimension(2) :: test_dim, test_max
    character(len=255)              :: f_name
    character(len=255)              :: dset_name
    integer                        :: hdferr

    ! cast dims to the appropriate integer size
    dset_dims = dims
    f_name    = trim(adjustl(file_name))
    dset_name = trim(adjustl(data_name))

    ! assume file and data_set don't already exist
    file_exists = .false.
    dset_exists = .false.

    ! initialize the fortran interface
    call h5open_f(hdferr)

    ! if file exists, append, else create a new one
    inquire(file=f_name, exist=file_exists)

    ! open and change geom dataset if file exists, else create
    if(file_exists) then
      call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, hdferr)
    else
      call h5fcreate_f(f_name, H5F_ACC_TRUNC_F, file_id, hdferr)
    endif

    ! check if dataset exists, if not create it:
    if (file_exists) call h5lexists_f(file_id, dset_name, dset_exists, hdferr)

    ! if data set exists, make sure it's big enough. If not, resize.
    ! It would seem the simpler thing would be to delete the dataset and recreate
    ! if it already exists, but HDF5 doesn't seem to make that straightforward.
    ! Specifically, you can't reclaim the space, so if done repeatedly, this 
    ! would lead to file size bloat.
    ! It may be worth revisiting this in the future.
    if (dset_exists) then
        call h5dopen_f(file_id, dset_name, dset_id, hdferr)
        call h5dget_space_f(dset_id, space_id, hdferr)
        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (rank /= 2) then
          print *,'error in write_dataset: existing rank /= 2'
          call exit(1)
        endif
        ! data set is too small, resize    
        call h5sget_simple_extent_dims_f(space_id, test_dim, test_max, hdferr)
        if (any(dset_dims > test_max)) call h5sset_extent_simple_f(space_id, 2, dset_dims, dset_dims, hdferr)
    else
        ! create the datasapce
        call h5screate_simple_f(2, dset_dims, space_id, hdferr)
        call h5dcreate_f(file_id, dset_name, h5kind_to_type(ib,H5_INTEGER_KIND), space_id, dset_id, hdferr)
    endif

    ! write the dataset
    call h5dwrite_f(dset_id, h5kind_to_type(ib,H5_INTEGER_KIND), data_set, dset_dims, hdferr)

    ! Close and release resources.
    call h5dclose_f(dset_id , hdferr)
    call h5sclose_f(space_id, hdferr)
    call h5fclose_f(file_id , hdferr)
    call h5close_f(hdferr)

  end subroutine write_dataset_int64

  !
  !
  !
  subroutine read_dataset_dble(file_name, data_name, dims, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    integer(is), intent(in)        :: dims(2)
    real(dp), intent(out)          :: data_set(dims(1), dims(2))
    integer(is), intent(out)       :: dim_read(2)
    
    logical                        :: file_exists
    logical                        :: dset_exists
    integer(is)                    :: i, rank

    integer(hid_t)                 :: file_id, space_id, dset_id
    integer(hsize_t), dimension(2) :: test_dim, test_max
    character(len=255)              :: f_name
    character(len=255)              :: dset_name
    integer                        :: hdferr

    ! initialize to zero
    data_set = 0.
    f_name    = trim(adjustl(file_name))
    dset_name = trim(adjustl(data_name))

    ! assume file and data_set don't already exist
    file_exists = .false.
    dset_exists = .false.

    ! initialize the fortran interface
    call h5open_f(hdferr)

    ! if file exists, append, else create a new one
    inquire(file=f_name, exist=file_exists)

    ! open and change geom dataset if file exists, else create
    if(file_exists) then
      call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, hdferr)

      ! check if dataset exists, if not create it:
      call h5lexists_f(file_id, dset_name, dset_exists, hdferr)

      if (dset_exists) then
        call h5dopen_f(file_id, dset_name, dset_id, hdferr)
        call h5dget_space_f(dset_id, space_id, hdferr)
        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (rank /= 2) then
          print *,'error in write_dataset: existing rank /= 2'
          call exit(1)
        endif

        ! data set is too small, resize    
        call h5sget_simple_extent_dims_f(space_id, test_dim, test_max, hdferr)

        do i = 1,2
          test_dim(i) = min(test_dim(i), dims(i))
          dim_read(i) = test_dim(i)
        enddo

        ! read the dataset
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data_set(1:dim_read(1),1:dim_read(2)), test_dim, hdferr)

        ! Close and release resources.
        call h5dclose_f(dset_id , hdferr)
        call h5sclose_f(space_id, hdferr)

      endif

      call h5fclose_f(file_id , hdferr)
    endif

    call h5close_f(hdferr)

  end subroutine read_dataset_dble

  !
  !
  !
  subroutine read_dataset_int64(file_name, data_name, dims, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    integer(is), intent(in)        :: dims(2)
    integer(ib), intent(out)       :: data_set(dims(1), dims(2))
    integer(is), intent(out)       :: dim_read(2)

    logical                        :: file_exists
    logical                        :: dset_exists
    integer(is)                    :: i, rank

    integer(hid_t)                 :: file_id, space_id, dset_id
    integer(hsize_t), dimension(2) :: test_dim, test_max
    character(len=255)              :: f_name
    character(len=255)              :: dset_name
    integer                        :: hdferr

    ! initialize to zero
    data_set = 0.
    f_name    = trim(adjustl(file_name))
    dset_name = trim(adjustl(data_name))

    ! assume file and data_set don't already exist
    file_exists = .false.
    dset_exists = .false.

    ! initialize the fortran interface
    call h5open_f(hdferr)

    ! if file exists, append, else create a new one
    inquire(file=f_name, exist=file_exists)

    ! open and change geom dataset if file exists, else create
    if(file_exists) then

      ! initialize the fortran interface
      call h5open_f(hdferr)

      call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, hdferr)

      ! check if dataset exists, if not create it:
      call h5lexists_f(file_id, dset_name, dset_exists, hdferr)

      if (dset_exists) then
        call h5dopen_f(file_id, dset_name, dset_id, hdferr)
        call h5dget_space_f(dset_id, space_id, hdferr)
        call h5sget_simple_extent_ndims_f(space_id, rank, hdferr)
        if (rank /= 2) then
          print *,'error in read_dataset: existing rank /= 2'
          call exit(1)
        endif

        ! data set is too small, resize
        call h5sget_simple_extent_dims_f(space_id, test_dim, test_max, hdferr)

        do i = 1,2
          test_dim(i) = min(test_dim(i), dims(i))
        enddo
        dim_read = test_dim

        ! read the dataset
        call h5dread_f(dset_id, h5kind_to_type(ib,H5_INTEGER_KIND), data_set(1:dim_read(1),1:dim_read(2)), test_dim, hdferr)

        call h5dclose_f(dset_id , hdferr)
        call h5sclose_f(space_id, hdferr)
      endif

      call h5fclose_f(file_id , hdferr)
    endif

    call h5close_f(hdferr)
  end subroutine read_dataset_int64

  !
  !
  !
  subroutine write_attribute_int(file_name, data_name, attr_name, attr_val)
    use iso_c_binding

    character(len=255), intent(in)    :: file_name
    character(len=255), intent(in)    :: data_name
    character(len=10), intent(in)     :: attr_name
    integer(is), intent(in)          :: attr_val

    logical                          :: file_exists
    logical                          :: dset_exists
    character(len=255)                :: f_name
    character(len=255)                :: dset_name
    character(len=10)                 :: a_name

    integer(hsize_t)                 :: dims(1)
    integer(hid_t)                   :: file_id, dset_id
    integer(hid_t)                   :: space_id, attr_id 
    integer(is)                      :: hdferr
    integer(is),target               :: attr(1)
    type(c_ptr)                      :: f_ptr


    dims      = (/ 1 /)
    attr      = (/ attr_val /)

    f_name    = trim(adjustl(file_name))
    dset_name = trim(adjustl(data_name))
    a_name    = trim(adjustl(attr_name))

    ! assume file and data_set don't already exist
    file_exists = .false.
    dset_exists = .false.

    ! if file exists, append, else create a new one
    inquire(file=f_name, exist=file_exists)

    ! open and change geom dataset if file exists, else create
    if(file_exists) then

      ! initialize the fortran interface
      call h5open_f(hdferr)

      call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, hdferr)

      ! check if dataset exists
      call h5lexists_f(file_id, dset_name, dset_exists, hdferr)

      ! if data set exits, write attribute to it
      if (dset_exists) then

        ! open dataset to get data_set id
        !
        call h5dopen_f(file_id, dset_name, dset_id, hdferr)
  
        ! Create dataspace.  Setting maximum size to be the current size.
        !
        call h5screate_simple_f(1, dims, space_id, hdferr)

        ! Create the attribute and write the array data to it.
        !
        call H5Acreate_f(dset_id, a_name, H5T_NATIVE_INTEGER, space_id, attr_id, hdferr)

        f_ptr = c_loc(attr(1))
        call H5Awrite_f(attr_id, H5T_NATIVE_INTEGER, f_ptr, hdferr)
        !
        ! Close and release resources.
        !
        call H5Aclose_f(attr_id, hdferr)
        call H5Dclose_f(dset_id, hdferr)
      endif

      call H5Fclose_f(file_id, hdferr)
      call h5close_f(hdferr)

    endif

  end subroutine write_attribute_int

end module chkpt
