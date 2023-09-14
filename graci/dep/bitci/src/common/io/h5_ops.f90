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
module h5_ops 

  use hdf5 
  use constants
  implicit none

  interface read_dataset
     module procedure read_dataset_int64
     !module procedure read_dataset_float16
     module procedure read_dataset_float
     module procedure read_dataset_dble
  end interface

  interface read_dataset_t
     !module procedure read_dataset_t_float16
     module procedure read_dataset_t_float
     module procedure read_dataset_t_dble
  end interface

  contains

!#########################################################################3

  !
  !
  !
  function dataset_exists(file_name, dset_name) result(exists)
    character(len=255), intent(in) :: file_name
    character(len=255), intent(in) :: dset_name

    logical                         :: exists
    integer(hid_t)                  :: file_id, dset_id, space_id
    integer(hid_t)                  :: data_type
    integer(hid_t)                  :: h5dims(10)
    integer(is)                     :: rank
    integer(is)                     :: error

    call open_dataset(file_name, dset_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    exists = error.eq.0

    call close_dataset(file_id, dset_id, space_id)

    return
  end function dataset_exists

  !
  !
  !
  subroutine dataset_dims(file_name, dset_name, rank, dims)

    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: dset_name
    integer(is), intent(out)        :: rank
    integer(is), intent(inout)      :: dims(:)

    integer(hid_t)                  :: file_id, space_id, dset_id
    integer(hid_t)                  :: data_type
    integer(hid_t)                  :: h5dims(10)
    integer(is)                     :: error

    call open_dataset(file_name, dset_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    if error.eq.0 then
      rank = rank
      dims = h5dims(1:rank)
    else
      rank = -1
      dims = -1
    endif

    call close_dataset(file_id, dset_id, space_id)

  end subroutine dataset_dims

  !
  ! open a dataset and return all related IDs. If it doesn't exist, create it.
  !
  subroutine open_dataset(file_name, dset_name, data_type, file_id, dset_id, space_id, rank, dims, oerror, wt_dset)
    implicit none
    character(len=255), intent(in)      :: file_name   ! name of integral file
    character(len=255),intent(in)       :: dset_name   ! name of data set (i.e. h_core, eri) 
    integer(hid_t),intent(inout)        :: data_type   ! data_type (i.e. NATIVE_DOUBLE)
    integer(hid_t),intent(out)          :: file_id     ! file identifier
    integer(hid_t),intent(out)          :: dset_id     ! geeric dataset id
    integer(hid_t),intent(out)          :: space_id    ! data_space id
    integer(is), intent(inout)          :: rank        ! rank of dataset, will set to this value if it doesn't exist
    integer(hid_t), intent(inout)       :: dims(2)     ! dims of dataset, will set to this value if it doesn't exist
    integer(is), intent(out)            :: oerror
    logical, intent(in), optional       :: wt_dset

    character(len=255)                  :: f_name
    character(len=255)                  :: d_name
    logical                             :: file_exists
    logical                             :: dset_exists
    logical                             :: dset_wt
    integer(is)                         :: error
    integer(hsize_t)                    :: chk_sze(2)
    integer(hsize_t)                    :: chk_max(2)

    oerror = 0
    if(present(wt_dset)) then
      dset_wt = wt_dset
    else
      dset_wt = .false.
    endif

    f_name  = trim(adjustl(file_name))
    d_name  = trim(adjustl(dset_name))

    ! 
    ! initialize h5 system
    !
    call h5open_f(error)

    ! if file exists, append, else create a new one
    inquire(file=f_name, exist=file_exists)

    ! open and change geom dataset if file exists, else create
    if(file_exists) then
      call h5fopen_f(f_name, H5F_ACC_RDWR_F, file_id, error)
    elseif(dset_wt) then
      call h5fcreate_f(f_name, H5F_ACC_TRUNC_F, file_id, error)
    else
      oerror = 1
      return
    endif

    ! check if dataset exists, if not create it:
    if (file_exists) call h5lexists_f(file_id, d_name, dset_exists, error)

    ! if data set exists, make sure it's big enough. If not, resize.
    ! It would seem the simpler thing would be to delete the dataset and recreate
    ! if it already exists, but HDF5 doesn't seem to make that straightforward.
    ! Specifically, you can't reclaim the space, so if done repeatedly, this 
    ! would lead to file size bloat.
    ! It may be worth revisiting this in the future.
    if (dset_exists) then
        call h5dopen_f(file_id, d_name, dset_id, error)
        call h5dget_type_f(dset_id, data_type, error)
        call h5dget_space_f(dset_id, space_id, error)
        call h5sget_simple_extent_ndims_f(space_id, rank, error)
        call h5sget_simple_extent_dims_f(space_id, chk_sze, chk_max, error)
        ! if dataset too small and we want to resize, do so
        if (any(dims > chk_max).and.dset_wt) then
          call h5sset_extent_simple_f(space_id, rank, dims, dims, error)
        ! else, set dims to the current size of the dataset
        else
          dims = chk_sze
        endif

    elseif(dset_wt) then
        ! create the datasapce
        call h5screate_simple_f(rank, dims, space_id, error)
        call h5dcreate_f(file_id, d_name, data_type, space_id, dset_id, error)
    else
        oerror = 1
    endif

    return
  end subroutine open_dataset

  !
  !
  !
  subroutine close_dataset(file_id, dset_id, space_id)
    integer(hid_t), intent(in)     :: file_id
    integer(hid_t), intent(in)     :: dset_id
    integer(hid_t), intent(in)     :: space_id
  
    integer(is)                    :: error

    ! Close and release resources.
    call h5dclose_f(dset_id , error)
    call h5sclose_f(space_id, error)
    call h5fclose_f(file_id , error)
    call h5close_f(error)

  end subroutine close_dataset

  !
  !
  !
  subroutine read_dataset_int64(file_name, data_name, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    integer(ib), intent(inout)      :: data_set(:, :)
    integer(is), intent(out)        :: dim_read(2)

    integer(hid_t)                  :: file_id, space_id, dset_id
    integer(hid_t)                  :: h5dims(2)
    integer(hid_t)                  :: data_type
    integer(is)                     :: error
    integer(is)                     :: rank

    ! initialize to zero
    data_set = 0.
    rank     = 2
    h5dims   = shape(data_set)

    !data_type = h5kind_to_type(ib,H5_INTEGER_KIND)
    call open_dataset(file_name, data_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    if(any(h5dims > shape(data_set))) stop 'read_dataset_int64 - array too small to hold data set'//data_name
    dim_read = h5dims

    data_type = h5kind_to_type(ib,H5_INTEGER_KIND)
    call h5dread_f(dset_id, data_type, data_set(1:dim_read(1),1:dim_read(2)), h5dims, error)

    call close_dataset(file_id, dset_id, space_id)

  end subroutine read_dataset_int64

  !
  !
  !
  subroutine read_dataset_float16(file_name, data_name, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    real(hp), intent(inout)         :: data_set(:,:)
    integer(is), intent(out)        :: dim_read(2)
  
    integer(hid_t)                  :: data_type
    integer(hid_t)                  :: file_id, space_id, dset_id
    integer(hid_t)                  :: h5dims(2)
    integer(is)                     :: error
    integer(is)                     :: rank

    ! initialize to zero
    data_set = 0.
    rank     = 2

    !data_type = H5T_NATIVE_DOUBLE
    call open_dataset(file_name, data_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    ! if we're wrong about size of data set, stop with error
    if(any(h5dims > shape(data_set))) stop 'read_dataset_float - array too small to hold data set'//data_name
    dim_read = h5dims

    data_type = h5kind_to_type(hp,H5_REAL_KIND)
    call h5dread_f(dset_id, data_type, data_set(1:dim_read(1),1:dim_read(2)), h5dims, error)

    call close_dataset(file_id, dset_id, space_id)

  end subroutine read_dataset_float16

  !
  !
  !
  subroutine read_dataset_float(file_name, data_name, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    real(sp), intent(inout)         :: data_set(:,:)
    integer(is), intent(out)        :: dim_read(2)

    integer(hid_t)                  :: data_type
    integer(hid_t)                  :: file_id, space_id, dset_id
    integer(hid_t)                  :: h5dims(2)
    integer(is)                     :: error
    integer(is)                     :: rank

    ! initialize to zero
    data_set = 0.
    rank     = 2

    !data_type = H5T_NATIVE_DOUBLE
    call open_dataset(file_name, data_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    ! if we're wrong about size of data set, stop with error
    if(any(h5dims > shape(data_set))) stop 'read_dataset_float - array too small to hold data set'//data_name
    dim_read = h5dims

    data_type = h5kind_to_type(sp,H5_REAL_KIND)
    call h5dread_f(dset_id, data_type, data_set(1:dim_read(1),1:dim_read(2)), h5dims, error)

    call close_dataset(file_id, dset_id, space_id)
 
  end subroutine read_dataset_float

  !
  !
  !
  subroutine read_dataset_dble(file_name, data_name, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    real(dp), intent(inout)         :: data_set(:,:)
    integer(is), intent(out)        :: dim_read(2)

    integer(hid_t)                  :: data_type
    integer(hid_t)                  :: file_id, space_id, dset_id
    integer(hid_t)                  :: h5dims(2)
    integer(is)                     :: error
    integer(is)                     :: rank 

    ! initialize to zero
    data_set = 0.
    rank     = 2

    !data_type = H5T_NATIVE_DOUBLE
    call open_dataset(file_name, data_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    ! if we're wrong about size of data set, stop with error
    if(any(h5dims > shape(data_set))) stop 'read_dataset_dble - array too small to hold data set'//data_name
    dim_read = h5dims

    data_type = h5kind_to_type(dp,H5_REAL_KIND)
    call h5dread_f(dset_id, data_type, data_set(1:dim_read(1),1:dim_read(2)), h5dims, error)

    call close_dataset(file_id, dset_id, space_id)

  end subroutine read_dataset_dble

  !
  !
  !
  subroutine read_dataset_t_float16(file_name, data_name, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    real(hp), intent(inout)         :: data_set(:, :)
    integer(is), intent(out)        :: dim_read(2)

    integer(is)                     :: buf_sze
    integer(is)                     :: n_remain
    integer(is)                     :: n_read
    integer(is)                     :: n_read_tot
    integer(hsize_t)                :: offset(2)
    integer(hsize_t)                :: offmem(2)
    integer(hsize_t)                :: buf_read(2)
    real(sp),allocatable            :: int_buf(:,:)

    integer(hid_t)                  :: file_id, space_id, dset_id, buf_space
    integer(hid_t)                  :: h5dims(2)
    integer(hid_t)                  :: data_type
    integer(is)                     :: error
    integer(is)                     :: rank

    ! initialize to zero
    data_set = 0.
    rank     = 2

    !data_type = H5T_NATIVE_DOUBLE
    call open_dataset(file_name, data_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    ! make sure we can hold the transposed data set
    ! if we're wrong about size of data set, stop with error
    if(any(h5dims(2:1:-1) > shape(data_set))) stop 'read_dataset_t_float - array too small to hold data set'//data_name
    ! dim_read is the dimensions of the final transposed array
    dim_read = h5dims(2:1:-1)

    ! for now make buf_sze 1/10 of full integral tensor. Clearly
    ! this can be optimized if need be. 
    buf_sze    = int(0.1 * dim_read(2))
    allocate(int_buf(buf_sze, dim_read(1)))

    n_remain   = h5dims(1)
    n_read_tot = 0
    offmem     = (/0, 0/)
    data_type  = h5kind_to_type(hp,H5_REAL_KIND)
    do while(n_remain > 0)

      n_read    = min(n_remain, buf_sze)
      offset    = (/ n_read_tot, 0 /)
      buf_read  = (/ n_read, dim_read(1) /)

      !
      ! Create memory dataspace.
      !
      call h5screate_simple_f(2, buf_read, buf_space, error)

      !
      ! Select hyperslab in the data set 
      !
      CALL h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, buf_read, error)

      !
      ! Select hyperslab in memory.
      !
      CALL h5sselect_hyperslab_f(buf_space, H5S_SELECT_SET_F, offmem, buf_read, error)

      !
      ! Read the dataset.
      !
      call h5dread_f(dset_id, data_type, int_buf, buf_read, error,  &
                 mem_space_id=buf_space, file_space_id=space_id)

      data_set(:dim_read(1), n_read_tot+1:n_read_tot+n_read) = transpose(int_buf(:n_read, :dim_read(1)))

      ! 
      ! close the memory space
      !
      call h5sclose_f(buf_space, error)

      ! update n_remain
      n_remain   = n_remain - n_read
      n_read_tot = n_read_tot + n_read

    enddo

    call close_dataset(file_id, dset_id, space_id)
    deallocate(int_buf)

  end subroutine read_dataset_t_float16

  !
  !
  !
  subroutine read_dataset_t_float(file_name, data_name, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    real(sp), intent(inout)         :: data_set(:, :)
    integer(is), intent(out)        :: dim_read(2)

    integer(is)                     :: buf_sze
    integer(is)                     :: n_remain
    integer(is)                     :: n_read
    integer(is)                     :: n_read_tot
    integer(hsize_t)                :: offset(2)
    integer(hsize_t)                :: offmem(2)
    integer(hsize_t)                :: buf_read(2)
    real(sp),allocatable            :: int_buf(:,:)

    integer(hid_t)                  :: file_id, space_id, dset_id, buf_space
    integer(hid_t)                  :: h5dims(2)
    integer(hid_t)                  :: data_type
    integer(is)                     :: error
    integer(is)                     :: rank

    ! initialize to zero
    data_set = 0.
    rank     = 2

    !data_type = H5T_NATIVE_DOUBLE
    call open_dataset(file_name, data_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    ! make sure we can hold the transposed data set
    ! if we're wrong about size of data set, stop with error
    if(any(h5dims(2:1:-1) > shape(data_set))) stop 'read_dataset_t_float - array too small to hold data set'//data_name
    ! dim_read is the dimensions of the final transposed array
    dim_read = h5dims(2:1:-1)

    ! for now make buf_sze 1/10 of full integral tensor. Clearly
    ! this can be optimized if need be. 
    buf_sze    = int(0.1 * dim_read(2))
    allocate(int_buf(buf_sze, dim_read(1)))

    n_remain   = h5dims(1)
    n_read_tot = 0
    offmem     = (/0, 0/)
    data_type  = h5kind_to_type(sp,H5_REAL_KIND)
    do while(n_remain > 0)

      n_read    = min(n_remain, buf_sze)
      offset    = (/ n_read_tot, 0 /)
      buf_read  = (/ n_read, dim_read(1) /)

      !
      ! Create memory dataspace.
      !
      call h5screate_simple_f(2, buf_read, buf_space, error)

      !
      ! Select hyperslab in the data set 
      !
      CALL h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, buf_read, error)

      !
      ! Select hyperslab in memory.
      !
      CALL h5sselect_hyperslab_f(buf_space, H5S_SELECT_SET_F, offmem, buf_read, error)

      !
      ! Read the dataset.
      !
      call h5dread_f(dset_id, data_type, int_buf, buf_read, error,  &
                 mem_space_id=buf_space, file_space_id=space_id)

      data_set(:dim_read(1), n_read_tot+1:n_read_tot+n_read) = transpose(int_buf(:n_read, :dim_read(1)))

      ! 
      ! close the memory space
      !
      call h5sclose_f(buf_space, error)

      ! update n_remain
      n_remain   = n_remain - n_read
      n_read_tot = n_read_tot + n_read

    enddo

    call close_dataset(file_id, dset_id, space_id)
    deallocate(int_buf)

  end subroutine read_dataset_t_float

  !
  !
  !
  subroutine read_dataset_t_dble(file_name, data_name, data_set, dim_read)
    character(len=255), intent(in)  :: file_name
    character(len=255), intent(in)  :: data_name
    real(dp), intent(inout)         :: data_set(:, :)
    integer(is), intent(out)        :: dim_read(2)

    integer(is)                     :: buf_sze
    integer(is)                     :: n_remain
    integer(is)                     :: n_read
    integer(is)                     :: n_read_tot
    integer(hsize_t)                :: offset(2)
    integer(hsize_t)                :: offmem(2)
    integer(hsize_t)                :: buf_read(2)
    real(dp),allocatable            :: int_buf(:,:)

    integer(hid_t)                  :: file_id, space_id, dset_id, buf_space
    integer(hid_t)                  :: h5dims(2)
    integer(hid_t)                  :: data_type
    integer(is)                     :: error
    integer(is)                     :: rank

    ! initialize to zero
    data_set = 0.
    rank     = 2

    !data_type = H5T_NATIVE_DOUBLE
    call open_dataset(file_name, data_name, data_type, file_id, dset_id, space_id, rank, h5dims, error)

    ! make sure we can hold the transposed data set
    ! if we're wrong about size of data set, stop with error
    if(any(h5dims(2:1:-1) > shape(data_set))) stop 'read_dataset_dble_t - array too small to hold data set'//data_name
    ! dim_read is the dimensions of the final transposed array
    dim_read = h5dims(2:1:-1)

    ! for now make buf_sze 1/10 of full integral tensor. Clearly
    ! this can be optimized if need be. 
    buf_sze    = int(0.1 * dim_read(2))
    allocate(int_buf(buf_sze, dim_read(1)))

    n_remain   = h5dims(1)
    n_read_tot = 0
    offmem     = (/0, 0/)
    data_type = h5kind_to_type(dp,H5_REAL_KIND)
    do while(n_remain > 0)

      n_read    = min(n_remain, buf_sze)
      offset    = (/ n_read_tot, 0 /)
      buf_read  = (/ n_read, dim_read(1) /)

      !
      ! Create memory dataspace.
      !
      call h5screate_simple_f(2, buf_read, buf_space, error)

      !
      ! Select hyperslab in the data set 
      !
      CALL h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, buf_read, error)

      !
      ! Select hyperslab in memory.
      !
      CALL h5sselect_hyperslab_f(buf_space, H5S_SELECT_SET_F, offmem, buf_read, error)

      !
      ! Read the dataset.
      !
      call h5dread_f(dset_id, data_type, int_buf, buf_read, error,  &
                 mem_space_id=buf_space, file_space_id=space_id)
 
      data_set(:dim_read(1), n_read_tot+1:n_read_tot+n_read) = transpose(int_buf(:n_read, :dim_read(1)))        

      ! 
      ! close the memory space
      !
      call h5sclose_f(buf_space, error)

      ! update n_remain
      n_remain   = n_remain - n_read
      n_read_tot = n_read_tot + n_read 

    enddo

    call close_dataset(file_id, dset_id, space_id)
    deallocate(int_buf)

  end subroutine read_dataset_t_dble

end module h5_ops 
