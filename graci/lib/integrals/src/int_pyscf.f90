!------------------------------------------------------------------
! int_pyscf
!
! module for accessing integrals written to h5 files, including
! either canonical 4-index electron repulsion integrals, or 
! 3-index quantities obtained from density-fitting procedure
!
!
!

module int_pyscf
  use hdf5
  use accuracy
  use map_module
  implicit none

  public intpyscf_initialise 

  interface mo_integral
    module procedure mo_integral_ijkl
    module procedure mo_integral_jkl
    module procedure mo_integral_kl
  end interface mo_integral

  ! number of molecular orbitals
  integer(ik)                                 :: n_mo
  ! lower triangle of n_mo x n_mo
  integer(ik)                                 :: n_ij
  ! number of auxiliary basis functions if density fitting employed.
  ! If no density fitting, this value is equal to n_mo
  integer(ik)                                 :: n_aux
  ! the symmetries of each MO
  integer(ik), allocatable                    :: mo_sym(:)

  ! core hamiltonian
  real(drk), allocatable                      :: h_core(:, :)

  ! temporary integral buffer -- 3 index quantities
  real(drk), allocatable                      :: integrals(:,:) 
  ! temporary  -- this should not be a global variable
  logical                                     :: ri_int
  logical                                     :: use_hash

  type(map_type)                              :: int_table
  
  contains

  !--------------------------------------------------------------------------------------
  !  Module Public
  !

  !
  ! Routine to load MO integrals from HDF5 file format.
  ! Arguments:
  !     f1e:     Name of one-electron integral file (character string)
  !     f2eri:   name of two-electron ERI integral file (character string)
  !     n_bra:   If full 4 index integral transformation is 
  !              performed, this will be the lower triangle 
  !              of the ij indicies, i.e. n_mo * (n_mo+1)/2
  !              if density fitting is employed, this will 
  !              be the number of auxility basis functions
  !     n_ket:   This is the lower triangle of the kl indices
  !              i.e. n_mo*(n_mo+1)/2
  !     use_ri:  Boolean that is true if DF is employed, else
  !              False
  !     use_rr:  Boolean that is true if the rank reduction of the
  !              DF 3-index integral tensor is to be performed
  !     rr_fac:  Rank reduced DF compression factor  
  !
#ifdef CBINDING
    subroutine intpyscf_initialise(f1e, f2eri, n_mo1, n_aux1, use_ri, use_rr, rr_fac, thresh, max_memory) &
         bind(c,name="intpyscf_initialise")
#else
      subroutine intpyscf_initialise(f1e, f2eri, n_mo, n_aux, use_ri, use_rr, rr_fac, thresh, max_memory)
#endif
    use iomod
    use rrdf
    use iso_c_binding, only: C_CHAR

    integer(ik), intent(in)                   :: n_mo1           ! number of MOs
    integer(ik), intent(in)                   :: n_aux1          ! number of auxility basis functions, if using DF, else
                                                                 ! it's the dimension of the AO space
    logical, intent(in)                       :: use_ri          ! whether or not to use density fitting
    logical, intent(in)                       :: use_rr          ! whether or not to use rank density fitting
    integer(ik), intent(in)                   :: rr_fac          ! rank reduced DF compression factor
    real(drk), intent(in)                     :: thresh          ! threshold for zero integrals
    integer(ik), intent(in), optional         :: max_memory      ! maximum memory that can be used by integral arrays     
                                                                 ! (in double precision numbers). If < 0, no limit is enforced
#ifdef CBINDING
    character(kind=C_CHAR), intent(in)       :: f1e(*)          ! file name, one electron integrals
    character(kind=C_CHAR), intent(in)       :: f2eri(*)        ! flie name, two-e ERI integrals
    character(len=255)                       :: f1name
    character(len=255)                       :: f2name
    integer(ik)                              :: length
#else
    character(len=*), intent(in)             :: f1e
    character(len=*), intent(in)             :: f2eri 
    character(len=255)                       :: f1name
    character(len=255)                       :: f2name
#endif

    ! H5 file variables
    integer(HID_T)                            :: data_type       ! data_type (i.e. NATIVE_DOUBLE)
    integer(HID_T)                            :: file_id         ! file_id for integral file
    integer(HID_T)                            :: dset_id         ! geeric dataset id
    integer(HID_T)                            :: data_space        ! data_space id
    integer(HID_T)                            :: mem_space       ! memspace identifier 
    integer(HSIZE_T)                          :: int_dims(2)     ! dimensions of the integral dataset (expected)
    integer(HSIZE_T)                          :: dims(2)         ! dimensions of dataset upon query
    integer(HSIZE_T)                          :: offset(2)       ! offset for the hyperslab 
    integer(HSIZE_T)                          :: offset_out(2)   ! offset for the array buffer
    integer(HSIZE_T)                          :: buf_read(2)     ! dimensions of the read buffer
    integer(ik)                               :: error
    integer(ik)                               :: rank

    ! buffers and counting variables
    character(len=8)                          :: dset_name
    integer(ik)                               :: n_read
    integer(ik)                               :: n_read_tot
    integer(ik)                               :: n_remain
    integer(ik)                               :: buf_sze
    integer(ik)                               :: bra, ket
    integer(ik)                               :: ket_off
    integer(ik)                               :: nd_ket, nd_bra
    real(drk), allocatable                    :: int_buffer(:,:)

    !----------------------------------------------------------------------
    ! If C bindings are on, then convert the calculation label from the
    ! C char type to the Fortran character type
    !----------------------------------------------------------------------
#ifdef CBINDING
    length=cstrlen(f1e)
    call c2fstr(f1e, f1name, length)
    length=cstrlen(f2eri)
    call c2fstr(f2eri, f2name, length)
#else
    f1name = adjustl(trim(f1e))
    f2name = adjustl(trim(f2eri))
#endif

    ! set module variables
    ri_int = use_ri
    n_aux  = n_aux1
    n_mo   = n_mo1

    !*********************************************************
    ! one-electron integral file(s)
    !
    int_dims = (/ n_mo, n_mo /)
    allocate(h_core(n_mo, n_mo))

    dset_name = 'hcore_mo'
    call open_dataset(f1name, dset_name, data_type, file_id, dset_id,&
         data_space, rank, dims)

    if(rank /= 2 .or. any(dims /= int_dims))then
      print *,'Error in hcore_mo file: read status=',error,' rank=',rank, &
              ' dims found=',dims, ' dims expected=',int_dims
      stop 'Error reading core Hamiltonian integral file'
    endif
 
    !
    ! Read the dataset.
    !
    call h5dread_f(dset_id, data_type, h_core, int_dims, error)
   
    ! 
    ! close the dataset 
    !
    call close_dataset(file_id, dset_id)
 
    !********************************************************
    ! determine dimensions of integral files
    !
    n_ij   = n_mo * (n_mo + 1) / 2

    nd_ket = n_ij
    if(use_ri) then
      nd_bra = n_aux
    else
      nd_bra = n_ij
    endif

    !*********************************************************
    ! set up integral storage scheme
    !
    ! allocate the buffer to hold all 2 ERI, or, 3-indx DF quantities
    call set_int_buffer_sizes(n_mo, nd_bra, nd_ket, use_ri, max_memory, buf_sze)

    !********************************************************
    ! if using hash storage, using initialize hash table
    !
    if(use_hash) then
      call map_init(int_table, mo_eri_integral_index(n_mo, n_mo, n_mo, n_mo))
    else
      allocate(integrals(nd_bra, nd_ket))
    endif

    !********************************************************
    ! two-electron electron repulsion integral file(s)
    !

    ! allocate the buffer to hold all 2 ERI, or, 3-indx DF quantities
    int_dims = (/ nd_ket, nd_bra /)
    allocate(int_buffer(buf_sze, nd_bra))

    !
    ! open the dataset
    !
    dset_name = 'eri_mo'
    call open_dataset(f2name, dset_name, data_type, file_id, dset_id, data_space, rank, dims)
    
    if(rank /= 2 .or. any(dims /= int_dims))then
      print *,'Error in 2ERI file: read status=',error,' rank=',rank, &
              ' dims found=',dims, ' dims expected=',int_dims
      stop 'Error reading 2e-ERI integral file'
    endif

    n_remain   = n_ij
    n_read_tot = 0
    offset_out = (/0, 0/)
    do while(n_remain > 0)

      n_read    = min(n_remain, buf_sze)
      offset    = (/ n_read_tot, 0 /)
      buf_read  = (/ n_read, nd_bra /)

      !
      ! Create memory dataspace.
      !
      call h5screate_simple_f(2, buf_read, mem_space, error)

      !
      ! Select hyperslab in the data set 
      !
      CALL h5sselect_hyperslab_f(data_space, H5S_SELECT_SET_F, offset, buf_read, error)

      !
      ! Select hyperslab in memory.
      !
      CALL h5sselect_hyperslab_f(mem_space, H5S_SELECT_SET_F, offset_out, buf_read, error)

      !
      ! Read the dataset.
      !
      call h5dread_f(dset_id, data_type, int_buffer, buf_read, error,  &
                     mem_space_id=mem_space, file_space_id=data_space)
    
      ket_off = offset(1) 
      if(use_hash) then
        call load_hash_batch(nd_bra, n_read, ket_off, int_buffer, real(thresh,kind=integral_kind))
      else
        integrals(:nd_bra, ket_off+1:n_read) = transpose(int_buffer(:n_read, :nd_bra))        
      endif

      ! 
      ! close the memory space
      !
      call h5sclose_f(mem_space, error)

      ! update n_remain
      n_remain   = n_remain - n_read
      n_read_tot = n_read_tot + n_read 

    enddo

    ! 
    ! deallocate the integral buffer
    !
    deallocate(int_buffer)

    !
    ! close the dataset
    !
    call close_dataset(file_id, dset_id)  

    !
    ! Close FORTRAN interface.
    !
    call h5close_f(error)

    !
    ! Rank reduction of the density fitting 3-index integrals
    !    
    if (use_ri .and. use_rr) call rank_reduce(n_aux, n_ij, integrals, rr_fac)
    
    return
  end subroutine intpyscf_initialise 

  !
  ! deallocates all the memory structures associated with integral
  ! integral storage
  !
  subroutine intpyscf_finalise()
    implicit none

    if(allocated(mo_sym))deallocate(mo_sym)
    deallocate(h_core)
    deallocate(integrals)

    return
  end subroutine intpyscf_finalise


  !-------------------------------------------------------------------------------------
  ! Module private
  !

  !
  ! load a batch of integrals into the hash table
  !
  subroutine load_hash_batch(n_bra, n_ket, ket_offset, int_buffer, thresh)
    implicit none
    integer(ik), intent(in)           :: n_bra                    ! number of bra values
    integer(ik), intent(in)           :: n_ket                    ! number of ket values
    integer(ik), intent(in)           :: ket_offset               ! ket values read in batches, this is offset
                                                                  ! of current batch
    real(rk), intent(in)              :: int_buffer(n_ket, n_bra) ! integrals read from buffer transposed
    real(integral_kind), intent(in)   :: thresh

    integer(key_kind),allocatable     :: keys(:)      
    real(integral_kind),allocatable   :: integrals(:)
    integer(ik)                       :: i,j,k,l
    integer(ik)                       :: bra, ket
    integer(ik)                       :: n_append, n_total

    ! read in terms of columns of the integral buffer
    allocate(keys(n_ket))
    allocate(integrals(n_ket))

    n_total = 0
    do bra = 1,n_bra
      call indx_pair(bra, i, j)

      n_append = 0
      do ket = 1,n_ket
        if(abs(int_buffer(ket, bra)) >= thresh) then
          call indx_pair(ket+ket_offset, k, l)
          n_append            = n_append + 1
          keys(n_append)      = mo_eri_integral_index(i, j, k, l)
          integrals(n_append) = int_buffer(ket, bra)
        endif
      enddo

      call map_update(int_table, keys(:n_append), integrals(:n_append), n_append, thresh, .false.)
      n_total = n_total + n_append

    enddo

    call map_merge(int_table)

    write(6,1000)n_total, n_ket*n_bra

    deallocate(keys)
    deallocate(integrals)

    return
1000 format('   --> Appending ',i10,' out of ',i10,' integrals in batch to hash table')
  end subroutine load_hash_batch

  ! 
  ! 
  ! 
  subroutine open_dataset(file_name, dset_name, data_type, file_id, dset_id, data_space, rank, dims)
    implicit none
    character(len=255), intent(in)      :: file_name   ! name of integral file
    character(len=8),intent(in)         :: dset_name   ! name of data set (i.e. h_core, eri) 
    integer(HID_T),intent(out)          :: data_type   ! data_type (i.e. NATIVE_DOUBLE)
    integer(HID_T),intent(out)          :: file_id     ! file identifier
    integer(HID_T),intent(out)          :: dset_id     ! geeric dataset id
    integer(HID_T),intent(out)          :: data_space  ! data_space id
    integer(ik), intent(out)            :: rank        ! rank of dataset (should be 2)
    integer(HSIZE_T),intent(out)        :: dims(:)     ! dimensions of dataset upon query

    logical                             :: exists
    integer(ik)                         :: error

    integer(HSIZE_T)                    :: max_dims(4)


    ! make sure integral files exist
    inquire(file=file_name, exist=exists)
    if(.not.exists) then
      print *, 'Cannot file file: '//file_name 
      stop 'Cannot find integral file.'
    endif

    ! 
    ! initialize h5 system
    !
    call h5open_f(error)

    !
    ! read in the one electron integrals
    !
    call h5fopen_f (file_name, H5F_ACC_RDONLY_F, file_id, error)

    !
    ! open the integral data set, will only accept the name
    !
    call h5dopen_f(file_id, dset_name, dset_id, error)

    !
    ! Get dataset's data type.
    !
    call h5dget_type_f(dset_id, data_type, error)
    call h5dget_space_f(dset_id, data_space, error)

    ! 
    ! Get the rank of the dataset and sanity check against 
    ! what is expected 
    !
    call h5sget_simple_extent_ndims_f(data_space, rank, error)
    call h5sget_simple_extent_dims_f(data_space, dims, max_dims, error)

    return
  end subroutine open_dataset

  ! 
  !
  ! 
  subroutine close_dataset(file_id, dset_id)
    implicit none
    integer(HID_T),intent(in)          :: dset_id    ! data_space id
    integer(HID_T),intent(in)          :: file_id    ! file_id for integral file

    integer(ik)                        :: error

    !
    ! Close the dataset
    !
    call h5dclose_f(dset_id, error)

    !
    ! Close the integral file
    !
    call h5fclose_f(file_id, error)

    return
  end subroutine close_dataset

  !
  !
  !
  subroutine set_int_buffer_sizes(n_mo, nd_bra, nd_ket, use_ri, max_memory, buf_sze)
    implicit none
    integer(ik), intent(in)                   :: n_mo
    integer(ik), intent(in)                   :: nd_bra
    integer(ik), intent(in)                   :: nd_ket
    logical, intent(in)                       :: use_ri
    integer(ik), intent(in)                   :: max_memory
    integer(ik), intent(out)                  :: buf_sze

    use_hash = .false.

    if(max_memory < 0) then
      buf_sze      = nd_ket ! can hold entire integral file in memory
      write(6,1000)
    else

      ! remove space required to store one and two electron integrals 
      buf_sze = max_memory - nd_bra * nd_ket - n_mo * n_mo
      buf_sze = int(buf_sze / nd_bra)+1

      ! if buf_sze is greater than nd_ket, we can slurp in all the 
      ! integrals at once
      if(buf_sze > 0) then
        buf_sze = min(nd_ket, buf_sze)
        write(6,1001)ceiling(1.*nd_ket/buf_sze)
      ! if negative, we don't have memory to store all integrals
      ! in memory and will have to a hash. However, we don't know
      ! how big hash table will be. A sensible buffer size assumes
      ! we can use at most 10% of the maximum memory limit
      elseif(buf_sze < 0) then

        if(use_ri) stop 'We do not currently support hashing of DF computed integrals'
        buf_sze = min(nd_ket, int(0.1*max_memory))
        use_hash = .true.
        write(6,1001)ceiling(1.*nd_ket/buf_sze)
      endif

    endif

    return
1000 format(/,' -- Storing all integrals in memory, reading in single batch --')
1001 format(/,' -- Storing integrals in hash table, buffered read will use ',i4,' batches --')
  end subroutine set_int_buffer_sizes

  !
  ! return the core Hamiltonian
  !
  function h_1e() result(h)
    real(drk)                                 :: h(n_mo, n_mo)

    h = h_core
  end function h_1e

  !
  ! return an element of the core Hamiltonian
  !
  function h_1e_ij(i, j) result(hij)
    integer(ik), intent(in)                   :: i, j
    real(drk)                                 :: hij

    hij = h_core(i,j)

  end function h_1e_ij
  
  !
  ! return the single integral, (ij|kl)
  !
  function mo_integral_ijkl(i, j, k, l) result(ijkl)
    integer(ik), intent(in)                   :: i
    integer(ik), intent(in)                   :: j
    integer(ik), intent(in)                   :: k
    integer(ik), intent(in)                   :: l

    integer(ik)                               :: ij, kl
    real(drk)                                 :: ijkl
    real(integral_kind)                       :: ijkl_tmp

    if (use_hash) then
       call map_get(int_table, mo_eri_integral_index(i,j,k,l), ijkl_tmp)
       ijkl = real(ijkl_tmp, kind=drk)
       return
    else
       ij = indx_ut(i,j)
       kl = indx_ut(k,l)
       if(ri_int) then
          ijkl = dot_product(integrals(:,ij),integrals(:,kl))
       else
          ijkl = integrals(ij, kl)
       endif
       return
    endif
    
  end function mo_integral_ijkl

  !
  ! with "jkl" fixed, return integrals running over "i"
  !
  function mo_integral_jkl(j, k, l) result(jkl)
    integer(ik), intent(in)                   :: j
    integer(ik), intent(in)                   :: k
    integer(ik), intent(in)                   :: l

    real(drk)                                 :: jkl(n_mo)
    integer(ik)                               :: i_mo

    do i_mo = 1,n_mo
      jkl(i_mo) = mo_integral_ijkl(i_mo, j, k, l)
    enddo

    return
  end function mo_integral_jkl

  !
  ! with "kl" fixed, return integrals running over i
  !
  function mo_integral_kl(k, l) result(kl)
    integer(ik), intent(in)                   :: k
    integer(ik), intent(in)                   :: l

    integer(key_kind)                         :: key_array(n_mo*n_mo)
    real(integral_kind)                       :: int_array(n_mo*n_mo)
    real(drk)                                 :: kl(n_mo, n_mo)
    integer(ik)                               :: i, j, sze

    if(use_hash) then

      sze = 0
      do i = 1,n_mo
        do j = 1,n_mo
            sze = sze + 1
            key_array(sze) = mo_eri_integral_index(i,j,k,l)
        enddo
      enddo

      call map_get_many(int_table, key_array, int_array, sze)

      sze = 0
      do i = 1,n_mo
        do j = 1,n_mo
          sze     = sze + 1
          kl(i,j) = real(int_array(sze), kind=drk)
        enddo
      enddo

    else

      do i = 1,n_mo
        do j = 1,n_mo
          kl(i,j) = mo_integral_ijkl(i, j, k, l)
        enddo
      enddo

    endif

    return
  end function mo_integral_kl

  !
  ! unpacks indices i,j from upper triangle index 'ij'
  !
  subroutine indx_pair(ij, i, j)
    integer(ik), intent(in)          :: ij
    integer(ik), intent(out)         :: i
    integer(ik), intent(out)         :: j

    ! solve quadratic equation for i: let j = 1, solve for i,
    ! then determine j
    i = int(0.5*(1 + sqrt(1. + 8.*(ij-1.))))
    !i = max(i,1)
    j = ij - i*(i-1)/2

    if(i*(i-1)/2+j /= ij) stop 'indx_pair error'


    return
  end subroutine indx_pair

  !
  ! returns upper-triangle index
  !
  function indx_ut(i, j) result(ij)
    integer(ik), intent(in)              :: i
    integer(ik), intent(in)              :: j
    integer(ik)                          :: ij

    ij = int(max(i,j)*(max(i,j) - 1)/2 + min(i,j))

    return
  end function indx_ut

  !
  ! Computes an unique index for i,j,k,l integrals
  !
  function mo_eri_integral_index(i, k, j, l) result(key)
    integer, intent(in)              :: i,j,k,l
    integer(key_kind)                    :: key
    integer(key_kind)                    :: p,q,r,s,i2
    p   = min(i,k)
    r   = max(i,k)
    p   = p + shiftr(r*r-r,1)
    q   = min(j,l)
    s   = max(j,l)
    q   = q + shiftr(s*s-s,1)
    key = min(p,q)
    i2  = max(p,q)
    key = key + shiftr(i2*i2-i2,1)
  end function mo_eri_integral_index


end module int_pyscf

