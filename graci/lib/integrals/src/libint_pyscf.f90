!------------------------------------------------------------------
! libint_pyscf
!
! module for accessing integrals written to h5 files, including
! either canonical 4-index electron repulsion integrals, or 
! 3-index quantities obtained from density-fitting procedure
!
!
!

module libint_pyscf
  implicit none
  use hdf5
  use accuracy
  use map_module

  module procedure load_mo_integrals

  interface mo_integral
    public mo_integral_ijkl
    public mo_integral_jkl
    public mo_integral_kl
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

  type(map_type)                              :: int_table
  
  contains

  !
  ! Routine to load MO integrals from HDF5 file format.
  ! Arguments:
  !     h5_name: Name of integral file (character string)
  !     n_bra:   If full 4 index integral transformation is 
  !              performed, this will be the lower triangle 
  !              of the ij indicies, i.e. n_mo * (n_mo+1)/2
  !              if density fitting is employed, this will 
  !              be the number of auxility basis functions
  !     n_ket:   This is the lower triangle of the kl indices
  !              i.e. n_mo*(n_mo+1)/2
  !     use_ri:  Boolean that is true if DF is employed, else
  !              False
  !     
  subroutine load_mo_integrals(h5_name, n_bra, n_ket, use_ri)
    character(len=100), intent(in)            :: h5_name
    integer(ik), intent(in)                   :: n_bra
    integer(ik), intent(in)                   :: n_ket
    logical, intent(in)                       :: use_ri

    logical                                   :: exists
    integer(HID_T)                            :: int_file_id
    integer(HID_T)                            :: dset_id
    integer(HID_T)                            :: dtype_id
    integer(ik)                               :: error
    integer(ik)                               :: nd 
    integer(ik)                               :: int_dims(2)

    real(drk), allocatable                    :: int_buffer(:, :)
    character(len=7)                          :: dset_name

    n_mo     = n_ket
    n_aux    = n_bra
    n_ij     = n_mo * (n_mo + 1) / 2

    nd_bra   = n_bra * (n_bra + 1) / 2
    nd_ket   = n_ket * (n_ket + 1) / 2
    int_dims = (/nd_bra, nd_ket/)
 
    ! allocate the buffer to hold all 2 ERI, or, 3-indx DF quantities
    allocate(int_bufer(nd_bra, nd_ket))

    ! name the data set to retrieve
    if(use_ri) then
      dset_name = adjustl('/j3c')
    else
      dset_name = '/eri_mo'
    endif

    !
    ! make sure the integral file exists
    !
    inquire(file=h5_name, exist=exists)

    if(.not.exists) stop 'cannot find file: '//trim(h5_name)
    
    !
    ! Open hdf5 integral file
    !
    call h5fopen_f (h5_name, H5F_ACC_RDONLY_F, int_file_id, error)

    !
    ! open the integral data set, will only accept the name "ERI"
    !
    call h5dopen_f(int_file_id, dset_name, dset_id, error)

    !
    ! Get dataset's data type.
    !
    call h5dget_type_f(dset_id, dtype_id, error)

    !
    ! Read the dataset.
    !
    call h5dread_f(dset_id, dtype_id, int_buffer, int_dims, error)

    !
    ! Close the dataset
    !
    !
    call h5dclose_f(dset_id, error)
    call h5tclose_f(dtype_id, error)

    !
    ! Close the integral file
    !
    call h5fclose_f(int_file_id, error)

    !
    ! Close FORTRAN interface.
    !
    call h5close_f(error)

    !
    ! Load the integral buffer into the hash
    !
    

    deallocate(int_buffer)
    return
  end subroutine load_mo_integrals

  !
  ! return the single integral, (ij|kl)
  !
  function mo_integral_ijkl(i, j, k, l) result(ijkl)
    integer(ik), intent(in)                   :: i
    integer(ik), intent(in)                   :: j
    integer(ik), intent(in)                   :: k
    integer(ik), intent(in)                   :: l
 
    real(drk)                                 :: ijkl


    return
  end function mo_integral_ijkl

  !
  ! with "jkl" fixed, return integrals running over "i"
  !
  function mo_integral_jkl(j, k, l) result(jkl)
    integer(ik), intent(in)                   :: j
    integer(ik), intent(in)                   :: k
    integer(ik), intent(in)                   :: l

    real(drk)                                 :: jkl(n_mo)


    return
  end function mo_integral_jkl

  !
  ! with "kl" fixed, return integrals running over i
  !
  function mo_integral_kl(k, l) result(kl)
    integer(ik), intent(in)                   :: j
    integer(ik), intent(in)                   :: k

    real(drk)                                 :: kl(n_ij)

    return
  end function mo_integral_kl

  !
  ! Computes an unique index for i,j,k,l integrals
  !
  function mo_eri_integral_index(i, j, k, l) result(key)
    integer(ik), intent(in)              :: i,j,k,l
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
  end subroutine mo_eri_integral_index


end module libint_pyscf

