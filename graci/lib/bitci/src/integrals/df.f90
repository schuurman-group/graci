module dfmod

  use h5_ops
  use constants
  use integrals
  
  implicit none

  !
  ! Exact integrals type
  !
  type, extends(eri) :: df 
   integer(is)                :: n_aux
   real(dp), allocatable      :: bra_ket(:,:)   

   contains
     procedure   :: init_pyscf => init_pyscf_df
     procedure   :: mo_ints    => mo_ints_df
     procedure   :: mo_int     => mo_int_df
     procedure   :: finalize   => finalize_df 
 end type df 

contains

  subroutine init_pyscf_df(ints, core_file, eri_file)

    class(df)               :: ints
    character(len=255)      :: core_file
    character(len=255)      :: eri_file

    character(len=255)      :: f_name
    character(len=255)      :: dset_name
    logical                 :: exists
    integer(is)             :: n_ij
    integer(is)             :: rank
    integer(is)             :: dims(2)
    integer(is)             :: dims_read(2)
    integer(is)             :: n_bra_ket

    ! load the one hamiltonian
    !------------------------------------------------------------
    f_name    = trim(adjustl(core_file))
    dset_name = 'hcore_mo'
    exists    = dataset_exists(f_name, dset_name)

    if(.not.exists) stop 'cannot find hcore_mo in file='//f_name

    call dataset_dims(f_name, dset_name, rank, dims)
    ints%nmo = dims(1) 

    if(allocated(ints%h_core))deallocate(ints%h_core)
    allocate(ints%h_core(ints%nmo, ints%nmo))

    call read_dataset_dble(f_name, dset_name, ints%h_core, dims_read)

    ! load ERI
    !--------------------------------------------------------------
    f_name    = trim(adjustl(eri_file))
    dset_name = 'eri_mo'
    exists    = dataset_exists(f_name, dset_name)

    if(.not.exists) stop 'cannot find eri_mo in file='//f_name

    call dataset_dims(f_name, dset_name, rank, dims)
    n_ij       = ints%nmo * (ints%nmo + 1)/2
    ints%n_aux = dims(2)    
    if(dims(1) /= n_ij) stop 'ERI data set wrong size, file='//f_name

    ! dataset written as (n_ij, n_aux) -- we want the transpose
    dims = (/ints%n_aux, n_ij/)

    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)
    allocate(ints%bra_ket(ints%n_aux, n_ij))

    ! this method performs buffered read of dataset and transposes on the fly
    ! since the DF tensor from PySCF is stored (n_ij, n_aux)
    call read_dataset_dble_t(f_name, dset_name, ints%bra_ket, dims_read)

    return

  end subroutine init_pyscf_df

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_df(ints, indices, int_vec)

    class(df)              :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp), intent(out)  :: int_vec(:)

    integer(is)            :: i, ij, kl 
    integer(is)            ::  nints

    if(size(indices, dim=1) /= 4) stop 'mo_ints: indices dim=1 must equal 4'
    nints = size(int_vec)

    do i = 1,nints
        ij = ints%indx_ut(indices(1,i), indices(2,i))
        kl = ints%indx_ut(indices(3,i), indices(4,i))
        int_vec(i) = dot_product(ints%bra_ket(:,ij), ints%bra_ket(:,kl))
    enddo

  end subroutine mo_ints_df

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_df(ints, i, j, k, l) result(int_val)

    class(df)              :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = dot_product(ints%bra_ket(:,ints%indx_ut(i,j)), ints%bra_ket(:,ints%indx_ut(k,l)))

    return
  end function mo_int_df

  !
  !
  !
  subroutine finalize_df(ints)

    class(df)              :: ints  

    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)

  end subroutine finalize_df

end module dfmod
