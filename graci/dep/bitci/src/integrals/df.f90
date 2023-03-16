module dfmod

  use h5_ops
  use constants
  use integrals
  
  implicit none

  !
  ! density fitting integrals type, double precision
  !
  type, extends(eri) :: df_dp 
   integer(is)                :: n_aux
   real(dp), allocatable      :: bra_ket(:,:)   

   contains
     procedure   :: init_pyscf => init_pyscf_df_dp
     procedure   :: mo_ints    => mo_ints_df_dp
     procedure   :: mo_int     => mo_int_df_dp
     procedure   :: finalize   => finalize_df_dp
  end type df_dp 

  !
  ! density fitting integrals type, single precision
  !
  type, extends(eri) :: df_sp
   integer(is)                :: n_aux
   real(sp), allocatable      :: bra_ket(:,:)

   contains
     procedure   :: init_pyscf => init_pyscf_df_sp
     procedure   :: mo_ints    => mo_ints_df_sp
     procedure   :: mo_int     => mo_int_df_sp
     procedure   :: finalize   => finalize_df_sp
  end type df_sp

  !
  ! density fitting integrals type, half precision
  !
  type, extends(eri) :: df_hp
   integer(is)                :: n_aux
   real(hp), allocatable      :: bra_ket(:,:)

   contains
     procedure   :: init_pyscf => init_pyscf_df_hp
     procedure   :: mo_ints    => mo_ints_df_hp
     procedure   :: mo_int     => mo_int_df_hp
     procedure   :: finalize   => finalize_df_hp
  end type df_hp

contains

 !-------------------------------------------------------------
 ! Double precision routines
 !

  !
  !
  !
  subroutine init_pyscf_df_dp(ints, core_file, eri_file)

    class(df_dp)            :: ints
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

    call read_dataset(f_name, dset_name, ints%h_core, dims_read)

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
    call read_dataset_t(f_name, dset_name, ints%bra_ket, dims_read)

    return

  end subroutine init_pyscf_df_dp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_df_dp(ints, indices, int_vec)

    class(df_dp)           :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp),intent(out)   :: int_vec(:)

    integer(is)            :: i, ij, kl 
    integer(is)            :: nints

    if(size(indices, dim=1) /= 4) stop 'mo_ints: indices dim=1 must equal 4'
    nints = size(int_vec)

    do i = 1,nints
        ij = ints%indx_ut(indices(1,i), indices(2,i))
        kl = ints%indx_ut(indices(3,i), indices(4,i))
        int_vec(i) = dot_product(ints%bra_ket(:,ij), ints%bra_ket(:,kl))
    enddo

  end subroutine mo_ints_df_dp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_df_dp(ints, i, j, k, l) result(int_val)

    class(df_dp)           :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = dot_product(ints%bra_ket(:,ints%indx_ut(i,j)), & 
                          ints%bra_ket(:,ints%indx_ut(k,l)))

    return
  end function mo_int_df_dp

  !
  !
  !
  subroutine finalize_df_dp(ints)
  
    class(df_dp)           :: ints  
  
    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)
  
  end subroutine finalize_df_dp

  !---------------------------------------------------------------------------------
  ! Single precision routines
  !
 
  subroutine init_pyscf_df_sp(ints, core_file, eri_file)

    class(df_sp)            :: ints
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

    call read_dataset(f_name, dset_name, ints%h_core, dims_read)

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
    call read_dataset_t(f_name, dset_name, ints%bra_ket, dims_read)

    return

  end subroutine init_pyscf_df_sp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_df_sp(ints, indices, int_vec)

    class(df_sp)           :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp), intent(out)  :: int_vec(:)

    integer(is)            :: i, ij, kl
    integer(is)            :: nints

    if(size(indices, dim=1) /= 4) stop 'mo_ints: indices dim=1 must equal 4'
    nints = size(int_vec)

    do i = 1,nints
        ij = ints%indx_ut(indices(1,i), indices(2,i))
        kl = ints%indx_ut(indices(3,i), indices(4,i))
        int_vec(i) = dot_product(ints%bra_ket(:,ij), ints%bra_ket(:,kl))
    enddo
 
  end subroutine mo_ints_df_sp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_df_sp(ints, i, j, k, l) result(int_val)

    class(df_sp)           :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = dot_product(ints%bra_ket(:,ints%indx_ut(i,j)), \
                          ints%bra_ket(:,ints%indx_ut(k,l)))

    return
  end function mo_int_df_sp

  !
  !
  !
  subroutine finalize_df_sp(ints)

    class(df_sp)              :: ints

    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)

  end subroutine finalize_df_sp

  !---------------------------------------------------------------------------------
  ! half precision routines
  !

  subroutine init_pyscf_df_hp(ints, core_file, eri_file)

    class(df_hp)            :: ints
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

    call read_dataset(f_name, dset_name, ints%h_core, dims_read)

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
    call read_dataset_t(f_name, dset_name, ints%bra_ket, dims_read)

    return
  end subroutine init_pyscf_df_hp 

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_df_hp(ints, indices, int_vec)

    class(df_hp)           :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp), intent(out)  :: int_vec(:)

    integer(is)            :: i, ij, kl
    integer(is)            :: nints

    if(size(indices, dim=1) /= 4) stop 'mo_ints: indices dim=1 must equal 4'
    nints = size(int_vec)

    do i = 1,nints
        ij = ints%indx_ut(indices(1,i), indices(2,i))
        kl = ints%indx_ut(indices(3,i), indices(4,i))
        int_vec(i) = dot_product(ints%bra_ket(:,ij), ints%bra_ket(:,kl))
    enddo

  end subroutine mo_ints_df_hp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_df_hp(ints, i, j, k, l) result(int_val)

    class(df_hp)           :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = dot_product(ints%bra_ket(:,ints%indx_ut(i,j)), \
                          ints%bra_ket(:,ints%indx_ut(k,l)))

    return
  end function mo_int_df_hp

  !
  !
  !
  subroutine finalize_df_hp(ints)

    class(df_hp)              :: ints

    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)

  end subroutine finalize_df_hp


end module dfmod
