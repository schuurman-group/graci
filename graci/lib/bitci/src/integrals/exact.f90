module exactmod

  use h5_ops
  use constants
  use integrals
  
  implicit none

  !
  ! Exact integrals type
  !
  type, extends(eri) :: exact
   real(dp), allocatable      :: bra_ket(:,:)   

   contains
     procedure   :: init_pyscf => init_pyscf_exact
     procedure   :: mo_ints    => mo_ints_exact
     procedure   :: mo_int     => mo_int_exact
     procedure   :: finalize   => finalize_exact
  end type exact

contains

  subroutine init_pyscf_exact(ints, core_file, eri_file)

    class(exact)            :: ints
    character(len=255)      :: core_file
    character(len=255)      :: eri_file

    character(len=255)      :: f_name
    character(len=255)      :: dset_name
    logical                 :: exists
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
    n_bra_ket = ints%nmo * (ints%nmo + 1)/2
    if( any(dims /= (/n_bra_ket, n_bra_ket /))) stop 'ERI data set wrong size, file='//f_name

    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)
    allocate(ints%bra_ket(n_bra_ket, n_bra_ket))

    call read_dataset_dble(f_name, dset_name, ints%bra_ket, dims_read)

    return

  end subroutine init_pyscf_exact

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_exact(ints, indices, int_vec)

    class(exact)           :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp), intent(out)  :: int_vec(:)

    integer(is)            :: i, ij, kl 
    integer(is)            :: nints

    if(size(indices, dim=1) /= 4) stop 'mo_ints: indices dim=1 must equal 4'
    nints = size(int_vec)

    do i = 1,nints
        ij = ints%indx_ut(indices(1,i), indices(2,i))
        kl = ints%indx_ut(indices(3,i), indices(4,i))
        int_vec(i) = ints%bra_ket(ij, kl)
    enddo

  end subroutine mo_ints_exact

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_exact(ints, i, j, k, l) result(int_val)

    class(exact)           :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = ints%bra_ket(ints%indx_ut(i,j), ints%indx_ut(k,l))

    return
  end function mo_int_exact

  !
  !
  !
  subroutine finalize_exact(ints)

    class(exact)           :: ints

    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)

  end subroutine finalize_exact


end module exactmod
