module exactmod

  !use h5_ops
  use constants
  use integrals
  use iomod

  implicit none

  !
  ! Exact integrals type
  !
  type, extends(eri) :: exact_dp
   real(dp), allocatable      :: bra_ket(:,:)   

   contains
     procedure   :: init_pyscf => init_pyscf_exact_dp
     procedure   :: mo_ints    => mo_ints_exact_dp
     procedure   :: mo_int     => mo_int_exact_dp
     procedure   :: finalize   => finalize_exact_dp
  end type exact_dp

  !
  ! Exact integrals type
  !
  type, extends(eri) :: exact_sp
   real(sp), allocatable      :: bra_ket(:,:)

   contains
     procedure   :: init_pyscf => init_pyscf_exact_sp
     procedure   :: mo_ints    => mo_ints_exact_sp
     procedure   :: mo_int     => mo_int_exact_sp
     procedure   :: finalize   => finalize_exact_sp
  end type exact_sp

  !
  ! Exact integrals type
  !
  type, extends(eri) :: exact_hp
   real(hp), allocatable      :: bra_ket(:,:)

   contains
     procedure   :: init_pyscf => init_pyscf_exact_hp
     procedure   :: mo_ints    => mo_ints_exact_hp
     procedure   :: mo_int     => mo_int_exact_hp
     procedure   :: finalize   => finalize_exact_hp
  end type exact_hp


contains

  !-----------------------------------------------------------
  ! Double precision routines
  !
  !

  !
  !
  !
  subroutine init_pyscf_exact_dp(ints, core_file, eri_file)

    class(exact_dp)         :: ints
    character(len=255)      :: core_file
    character(len=255)      :: eri_file

    character(len=255)      :: f_name
    character(len=255)      :: dset_name
    real(dp),allocatable    :: dp_record(:)
    logical                 :: exists
    integer(is)             :: i
    integer(is)             :: unit
    integer(is)             :: rank
    integer(is)             :: dims(2)
    integer(is)             :: dims_read(2)
    integer(is)             :: n_bra_ket

    ! load the one hamiltonian
    !------------------------------------------------------------
    f_name    = trim(adjustl(core_file))
    dset_name = 'hcore_mo'
    inquire(file=f_name, exist=exists)

    if(.not.exists) stop 'cannot find hcore_mo file='//f_name

    ! read in the tensor dimensions
    ! -----------------------------
    call freeunit(unit)
    open(unit, file=f_name, form='unformatted')
    
    do i = 1,2
      read(unit) dims(i)
    enddo
    ints%nmo = dims(1)

    if(allocated(ints%h_core))deallocate(ints%h_core)
    allocate(ints%h_core(ints%nmo, ints%nmo))

    do i = 1,ints%nmo
      read(unit)ints%h_core(:, i)
    enddo

    close(unit)

    ! load ERI
    !--------------------------------------------------------------
    f_name    = trim(adjustl(eri_file))
    dset_name = 'eri_mo'
    inquire(file=f_name, exist=exists)

    if(.not.exists) stop 'cannot find eri_mo in file='//f_name

    ! read in the tensor dimensions
    ! -----------------------------
    call freeunit(unit)
    open(unit, file=f_name, form='unformatted')
    
    do i = 1,2
      read(unit) dims(i)
    enddo

    n_bra_ket  = ints%nmo * (ints%nmo + 1)/2
    if( any(dims /= (/n_bra_ket, n_bra_ket /))) stop 'ERI data set wrong size, file='//f_name

    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)
    allocate(ints%bra_ket(n_bra_ket, n_bra_ket))

    ! reads in the 2e ERIs column by column
    do i = 1, n_bra_ket
      read(unit)ints%bra_ket(:,i)
    enddo
    close(unit)

    return

  end subroutine init_pyscf_exact_dp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_exact_dp(ints, indices, int_vec)

    class(exact_dp)        :: ints
    integer(is),intent(in) :: indices(:,:)
    real(dp),intent(out)   :: int_vec(:)

    integer(is)            :: i, ij, kl 
    integer(is)            :: nints

    if(size(indices, dim=1) /= 4) stop 'mo_ints: indices dim=1 must equal 4'
    nints = size(int_vec)

    do i = 1,nints
        ij = ints%indx_ut(indices(1,i), indices(2,i))
        kl = ints%indx_ut(indices(3,i), indices(4,i))
        int_vec(i) = ints%bra_ket(ij, kl)
    enddo

  end subroutine mo_ints_exact_dp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_exact_dp(ints, i, j, k, l) result(int_val)

    class(exact_dp)        :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = ints%bra_ket(ints%indx_ut(i,j), ints%indx_ut(k,l))

    return
  end function mo_int_exact_dp

  !
  !
  !
  subroutine finalize_exact_dp(ints)
  
    class(exact_dp)        :: ints
  
    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)
  
  end subroutine finalize_exact_dp

 !---------------------------------------------------------------------
 ! Single precision routines
 !
 !

  !
  !
  !
  subroutine init_pyscf_exact_sp(ints, core_file, eri_file)

    class(exact_sp)         :: ints
    character(len=255)      :: core_file
    character(len=255)      :: eri_file

    character(len=255)      :: f_name
    character(len=255)      :: dset_name
    real(dp),allocatable    :: dp_record(:)
    logical                 :: exists
    integer(is)             :: i
    integer(is)             :: unit
    integer(is)             :: rank
    integer(is)             :: dims(2)
    integer(is)             :: dims_read(2)
    integer(is)             :: n_bra_ket

    ! load the one hamiltonian
    !------------------------------------------------------------
    f_name    = trim(adjustl(core_file))
    dset_name = 'hcore_mo'
    inquire(file=f_name, exist=exists)

    if(.not.exists) stop 'cannot find hcore_mo file='//f_name

    ! read in the tensor dimensions
    ! -----------------------------
    call freeunit(unit)
    open(unit, file=f_name, form='unformatted')
    
    do i = 1,2
      read(unit) dims(i)
    enddo
    ints%nmo = dims(1)

    if(allocated(ints%h_core))deallocate(ints%h_core)
    allocate(ints%h_core(ints%nmo, ints%nmo))

    do i = 1,ints%nmo
      read(unit)ints%h_core(:, i)
    enddo

    close(unit)

    ! load ERI
    !--------------------------------------------------------------
    f_name    = trim(adjustl(eri_file))
    dset_name = 'eri_mo'
    inquire(file=f_name, exist=exists)

    if(.not.exists) stop 'cannot find eri_mo in file='//f_name

     ! read in the tensor dimensions
    ! -----------------------------
    call freeunit(unit)
    open(unit, file=f_name, form='unformatted')
    
    do i = 1,2
      read(unit) dims(i)
    enddo

    n_bra_ket       = ints%nmo * (ints%nmo + 1)/2
    if( any(dims /= (/n_bra_ket, n_bra_ket /))) stop 'ERI data set wrong size,file='//f_name

    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)
    allocate(ints%bra_ket(n_bra_ket, n_bra_ket))
    allocate(dp_record(n_bra_ket))

    do i = 1, n_bra_ket
      read(unit)dp_record
      ints%bra_ket(:,i) = dp_record
    enddo
    close(unit)

    deallocate(dp_record)

    return

  end subroutine init_pyscf_exact_sp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_exact_sp(ints, indices, int_vec)

    class(exact_sp)        :: ints
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

  end subroutine mo_ints_exact_sp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_exact_sp(ints, i, j, k, l) result(int_val)

    class(exact_sp)        :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = ints%bra_ket(ints%indx_ut(i,j), ints%indx_ut(k,l))

    return
  end function mo_int_exact_sp

  !
  !
  !
  subroutine finalize_exact_sp(ints)

    class(exact_sp)        :: ints

    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)

  end subroutine finalize_exact_sp

 !---------------------------------------------------------------------
 ! Half precision routines
 !
 !

  !
  !
  !
  subroutine init_pyscf_exact_hp(ints, core_file, eri_file)

    class(exact_hp)         :: ints
    character(len=255)      :: core_file
    character(len=255)      :: eri_file

    character(len=255)      :: f_name
    character(len=255)      :: dset_name
    real(dp),allocatable    :: dp_record(:)
    logical                 :: exists
    integer(is)             :: i
    integer(is)             :: unit
    integer(is)             :: rank
    integer(is)             :: dims(2)
    integer(is)             :: dims_read(2)
    integer(is)             :: n_bra_ket

    ! load the one hamiltonian
    !------------------------------------------------------------
    f_name    = trim(adjustl(core_file))
    dset_name = 'hcore_mo'
    inquire(file=f_name, exist=exists)

    if(.not.exists) stop 'cannot find hcore_mo file='//f_name

    ! read in the tensor dimensions
    ! -----------------------------
    call freeunit(unit)
    open(unit, file=f_name, form='unformatted')
    
    do i = 1,2
      read(unit) dims(i)
    enddo
    ints%nmo = dims(1)

    if(allocated(ints%h_core))deallocate(ints%h_core)
    allocate(ints%h_core(ints%nmo, ints%nmo))

    do i = 1,ints%nmo
      read(unit)ints%h_core(:, i)
    enddo

    close(unit)

    ! load ERI
    !--------------------------------------------------------------
    f_name    = trim(adjustl(eri_file))
    dset_name = 'eri_mo'
    inquire(file=f_name, exist=exists)

    if(.not.exists) stop 'cannot find eri_mo in file='//f_name

    ! read in the tensor dimensions
    ! -----------------------------
    call freeunit(unit)
    open(unit, file=f_name, form='unformatted')
    
    do i = 1,2
      read(unit) dims(i)
    enddo

    n_bra_ket       = ints%nmo * (ints%nmo + 1)/2
    if( any(dims /= (/n_bra_ket, n_bra_ket /))) stop 'ERI data set wrong size,file='//f_name

    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)
    allocate(ints%bra_ket(n_bra_ket, n_bra_ket))
    allocate(dp_record(n_bra_ket))

    do i = 1, n_bra_ket 
      read(unit)dp_record
      ints%bra_ket(:,i) = dp_record
    enddo
    close(unit)

    deallocate(dp_record)

    return
   end subroutine init_pyscf_exact_hp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  subroutine mo_ints_exact_hp(ints, indices, int_vec)

    class(exact_hp)        :: ints
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

  end subroutine mo_ints_exact_hp

  !
  ! Function that takes a list of mo indices and returns list
  ! of integrals 
  !
  function mo_int_exact_hp(ints, i, j, k, l) result(int_val)

    class(exact_hp)        :: ints
    integer(is),intent(in) :: i, j, k, l

    real(dp)               :: int_val

    int_val = ints%bra_ket(ints%indx_ut(i,j), ints%indx_ut(k,l))

    return
  end function mo_int_exact_hp

  !
  !
  !
  subroutine finalize_exact_hp(ints)

    class(exact_hp)        :: ints

    if(allocated(ints%h_core))deallocate(ints%h_core)
    if(allocated(ints%bra_ket))deallocate(ints%bra_ket)

  end subroutine finalize_exact_hp


end module exactmod
