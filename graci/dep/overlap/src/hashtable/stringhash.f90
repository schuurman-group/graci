!######################################################################
! stringhash: Hash table data structure with:
!
!             Keys: alpha/beta hole strings
!
!             Values: tuples of (i)  parent string indices
!                               (ii) annihilated MO indices
!
!----------------------------------------------------------------------
!             Uses open addressing with linear probing and the 32-bit
!             MurmurHash3 hash function.
!----------------------------------------------------------------------
!             Adapted from the bitci dethash module, which in turn
!             borrows some ideas and code snippets from the FFHASH
!             code of Jannis Teunissen
!             [https://github.com/jannisteunissen/ffhash].
!######################################################################
module stringhash

  use constants

  private

  ! Derived data type for a hash table storing alpha/beta hole strings
  type, public :: shtbl

     ! Number of (n_bits)-bit integers required to represent
     ! each hole string
     integer(is)                :: n_int
     
     ! Number of buckets in the hash table
     integer(ib)                :: n_buckets=0

     ! Number of keys stored in the hash table
     integer(ib)                :: n_keys_stored=0

     ! Number of collisions
     integer(ib)                :: n_collisions=0
     
     ! Keys (hole strings) stored in the hash table
     integer(ib), allocatable   :: keys(:,:)

     ! Values stored in the hash table
     integer(is), allocatable   :: values(:,:)

     ! Hashes of the keys stored in the hash table
     integer(ib), allocatable   :: hashes(:)

     ! Empty/full bucket flag
     integer(is), allocatable   :: full(:)

     ! Deleted bucket flag
     integer(is), allocatable   :: deleted(:)

     ! Maximum load factor for the hash table
     real(dp)                   :: max_load_factor=0.5d0

     ! Hash table resizing factor
     real(dp)                   :: resize_factor=2.0d0
     
   contains

     ! Initialise the hash table
     procedure, non_overridable :: initialise_table

     ! Delete the hash table
     procedure, non_overridable :: delete_table

     ! Insert the key (determinant) into the hash table
     procedure, non_overridable :: insert_key
     
  end type shtbl

contains

!######################################################################
! initialise_table: initialises the hash table to a user-specified
!                   size (no. buckets)
!######################################################################
  subroutine initialise_table(h,n_int,initial_size)

    use constants

    implicit none

    class(shtbl), intent(inout) :: h
    integer(is), intent(in)     :: initial_size
    integer(is), intent(in)     :: n_int

    ! Set n_int
    h%n_int=n_int
    
    ! Set the initial hash table size
    h%n_buckets=initial_size

    ! Round up the initial hash table size to a power of 2
    ! Note that this is essential for both the replacement of modulo
    ! operations by iand operations, and is also necessary if we
    ! ever switch to quadratic probing
    h%n_buckets=2**ceiling(log(dble(h%n_buckets))/log(2.0d0))

     ! Set the number of keys stored to zero
    h%n_keys_stored=0_ib
    
    ! Allocate arrays
    if (allocated(h%keys)) deallocate(h%keys)
    if (allocated(h%values)) deallocate(h%values)
    if (allocated(h%hashes)) deallocate(h%hashes)
    if (allocated(h%full)) deallocate(h%full)
    if (allocated(h%deleted)) deallocate(h%deleted)
    allocate(h%keys(h%n_int,h%n_buckets))
    allocate(h%values(2,h%n_buckets))
    allocate(h%hashes(h%n_buckets))
    allocate(h%full(h%n_buckets))
    allocate(h%deleted(h%n_buckets))

    ! Initialise arrays
    h%values=0
    h%full=0
    h%deleted=0
    
    return
    
  end subroutine initialise_table

!######################################################################
! delete_table: deletes the hash table
!######################################################################
  subroutine delete_table(h)

    use constants

    implicit none

    class(shtbl), intent(inout) :: h

    ! Deallocate arrays
    deallocate(h%keys)
    deallocate(h%values)
    deallocate(h%hashes)
    deallocate(h%full)
    deallocate(h%deleted)
    
    ! Zero the number of buckets, keys stored, and collisions
    h%n_buckets=0
    h%n_keys_stored=0
    h%n_collisions=0
    
    return
    
  end subroutine delete_table

!######################################################################
! insert_key: stores a single key and a value tuple into the hash
!######################################################################
  subroutine insert_key(h,key,val)

    use constants

    implicit none

    class(shtbl), intent(inout)       :: h
    integer(ib), intent(in)           :: key(h%n_int)
    integer(is), optional, intent(in) :: val(2)

!----------------------------------------------------------------------
! Get the key index
!----------------------------------------------------------------------
    
    
    STOP
    
    return
    
  end subroutine insert_key
    
!######################################################################
  
end module stringhash
