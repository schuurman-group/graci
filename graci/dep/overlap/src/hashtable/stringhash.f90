!######################################################################
! stringhash: Hash table data structure for the storage of alpha/beta
!             (hole) strings
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
     integer(is), allocatable   :: values(:)

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

     ! Resize the hash table
     procedure, non_overridable :: resize_table
     
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
    allocate(h%values(h%n_buckets))
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
    integer(is), optional, intent(in) :: val
    integer(ib)                       :: hash,indx
    integer(ib)                       :: nb,i1,step
    logical                           :: collision

!----------------------------------------------------------------------
! Get the key index
!----------------------------------------------------------------------
        !
    ! Initialise the key index to -1 (a value that cannot correspond to a
    ! hash index for a determinant key).
    ! If the key is found to pre-exist in the table, then this is the
    ! value that will be returned
    !
    indx=-1

    !
    ! Get the hash of the key
    !
    hash=string_hash(h%n_int,key)

    !
    ! Get the hash index for the determinant key
    !
    nb=h%n_buckets
    i1=iand(hash,nb-1)+1

    !
    ! Loop over buckets, and check to see if the key already exists in
    ! the table or if another key with the same index is already present
    !
    collision=.false.
    ! Loop over buckets
    do step=1,h%n_buckets
       
       ! Exit the loop when an empty bucket or the key is found
       if (h%full(i1)==0) then
          ! Empty bucket: the key is not yet in the hash table
          indx=i1
          exit
       else if (keys_equal(h%n_int,h%keys(:,i1),key)) then
          ! The key is already stored in the hash table
          exit
       else
          collision=.true.
          ! Quadratic probing
          !i1=iand(hash+int(0.5*step+0.5*step**2),nb-1)+1
          ! Linear probing
          i1=iand(hash+step,nb-1)+1
       endif

    enddo 

!----------------------------------------------------------------------
! Return if the key is already stored in the hash table
!----------------------------------------------------------------------
    if (indx == -1) return

!----------------------------------------------------------------------
! Update the number of collisions
!----------------------------------------------------------------------
    if (collision) h%n_collisions=h%n_collisions+1

!----------------------------------------------------------------------
! Update the number of filled buckets
!----------------------------------------------------------------------
    h%n_keys_stored=h%n_keys_stored+1

!----------------------------------------------------------------------
! Insert the key into the hash table
!----------------------------------------------------------------------
    h%keys(:,indx)=key

!----------------------------------------------------------------------
! Optional: insert the value into the hash table
!----------------------------------------------------------------------
    if (present(val)) h%values(indx)=val

!----------------------------------------------------------------------
! Insert the hash function value into the auxiliary hashes array
!----------------------------------------------------------------------
    h%hashes(indx)=hash

!----------------------------------------------------------------------
! Switch the full flag
!----------------------------------------------------------------------
    h%full(indx)=1

!----------------------------------------------------------------------
! Increase the size of the hash table if the no. filled buckets
! exceeds the maximum load factor
!----------------------------------------------------------------------
    if (h%n_keys_stored >= h%n_buckets*h%max_load_factor) &
         call h%resize_table
    
    return
    
  end subroutine insert_key

!######################################################################
! string_hash: hash of the bit string representation of an alpha/beta
!              (hole) string using the MurmerHash3 32-bit hash function
!######################################################################
  function string_hash(n_int,string) result(hash)

    use constants

    implicit none

    ! Function result
    integer(ib)             :: hash

    ! Input string
    integer(is), intent(in) :: n_int
    integer(ib), intent(in) :: string(n_int)

    ! Working variables
    integer(is)             :: chunks(2*n_int)
    integer(is)             :: len
    integer(is)             :: seed
    integer(is)             :: i,j,k,t

    ! Bit masks
    integer(ib), parameter  :: mask1=4294967295_ib  ! 11...100...0
    integer(ib), parameter  :: mask2=-4294967296_ib ! 00...011...1

    ! MurmurHash3 32-bit hash function parameters
    integer(is), parameter  :: c1=-862048943 ! 0xcc9e2d51
    integer(is), parameter  :: c2=461845907  ! 0x1b873593
    integer(is), parameter  :: m=5
    integer(is), parameter  :: n=430675100   ! 0xe6546b64
    integer(is), parameter  :: p=-2048144789 ! 0x85ebca6b
    integer(is), parameter  :: q=-1028477387 ! 0xc2b2ae35
    integer(is), parameter  :: r1=15
    integer(is), parameter  :: r2=13

!----------------------------------------------------------------------
! Split the 64-bit bit strings into 32-bit chunks
!----------------------------------------------------------------------
    j=0
    do i=1,n_int
       j=j+1
       chunks(j)=iand(mask1,string(i))
       j=j+1
       chunks(j)=ishft(iand(mask2,string(i)),-32)
    enddo

!----------------------------------------------------------------------    
! Compute the hash function value
!----------------------------------------------------------------------
    seed=0
    len=8*n_int

    t=seed

    do i=1,2*n_int
       k=chunks(i)
       k=k*c1
       k=ishftc(k,r1)
       k=k*c2
       t=ieor(t,k)
       t=ishftc(t,r2)
       t=t*m-n
    enddo

    t=ieor(t,len)

    t=ieor(t,ishft(t,-16))
    t=t*p
    t=ieor(t,ishft(t,-13))
    t=t*q
    t=ieor(t,ishft(t,-16))

    hash=t
    
    return
    
  end function string_hash

!######################################################################
! keys_equal: returns .true. if the two keys are equal, .false.
!             otherwise
!######################################################################
  function keys_equal(n_int,key1,key2)

    use constants
    
    implicit none

    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: key1(n_int),key2(n_int)
    logical                  :: keys_equal
    integer(is)              :: k

    keys_equal=.true.
    
    do k=1,n_int
       if (key1(k) /= key2(k)) then
          keys_equal=.false.
          exit
       endif
    enddo
    
    return
    
  end function keys_equal

!######################################################################
! resize_table: increases the size of the hash table and re-maps all
!               existing keys stored in it
!######################################################################
  subroutine resize_table(h)

    use constants

    implicit none

    class(shtbl), intent(inout) :: h

    integer(ib), allocatable    :: keys_save(:,:)
    integer(is), allocatable    :: values_save(:)
    integer(ib), allocatable    :: hashes_save(:)
    integer(is), allocatable    :: full_save(:)
    integer(is), allocatable    :: deleted_save(:)
    integer(ib), allocatable    :: key(:)
    integer(is)                 :: val
    integer(ib)                 :: hash
    integer(is)                 :: k,bucket,nkeys,n_buckets_old
    integer(ib)                 :: kindx,i1,indx,nb,step

!----------------------------------------------------------------------
! Make a copy of the keys and their hashes
!----------------------------------------------------------------------
    ! Save the old hash table size
    n_buckets_old=h%n_buckets

    ! Allocate the temporary arrays
    allocate(keys_save(h%n_int,n_buckets_old))
    allocate(values_save(n_buckets_old))
    allocate(hashes_save(n_buckets_old))
    allocate(full_save(n_buckets_old))
    allocate(deleted_save(n_buckets_old))
    allocate(key(h%n_int))

    ! Save the old key, value and hash values
    keys_save=h%keys
    values_save=h%values
    hashes_save=h%hashes
    full_save=h%full
    deleted_save=h%deleted

!----------------------------------------------------------------------
! Re-allocate the arrays holding the keys and their hashes
!----------------------------------------------------------------------
    ! Deallocate arrays
    deallocate(h%keys)
    deallocate(h%values)
    deallocate(h%hashes)
    deallocate(h%full)
    deallocate(h%deleted)
    
    ! New hash table size
    h%n_buckets=h%resize_factor*h%n_buckets
    
    ! Allocate arrays
    allocate(h%keys(h%n_int,h%n_buckets))
    allocate(h%values(h%n_buckets))
    allocate(h%hashes(h%n_buckets))
    allocate(h%full(h%n_buckets))
    allocate(h%deleted(h%n_buckets))
    
    ! Initialise the full and deleted flag arrays to 0
    h%full=0
    h%deleted=0

!----------------------------------------------------------------------
! Re-mapping of the keys and their hashes
!----------------------------------------------------------------------
    nkeys=0

    nb=h%n_buckets
    
    ! Loop over old buckets
    do bucket=1,n_buckets_old

       ! Current key
       key=keys_save(:,bucket)

       ! Current value
       val=values_save(bucket)
       
       ! Current hashed key value
       hash=hashes_save(bucket)
       
       ! Cycle if the old bucket is empty...
       if (full_save(bucket)==0) cycle

       ! Cycle if the old bucket was deleted
       if (deleted_save(bucket)==1) cycle
       
       ! ...else we are at the next stored key
       nkeys=nkeys+1
       
       ! Initialise the new key index to -1 (a value that cannot
       ! correspond to a hashed key value)
       indx=-1

       ! New hash index
       i1=iand(hash,nb-1)+1
       
       !
       ! Insert the determinant key into the new hash table
       !
       ! Loop over new buckets
       do step=1,nb

          ! Exit the loop when an empty bucket or the key is found
          if (h%full(i1)==0) then
             ! Empty bucket: the key is not yet in the hash table
             indx=i1
             exit
          else if (keys_equal(h%n_int,h%keys(:,i1),key(:))) then
             ! The key is already stored in the hash table
             print*,'Something terrible has happend...'
             stop
             exit
          else
             ! Quadratic probing
             !i1=iand(hash+int(0.5*step+0.5*step**2),nb-1)+1
             ! Linear probing
             i1=iand(hash+step,nb-1)+1
          endif
          
       enddo

       ! Insert the determinant key
       h%keys(:,indx)=key

       ! Insert the value
       h%values(indx)=val
       
       ! Insert the hash function value into the auxiliary hashes array
       h%hashes(indx)=hash

       ! Switch the full flag
       h%full(indx)=1
       
    enddo

!----------------------------------------------------------------------
! Deallocate the temporary arrays
!----------------------------------------------------------------------
    deallocate(keys_save)
    deallocate(values_save)
    deallocate(hashes_save)
    deallocate(full_save)
    deallocate(deleted_save)
    deallocate(key)
    
    return
    
  end subroutine resize_table
  
!######################################################################
  
end module stringhash
