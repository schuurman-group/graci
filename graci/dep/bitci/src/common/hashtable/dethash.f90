!######################################################################
! dethash: Routines for the creation of a hash table to store the
!          unique determinants, configurations and SOPs. Also allows
!          for the optional storage of 4-byte integer values, which
!          can be used for simple labeling of the stored det/conf/SOP
!          keys.
!----------------------------------------------------------------------
!          Uses open addressing with linear probing and the 32-bit
!          MurmurHash3 hash function.
!----------------------------------------------------------------------
!          The resizing algorithm is quite crude: upon increasing the
!          number of buckets, all the previous keys are re-mapped.
!          To make this fast, the hashes of the keys are also stored.
!          Then, each re-mapping only requires a single modulo
!          operation. Because each key is a determinant represented
!          by 2*N_int (N_bits)-bit integers, the memory overhead
!          introduced by also storing the hashed keys is relatively
!          small - only one extra (N_bits)-bit integer per key.
!----------------------------------------------------------------------
!          This borrows some ideas and code snippets from the FFHASH
!          code of Jannis Teunissen
!          [https://github.com/jannisteunissen/ffhash].
!######################################################################
module dethash

  use constants
  use detutils, only: det_hash
  
  private

  ! Derived data type for a hash table storing determinants,
  ! configurations and SOPs
  type, public :: dhtbl

     ! Number of buckets in the hash table
     integer(ib)                :: n_buckets=0

     ! Number of keys stored in the hash table
     integer(ib)                :: n_keys_stored=0

     ! Number of collisions
     integer(ib)                :: n_collisions=0
     
     ! Keys (determinants) stored in the hash table
     integer(ib), allocatable   :: keys(:,:,:)

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
     
     ! Retrieve the index of a given key (determinant)
     procedure, non_overridable :: get_index

     ! Retrieve the value for a given key
     procedure, non_overridable :: get_value
     
     ! Insert the key (determinant) into the hash table
     procedure, non_overridable :: insert_key

     ! Check if the key is already stored in the hash
     ! table
     procedure, non_overridable :: key_exists
     
     ! Resize the hash table
     procedure, non_overridable :: resize_table

     ! Dump the stored keys to disk
     procedure, non_overridable :: dump_keys

     ! Return the keys in a single array
     procedure, non_overridable :: retrieve_keys

     ! Delete the key (determinant) from the hash table
     procedure, non_overridable :: delete_key
     
  end type dhtbl

contains

!######################################################################
! initialise_table: initialises the hash table to a user-specified
!                   size (no. buckets)
!######################################################################
  subroutine initialise_table(h,initial_size)

    use constants

    implicit none

    class(dhtbl), intent(inout) :: h
    integer(is), intent(in)     :: initial_size
    
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
    allocate(h%keys(n_int,2,h%n_buckets))
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

    class(dhtbl), intent(inout) :: h

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
! get_index: retrieve the index of a given key (determinant)
!######################################################################
  function get_index(h,key) result(kindx)

    use constants
    
    implicit none

    class(dhtbl), intent(in) :: h
    integer(ib), intent(in)  :: key(n_int,2)
    integer(ib)              :: kindx
    integer(ib)              :: hash,i1,step,nb

!----------------------------------------------------------------------
! Initialise the key index to -1 (a value that cannot correspond to a
! hash index for a determinant key).
! If the key is found to pre-exist in the table, then this is the
! value that will be returned
!----------------------------------------------------------------------
    kindx=-1
    
!----------------------------------------------------------------------
! Get the hash index for the determinant key
!----------------------------------------------------------------------
    hash=det_hash(key)
    nb=h%n_buckets
    i1=iand(hash,nb-1)+1

!----------------------------------------------------------------------
! Loop over buckets, and check to see if the key already exists in
! the table or if another key with the same index is already present
!----------------------------------------------------------------------
    ! Loop over buckets
    do step=1,h%n_buckets

       ! Exit the loop when an empty bucket or the key is found
       if (h%full(i1)==0) then
          ! Empty bucket: the key is not yet in the hash table
          exit
       else if (keys_equal(h%keys(:,:,i1),key)) then
          ! The key is already stored in the hash table
          kindx=i1
          exit
       else
          ! Quadratic probing
          !i1=iand(hash+int(0.5*step+0.5*step**2),nb-1)+1
          ! Linear probing
          i1=iand(hash+step,nb-1)+1
       endif
       
    enddo
    
    return
    
  end function get_index

!######################################################################
! get_value: retrieve the value associated with a given key
!######################################################################
  function get_value(h,key) result(val)

    use constants

    implicit none

    class(dhtbl), intent(in) :: h
    integer(ib), intent(in)  :: key(n_int,2)
    integer(is)              :: val
    integer(ib)              :: indx

    !
    ! Retrieve the index of the key
    !
    indx=h%get_index(key)

    !
    ! Return a value of -1 if the key is not in the table
    !
    if (indx == -1) then
       val=-1
       return
    endif
    
    !
    ! Value
    !
    val=h%values(indx)
    
    return
    
  end function get_value
    
!######################################################################
! keys_equal: returns .true. if the two determinant keys are equal,
!             .false. otherwise
!######################################################################
  function keys_equal(key1,key2)

    use constants
    
    implicit none

    integer(ib), intent(in)  :: key1(n_int,2),key2(n_int,2)
    logical                  :: keys_equal
    integer(is)              :: i,k

    keys_equal=.true.
    
    do i=1,2
       do k=1,n_int
          if (key1(k,i).ne.key2(k,i)) then
             keys_equal=.false.
             exit
          endif
       enddo
    enddo
    
    return
    
  end function keys_equal

!######################################################################
! insert_key: stores a single key and, optionally, value in the hash
!             table if the key is not already stored
!######################################################################
  subroutine insert_key(h,key,val)
    
    use constants
    
    implicit none

    class(dhtbl), intent(inout)       :: h
    integer(ib), intent(in)           :: key(n_int,2)
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
    hash=det_hash(key)
    
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
       else if (keys_equal(h%keys(:,:,i1),key)) then
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
! Insert the key (determinant) into the hash table
!----------------------------------------------------------------------
    h%keys(:,:,indx)=key

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
    if (h%n_keys_stored>=h%n_buckets*h%max_load_factor) &
         call h%resize_table

    return
    
  end subroutine insert_key

!######################################################################
! key_exists: returns .true. if a given key already exists, else
!             returns .false.
!######################################################################
  function key_exists(h,key)

    use constants

    implicit none

    logical                     :: key_exists
    
    class(dhtbl), intent(inout) :: h
    integer(ib), intent(in)     :: key(n_int,2)
    integer(ib)                 :: hash
    integer(ib)                 :: nb,i1,step

    !
    ! Get the hash of the key
    !
    hash=det_hash(key)
    
    !
    ! Get the hash index for the key
    !
    nb=h%n_buckets
    i1=iand(hash,nb-1)+1

    !
    ! Loop over buckets, and check to see if the key already exists in
    ! the table
    !
    ! Loop over buckets
    do step=1,h%n_buckets
       
       ! Exit the loop when an empty bucket or the key is found
       if (h%full(i1)==0) then
          ! Empty bucket: the key is not yet in the hash table
          key_exists=.false.
          return
       else if (keys_equal(h%keys(:,:,i1),key)) then
          ! The key is already stored in the hash table
          key_exists=.true.
          return
       else
          ! Collision: the key could still be in the hash table,
          ! move to the next bucket
          ! Quadratic probing
          !i1=iand(hash+int(0.5*step+0.5*step**2),nb-1)+1
          ! Linear probing
          i1=iand(hash+step,nb-1)+1
       endif

    enddo

    !
    ! If we are here, then the key is not stored in the hash table
    !
    key_exists=.false.
    
    return
    
  end function key_exists
    
!######################################################################
! resize_table: increases the size of the hash table and re-maps all
!               existing keys stored in it
!######################################################################
  subroutine resize_table(h)

    use constants
    
    implicit none

    class(dhtbl), intent(inout) :: h

    integer(ib), allocatable    :: keys_save(:,:,:)
    integer(is), allocatable    :: values_save(:)
    integer(ib), allocatable    :: hashes_save(:)
    integer(is), allocatable    :: full_save(:)
    integer(is), allocatable    :: deleted_save(:)
    integer(ib), allocatable    :: key(:,:)
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
    allocate(keys_save(n_int,2,n_buckets_old))
    allocate(values_save(n_buckets_old))
    allocate(hashes_save(n_buckets_old))
    allocate(full_save(n_buckets_old))
    allocate(deleted_save(n_buckets_old))
    allocate(key(n_int,2))

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
    allocate(h%keys(n_int,2,h%n_buckets))
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
       key=keys_save(:,:,bucket)

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
          else if (keys_equal(h%keys(:,:,i1),key(:,:))) then
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
       h%keys(:,:,indx)=key

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
! dump_keys: buffered write of the stored keys to disk
!######################################################################
  subroutine dump_keys(h,filename)

    use constants
    use iomod, only: freeunit
    
    implicit none
    
    class(dhtbl), intent(in)     :: h
    character(len=*), intent(in) :: filename

    integer(is)                  :: unit,bucket,counter
    integer(is)                  :: buffer_size,n_records
    integer(ib), allocatable     :: buffer(:,:,:)

!----------------------------------------------------------------------
! Set the buffer size and allocate the buffer
!----------------------------------------------------------------------
    ! Buffer size: maximum number of determinant keys that fit into
    ! 10 MB
    buffer_size=10*1024**2/(8*2*n_int)

    ! Allocate the buffer
    allocate(buffer(n_int,2,buffer_size))
    
    ! Number of records
    n_records=ceiling(dble(h%n_keys_stored)/dble(buffer_size))
    
!----------------------------------------------------------------------
! Open the determinant file
!----------------------------------------------------------------------
    call freeunit(unit)
    open(unit,file=filename,form='unformatted',status='unknown')

!----------------------------------------------------------------------
! File header
!----------------------------------------------------------------------
    ! Number of determinants
    write(unit) h%n_keys_stored

    ! Buffer size
    write(unit) buffer_size

    ! Number of records
    write(unit) n_records
    
!----------------------------------------------------------------------
! Write the determinants to disk
!----------------------------------------------------------------------
    ! Initialise the counter
    counter=0
    
    ! Loop over buckets
    do bucket=1,h%n_buckets

       ! If the bucket is filled, put the key into the buffer
       if (h%full(bucket)==1 .and. h%deleted(bucket)==0) then
          counter=counter+1
          buffer(:,:,counter)=h%keys(:,:,bucket)
       endif

       ! If the buffer is full, write it to disk
       if (counter == buffer_size) then
          write(unit) counter,buffer
          counter=0
       endif
       
    enddo

    ! Last record
    if (counter > 0) write(unit) counter,buffer

!----------------------------------------------------------------------
! Close the determinant file
!----------------------------------------------------------------------
    close(unit)

!----------------------------------------------------------------------
! Deallocate the buffer
!----------------------------------------------------------------------
    deallocate(buffer)
    
    return
    
  end subroutine dump_keys

!######################################################################
! retrieve_keys: returns the stored keys in the array key_array
!######################################################################
  subroutine retrieve_keys(h,key_array)

    use constants

    implicit none

    class(dhtbl), intent(in) :: h
    integer(ib), intent(out) :: key_array(n_int,2,h%n_keys_stored)

    integer(is) :: counter,bucket

!----------------------------------------------------------------------
! Extract the stored keys from the hash table
!----------------------------------------------------------------------
    ! Initialise the counter
    counter=0
    
    ! Loop over buckets
    do bucket=1,h%n_buckets

       ! If the bucket is filled, put the key into the array
       if (h%full(bucket)==1 .and. h%deleted(bucket)==0) then
          counter=counter+1
          key_array(:,:,counter)=h%keys(:,:,bucket)
       endif

    enddo

    return
    
  end subroutine retrieve_keys

!######################################################################
! delete_key: deletes a single key that is already stored in the hash
!             table
!######################################################################
  subroutine delete_key(h,key)

    use constants

    class(dhtbl), intent(inout) :: h
    integer(ib), intent(in)     :: key(n_int,2)
    integer(ib)                 :: hash,indx
    integer(ib)                 :: nb,i1,step
    
!----------------------------------------------------------------------
! Get the key index
!----------------------------------------------------------------------
    !
    ! Initialise the key index to -1 (a value that cannot correspond to a
    ! hash index for a determinant key).
    ! If the key is not found to pre-exist in the table, then this is the
    ! value that will be returned
    !
    indx=-1

    !
    ! Get the hash of the key
    !
    !hash=key_hash(h,key)
    hash=det_hash(key)
    
    !
    ! Get the hash index for the determinant key
    !
    nb=h%n_buckets
    i1=iand(hash,nb-1)+1

    !
    ! Loop over buckets, and check to see if the key already exists in
    ! the table or if another key with the same index is already present
    !
    ! Loop over buckets
    do step=1,h%n_buckets
       
       ! Exit the loop when an empty bucket or the key is found
       if (h%full(i1)==0) then
          ! Empty bucket: the key is not yet in the hash table
          exit
       else if (keys_equal(h%keys(:,:,i1),key)) then
          ! The key is already stored in the hash table
          indx=i1
          exit
       else
          ! Quadratic probing
          !i1=iand(hash+int(0.5*step+0.5*step**2),nb-1)+1
          ! Linear probing
          i1=iand(hash+step,nb-1)+1
       endif

    enddo

!----------------------------------------------------------------------
! Return if the key is not stored in the hash table
!----------------------------------------------------------------------
    if (indx == -1) return
    
!----------------------------------------------------------------------
! Update the number of filled buckets
!----------------------------------------------------------------------
    h%n_keys_stored=h%n_keys_stored-1

!----------------------------------------------------------------------
! Delete the key (determinant) from the hash table
!----------------------------------------------------------------------
    h%keys(:,:,indx)=0_ib

!----------------------------------------------------------------------
! Delete the value from the hash table
!----------------------------------------------------------------------
    h%values(indx)=0
    
!----------------------------------------------------------------------
! Delete the hash function value from the auxiliary hashes array
!----------------------------------------------------------------------
    h%hashes(indx)=0_ib

!----------------------------------------------------------------------
! Flag the bucket as deleted
!----------------------------------------------------------------------
    h%deleted(indx)=1
    
    return
    
  end subroutine delete_key
  
!######################################################################
  
end module dethash
