!######################################################################
! dethash: Routines for the creation of a hash table to store the
!          unique determinants generated from single and double
!          excitations out of the reference space.
!----------------------------------------------------------------------
!          Uses open addressing with linear probing.
!          The hash algorithm used is the boost hash_combine function 
!          modified for use with 64-bit integer array keys.
!----------------------------------------------------------------------
!          Note that the sole purpose of this hash table is to be able
!          to efficiently filter out duplicate determinants as they
!          are generated. As such, it is a simple and somewhat
!          incomplete data structure (e.g., no deletion function).
!----------------------------------------------------------------------
!          The resizing algorithm is quite crude: upon increasing the
!          number of buckets, all the previous keys are re-mapped.
!          To make this fast, the hashes of the keys are also stored.
!          Then, each re-mapping only requires a single modulo
!          operation. Because each key is a determinant represented
!          by 2*N_int 64-bit integers, the memory overhead introduced
!          by also storing the hashed keys is relatively small - only
!          one extra 64-bit integer per key.
!----------------------------------------------------------------------
!          This borrows some ideas and code snippets from the FFHASH
!          code of Jannis Teunissen
!          [https://github.com/jannisteunissen/ffhash].
!######################################################################
module dethash

  use constants
  
  private

  ! Derived data type for a hash table storing determinants
  type, public :: dhtbl

     ! Number of buckets in the hash table
     integer(ib)                :: n_buckets=0

     ! Number of keys stored in the hash table
     integer(ib)                :: n_keys_stored=0

     ! Keys (determinants) stored in the hash table
     integer(ib), allocatable   :: keys(:,:,:)

     ! Hashes of the keys stored in the hash table
     integer(ib), allocatable   :: hashes(:)

     ! Empty/full bucket flag
     integer(is), allocatable   :: full(:)
     
     ! Maximum load factor for the hash table
     real(dp)                   :: max_load_factor=0.7d0

     ! Hash table resizing factor
     real(dp)                   :: resize_factor=1.5d0
     
   contains

     ! Initialise the hash table
     procedure, non_overridable :: initialise_table

     ! Delete the hash table
     procedure, non_overridable :: delete_table
     
     ! Retrieve the index of a given key (determinant)
     procedure, non_overridable :: get_index

     ! Insert the key (determinant) into the hash table
     procedure, non_overridable :: insert_key

     ! Resize the hash table
     procedure, non_overridable :: resize_table

     ! Dump the stored keys to disk
     procedure, non_overridable :: dump_keys
     
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

    ! Set the number of keys stored to zero
    h%n_keys_stored=0_ib
    
    ! Allocate arrays
    if (allocated(h%keys)) deallocate(h%keys)
    if (allocated(h%hashes)) deallocate(h%hashes)
    if (allocated(h%full)) deallocate(h%full)
    allocate(h%keys(n_int,2,initial_size))
    allocate(h%hashes(initial_size))
    allocate(h%full(initial_size))

    ! Initialise the full flag array to 0
    h%full=0
    
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
    deallocate(h%hashes)
    deallocate(h%full)

    ! Zero the number of buckets and keys stored
    h%n_buckets=0
    h%n_keys_stored=0
    
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
!    i1=hash_index(h,key)

    hash=key_hash(h,key)
    nb=h%n_buckets
    i1=modulo(hash,nb)+1
    
!----------------------------------------------------------------------
! Loop over buckets, and check to see if the key already exists in
! the table or if another key with the same index is already present
!----------------------------------------------------------------------
    ! Loop over buckets
    do step=1,h%n_buckets

       ! Exit the loop when an empty bucket or the key is found
       if (h%full(i1)==0) then
          ! Empty bucket: the key is not yet in the hash table
          kindx=i1
          exit
       else if (keys_equal(h%keys(:,:,i1),key)) then
          ! The key is already stored in the hash table
          exit
       else
          i1=modulo(hash+step,nb)+1
       endif
          
    enddo
    
    return
    
  end function get_index

!######################################################################
! key_hash: returns the hash of a determinant key. Uses a modified
!            version of the boost hash_combine function.
!######################################################################
    function key_hash(h,key) result(hash)

    use constants

    implicit none

    class(dhtbl), intent(in) :: h
    integer(ib), intent(in)  :: key(n_int,2)
    integer(ib)              :: hash
    integer(is)              :: i,k
    integer(ib)              :: nb

!----------------------------------------------------------------------
! Hash of the determinant key
!----------------------------------------------------------------------
    hash=0_il
    do i=1,2
       do k=1,n_int
          hash=ieor(hash,key(k,i)+7046029254386353130_il+ishft(hash,6)&
               +ishft(hash,-2))
       enddo
    enddo
    
    return
    
  end function key_hash
  
!######################################################################
! hash_index: returns the hash of a key modulo the size of the hash
!             table. Uses a modified version of the boost hash_combine
!             function
!######################################################################
  function hash_index(h,key) result(hindx)

    use constants

    implicit none

    class(dhtbl), intent(in) :: h
    integer(ib), intent(in)  :: key(n_int,2)
    integer(ib)              :: hindx
    integer(ib)              :: hash
    integer(is)              :: i,k
    integer(ib)              :: nb

!----------------------------------------------------------------------
! Hash of the determinant key
!----------------------------------------------------------------------
    hash=0_il
    do i=1,2
       do k=1,n_int
          hash=ieor(hash,key(k,i)+7046029254386353130_il+ishft(hash,6)&
               +ishft(hash,-2))
       enddo
    enddo

!----------------------------------------------------------------------
! Hash of the determinant key modulo the hash table size
!----------------------------------------------------------------------
    nb=h%n_buckets
    hindx=modulo(hash,nb)+1
    
    return
    
  end function hash_index
    
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
! insert_key: stores a single key in the hash table if it is not
!             already in the hash table
!######################################################################
  subroutine insert_key(h,key)
    
    use constants
    
    implicit none

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
    ! If the key is found to pre-exist in the table, then this is the
    ! value that will be returned
    !
    indx=-1

    !
    ! Get the hash of the key
    !
    hash=key_hash(h,key)
    
    !
    ! Get the hash index for the determinant key
    !
    nb=h%n_buckets
    i1=modulo(hash,nb)+1
    
    !
    ! Loop over buckets, and check to see if the key already exists in
    ! the table or if another key with the same index is already present
    !
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
          i1=modulo(hash+step,nb)+1
       endif

    enddo
    
!----------------------------------------------------------------------
! Return if the key is already stored in the hash table
!----------------------------------------------------------------------
    if (indx.eq.-1) return
       
!----------------------------------------------------------------------
! Update the number of filled buckets
!----------------------------------------------------------------------
    h%n_keys_stored=h%n_keys_stored+1

!----------------------------------------------------------------------
! Insert the key (determinant) into the hash table
!----------------------------------------------------------------------
    h%keys(:,:,indx)=key

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
! resize_table: increases the size of the hash table and re-maps all
!               existing keys stored in it
!######################################################################
  subroutine resize_table(h)

    use constants
    
    implicit none

    class(dhtbl), intent(inout) :: h

    integer(ib), allocatable    :: keys_save(:,:,:)
    integer(ib), allocatable    :: hashes_save(:)
    integer(is), allocatable    :: full_save(:)
    integer(ib), allocatable    :: key(:,:)
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
    allocate(hashes_save(n_buckets_old))
    allocate(full_save(n_buckets_old))
    allocate(key(n_int,2))

    ! Save the old key and hash values
    keys_save=h%keys
    hashes_save=h%hashes
    full_save=h%full
    
!----------------------------------------------------------------------
! Re-allocate the arrays holding the keys and their hashes
!----------------------------------------------------------------------
    ! Deallocate arrays
    deallocate(h%keys)
    deallocate(h%hashes)
    deallocate(h%full)
    
    ! New hash table size
    h%n_buckets=h%resize_factor*h%n_buckets
    
    ! Allocate arrays
    allocate(h%keys(n_int,2,h%n_buckets))
    allocate(h%hashes(h%n_buckets))
    allocate(h%full(h%n_buckets))

    ! Initialise the full flag array to 0
    h%full=0
    
!----------------------------------------------------------------------
! Re-mapping of the keys and their hashes
!----------------------------------------------------------------------
    nkeys=0

    nb=h%n_buckets
    
    ! Loop over old buckets
    do bucket=1,n_buckets_old

       ! Current key
       key=keys_save(:,:,bucket)

       ! Current hashed key value
       hash=hashes_save(bucket)
       
       ! Cycle if the old bucket is empty...
       if (full_save(bucket)==0) cycle
       
       ! ...else we are at the next stored key
       nkeys=nkeys+1
       
       ! Initialise the new key index to -1 (a value that cannot
       ! correspond to a hashed key value)
       indx=-1

       ! New hash index
       i1=modulo(hash,nb)+1

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
             i1=modulo(hash+step,nb)+1
          endif
          
       enddo

       ! Insert the determinant key
       h%keys(:,:,indx)=key

       ! Insert the hash function value into the auxiliary hashes array
       h%hashes(indx)=hash

       ! Switch the full flag
       h%full(indx)=1
       
    enddo

!----------------------------------------------------------------------
! Deallocate the temporary arrays
!----------------------------------------------------------------------
    deallocate(keys_save)
    deallocate(hashes_save)
    deallocate(full_save)
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
       if (h%full(bucket)==1) then
          counter=counter+1
          buffer(:,:,counter)=h%keys(:,:,bucket)
       endif

       ! If the buffer is full, write it to disk
       if (counter.eq.buffer_size) then
          write(unit) counter,buffer
          counter=0
       endif
       
    enddo

    ! Last record
    if (counter.gt.0) write(unit) counter,buffer

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
  
end module dethash
