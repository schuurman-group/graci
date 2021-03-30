!**********************************************************************
! Routines for the generation of guess vectors for the iterative
! diagonalisation of the Hamiltonian matrix
!**********************************************************************

module guessvecs

contains

!######################################################################
! guess_space: Determines the subset of determinants to use in the
!              guess vector generation
!######################################################################
  subroutine guess_space(hdiag,nsym,irrep,guessdim,d,ndet,offset,&
       offdim)

    use constants
    use bitglobal
    use utils
    use dethash
  
    implicit none

    ! Irrep and numbers of determinants
    integer(is), intent(in)    :: irrep
    integer(is), intent(in)    :: nsym(0:nirrep-1)

    ! On-diagonal elements of the Hamiltonian matrix
    real(dp), intent(in)       :: hdiag(nsym(irrep))

    ! Requested guess vector space dimension
    integer(is), intent(inout) :: guessdim

    ! Determinant arrays
    integer(is), intent(in)    :: ndet,offdim
    integer(ib), intent(in)    :: d(n_int,2,ndet)
    integer(is), intent(in)    :: offset(offdim,0:nirrep-1)
  
    ! Sorting arrays
    integer(is), allocatable   :: indx(:)
    
    ! SOPs
    integer(is)                :: nsop
    integer(ib), allocatable   :: sop(:,:,:)
    
    ! Everything else
    integer(is)                :: i

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(indx(nsym(irrep)))
    indx=0
  
!----------------------------------------------------------------------
! Sort the on-diagonal matrix elements in order of increasing value
!----------------------------------------------------------------------  
    call dsortindxa1('A',nsym(irrep),hdiag,indx)

!----------------------------------------------------------------------  
! Get the unique SOPs corresponding to the determinants with the
! lowest value on-diagonal Hamiltonian matrix elements
!----------------------------------------------------------------------  
    call guess_space_sops(sop,nsop,nsym,irrep,d,ndet,offset,offdim,&
         guessdim,indx(1:guessdim))

!----------------------------------------------------------------------
! Generate all possible determinants corresponding to the SOPs
!----------------------------------------------------------------------
    call guess_space_dets(sop,nsop)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(indx)
    deallocate(sop)
  
    return
  
  end subroutine guess_space

!######################################################################
! guess_space_sops: Fills in the hash table h with the unique SOPs
!                   corresponding the the guess space determinants
!######################################################################
  subroutine guess_space_sops(sop,nsop,nsym,irrep,d,ndet,offset,&
       offdim,guessdim,indx)

    use constants
    use bitglobal
    use detutils
    use dethash

    implicit none

    ! SOPs
    integer(is)                :: nsop
    integer(ib), allocatable   :: sop(:,:,:)
  
    ! Irrep and numbers of determinants
    integer(is), intent(in)    :: irrep
    integer(is), intent(in)    :: nsym(0:nirrep-1)
    
    ! Requested guess space determinants
    integer(is), intent(inout) :: guessdim
    integer(is), intent(in)    :: indx(guessdim)
    
    ! Determinant arrays
    integer(is), intent(in)    :: ndet,offdim
    integer(ib), intent(in)    :: d(n_int,2,ndet)
    integer(is), intent(in)    :: offset(offdim,0:nirrep-1)

    ! SOP hash table
    type(dhtbl)                :: h
    integer(is)                :: initial_size
    integer(ib)                :: key(n_int,2)
  
    ! Everything else
    integer(is)                :: i,i1
    integer(is)                :: nsym_sum(0:nirrep-1)

!----------------------------------------------------------------------
! Cumulative numbers of determinants in the symmetry blocks
!----------------------------------------------------------------------
    nsym_sum=0
    nsym_sum(1:nirrep-1)=nsym(0:nirrep-2)
    do i=1,nirrep-1
       nsym_sum(i)=nsym_sum(i)+nsym_sum(i-1)
    enddo
    
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------  
    initial_size=10
    call h%initialise_table(initial_size)

!----------------------------------------------------------------------
! Insert the guess space SOPs into the hash table
!----------------------------------------------------------------------
    ! Loop over guess space determinants
    do i1=1,guessdim
       
       ! Index of the current determinant
       i=indx(i1)+nsym_sum(irrep)
       
       ! Calculate the SOP
       call get_sop(d(:,:,i),key)

       ! Insert the SOP into the hash table
       call h%insert_key(key)
     
    enddo

!----------------------------------------------------------------------
! Fill in the SOP array
!----------------------------------------------------------------------
    nsop=h%n_keys_stored
    allocate(sop(n_int,2,nsop))
    call h%retrieve_keys(sop)
  
!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
    call h%delete_table
  
    return
  
  end subroutine guess_space_sops
  
!######################################################################
! guess_space_dets: From a set of guess space SOPs, generates all
!                   possible determinants. The resulting set of
!                   determinants can be used to construct spin-pure
!                   guess vectors
!######################################################################
  subroutine guess_space_dets(sop,nsop)

    use constants
    use bitglobal
    use detutils
    
    implicit none

    ! SOPs
    integer(is), intent(in) :: nsop
    integer(ib), intent(in) :: sop(n_int,2,nsop)

    print*,"HERE"
    stop
    
    return
    
  end subroutine guess_space_dets

!######################################################################
  
end module guessvecs
