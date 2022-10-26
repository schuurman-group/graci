!**********************************************************************
! Routines for the generation of unique sigma-hole strings, i.e., the
! result of the operation of all allowed annihilation operators on the
! set of unique sigma strings
!**********************************************************************
module hole_strings

  implicit none
  
contains

!######################################################################
! get_sigma_holes: given an array of sigma strings, generates all
!                  unique sigma-hole strings
!######################################################################
  subroutine get_sigma_holes(irrepB,irrepK,isigma,nmo,mosym,n_int,&
       nsigma,sigma)

    use constants
    use detfuncs

    implicit none

    ! Irrep generated by the wave functions
    integer(is), intent(in)  :: irrepB,irrepK
    
    ! Sigma spin index
    integer(is), intent(in)  :: isigma

    ! Irreps generated by the MOs
    integer(is), intent(in)  :: nmo
    integer(ib), intent(in)  :: mosym(nmo)
    
    ! Number of unique sigma strings
    integer(is), intent(in)  :: nsigma

    ! Unique sigma strings
    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: sigma(n_int,nsigma)

    ! Hole strings
    integer(ib)              :: hole_string(n_int)
    
    ! Occpied MOs
    integer(is)              :: nocc
    integer(is), allocatable :: occ(:)
    
    ! Everything else
    integer(is)              :: i,j
    integer(is)              :: n_ap
    integer(is), allocatable :: ap(:)
    integer(ib)              :: iBra,iKet,irrepBK

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(occ(nmo))
    allocate(ap(nmo))
    
!----------------------------------------------------------------------
! Direct product of the bra and ket wave function irreps
!----------------------------------------------------------------------
    ! Convert the irreps to 64-bit integers
    iBra=irrepB
    iKet=irrepK

    ! Direct product of the irreps
    irrepBK=ieor(iBra,iKet)
    
!----------------------------------------------------------------------
! Determine the unique sigma-hole strings along with the corresponding
! annihilation operator indices, phases and parent sigma-string
! indices
!----------------------------------------------------------------------
    ! Loop over unique sigma strings
    do i=1,nsigma

       ! Get the list of occupied MOs
       call mo_occ_string(n_int,sigma(:,i),nmo,nocc,occ)

       ! Reduce the list of MOs to those of symmetry
       ! irrepB \otimes irrepK
       n_ap=0
       do j=1,nocc
          if (mosym(occ(j)) == irrepBK) then
             n_ap=n_ap+1
             ap(n_ap)=occ(j)
          endif
       enddo

       ! Generate the sigma-hole strings and insert into the hash
       ! table
       do j=1,n_ap

          hole_string=annihilate_electron_string(n_int,sigma(:,i),&
               ap(j))

       enddo
       
    enddo

    
    STOP

    
    return
    
  end subroutine get_sigma_holes

!######################################################################
  
end module hole_strings
