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
  subroutine get_sigma_holes(irrepB,irrepK,nel_before,nmo,mosym,n_int,&
       nsigma,sigma,nsigmaH,sigmaH,Hinfo_dim,Hinfo,offHinfo)

    use constants
    
    implicit none

    ! Irrep generated by the wave functions
    integer(is), intent(in)  :: irrepB,irrepK
    
    ! Number of electrons before the sigma strings in the
    ! determinants
    integer(is), intent(in)  :: nel_before

    ! Irreps generated by the MOs
    integer(is), intent(in)  :: nmo
    integer(ib), intent(in)  :: mosym(nmo)
    
    ! Number of unique sigma strings
    integer(is), intent(in)  :: nsigma

    ! Unique sigma strings
    integer(is), intent(in)  :: n_int
    integer(ib), intent(in)  :: sigma(n_int,nsigma)

    ! Unique sigma-hole strings
    integer(is), intent(out) :: nsigmaH
    integer(ib), allocatable :: sigmaH(:,:)

    ! Sigma-hole string info
    integer(is), intent(out) :: Hinfo_dim
    integer(is), allocatable :: Hinfo(:,:)

    ! Sigma-hole string info offsets
    integer(is), allocatable :: offHinfo(:)
    
    ! Everything else
    integer(ib)              :: iBra,iKet,irrepBK
    
!----------------------------------------------------------------------
! Direct product of the bra and ket wave function irreps
!----------------------------------------------------------------------
    ! Convert the irreps to 64-bit integers
    iBra=irrepB
    iKet=irrepK

    ! Direct product of the irreps
    irrepBK=ieor(iBra,iKet)

!----------------------------------------------------------------------
! (1) Determine the number of unique sigma-hole strings
!----------------------------------------------------------------------
    call get_nsigmaH(n_int,nsigma,sigma,irrepBK,nmo,mosym,nsigmaH,&
         Hinfo_dim)
    
!----------------------------------------------------------------------
! (2) Save the unique sigma-hole strings along with the corresponding
!     lists of annihilation operator indices, parent sigma
!     string indices and phases
!----------------------------------------------------------------------
    allocate(sigmaH(n_int,nsigmaH))
    sigmaH=0_ib
    
    allocate(Hinfo(3,Hinfo_dim))
    Hinfo=0

    allocate(offHinfo(nsigmaH+1))
    offHinfo=0
    
    call get_hole_strings(n_int,nsigma,sigma,nel_before,irrepBK,nmo,&
         mosym,nsigmaH,sigmaH,Hinfo_dim,Hinfo,offHinfo)

    return
    
  end subroutine get_sigma_holes

!######################################################################
! get_nsigmaH: given a set of unique sigma strings + the irreps
!              generated by the bra and ket WFs and MOs, determines
!              the number of unique sigma-hole strings
!######################################################################
  subroutine get_nsigmaH(n_int,nsigma,sigma,irrepBK,nmo,mosym,nsigmaH,&
       Hinfo_dim)

    use constants
    use detfuncs
    use stringhash
    
    implicit none

    ! Unique sigma strings
    integer(is), intent(in)  :: n_int,nsigma
    integer(ib), intent(in)  :: sigma(n_int,nsigma)

    ! Irrep generatated by the product of bra and ket wave functions
    integer(ib), intent(in)  :: irrepBK

    ! Irreps generated by the ket MOs
    integer(is), intent(in)  :: nmo
    integer(ib), intent(in)  :: mosym(nmo)
    
    ! Number of unique sigma-hole strings
    integer(is), intent(out) :: nsigmaH

    ! Number of (hole index, parent string, phase) tuples
    ! needed to characterise the sigma-hole strings
    integer(is), intent(out) :: Hinfo_dim
    
    ! Hash table
    type(shtbl)              :: h
    integer(is)              :: initial_size

    ! Hole strings
    integer(ib)              :: key(n_int)
    
    ! Occpied MOs
    integer(is)              :: nocc
    integer(is), allocatable :: occ(:)
    
    ! Everything else
    integer(is)              :: i,j
    integer(is)              :: n_ap
    integer(is), allocatable :: ap(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(occ(nmo))
    allocate(ap(nmo))
    
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------
    initial_size=nsigma
    call h%initialise_table(n_int,initial_size)

!----------------------------------------------------------------------
! Determine the number of unique sigma-hole strings
!----------------------------------------------------------------------
    Hinfo_dim=0

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

          ! Update the no. (hole index, parent string, phase) tuples
          Hinfo_dim=Hinfo_dim+1
          
          ! Key: sigma-hole string
          key=annihilate_electron_string(n_int,sigma(:,i),&
               ap(j))

          ! Insert the key
          call h%insert_key(key)
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! The number of unique sigma-hole strings is the number of keys
! stored in the hash table
!----------------------------------------------------------------------
    nsigmaH=h%n_keys_stored
    
    return
    
  end subroutine get_nsigmaH

!######################################################################
! get_hole_strings: returns the unique sigma-hole strings along with
!                   the generating annihilation operator indices,
!                   parent sigma string indices and phase factors
!######################################################################
  subroutine get_hole_strings(n_int,nsigma,sigma,nel_before,irrepBK,&
       nmo,mosym,nsigmaH,sigmaH,Hinfo_dim,Hinfo,offHinfo)

    use constants
    use detfuncs
    use stringhash
    
    implicit none

    ! Unique sigma strings
    integer(is), intent(in)  :: n_int,nsigma
    integer(ib), intent(in)  :: sigma(n_int,nsigma)

    ! Number of electrons before the sigma strings in the
    ! determinants
    integer(is), intent(in)  :: nel_before
    
    ! Irrep generatated by the product of bra and ket wave functions
    integer(ib), intent(in)  :: irrepBK

    ! Irreps generated by the ket MOs
    integer(is), intent(in)  :: nmo
    integer(ib), intent(in)  :: mosym(nmo)
    
    ! Unique sigma-hole strings
    integer(is), intent(in)  :: nsigmaH
    integer(ib), intent(out) :: sigmaH(n_int,nsigmaH)
    
    ! (hole index, parent string, phase) tuples
    ! needed to characterise the sigma-hole strings
    integer(is), intent(in)  :: Hinfo_dim
    integer(is), intent(out) :: Hinfo(3,Hinfo_dim)

    ! Sigma-hole string info offsets
    integer(is), intent(out) :: offHinfo(nsigmaH+1)
    
    ! Number of sigma strings generated by each sigma-hole string
    integer(is), allocatable :: counter(:)

    ! Hash table
    type(shtbl)              :: h
    integer(is)              :: initial_size

    ! Hole strings
    integer(ib)              :: key(n_int)
    integer(is)              :: values(2)
    
    ! Occpied MOs
    integer(is)              :: nocc
    integer(is), allocatable :: occ(:)
    
    ! Everything else
    integer(is)              :: i,j,nold,indx,iloc
    integer(is)              :: n_ap
    integer(is), allocatable :: ap(:)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(counter(nsigmaH))
    counter=0

    allocate(occ(nmo))
    occ=0

    allocate(ap(nmo))
    ap=0
    
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------
    initial_size=nsigma
    call h%initialise_table(n_int,initial_size)

!----------------------------------------------------------------------
! Insert the sigma-hole strings into the hash table along with their
! indices
!----------------------------------------------------------------------
    nold=0
    indx=1

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

          ! Key: sigma-hole string
          key=annihilate_electron_string(n_int,sigma(:,i),&
               ap(j))

          ! Insert the (key, indx) pair
          call h%insert_key(key,indx)

          ! Update the unique sigma-hole index
          if (h%n_keys_stored > nold) then
             nold=h%n_keys_stored
             indx=indx+1
          endif
          
       enddo
       
    enddo
        
!----------------------------------------------------------------------
! Store the unique sigma-hole strings and determine the no. sigma
! strings generated by each
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

          ! Key: sigma-hole string
          key=annihilate_electron_string(n_int,sigma(:,i),&
               ap(j))

          ! Get the sigma-hole string index
          indx=h%get_value(key)

          ! Update the number of sigma strings generated by this
          ! sigma-hole string
          counter(indx)=counter(indx)+1
          
          ! Save the sigma-hole string
          if (counter(indx) == 1) sigmaH(:,indx)=key
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Fill in the sigma-hole string offset array
!----------------------------------------------------------------------
    offHinfo(1)=1
    do i=2,nsigmaH+1
       offHinfo(i)=offHinfo(i-1)+counter(i-1)
    enddo

!----------------------------------------------------------------------
! Save the (hole index, parent string, phase) tuple
!----------------------------------------------------------------------
    counter=0

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

          ! Key: sigma-hole string
          key=annihilate_electron_string(n_int,sigma(:,i),&
               ap(j))

          ! Get the sigma-hole string index
          indx=h%get_value(key)

          ! Update the counter for this sigma-hole string
          counter(indx)=counter(indx)+1

          ! Location in the Hinfo array
          iloc=offHinfo(indx)+counter(indx)-1

          ! Hole index
          Hinfo(1,iloc)=ap(j)

          ! Parent sigma string index
          Hinfo(2,iloc)=i

          ! Phase factor
          Hinfo(3,iloc)=phase_fac(n_int,sigma(:,i),ap(j),nel_before)
          
       enddo
       
    enddo
    
    return
    
  end subroutine get_hole_strings
    
!######################################################################
! phase_fac: given a sigma string, determines the phase factor
!            corresponding to annihilation of an electron in the imo'th
!            position
!            Here, nel_before is the no. electrons appearing before
!            the sigma-spin block
!######################################################################
  function phase_fac(n_int,string,imo,nel_before)

    use constants
    
    implicit none

    integer(is)             :: phase_fac
    integer(is), intent(in) :: n_int,imo,nel_before
    integer(ib), intent(in) :: string(n_int)
    integer(is)             :: kimo,k,i,imo1,nset

    !
    ! No. bits (MO's) set before the imo'th position
    !
    nset=0
    
    !
    ! Block corresponding to the index imo
    !
    k=(imo-1)/n_bits+1

    !
    ! Contributions from the blocks 1,...,kimo-1
    !
    do i=1,k-1
       nset=nset+popcnt(string(i))
    enddo

    !
    ! Contribution from the block that the bit indexed imo sits in
    !
    imo1=imo-(k-1)*n_bits-1
    nset=nset+popcnt(ishft(string(k),n_bits-imo1))
    
    !
    ! Phase factor
    !
    phase_fac=(-1)**(nset+nel_before)
    
    return
    
  end function phase_fac

!######################################################################
  
end module hole_strings
