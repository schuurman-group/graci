!**********************************************************************
! Utility routines useful for MRCI calculations
!**********************************************************************

module mrciutils

  implicit none

contains

!######################################################################
! sop_sym_mrci: Given an SOP in the MRCI ordering, returns the integer
!               label of the irrep thatit generates. Here, m2c is the
!               canonical MO to MRCI-ordered MO mapping array.
!######################################################################
  function sop_sym_mrci(sop,m2c)
    
    use constants
    use bitglobal
    
    implicit none

    integer(is)             :: sop_sym_mrci
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: m2c(nmo)
    integer(is)             :: k,ipos,imo
    integer(ib)             :: h,isym
    
    !
    ! Direct product of the irreps generated by the singly-occupied
    ! orbitals
    !
    ! Initialisation to the totally symmetric irrep (which is always
    ! encoded as 0)
    isym=0_ib
    
    ! Loop over blocks
    do k=1,n_int

       ! Initialise the work array
       h=sop(k,1)

       ! Loop over the singly-occupied orbital indices
       do while (h /= 0_ib)

          ! Number of trailing zeros left in h
          ipos=trailz(h)

          ! Index of the next singly-occupied orbital
          imo=1+ipos+(k-1)*64

          ! Direct product
          isym=ieor(mosym(m2c(imo)),isym)
          
          ! Clear the bits up to the current singly-occupied orbital
          h=ibclr(h,ipos)
          
       enddo
          
    enddo

    sop_sym_mrci=isym
    
    return
    
  end function sop_sym_mrci
  
!######################################################################
! exc_degree_conf: Computes the excitation degree between two
!                  configurations c1 and c2
!######################################################################
  function exc_degree_conf(c1,c2,ldc) result(nexci)

    use constants
        
    implicit none

    integer(is)              :: nexci
    integer(is), intent(in)  :: ldc
    integer(ib), intent(in)  :: c1(ldc,2),c2(ldc,2)
    integer(is)              :: k

    ! Compute the number of creation and annihilation operators
    ! connecting c1 and c2
    nexci=popcnt(ieor(c1(1,1),c2(1,1)))+popcnt(ieor(c1(1,2),c2(1,2)))
    do k=2,ldc
       nexci=nexci+popcnt(ieor(c1(k,1),c2(k,1))) &
            +popcnt(ieor(c1(k,2),c2(k,2)))
    enddo

    ! Divide by 2 to get the no. excitations connecting c1 and c2
    ! Note that this is equivalent to performing a right bit shift
    ! by 1 place
    nexci=shiftr(nexci,1)
    
    return
    
  end function exc_degree_conf

!######################################################################
! n_bits_set_before: Given an array of bit strings, b, returns the
!                    number of set bits before the position ipos.
!                    Here, the bit in the i'th position in b(k) has
!                    the index k*64+(i-1), where i=0,1,...,63.
!######################################################################
  function n_bits_set_before(b,bdim,ipos) result(nbits)

    use constants
    
    implicit none

    integer(is)             :: nbits
    integer(is), intent(in) :: bdim,ipos
    integer(ib), intent(in) :: b(bdim)
    integer(is)             :: kipos,k,ipos1

    !
    ! Initialisation
    !
    nbits=0
    
    !
    ! Block corresponding to the index ipos
    !
    kipos=(ipos-1)/64+1

    !
    ! Contributions from the blocks 1,...,kipos-1
    !
    do k=1,kipos-1
       nbits=nbits+popcnt(b(k))
    enddo

    !
    ! Contribution from the block that the bit indexed ipos sits in
    !
    ipos1=ipos-(kipos-1)*64-1
    nbits=nbits+popcnt(ishft(b(kipos),64-ipos1))
    
    return
    
  end function n_bits_set_before

!######################################################################
! conf_occ_list: Given a configuration c, returns the list of occupied
!                MOs. nmo1 is the number of MOs represented by c.
!######################################################################
  subroutine conf_occ_list(c,ldc,list,nmo1,n)

    use constants
    
    implicit none

    integer(is), intent(in)  :: ldc,nmo1
    integer(ib), intent(in)  :: c(ldc,2)
    integer(is), intent(out) :: list(nmo1)
    integer(is), intent(out) :: n
    integer(ib)              :: iocc(ldc)
    integer(is)              :: k,e
    
    !
    ! Bit string encoding of the occupied MOs
    !
    do k=1,ldc
       iocc(k)=ior(c(k,1),c(k,2))
    enddo
    
    !
    ! Indices of the occupied MOs
    !
    ! Initialisation
    list=0
    n=1
    ! Loop over blocks
    do k=1,ldc
       
       ! Determine the indices of any set bits
       do while (iocc(k) /= 0_ib)
          e=trailz(iocc(k))
          iocc(k)=ibclr(iocc(k),e)
          list(n)=e+(k-1)*64+1
          n=n+1
       enddo
       
    enddo

    !
    ! Number of occupied spatial MOs
    !
    n=n-1
    
    return
    
  end subroutine conf_occ_list

!######################################################################
! conf_unocc_list: Given a configuration c, returns the list of
!                  unoccupied MOs. nmo1 is the number of MOs
!                  represented by c.
!######################################################################
  subroutine conf_unocc_list(c,ldc,list,nmo1,n)

    use constants
    
    implicit none

    integer(is), intent(in)  :: ldc,nmo1
    integer(ib), intent(in)  :: c(ldc,2)
    integer(is), intent(out) :: list(nmo1)
    integer(is), intent(out) :: n
    integer(ib)              :: iunocc(ldc)
    integer(is)              :: k,e,last

    !
    ! Bit string encoding of the unoccupied MOs
    !
    do k=1,ldc
       iunocc(k)=iand(not(c(k,1)),not(c(k,2)))
    enddo

    !
    ! Clear the unused bits at the end of the bitstring
    !
    last=nmo1-(ldc-1)*64
    iunocc(ldc)=ibits(iunocc(ldc),0,last)
    
    !
    ! Indices of the unoccupied MOs
    !
    ! Initialisation
    list=0
    n=1
    ! Loop over blocks
    do k=1,ldc
       
       ! Determine the indices of any set bits
       do while (iunocc(k) /= 0_ib)
          e=trailz(iunocc(k))
          iunocc(k)=ibclr(iunocc(k),e)
          list(n)=e+(k-1)*64+1
          n=n+1
       enddo
       
    enddo

    !
    ! Number of unoccupied spatial MOs
    !
    n=n-1
    
    return
    
  end subroutine conf_unocc_list

!######################################################################
! sop_socc_list: Given an SOP, returns the list of singly-occupied MOs.
!                nmo1 is the number of MOs represented by the SOP.
!######################################################################
  subroutine sop_socc_list(sop,ldsop,list,nmo1,n)

    use constants
    
    implicit none

    integer(is), intent(in)  :: ldsop,nmo1
    integer(ib), intent(in)  :: sop(ldsop,2)
    integer(is), intent(out) :: list(nmo1)
    integer(is), intent(out) :: n
    integer(is)              :: k,e
    integer(ib)              :: sop1(ldsop)

    !
    ! Initialisation
    !
    list=0
    n=1

    !
    ! Working array that we can destroy
    !
    sop1=sop(:,1)

    !
    ! Indices of the singly-occupied MOs
    !
    ! Loop over blocks
    do k=1,ldsop
       
       ! Determine the indices of any set bits
       do while (sop1(k) /= 0_ib)
          e=trailz(sop1(k))
          sop1(k)=ibclr(sop1(k),e)
          list(n)=e+(k-1)*64+1
          n=n+1
       enddo
       
    enddo

    !
    ! Number of singly-occupied MOs
    !
    n=n-1

    return
    
  end subroutine sop_socc_list

!######################################################################
! sop_docc_list: Given an SOP, returns the list of doubly-occupied MOs.
!                nmo1 is the number of MOs represented by the SOP.
!######################################################################
  subroutine sop_docc_list(sop,ldsop,list,nmo1,n)

    use constants
    
    implicit none

    integer(is), intent(in)  :: ldsop,nmo1
    integer(ib), intent(in)  :: sop(ldsop,2)
    integer(is), intent(out) :: list(nmo1)
    integer(is), intent(out) :: n
    integer(is)              :: k,e
    integer(ib)              :: sop2(ldsop)

    !
    ! Initialisation
    !
    list=0
    n=1

    !
    ! Working array that we can destroy
    !
    sop2=sop(:,2)

    !
    ! Indices of the doubly-occupied MOs
    !
    ! Loop over blocks
    do k=1,ldsop
       
       ! Determine the indices of any set bits
       do while (sop2(k) /= 0_ib)
          e=trailz(sop2(k))
          sop2(k)=ibclr(sop2(k),e)
          list(n)=e+(k-1)*64+1
          n=n+1
       enddo
       
    enddo

    !
    ! Number of doubly-occupied MOs
    !
    n=n-1
    
    return
    
  end subroutine sop_docc_list

!######################################################################
! sop_unocc_list: Given an SOP, returns the list of unoccupied MOs.
!                 nmo1 is the number of MOs represented by the SOP.
!######################################################################
  subroutine sop_unocc_list(sop,ldsop,list,nmo1,n)

     use constants
    
    implicit none

    integer(is), intent(in)  :: ldsop,nmo1
    integer(ib), intent(in)  :: sop(ldsop,2)
    integer(is), intent(out) :: list(nmo1)
    integer(is), intent(out) :: n
    integer(is)              :: k,e,last
    integer(ib)              :: iunocc(ldsop)
    
    !
    ! Bit string encoding of the unoccupied MOs
    !
    do k=1,ldsop
       iunocc(k)=iand(not(sop(k,1)),not(sop(k,2)))
    enddo

    !
    ! Clear the unused bits at the end of the bitstring
    !
    last=nmo1-(ldsop-1)*64
    iunocc(ldsop)=ibits(iunocc(ldsop),0,last)

    !
    ! Indices of the unoccupied MOs
    !
    ! Initialisation
    list=0
    n=1
    ! Loop over blocks
    do k=1,ldsop
       
       ! Determine the indices of any set bits
       do while (iunocc(k) /= 0_ib)
          e=trailz(iunocc(k))
          iunocc(k)=ibclr(iunocc(k),e)
          list(n)=e+(k-1)*64+1
          n=n+1
       enddo
       
    enddo

    !
    ! Number of unoccupied spatial MOs
    !
    n=n-1
    
    return
    
  end subroutine sop_unocc_list
    
!######################################################################
! sop_nopen: Given a SOP, returns the number of open shells
!######################################################################
  function sop_nopen(sop,ldsop) result(nopen)

    use constants
    
    implicit none

    integer(is)             :: nopen
    integer(is), intent(in) :: ldsop
    integer(ib), intent(in) :: sop(ldsop,2)
    integer(is)             :: k
    
    !
    ! Determine the number of singly occupied MOs
    !
    nopen=0
    do k=1,ldsop
       nopen=nopen+popcnt(sop(k,1))
    enddo
    
    return
    
  end function sop_nopen
    
!######################################################################
! get_exci_indices: Determines the indices of the MOs in the
!                   excitation(s) linking two configurations, c1 and
!                   c2. Here, c1 is the ket configuration and c2 is
!                   the bra configuration.
!######################################################################
  subroutine get_exci_indices(c1,c2,ldc,hlist,plist,nexci)

    use constants
    
    implicit none

    integer(is), intent(in)  :: ldc
    integer(ib), intent(in)  :: c1(ldc,2),c2(ldc,2)
    integer(is), intent(in)  :: nexci
    integer(is), intent(out) :: hlist(nexci),plist(nexci)
    integer(ib)              :: x(ldc,2)
    integer(ib)              :: h,p
    integer(is)              :: i,k,ipos,imo,n
    integer(is)              :: hlist1(nexci),plist1(nexci)
    
!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs involved in the
! excitation
!----------------------------------------------------------------------
    do k=1,ldc
       x(k,1)=ieor(c1(k,1),c2(k,1))
       x(k,2)=ieor(c1(k,2),c2(k,2))
    enddo

!----------------------------------------------------------------------
! Hole indices
!----------------------------------------------------------------------
    hlist1=0
    n=0
    do i=1,2
       do k=1,ldc
          
          h=iand(x(k,i),c1(k,i))

          do while (h /= 0_ib)
             ipos=trailz(h)
             imo=1+ipos+(k-1)*64
             n=n+1
             hlist1(n)=imo
             h=ibclr(h,ipos)
          enddo
    
       enddo
    enddo

!----------------------------------------------------------------------
! Particle indices
!----------------------------------------------------------------------
    plist1=0
    n=0
    do i=1,2
       do k=1,ldc
    
          p=iand(x(k,i),c2(k,i))
    
          do while (p /= 0_ib)
             ipos=trailz(p)
             imo=1+ipos+(k-1)*64
             n=n+1
             plist1(n)=imo
             p=ibclr(p,ipos)
          enddo
          
       enddo
    enddo

!----------------------------------------------------------------------
! Sort the particle and hole indices into ascending order
!
! This isn't necessary...
!----------------------------------------------------------------------
    !if (nexci == 1) then
    !   hlist=hlist1
    !   plist=plist1
    !else
    !   hlist(1)=minval(hlist1)
    !   hlist(2)=maxval(hlist1)
    !   plist(1)=minval(plist1)
    !   plist(2)=maxval(plist1)
    !endif

    hlist=hlist1
    plist=plist1
    
    return
    
  end subroutine get_exci_indices
  
!######################################################################
! diffconf: Determines the list of non-zero MO occupation differences
!           of the configuration c relative to the base configuration.
!           Only the blocks 1:ldc are considered.
!######################################################################
  subroutine diffconf(c,ldc,diff,dim,ndiff)
    
    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in)  :: ldc
    integer(ib), intent(in)  :: c(ldc,2)
    integer(is), intent(in)  :: dim
    integer(is), intent(out) :: diff(dim,2)
    integer(is), intent(out) :: ndiff
    
    integer(ib)              :: x(ldc,2)
    integer(ib)              :: h(ldc,2),p(ldc,2)
    integer(ib)              :: d(ldc,2),a(ldc,2)
    integer(ib)              :: g
    integer(is)              :: i,k,ipos,imo

!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs involved in the
! excitation
!----------------------------------------------------------------------
    do k=1,ldc
       x(k,1)=ieor(c(k,1),conf0(k,1))
       x(k,2)=ieor(c(k,2),conf0(k,2))
    enddo

!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs with Delta w_i < 0
!----------------------------------------------------------------------
    do k=1,ldc
       h(k,1)=iand(x(k,1),conf0(k,1))
       h(k,2)=iand(x(k,2),conf0(k,2))
    enddo

!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs with Delta w_i > 0
!----------------------------------------------------------------------
    do k=1,ldc
       p(k,1)=iand(x(k,1),c(k,1))
       p(k,2)=iand(x(k,2),c(k,2))
    enddo

!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs with Delta w_i = -1
!----------------------------------------------------------------------
    do k=1,ldc
       d(k,1)=ieor(h(k,1),h(k,2))
    enddo

!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs with Delta w_i = -2
!----------------------------------------------------------------------
    do k=1,ldc
       d(k,2)=iand(h(k,1),h(k,2))
    enddo
    
!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs with Delta w_i = +1
!----------------------------------------------------------------------
    do k=1,ldc
       a(k,1)=ieor(p(k,1),p(k,2))
    enddo

!----------------------------------------------------------------------
! Bit string encoding of the positions of the MOs with Delta w_i = +2
!----------------------------------------------------------------------
    do k=1,ldc
       a(k,2)=iand(p(k,1),p(k,2))
    enddo

!----------------------------------------------------------------------
! List of Delta w_i values
!----------------------------------------------------------------------
    ! Initialisation
    ndiff=0
    diff=0

    ! Negative Delta w_i values
    do i=1,2
       do k=1,ldc
          g=d(k,i)
          do while (g /= 0_ib)
             ipos=trailz(g)
             imo=1+ipos+(k-1)*64
             ndiff=ndiff+1
             diff(ndiff,1)=imo
             diff(ndiff,2)=-i
             g=ibclr(g,ipos)          
          enddo
       enddo
    enddo

    ! Positive Delta w_i values
    do i=1,2
       do k=1,ldc
          g=a(k,i)
          do while (g /= 0_ib)
             ipos=trailz(g)
             imo=1+ipos+(k-1)*64
             ndiff=ndiff+1
             diff(ndiff,1)=imo
             diff(ndiff,2)=i
             g=ibclr(g,ipos)          
          enddo
       enddo
    enddo

    return
    
  end subroutine diffconf
  
!######################################################################
! nobefore: Given an SOP, determines the the number of open shells
!           preceding each MO. Returns the result in the array 'list'.
!######################################################################
  subroutine nobefore(sop,list)

    use constants
    use bitglobal
    
    implicit none

    integer(ib), intent(in)  :: sop(n_int,2)
    integer(is), intent(out) :: list(nmo)
    integer(ib)              :: h
    integer(is)              :: k,ipos,imo,count,ilast

    integer(is) :: no
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    list=0

!----------------------------------------------------------------------
! Determine the numbers of open shells before each MO
!----------------------------------------------------------------------
    ! Counters
    count=0
    ilast=1
    
    ! Loop over blocks
    do k=1,n_int

       ! Initialise the work array
       h=sop(k,1)

       ! Loop over all singly-occupied orbital indices
       do while (h /= 0_ib)
          
          ! Number of trailing zeros left in h
          ipos=trailz(h)

          ! Index of the next singly-occupied orbital
          imo=1+ipos+(k-1)*64

          ! Fill in the array
          list(ilast:imo)=count
          
          ! Update the counters
          count=count+1
          ilast=imo+1
          
          ! Clear the bits up to the current singly-occupied orbital
          h=ibclr(h,ipos)
          
       enddo
       
    enddo

    ! Remaining MOs
    list(ilast:nmo)=count
    
    return
    
  end subroutine nobefore

!######################################################################
! conf_to_sop: given a configuration bitstring, returns the
!              SOP bitstring
!######################################################################
  function conf_to_sop(conf,ldconf) result(sop)

    use constants
    use bitglobal
    
    implicit none

    ! Configuration
    integer(is), intent(in) :: ldconf
    integer(ib), intent(in) :: conf(ldconf,2)
    
    ! Function result
    integer(ib)             :: sop(ldconf,2)
    
    integer(is)             :: k

    ! Loop over blocks
    do k=1,ldconf

       ! Encoding of the singly-occupied orbitals
       sop(k,1)=ieor(conf(k,1),conf(k,2))

       ! Encoding of the doubly-occupied orbitals
       sop(k,2)=iand(conf(k,1),conf(k,2))
       
    enddo
    
    return
    
  end function conf_to_sop

!######################################################################
! annihilate_electron: given a configuration bitstring, annihilates
!                      an electron in the MO indexed ia
!######################################################################
  function annihilate_electron(conf,ldconf,ia) result(confa)

    use constants
    use bitglobal
    
    implicit none

    ! Input configuration
    integer(is), intent(in) :: ldconf
    integer(ib), intent(in) :: conf(ldconf,2)

    ! Annihilation operator index
    integer(is), intent(in) :: ia

    ! Function result
    integer(ib)             :: confa(ldconf,2)
    
    ! Everything else
    integer(is)             :: k,i

    !
    ! Initialisation
    !
    confa=conf
    
    !
    ! Block index
    !
    k=(ia-1)/64+1

    !
    ! Position of the bit within the k'th block 
    !
    i=ia-(k-1)*64-1

    !
    ! Annihilate the electron
    !
    if (btest(conf(k,2),i)) then
       ! Doubly-occupied MO
       confa(k,2)=ibclr(confa(k,2),i)
    else
       ! Singly-occupied MO
       confa(k,1)=ibclr(confa(k,1),i)
    endif
    
    return
    
  end function annihilate_electron

!######################################################################
! create_electron: given a configuration bitstring, creates an
!                  an electron in the MO indexed ic
!######################################################################
  function create_electron(conf,ldconf,ic) result(confa)

    use constants
    use bitglobal
    
    implicit none

    ! Input configuration
    integer(is), intent(in) :: ldconf
    integer(ib), intent(in) :: conf(ldconf,2)

    ! Creation operator index
    integer(is), intent(in) :: ic

    ! Function result
    integer(ib)             :: confa(ldconf,2)
    
    ! Everything else
    integer(is)             :: k,i

    !
    ! Initialisation
    !
    confa=conf
    
    !
    ! Block index
    !
    k=(ic-1)/64+1

    !
    ! Position of the bit within the k'th block 
    !
    i=ic-(k-1)*64-1

    !
    ! Create the electron
    !
    if (btest(conf(k,1),i)) then
       ! Creation of a doubly-occupied MO
       confa(k,2)=ibset(confa(k,2),i)
    else
       ! Creation of a singly-occupied MO
       confa(k,1)=ibset(confa(k,1),i)
    endif
    
    return
    
  end function create_electron
  
!######################################################################
! n_create_annihilate: Given two configurations, computes the number
!                      of creation and annihilation operators linking
!                      the them
!######################################################################
  function n_create_annihilate(c1,c2,ldc) result(nac)

    use constants
    
    implicit none

    integer(is)              :: nac
    integer(is), intent(in)  :: ldc
    integer(ib), intent(in)  :: c1(ldc,2),c2(ldc,2)
    integer(is)              :: k

    ! Compute the number of creation and annihilation operators
    ! connecting c1 and c2
    nac=popcnt(ieor(c1(1,1),c2(1,1)))+popcnt(ieor(c1(1,2),c2(1,2)))
    do k=2,ldc
       nac=nac+popcnt(ieor(c1(k,1),c2(k,1))) &
            +popcnt(ieor(c1(k,2),c2(k,2)))
    enddo
    
  end function n_create_annihilate
  
!######################################################################
  
end module mrciutils
