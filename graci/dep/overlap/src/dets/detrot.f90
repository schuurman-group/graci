!**********************************************************************
! Routines for the rotation of orbitals s.t. the fixed-occupation
! orbitals come first in the alpha- and beta-strings
!**********************************************************************
module detrot

  implicit none

contains

!######################################################################
! rotate_orbitals: top level routine for the application of pi/2
!                  orbital rotations to ensure that the
!                  fixed-occupation spin orbitals appear first in all
!                  strings
!######################################################################
  subroutine rotate_orbitals

    use constants
    use global

    implicit none

    ! Phase factors
    integer(is), allocatable :: phase_alpha(:),phase_beta(:)

    ! Everything else
    integer(is)              :: n,iaX,ibX,id
    
!----------------------------------------------------------------------
! Bra determinants
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(phase_alpha(nalphaB))
    allocate(phase_beta(nbetaB))
    
    ! Rotate the alpha strings and comput the associated phase factors
    call rotate_string(n_intB,nalphaB,alphaB,phase_alpha)

    ! Rotate the beta strings and comput the associated phase factors
    call rotate_string(n_intB,nbetaB,betaB,phase_beta)

    ! Apply the phase factors to the eigenvectors
    !
    ! Loop over roots
    do n=1,nrootsB
       
       ! Loop over alpha strings
       do iaX=1,nalphaB

          ! Loop over the determinants corresponding to this alpha
          ! string
          do id=offsetB(iaX),offsetB(iaX+1)-1

             ! Beta string for this determinant
             ibX=det2betaB(id)
             
             ! Apply the phase factors
             vecB(id,n)=phase_alpha(iaX)*phase_beta(ibx)*vecB(id,n)
             
          enddo
       
       enddo

    enddo
       
    ! Deallocate arrays
    deallocate(phase_alpha)
    deallocate(phase_beta)

!----------------------------------------------------------------------
! Ket determinants
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(phase_alpha(nalphaK))
    allocate(phase_beta(nbetaK))
    
    ! Rotate the alpha strings and comput the associated phase factors
    call rotate_string(n_intK,nalphaK,alphaK,phase_alpha)

    ! Rotate the beta strings and comput the associated phase factors
    call rotate_string(n_intK,nbetaK,betaK,phase_beta)

    ! Apply the phase factors to the eigenvectors
    !
    ! Loop over roots
    do n=1,nrootsK
       
       ! Loop over alpha strings
       do iaX=1,nalphaK

          ! Loop over the determinants corresponding to this alpha
          ! string
          do id=offsetK(iaX),offsetK(iaX+1)-1

             ! Beta string for this determinant
             ibX=det2betaK(id)
             
             ! Apply the phase factors
             vecK(id,n)=phase_alpha(iaX)*phase_beta(ibx)*vecB(id,n)
             
          enddo
       
       enddo

    enddo
       
    ! Deallocate arrays
    deallocate(phase_alpha)
    deallocate(phase_beta)

!----------------------------------------------------------------------
! Rotate the MO overlap matrix
!----------------------------------------------------------------------
    ! We need the lists of orbital pair indices first...
    ! i.e., the plist and hlist arrays used in rotate_string
    
    STOP
    
    return
    
  end subroutine rotate_orbitals
  
!######################################################################
! rotate_string: performs pairwise pi/2 orbital rotations in order s.t.
!                         the fixed-occupation orbitals come first in
!                         all sigma-strings, sigma in {alpha, beta}
!                         also computes and returns the associated
!                         phase factors
!######################################################################
  subroutine rotate_string(n_int,nsigma,sigma,phase)

    use constants
    use detfuncs, only: list_from_bitstring, exc_degree_string, &
                        exc_string
    
    implicit none

    ! Bit strings
    integer(is), intent(in)    :: n_int,nsigma
    integer(ib), intent(inout) :: sigma(n_int,nsigma)

    ! Phase factors
    integer(is), intent(out)   :: phase(nsigma)
    
    ! Orbital rotation information
    integer(is)                :: nrot
    integer(ib)                :: p(n_int),h(n_int),p1(n_int),h1(n_int)
    integer(is), allocatable   :: plist(:),hlist(:),plist1(:),hlist1(:)
    
    ! Everything else
    integer(is)                :: n,j,ih,kh,ip,kp
    integer(ib)                :: string(n_int)
    integer(is)                :: nchange

!----------------------------------------------------------------------
! Get the lists of orbital pairs to be rotated
!----------------------------------------------------------------------
    call get_pairs(n_int,nsigma,sigma,nrot,plist,hlist)

    ! If there are no orbital rotations to perform, then return here
    if (nrot == 0) return
    
!----------------------------------------------------------------------
! Apply the pi/2 pairwise orbital rotations and compute the associated
! phase factors on the fly
!----------------------------------------------------------------------
    ! Allocate arrays
    allocate(plist1(nrot))
    allocate(hlist1(nrot))
    plist1=0
    hlist1=0

    ! Loop over sigma-strings
    do n=1,nsigma

       ! Number of occupation changes
       nchange=0
       
       !
       ! Create the orbital-rotated string
       !
       ! Initialisation
       string=sigma(:,n)

       ! Loop over orbital rotations
       do j=1,nrot

          ! If the hlist(j)'th is not set in the current sigma-string,
          ! then set it and unset the plist(j)'th bit
          kh=(hlist(j)-1)/n_bits+1
          ih=hlist(j)-1-(kh-1)*n_bits
          if (.not. btest(string(kh),ih)) then
             nchange=nchange+1
             string(kh)=ibset(string(kh),ih)
             kp=(plist(j)-1)/n_bits+1
             ip=plist(j)-1-(kp-1)*n_bits
             string(kp)=ibclr(string(kp),ip)
          endif
          
       enddo

       !
       ! Compute the associated phase factor
       !
       if (nchange > 0) then

          ! Get the lists of creation and annihilation
          ! operators
          call exc_string(n_int,sigma(:,n),string,p1,h1)
          call list_from_bitstring(n_int,p1,plist1,nrot)
          call list_from_bitstring(n_int,h1,hlist1,nrot)

          ! Compute the phase factor
          phase(n)=1
          do j=1,nchange

             ! Accumulate the phase factor for this pi/2 rotation
             phase(n)=phase(n)&
                  *phase_single(n_int,sigma(:,n),hlist1(j),plist1(j))

             ! Apply the annihilation operator to the string
             kh=(hlist1(j)-1)/n_bits+1
             ih=hlist1(j)-1-(kh-1)*n_bits
             sigma(kh,n)=ibclr(sigma(kh,n),ih)

             ! Apply the creation operator to the string
             kp=(plist1(j)-1)/n_bits+1
             ip=plist1(j)-1-(kp-1)*n_bits
             sigma(kp,n)=ibset(sigma(kp,n),ip)
             
          enddo

       endif
       
    enddo
    
    return
  
  end subroutine rotate_string

!######################################################################
! get_pairs: for a given set of sigma-strings, sigma in {alpha, beta},
!            determines the set of pairwise orbital rotations needed
!            to change the spin orbital ordering s.t. the fixed-
!            occupation orbitals come first
!######################################################################
  subroutine get_pairs(n_int,nsigma,sigma,nrot,plist,hlist)

    use constants
    use detfuncs, only: list_from_bitstring, exc_degree_string, &
                        exc_string

    implicit none

    ! Bit strings
    integer(is), intent(in)    :: n_int,nsigma
    integer(ib), intent(inout) :: sigma(n_int,nsigma)

    ! Orbital indices
    integer(is), intent(out)   :: nrot
    integer(is), allocatable   :: plist(:),hlist(:)
    
    ! Bit string encoding of the fixed-occupation spin orbitals
    integer(ib)                :: fstring(n_int)

    ! Number of fixed-occupation spin orbitals
    integer(is)                :: nfixed

    ! Orbital rotation information
    integer(ib)                :: p(n_int),h(n_int)
    
    ! Everything else
    integer(is)                :: n,i,k
    integer(ib)                :: string(n_int)
    integer(ib)                :: tstring(n_int)

!----------------------------------------------------------------------
! Construct the bit string encoding of the fixed-occupation spin
! orbitals
!----------------------------------------------------------------------
    ! Initialisation
    fstring=sigma(:,1)

    ! Loop over sigma-strings
    do n=1,nsigma

       ! Loop over blocks
       do k=1,n_int
          fstring(k)=iand(fstring(k),sigma(k,n))
       enddo
          
    enddo

!----------------------------------------------------------------------
! Number of fixed-occupation spin orbitals
!----------------------------------------------------------------------
    nfixed=0
    do k=1,n_int
       nfixed=nfixed+popcnt(fstring(k))
    enddo

!----------------------------------------------------------------------
! Target fixed-occupation orbital ordering: first nfixed bits set
!----------------------------------------------------------------------
    ! Initialisation
    tstring=0_ib

    ! Loop over fixed-occupation orbitals
    do n=1,nfixed

       ! Block index
       k=(n-1)/n_bits+1

       ! Bit index within the block
       i=n-1-(k-1)*n_bits

       ! Set the bit
       tstring(k)=ibset(tstring(k),i)
       
    enddo

!----------------------------------------------------------------------
! Determine the pairwise orbital rotation indices
!----------------------------------------------------------------------
    ! Number of orbital rotations
    nrot=exc_degree_string(n_int,tstring,fstring)

    ! Allocate arrays
    allocate(plist(nrot))
    allocate(hlist(nrot))

    ! If there are no orbital rotations to perform, then return here
    if (nrot == 0) return
    
    ! Get the bit string encodings of the orbital pairs to be
    ! rotated
    call exc_string(n_int,tstring,fstring,p,h)

    ! Get the lists of orbital pair indices
    call list_from_bitstring(n_int,p,plist,nrot)
    call list_from_bitstring(n_int,h,hlist,nrot)
    
    return
    
  end subroutine get_pairs
  
!######################################################################
! phase_single: for a sigma-string s, sigma in {alpha, beta}, and
!               hole and particle indices i and a, returns the phase
!               factor associated with operating on s with the single
!               excitation operator a^dagger j
!######################################################################
  function phase_single(n_int,s,j,a) result(phase)

    use constants
    
    implicit none

    ! Function result
    integer(is)             :: phase
    
    ! Bit string
    integer(is), intent(in) :: n_int
    integer(ib), intent(in) :: s(n_int)

    ! Annihilation and creation operator indices
    integer(is), intent(in) :: j,a

    ! Everything else
    integer(is)             :: high,low,kl,kh,k,l,h
    integer(ib)             :: mask(n_int)
    integer(ib)             :: ones,nperm,par
        
    high=max(j,a)-1
    low=min(j,a)-1

    kl=low/n_bits+1
    kh=high/n_bits+1
    
    l=mod(low,n_bits)
    h=mod(high,n_bits)

    ones=not(0_ib)

    mask=0_ib
    mask(kl:kh-1)=ones
    mask(kh)=ishft(ones,h+1)
    mask(kl)=ieor(mask(kl),ishft(ones,l))

    nperm=0

    do k=kl,kh
       nperm=nperm+popcnt(iand(s(k),mask(k)))
    enddo

    par=iand(nperm,1_ib)

    if (par == 0) then
       phase=1
    else if (par == 1) then
       phase=-1
    else
       write(6,'(/,2x,a,/)') 'Something terrible has happened...'
       stop
    endif
    
    return
    
  end function phase_single
    
!######################################################################
  
end module detrot
