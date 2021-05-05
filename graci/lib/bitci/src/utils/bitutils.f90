module bitutils

  implicit none

contains

!######################################################################
! swap_bits: given a bit string b, swaps the bits at positions p and q
!            *note that we are assuming that the two bits have
!             different values here*
!######################################################################
  subroutine swap_bits(b,p,q)

    use constants
    
    implicit none

    integer(ib), intent(inout) :: b
    integer(is), intent(in)    :: p,q
    integer(ib)                :: t,u

    u=ishft(1_ib,p)
    t=ieor(u,b)
    u=ishft(1_ib,q)
    b=ieor(u,t)
    
    return
    
  end subroutine swap_bits

!######################################################################
! get_permutations: for a bit string composed of n bits set to 1 and m
!                   bits set to 0, computes all possible permutations
!----------------------------------------------------------------------
!                   Uses Anderson's algorithm:
!                   graphics.stanford.edu/~seander/bithacks.html
!----------------------------------------------------------------------
!                   Also see this paper by Anthony Scemama:
!                   https://arxiv.org/abs/1812.06902
!######################################################################
  subroutine get_permutations(n,m,v,vdim)

    use constants
    
    implicit none

    integer(is), intent(in)    :: n,m,vdim
    integer(ib), intent(inout) :: v(vdim)
    integer(is)                :: k
    integer(ib)                :: u,t,tp,tpp,lim

    !
    ! Special case: n=0
    !
    if (n==0) then
       v(1)=0_ib
       return
    endif
    
    !
    ! Initialise the permutation counter
    !
    k=1

    !
    ! Set the first permutation to 00...011...1
    !
    u=ishft(1_ib,n)-1_ib
    
    !
    ! Remaining permutations
    !
    lim=ishft(1_ib,n+m)
    do while(u.lt.lim)
       v(k)=u
       k=k+1
       t=ior(u,u-1_ib)
       tp=t+1_ib
       tpp=ishft((iand(not(t),tp)-1_ib),-(trailz(u)+1_ib))
       u=ior(tp,tpp)
    enddo
    
    return
    
  end subroutine get_permutations
  
!######################################################################
! most_significant_bit: returns the most significant bit in the
!                       64 bit integer i
!######################################################################
  function most_significant_bit(i) result(msb)

    use constants

    implicit none

    integer(ib), intent(in)  :: i
    integer(is)              :: msb
    integer(ib)              :: ii
    
    msb=0
    ii=i
    do while (ii.ne.0_ib)
       msb=msb+1
       ii=ishft(ii,-1)
    enddo

    return
    
  end function most_significant_bit

!######################################################################
! list_from_bitstring: Given a bitstring I, fills in the array list
!                      with the indices of the set bits in I. Note
!                      that I is destroyed on output.
!######################################################################
  subroutine list_from_bitstring(I,list,listdim)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in)    :: listdim
    integer(ib), intent(inout) :: I(n_int)
    integer(is), intent(out)   :: list(listdim)
    integer(is)                :: ispin,k,n
    integer(is)                :: e
    
    ! Initialisation
    list=0
    n=1
    
    ! Loop over blocks
    do k=1,n_int
       
       ! Determine the indices of any set bits
       do while (I(k)/=0_ib)
          e=trailz(I(k))
          I(k)=ibclr(I(k),e)
          list(n)=e+(k-1)*64+1
          n=n+1
       enddo
       
    enddo

    return
    
  end subroutine list_from_bitstring
  
!######################################################################
  
end module bitutils
