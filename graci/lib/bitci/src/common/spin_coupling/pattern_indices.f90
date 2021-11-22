!**********************************************************************
! Routines for the computation of spin coupling coefficient pattern
! indices from SOPs and creation/annihilation operator indices
!**********************************************************************
module pattern_indices

  implicit none

contains
  
!######################################################################
! get_icase: Given an SOP and a pair of creation and annihilation
!            operators, returns the bitstring encoding of the
!            correspoinding spin-coupling coefficient case type
!######################################################################
  function get_icase(sop,ic,ja) result(indx)

    use constants
    use bitglobal
    
    implicit none

    ! SOP
    integer(ib) :: sop(n_int,2)

    ! Creation and annihilation operator indices
    integer(is) :: ic,ja

    ! Bitstring encoding of the sub-case
    integer(ib) :: indx

    ! Everything else
    integer(is) :: k,i
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    indx=0_ib
    
!----------------------------------------------------------------------
! Annihilation operator
!----------------------------------------------------------------------
    ! Block index
    k=(ja-1)/64+1

    ! Orbital index in the bit string
    i=ja-(k-1)*64-1
    
    ! Set the annihilation operator bits
    if (btest(sop(k,1),i)) then
       ! Singly-occupied MO
       indx=ibset(indx,0)
    else
       ! Doubly-occupied MO
       indx=ibset(indx,0)
       indx=ibset(indx,1)
    endif

!----------------------------------------------------------------------
! Creation operator
!----------------------------------------------------------------------
    ! Block index
    k=(ic-1)/64+1

    ! Orbital index in the bit string
    i=ic-(k-1)*64-1

    ! Set the creation operator bits
    if (btest(sop(k,1),i)) then
       ! Singly-occupied MO
       indx=ibset(indx,2)
    endif

    return
    
  end function get_icase
  
!######################################################################
! pattern_index:  Given an SOP with nopen open shells, creation and
!                 annihilation operator indices, ic and ja, and a
!                 sub-case bitstring encoding, returns a spin
!                 coupling coefficient pattern index.
!######################################################################
  function pattern_index(sop,ic,ja,nc,na,nopen,icase,transpose) &
       result(pattern)

    use constants
    use bitglobal
    use bitstrings
    use iomod

    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nc,na,nopen
    integer(ib), intent(in) :: icase
    logical, intent(in)     :: transpose
    integer(ib)             :: ip
    integer(is)             :: shift

    select case(icase)
    
    case(i1a) ! Sub-case 1a
       pattern=pattern_index_case1a(sop,ic,ja,nc,na,nopen,transpose)
    
    case(i1b) ! Sub-case 1b
       pattern=pattern_index_case1b(sop,ic,ja,nc,na,nopen,transpose)
    
    case(i2a) ! Sub-case 2a
       pattern=pattern_index_case2a(sop,ic,ja,nc,na,nopen,transpose)
    
    case(i2b) ! Sub-case 2b
       pattern=pattern_index_case2b(sop,ic,ja,nc,na,nopen,transpose)
    
    case default ! Unrecognised bitstring encoding
       errmsg='Error in pattern_index: unrecognised bitstring encoding'
       call error_control
       
    end select
    
    return
    
  end function pattern_index
    
!######################################################################
! pattern_index_case1a: Given an SOP with nopen open shells and
!                       creation and annihilation operator indices,
!                       ic and ja, returns a Case 1a spin coupling
!                       coefficient pattern index, i.e., for a
!                       singly-occupied -> unoccupied excitation
!######################################################################
  function pattern_index_case1a(sop,ic,ja,nc,na,nopen,transpose) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
    logical, intent(in)     :: transpose
    integer(ib)             :: ip
    integer(is)             :: shift

!----------------------------------------------------------------------
! Compute the pattern number
!----------------------------------------------------------------------
    ! Bit string with the first nopen+1 bits set
    ip=N1s(nopen+1)

    ! Clear the two bits corresponding to the positions of annihilated
    ! and created MOs in the simplified spatial occupation vector.
    ! This will directly yield the value of the pattern number
    if (ja < ic) then
       ! ws  ...1...0...
       ! ws' ...0...1...
       ip=ibclr(ip,na)
       ip=ibclr(ip,nc)
       if (transpose) then
          shift=offspincp(1)
       else
          shift=0
       endif
    else
       ! ws  ...0...1...
       ! ws' ...1...0...
       ip=ibclr(ip,nc)
       ip=ibclr(ip,na+1)
       if (transpose) then
          shift=0
       else
          shift=offspincp(1)
       endif
    endif

!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patternmap(ip)+shift
    
    return

  end function pattern_index_case1a

!######################################################################
! pattern_index_case1b: Given an SOP with nopen open shells and
!                       creation and annihilation operator indices,
!                       ic and ja, returns a Case 1b spin coupling
!                       coefficient pattern index, i.e., for a
!                       doubly-occupied -> singly-occupied excitation
!######################################################################
  function pattern_index_case1b(sop,ic,ja,nc,na,nopen,transpose) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
    logical, intent(in)     :: transpose
    integer(ib)             :: ip
    integer(is)             :: shift
    
!----------------------------------------------------------------------
! Compute the pattern number
!----------------------------------------------------------------------
    ! Bit string with the first nopen+1 bits set
    ip=N1s(nopen+1)

    ! Clear the two bits corresponding to the positions of annihilated
    ! and created MOs in the simplified spatial occupation vector.
    ! This will directly yield the value of the pattern number
    if (ja < ic) then
       ! ws  ...2...1...
       ! ws' ...1...2...
       ip=ibclr(ip,na)
       ip=ibclr(ip,nc+1)
       if (transpose) then
          shift=offspincp(3)
       else
          shift=offspincp(2)
       endif
    else
       ! ws  ...1...2...
       ! ws' ...2...2...
       ip=ibclr(ip,nc)
       ip=ibclr(ip,na)
       if (transpose) then
          shift=offspincp(2)
       else
          shift=offspincp(3)
       endif
    endif

!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patternmap(ip)+shift
    
    return
    
  end function pattern_index_case1b

!######################################################################
! pattern_index_case2a: Given an SOP with nopen open shells and
!                       creation and annihilation operator indices,
!                       ic and ja, returns a Case 2a spin coupling
!                       coefficient pattern index, i.e., for a
!                       doubly-occupied -> unoccupied excitation
!######################################################################
  function pattern_index_case2a(sop,ic,ja,nc,na,nopen,transpose) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
    logical, intent(in)     :: transpose
    integer(ib)             :: ip
    integer(is)             :: shift

!----------------------------------------------------------------------
! Compute the pattern number
!----------------------------------------------------------------------
    ! Bit string with the first nopen+2 bits set
    ip=N1s(nopen+2)
    
    ! Clear the two bits corresponding to the positions of annihilated
    ! and created MOs in the simplified spatial occupation vector.
    ! This will directly yield the value of the pattern number
    if (ja < ic) then
       ! ws  ...2...0...
       ! ws' ...1...1...
       ip=ibclr(ip,na)
       ip=ibclr(ip,nc+1)
    else
       ! ws  ...0...2...
       ! ws' ...1...1...
       ip=ibclr(ip,nc)
       ip=ibclr(ip,na+1)
    endif

    ! Transposition shift
    if (transpose) then
       shift=offspincp(4)
    else
       shift=0
    endif
        
!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patternmap(ip)+shift
    
    return
    
  end function pattern_index_case2a

!######################################################################
! pattern_index_case2b: Given an SOP with nopen open shells and
!                       creation and annihilation operator indices,
!                       ic and ja, returns a Case 2b spin coupling
!                       coefficient pattern index, i.e., for a
!                       singly-occupied -> singly-occupied excitation
!######################################################################
  function pattern_index_case2b(sop,ic,ja,nc,na,nopen,transpose) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
    logical, intent(in)     :: transpose
    integer(ib)             :: ip
    integer(is)             :: shift

!----------------------------------------------------------------------
! Compute the pattern number
!----------------------------------------------------------------------
    ! Bit string with the first nopen bits set
    ip=N1s(nopen)

    ! Clear the two bits corresponding to the positions of annihilated
    ! and created MOs in the simplified spatial occupation vector.
    ! This will directly yield the value of the pattern number
    !
    ! ws  ...1...1... or ...1...1... (this makes no difference because
    ! ws' ...0...2...    ...2...0... both MOs are singly-occupied in ws)
    ip=ibclr(ip,na)
    ip=ibclr(ip,nc)

    ! Transposition shift
    if (transpose) then
       shift=0
    else
       shift=offspincp(4)
    endif
    
!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patternmap(ip)+shift

    return
    
  end function pattern_index_case2b
  
!######################################################################
  
end module pattern_indices
