!**********************************************************************
! Routines for the calculation of pattern indices for the k=+1
! component of the triplet spin tensor operator
!**********************************************************************
module pattern_indices_kp1

  implicit none

contains

!######################################################################
! pattern_index_kp1: Given an SOP with nopen open shells, creation and
!                    annihilation operator indices, ic and ja, and a
!                    sub-case bitstring encoding, returns a spin
!                    coupling coefficient pattern index.
!######################################################################
  function pattern_index_kp1(sop,ic,ja,nc,na,nopen,icase) &
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

    select case(icase)
    
    case(i1a) ! Sub-case 1a
       pattern=pattern_index_case1a_kp1(sop,ic,ja,nc,na,nopen)
    
    case(i1b) ! Sub-case 1b
       pattern=pattern_index_case1b_kp1(sop,ic,ja,nc,na,nopen)
    
    case(i2a) ! Sub-case 2a
       pattern=pattern_index_case2a_kp1(sop,ic,ja,nc,na,nopen)
    
    case(i2b) ! Sub-case 2b
       pattern=pattern_index_case2b_kp1(sop,ic,ja,nc,na,nopen)
    
    case default ! Unrecognised bitstring encoding
       errmsg='Error in pattern_index_kp1: unrecognised bitstring ' &
            //'encoding'
       call error_control
       
    end select
    
    return

  end function pattern_index_kp1

!######################################################################
! pattern_index_case1a_kp1: Given an SOP with nopen open shells and
!                           creation and annihilation operator indices,
!                           ic and ja, returns a k=+1 Case 1a triplet
!                           spin coupling coefficient pattern index,
!                           i.e., for a singly-occupied -> unoccupied
!                           excitation
!######################################################################
  function pattern_index_case1a_kp1(sop,ic,ja,nc,na,nopen) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
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
       shift=0
    else
       ! ws  ...0...1...
       ! ws' ...1...0...
       ip=ibclr(ip,nc)
       ip=ibclr(ip,na+1)
       shift=offplus(1)
    endif

!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patmap1(ip)+shift
    
    return

  end function pattern_index_case1a_kp1

!######################################################################
! pattern_index_case1b_kp1: Given an SOP with nopen open shells and
!                           creation and annihilation operator indices,
!                           ic and ja, returns a k=+1 Case 1b triplet
!                           spin coupling coefficient pattern index,
!                           i.e., doubly-occupied -> singly-occupied
!                           excitation
!######################################################################
  function pattern_index_case1b_kp1(sop,ic,ja,nc,na,nopen) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
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
       shift=offplus(2)
    else
       ! ws  ...1...2...
       ! ws' ...2...2...
       ip=ibclr(ip,nc)
       ip=ibclr(ip,na)
       shift=offplus(3)
    endif

!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patmap1(ip)+shift
    
    return
    
  end function pattern_index_case1b_kp1

!######################################################################
! pattern_index_case2a_kp1: Given an SOP with nopen open shells and
!                           creation and annihilation operator indices,
!                           ic and ja, returns a k=+1 Case 2a triplet
!                           spin coupling coefficient pattern index,
!                           i.e, doubly-occupied -> unoccupied
!                           excitation
!######################################################################
  function pattern_index_case2a_kp1(sop,ic,ja,nc,na,nopen) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
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
       shift=0
    else
       ! ws  ...0...2...
       ! ws' ...1...1...
       ip=ibclr(ip,nc)
       ip=ibclr(ip,na+1)
       shift=offplus(4)
    endif

!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patmap2a_plus(ip)+shift
    
    return
    
  end function pattern_index_case2a_kp1

!######################################################################
! pattern_index_case2b_kp1: Given an SOP with nopen open shells and
!                           creation and annihilation operator indices,
!                           ic and ja, returns a k=+1 Case 2b triplet
!                           spin coupling coefficient pattern index,
!                           i.e, singly-occupied -> singly-occupied
!                           excitation
!######################################################################
  function pattern_index_case2b_kp1(sop,ic,ja,nc,na,nopen) &
       result(pattern)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is)             :: pattern
    integer(ib), intent(in) :: sop(n_int,2)
    integer(is), intent(in) :: ic,ja,nopen
    integer(is), intent(in) :: nc,na
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

    ! Pattern map shift
    if (ja < ic) then
       shift=0
    else
       shift=offplus(5)
    endif
       
!----------------------------------------------------------------------
! Pattern index
!----------------------------------------------------------------------
    pattern=patmap2b_plus(ip)+shift

    return
    
  end function pattern_index_case2b_kp1
  
!######################################################################
  
end module pattern_indices_kp1
