!**********************************************************************
! Routines for the merging the reference spaces of a pair of bra and
! ket configuration derived types
!**********************************************************************
module merge

  implicit none

contains

!######################################################################
! merge_ref_space: given a pair of bra and ket configuration derived
!                  types, re-structures the conf and SOP bit strings
!                  to correspond to the union of the two reference
!                  spaces
!######################################################################
  subroutine merge_ref_space(cfgB,cfgK)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(inout) :: cfgB,cfgK

!----------------------------------------------------------------------
! Do nothing if the two reference spaces are the same
!----------------------------------------------------------------------
    if (sameref(cfgB,cfgK)) return

!----------------------------------------------------------------------
! Get the indices of the MOs in the union of the bra and ket reference
! spaces
!----------------------------------------------------------------------
    errmsg='Different bra and ket ref spaces: '&
         //'merge_ref_space not yet written'
    call error_control
    
    return
    
  end subroutine merge_ref_space
  
!######################################################################

  function sameref(cfgB,cfgK)

    use constants
    use conftype
    
    implicit none

    ! Function result
    logical                 :: sameref

    ! MRCI configuration derived type
    type(mrcfg), intent(in) :: cfgB,cfgK

    ! Everything else
    integer                 :: n,k,i

!----------------------------------------------------------------------
! Different number ref space confs across all irreps <-> different ref
! spaces    
!----------------------------------------------------------------------
    if (cfgB%nR /= cfgK%nR) then
       sameref=.false.
       return
    endif

!----------------------------------------------------------------------
! Different values of n_int_I <-> different ref spaces
!----------------------------------------------------------------------
    if (cfgB%n_int_I /= cfgK%n_int_I) then
       sameref=.false.
       return
    endif
    
!----------------------------------------------------------------------
! Same number of ref space confs across all irreps: check for
! differences
!----------------------------------------------------------------------
    sameref=.true.

    ! Loop over all ref confs
    do n=1,cfgB%nR

       ! Different bra and ref confs?
       do i=1,2
          do k=1,cfgB%n_int_I
             if (cfgB%confR(k,i,n) /= cfgK%confR(k,i,n)) then
                sameref=.false.
                return
             endif
          enddo
       enddo
          
    enddo
    
    return
    
  end function sameref
  
!######################################################################
  
end module merge
