!**********************************************************************
! Routines for the extraction of the Slater determinant representation
! of bitci MRCI-type wave functions
!**********************************************************************

!######################################################################
! detwf: top level routine for the conversion from the CSF to
!        determinant basis of all states of a given irrep, either
!        for the bra or ket manifold (bkstr = 'bra' or 'ket')
!######################################################################
#ifdef CBINDING
subroutine detwf(conffile_in,vecfile_in,nroots,bkstr_in,wfscr) &
     bind(c,name='detwf')
#else
subroutine detwf(conffile_in,vecfile_in,nroots,bkstr_in,wfscr)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use iomod
  use conftype

  ! bitCI configuration and eigenvector file names
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: conffile_in(*),vecfile_in(*),&
                                        bkstr_in
  character(len=255)                 :: conffile,vecfile,bkstr
       
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: conffile_in,vecfile_in,bkstr_in
  character(len=255)                 :: conffile,vecfile,bkstr
#endif

  ! No. roots
  integer(is), intent(in) :: nroots

  ! Determinant wave function scratch file number
  integer(is), intent(out) :: wfscr

  ! Everything else
  
!----------------------------------------------------------------------
! If C bindings are on, then convert the bitCI configuration and
! eigenvector file names, and bra/ket string from the C char type
! to the Fortran character type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(conffile_in)
  call c2fstr(conffile_in,conffile,length)
  length=cstrlen(vecfile_in)
  call c2fstr(vecfile_in,vecfile,length)
  length=cstrlen(bkstr_in)
  call c2fstr(bkstr_in,bkstr,length)
#else
  conffile=adjustl(trim(conffile_in))
  vecfileB=adjustl(trim(vecfileB_in))
  bkstr=adjustl(trim(bkstr_in))
#endif
  
  return
  
end subroutine detwf
