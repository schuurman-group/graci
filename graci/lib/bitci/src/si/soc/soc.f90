!**********************************************************************
! Calculation of the spin-part of SOC integrals for MRCI wavefunctions
!**********************************************************************
! That is, the matrix elements < psi_I | T_ij^(1,k) | psi'_J >
! for the k=0,+1 components of the triplet excitation operators
! T_ij^(1,k)
!**********************************************************************
#ifdef CBINDING
subroutine soc_mrci(irrepB,irrepK,nrootsB,nrootsK,npairs,iroots,&
     Tij,conffileB_in,vecfileB_in,conffileK_in,vecfileK_in) &
     bind(c,name='soc_mrci')
#else
subroutine soc_mrci((irrepB,irrepK,nrootsB,nrootsK,npairs,iroots,&
     Tij,conffileB_in,vecfileB_in,conffileK_in,vecfileK_in)
#endif

  use iso_c_binding, only: C_CHAR
  use constants
  use bitglobal
  use iomod
  use conftype
  use merge
  use tdm

  implicit none

  ! MRCI configuration and eigenvector file names
#ifdef CBINDING
  character(kind=C_CHAR), intent(in) :: conffileB_in(*),vecfileB_in(*),&
                                        conffileK_in(*),vecfileK_in(*)
  character(len=255)                 :: conffileB,vecfileB,&
                                        conffileK,vecfileK
  integer(is)                        :: length
#else
  character(len=*), intent(in)       :: conffileB_in,vecfileB_in,&
                                        conffileK_in,vecfileK_in
  character(len=255)                 :: conffileB,vecfileB,&
                                        conffileK,vecfileK
#endif

  ! Irreps and no. roots
  integer(is), intent(in)  :: irrepB,irrepK,nrootsB,nrootsK

  ! Indices of the pairs of states for which 1-TDMs are requested
  integer(is), intent(in)  :: npairs
  integer(is), intent(in)  :: iroots(npairs,2)
  
  ! Triplet transition density matrix
  real(dp), intent(out)    :: Tij(nmo,nmo,npairs)

  ! MRCI configuration derived types
  type(mrcfg)              :: cfgB,cfgK

  ! Scratch file numbers
  integer(is)              :: confscrB,vecscrB,confscrK,vecscrK

  ! Eigenpairs
  real(dp), allocatable    :: vecB(:,:),enerB(:),vecK(:,:),enerK(:)
  
  ! Everything else
  integer(is)              :: i,j,k
  integer(is)              :: nvecB,nvecK
  integer(is), allocatable :: iBra(:),iKet(:)
  integer(is), allocatable :: Bmap(:),Kmap(:)
  integer(is), allocatable :: ireadB(:),ireadK(:)

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,52a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'Triplet transition density matrix calculation'
  write(6,'(52a)') ('-',i=1,52)

!----------------------------------------------------------------------
! If C bindings are on, then convert the MRCI configuration and
! eigenvector file names from the C char type to the Fortran character
! type
!----------------------------------------------------------------------
#ifdef CBINDING
  length=cstrlen(conffileB_in)
  call c2fstr(conffileB_in,conffileB,length)
  length=cstrlen(vecfileB_in)
  call c2fstr(vecfileB_in,vecfileB,length)
  length=cstrlen(conffileK_in)
  call c2fstr(conffileK_in,conffileK,length)
  length=cstrlen(vecfileK_in)
  call c2fstr(vecfileK_in,vecfileK,length)
#else
  conffileB=adjustl(trim(conffileB_in))
  vecfileB=adjustl(trim(conffileB_in))
  conffileK=adjustl(trim(conffileK_in))
  vecfileK=adjustl(trim(conffileK_in))
#endif

!----------------------------------------------------------------------
! Register the configuration and eigenvector scratch files
!----------------------------------------------------------------------
  call register_scratch_file(confscrB,conffileB)
  call register_scratch_file(vecscrB,vecfileB)

  call register_scratch_file(confscrK,conffileK)
  call register_scratch_file(vecscrK,vecfileK)
  
!----------------------------------------------------------------------
! Set up the bra and ket configuration derived types
!----------------------------------------------------------------------
  call cfgB%initialise(irrepB,confscrB)
  call cfgK%initialise(irrepK,confscrK)
  
  write(6,'(/,x,a,x,i0)') 'Bra CSF basis dimension:',cfgB%csfdim
  write(6,'(x,a,x,i0)') 'Ket CSF basis dimension:',cfgK%csfdim

  stop
  
!----------------------------------------------------------------------
! Merge the bra and ket reference spaces
!----------------------------------------------------------------------
  !call merge_ref_space(cfgB,cfgK,ncsfsB,ncsfsK)
  
  return
  
end subroutine soc_mrci
