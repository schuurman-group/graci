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

  print*,''
  print*,'here'
  stop
  
  return
  
end subroutine soc_mrci
