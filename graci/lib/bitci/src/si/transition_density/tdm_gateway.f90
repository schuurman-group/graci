!**********************************************************************
! Calculation of 1-TDMs for MRCI wavefunctions
!**********************************************************************

!######################################################################
! transition_density_mrci: Controls the calculation of MRCI 1-TDMs
!######################################################################
#ifdef CBINDING
subroutine transition_density_mrci(irrepB,irrepK,nrootsB,nrootsK,&
     npairs,iroots,rhoij,conffileB_in,vecfileB_in,conffileK_in,&
     vecfileK_in) bind(c,name='density_mrci')
#else
  subroutine transition_density_mrci(irrepB,irrepK,nrootsB,nrootsK,&
     npairs,iroots,rhoij,conffileB_in,vecfileB_in,conffileK_in,&
     vecfileK_in)
#endif

  use constants
  use bitglobal
  use iomod
  use iso_c_binding, only: C_CHAR
  use conftype
    
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
  integer(is), intent(in) :: irrepB,irrepK,nrootsB,nrootsK

  ! Indices of the pairs of states for which 1-TDMs are requested
  integer(is), intent(in) :: npairs
  integer(is), intent(in) :: iroots(npairs,2)
  
  ! Transition density matrix
  real(dp), intent(out)   :: rhoij(nmo,nmo,npairs)

  ! MRCI configuration derived type
  type(mrcfg)             :: cfg

  ! Scratch file numbers
  integer(is)             :: confscrB,vecscrB,confscrK,vecscrK

  ! Eigenpairs
  real(dp), allocatable   :: vecB(:,:),enerB(:),vecK(:,:),enerK(:)
  
  ! Everything else
  integer(is)             :: i

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  !write(6,'(/,52a)') ('-',i=1,52)
  !write(6,'(3(x,a))') 'Transition density matrix calculation for the',&
  !     trim(irreplbl(irrep,ipg)),'subspace'
  !write(6,'(52a)') ('-',i=1,52)
  
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
  
  print*,''
  print*,conffileB
  print*,conffileK
  print*,vecfileB
  print*,vecfileK
  
!----------------------------------------------------------------------
! Register the configuration and eigenvector scratch files
!----------------------------------------------------------------------
  !call register_scratch_file(confscrK,conffileK)
  !call register_scratch_file(vecscrK,vecfileK)

!----------------------------------------------------------------------
! Set up the bra and ket configuration derived types
!----------------------------------------------------------------------
!call cfg%initialise(irrep,confscr)
!
!write(6,'(/,x,a,x,i0)') 'CSF basis dimension:',cfg%csfdim
!
!---------------------------------------------------------------------
!Read in the eigenvectors
!---------------------------------------------------------------------
!! Allocate arrays
!allocate(vec(cfg%csfdim,nroots))
!allocate(ener(nroots))
!
!! Read in the eigenvectors
!call read_some_eigenpairs(vecscr,vec,ener,cfg%csfdim,nroots,iroots)
!
!---------------------------------------------------------------------
!Compute the 1-RDMs
!---------------------------------------------------------------------
!call rdm_mrci(cfg,cfg%csfdim,nroots,vec,dmat)
!
!---------------------------------------------------------------------
!Deallocate arrays
!---------------------------------------------------------------------
!call cfg%finalise
!deallocate(vec)
!deallocate(ener)
  
!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
    
end subroutine transition_density_mrci
