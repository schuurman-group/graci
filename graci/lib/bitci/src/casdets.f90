!######################################################################
! generate_cas_dets: generates the set of determinants corresponding
!                    to a CAS(n,m) active space and the base
!                    determinant
!######################################################################
#ifdef CBINDING
subroutine generate_cas_dets(icas,n,m,ndet) &
     bind(c,name="generate_cas_dets")
#else
subroutine generate_cas_dets(icas,n,m,ndet)
#endif
  
  use constants
  use bitglobal
  use bitutils
  use utils
  use math
  use iomod

  integer(is), intent(inout) :: icas(nmo)
  integer(is), intent(in)    :: n,m
  integer(is), intent(out)   :: ndet
  integer(is)                :: i,i1,imo,k,ic,ialpha,ibeta
  integer(is)                :: indx(m),aux(m)
  integer(is)                :: na,nb,ndeta,ndetb
  integer(ib), allocatable   :: detcas(:,:,:)
  integer(is)                :: vdim
  integer(ib), allocatable   :: v(:,:)
  integer(is)                :: iscratch
  character(len=60)          :: casfile
  
!----------------------------------------------------------------------
! Perform some sanity checks
!----------------------------------------------------------------------
  ! Exit if the base determinant has not been generated
  if (.not.allocated(det0)) then
     errmsg='Error in generate_cas_excitations: the base determinant'&
          //' has not yet been generated'
     call error_control
  endif

  ! Exit if the number of CAS orbitals is greater than 64 (this is both
  ! an insane number and would over-fill the permutation bit string)
  if (m.gt.64) then
     errmsg='Error in generate_cas_excitations: the no. CAS orbitals'&
          //' is greater than 64'
     call error_control
  endif

!----------------------------------------------------------------------
! Get the number of alpha and beta electrons in the active space
! orbitals
!----------------------------------------------------------------------
  na=0
  nb=0
  
  ! Loop over active space orbitals
  do i1=1,m

     ! Current active space orbital
     imo=icas(i1)

     ! Block index
     k=(imo-1)/64+1

     ! Orbital index in the bit string
     i=imo-(k-1)*64-1

     ! Cumulative numbers of alpha and beta electrons
     if (btest(det0(k,1),i)) na=na+1
     if (btest(det0(k,2),i)) nb=nb+1
     
  enddo

!----------------------------------------------------------------------
! One last sanity check: exit if the user-specified number of active
! space electrons doesn't match what is in the icas array
!----------------------------------------------------------------------
  if (na+nb/=n) then
     errmsg='Error in generate_cas_excitations: the number of CAS'&
          //' electrons is not consistent with the base orbital'&
          //' occupations'
     call error_control
  endif

!----------------------------------------------------------------------
! Sort the active space orbital indices
!----------------------------------------------------------------------
  call i4sortindx('A',m,icas(1:m),indx)
  do i=1,m
     aux(i)=icas(indx(i))
  enddo
  icas(1:m)=aux(1:m)

!----------------------------------------------------------------------
! Determine the number of CAS determinants and allocate the detcas
! array
!----------------------------------------------------------------------
  ! No. determinants resulting from alpha-spin excitations
  ndeta=bico(m,na)

  ! No. determinants resulting from beta-spin excitations
  ndetb=bico(m,nb)

  ! Total number of determinants
  ndet=ndeta*ndetb

  ! Allocate the CAS determinant array
  allocate(detcas(n_int,2,ndet))

!----------------------------------------------------------------------
! Compute the permutations of the alpha and beta electrons in the CAS
! orbitals
!----------------------------------------------------------------------
  ! Allocate the electron permutation array
  vdim=max(ndeta,ndetb)
  allocate(v(vdim,2))
  v=0_ib
  
  ! Alpha-spin excitations
  call get_permutations(na,m-na,v(1:ndeta,1),ndeta)

  ! Beta-spin excitations
  call get_permutations(nb,m-nb,v(1:ndetb,2),ndetb)

!----------------------------------------------------------------------
! Generate the CAS determinants
!----------------------------------------------------------------------
  ! Initialise all the CAS determinant to the base determinant
  ! We will subsequently alter the occupations of the CAS orbitals
  do i=1,ndet
     detcas(:,:,i)=det0
  enddo

  ! Determinant counter
  ic=0

  ! Loop over alpha strings
  do ialpha=1,ndeta

     ! Loop over beta strings
     do ibeta=1,ndetb
        
        ! Increment the determinant counter
        ic=ic+1
        
        ! Construct the current alpha string
        !
        ! Loop over CAS orbitals
        do i1=1,m
           ! Block index
           k=(icas(i1)-1)/64+1
           ! Index of the orbital in the bit string block
           i=icas(i1)-1-(k-1)*64
           if (btest(v(ialpha,1),i1-1)) then
              ! Bit set to 1
              detcas(k,1,ic)=ibset(detcas(k,1,ic),i)
           else
              ! Bit set to 0
              detcas(k,1,ic)=ibclr(detcas(k,1,ic),i)
           endif
        enddo
        
        ! Construct the current beta string
        !
        ! Loop over CAS orbitals
        do i1=1,m
           ! Block index
           k=(icas(i1)-1)/64+1
           ! Index of the orbital in the bit string block
           i=icas(i1)-1-(k-1)*64
           if (btest(v(ibeta,2),i1-1)) then
              ! Bit set to 1
              detcas(k,2,ic)=ibset(detcas(k,2,ic),i)
           else
              ! Bit set to 0
              detcas(k,2,ic)=ibclr(detcas(k,2,ic),i)
           endif
        enddo

     enddo
     
  enddo

!----------------------------------------------------------------------
! Write the CAS determinants to disk
!----------------------------------------------------------------------
  ! Open the scratch file
  call freeunit(iscratch)
  call scratch_name('detcas',casfile)
  open(iscratch,file=casfile,form='unformatted',status='unknown')

  ! Number of CAS determinants
  write(iscratch) ndet

  ! CAS determinants
  write(iscratch) detcas
  
  ! Close the scratch file
  close(iscratch)
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(detcas,v)
  
  return
  
end subroutine generate_cas_dets
