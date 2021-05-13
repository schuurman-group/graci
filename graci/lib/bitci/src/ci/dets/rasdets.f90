!**********************************************************************
! Routines for the generation of bit string representations of the
! Slater determinants of a RAS wavefunction
!**********************************************************************

!######################################################################
! rasperm: arrays holding the information about the possible
!          permutations of the electrons in the orbitals in the 
!          RAS1, RAS2 and RAS3 orbitals
!######################################################################
module rasperm

  use constants
  
  implicit none

  save

  ! Number of alpha and beta electrons in the RAS1, RAS2 and RAS3
  ! orbitals
  integer(is)              :: na_ras1,nb_ras1
  integer(is)              :: na_ras2,nb_ras2
  integer(is)              :: na_ras3,nb_ras3

  ! Number of permutations of the alpha and beta electrons in the
  ! RAS1, RAS2 and RAS3 orbitals
  integer(is)              :: nperma_ras1(-2:0),npermb_ras1(-2:0)
  integer(is)              :: nperma_ras2(-2:2),npermb_ras2(-2:2)
  integer(is)              :: nperma_ras3(0:2),npermb_ras3(0:2)

  ! Encoding of the permutations in 64-bit integers
  integer(is)              :: v_ras1_dim,v_ras2_dim,v_ras3_dim
  integer(ib), allocatable :: v_ras1(:,:,:),v_ras2(:,:,:),v_ras3(:,:,:)

  ! Logical flags indicating whether there are a non-zero number of
  ! RAS1, RAS2, or RAS3 orbitals
  logical                  :: okras1,okras2,okras3
  
end module rasperm
  
!######################################################################
! generate_ras_dets: generates the set of determinants corresponding
!                    to a RAS active space and the base determinant
!----------------------------------------------------------------------
! Input variables:   -iras1: indices of the RAS1 orbitals
!                    -iras2: indices of the RAS2 orbitals
!                    -iras3: indices of the RAS3 orbitals
!                    -nras1: maximum number of holes in the RAS1
!                            space
!                    -mras1: number of RAS1 orbitals
!                    -nras2: number of electrons in the RAS2 space
!                    -mras2: number of orbitals in the RAS2 space
!                    -nras3: maximum number of electrons in the RAS3
!                            space
!                    -mras3: number of orbitals in the RAS3 space
!                    -sort:  True if the determinants are to be double
!                            sorted into alpha- and beta-major order
!----------------------------------------------------------------------
! Output variables:  -ndet:    total number of RAS determinants
!                    -scrnum:  the scratch file number for the RAS
!                              determinants
!######################################################################
subroutine generate_ras_dets(iras1,iras2,iras3,nras1,mras1,nras2,&
     mras2,nras3,mras3,ndet,scrnum,sort)
  
  use constants
  use bitglobal
  use rasperm
  use detsort
  
  implicit none

  integer(is), intent(inout) :: iras1(nmo),iras2(nmo),iras3(nmo)
  integer(is), intent(in)    :: nras1,mras1,nras2,mras2,nras3,mras3
  integer(is), intent(out)   :: ndet
  integer(ib), allocatable   :: detras_a(:,:,:),detras_b(:,:,:)
  integer(is), intent(out)   :: scrnum
  logical, intent(in)        :: sort
  
  integer(is)                :: offdim_a,offdim_b
  integer(is), allocatable   :: mapab(:)
  integer(is), allocatable   :: nsym(:)
  integer(is), allocatable   :: offset_a(:,:),offset_b(:,:)
  integer(is), allocatable   :: nunique_a(:),nunique_b(:)

!----------------------------------------------------------------------
! Sanity checks
!----------------------------------------------------------------------
  call check_ras_input(iras1,iras2,iras3,nras1,mras1,nras2,mras2,&
       nras3,mras3)

!----------------------------------------------------------------------
! Generate all the information needed to construct the RAS determinants
!----------------------------------------------------------------------
  call get_common_ras_data(iras1,iras2,iras3,nras1,mras1,nras2,mras2,&
       nras3,mras3)

!----------------------------------------------------------------------
! Determine the total number of RAS determinants using a dumb loop
! over all allowed excitation classes
!----------------------------------------------------------------------
  call get_ndet_ras_dumb(ndet,nras1,nras2,nras3)
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(detras_a(n_int,2,ndet))
  allocate(detras_b(n_int,2,ndet))
  allocate(mapab(ndet))
  allocate(nsym(0:nirrep-1))
  allocate(nunique_a(0:nirrep-1))
  allocate(nunique_b(0:nirrep-1))
  
!----------------------------------------------------------------------
! Generate all the determinants of the RAS wavefunction
!----------------------------------------------------------------------
  call make_ras_dets(detras_a,ndet,iras1,iras2,iras3,nras1,mras1,nras2,&
       mras2,nras3,mras3)

!----------------------------------------------------------------------
! Sort the RAS determinants
!----------------------------------------------------------------------
  if (sort) call det_sort_all(detras_a,detras_b,ndet,nsym,mapab,&
       offset_a,offset_b,offdim_a,offdim_b,nunique_a,nunique_b)
  
!----------------------------------------------------------------------
! Write the sorted determinant arrays to disk
!----------------------------------------------------------------------
  if (sort) then
     call write_ras_dets_sorted(scrnum,detras_a,detras_b,ndet,nsym,&
          mapab,offset_a,offset_b,offdim_a,offdim_b,nunique_a,nunique_b)
  else
     call write_ras_dets_unsorted(scrnum,detras_a,ndet)
  endif
     
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  if (allocated(v_ras1)) deallocate(v_ras1)
  if (allocated(v_ras2)) deallocate(v_ras2)
  if (allocated(v_ras3)) deallocate(v_ras3)
  if (allocated(detras_a)) deallocate(detras_a)
  if (allocated(detras_b)) deallocate(detras_b)
  if (allocated(mapab)) deallocate(mapab)
  if (allocated(nsym)) deallocate(nsym)
  if (allocated(nunique_a)) deallocate(nunique_a)
  if (allocated(nunique_b)) deallocate(nunique_b)
  if (allocated(offset_a)) deallocate(offset_a)
  if (allocated(offset_b)) deallocate(offset_b)

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine generate_ras_dets       

!######################################################################
! check_ras_input: sanity checks on the user-specifed RAS parameters
!######################################################################
subroutine check_ras_input(iras1,iras2,iras3,nras1,mras1,nras2,mras2,&
     nras3,mras3)

  use constants
  use bitglobal
  use rasperm
  use iomod

  implicit none

  integer(is), intent(inout) :: iras1(nmo),iras2(nmo),iras3(nmo)
  integer(is), intent(in)    :: nras1,mras1,nras2,mras2,nras3,mras3
  integer(is)                :: i,i1,imo,k,n,na,nb

!----------------------------------------------------------------------
! Set the logical flags indicating whether there are a non-zero number
! of RAS1, RAS2, or RAS3 orbitals
!----------------------------------------------------------------------
  ! RAS1
  if (mras1>0) then
     okras1=.true.
  else
     okras1=.false.
  endif

  ! RAS2
  if (mras2>0) then
     okras2=.true.
  else
     okras2=.false.
  endif

  ! RAS3
  if (mras3>0) then
     okras3=.true.
  else
     okras3=.false.
  endif
  
!----------------------------------------------------------------------
! Exit if the base determinant has not been generated
!----------------------------------------------------------------------
  if (.not.allocated(det0)) then
     errmsg='Error in generate_cas_excitations: the base determinant'&
          //' has not yet been generated'
     call error_control
  endif

!----------------------------------------------------------------------  
! Exit if the number of CAS orbitals is greater than 64 (this is both
! an insane number and would over-fill the permutation bit string)
!----------------------------------------------------------------------
  if (mras2.gt.64) then
     errmsg='Error in generate_cas_excitations: the no. CAS orbitals'&
          //' is greater than 64'
     call error_control
  endif

!----------------------------------------------------------------------
! Exit if the user-specified number of active space electrons doesn't
! match what is in the icas array
!----------------------------------------------------------------------
  na=0
  nb=0
  
  ! Loop over RAS2 orbitals
  do i1=1,mras2

     ! Current active space orbital
     imo=iras2(i1)

     ! Block index
     k=(imo-1)/64+1

     ! Orbital index in the bit string
     i=imo-(k-1)*64-1

     ! Cumulative numbers of alpha and beta electrons
     if (btest(det0(k,1),i)) na=na+1
     if (btest(det0(k,2),i)) nb=nb+1
     
  enddo

  if (na+nb/=nras2) then
     errmsg='Error in check_ras_input: the number of CAS'&
          //' electrons is not consistent with the base orbital'&
          //' occupations'
     call error_control
  endif
  
!----------------------------------------------------------------------
! The maximum number of holes in RAS1 cannot be greater than 2
!----------------------------------------------------------------------
  if (okras1.and.nras1<1.or.nras1>2) then
     errmsg='Error in check_ras_input: the maximum no. RAS1 holes'&
          //' can only be 1 or 2'
     call error_control
  endif
  
!----------------------------------------------------------------------
! The maximum number of electrons in RAS3 cannot be greater than 2
!----------------------------------------------------------------------
  if (okras3.and.nras3<1.or.nras3>2) then
     errmsg='Error in check_ras_input: the maximum no. RAS3'&
          //' electrons can only be 1 or 2'
     call error_control
  endif

!----------------------------------------------------------------------
! Check the RAS1 orbitals: they have to be doubly-occupied in the
! current incarnation of the code
!----------------------------------------------------------------------
  n=0
  ! Loop over RAS1 orbitals
  do i1=1,mras1

     ! Current RAS1 orbital
     imo=iras1(i1)

     ! Block index
     k=(imo-1)/64+1

     ! Orbital index in the bit string
     i=imo-(k-1)*64-1

     ! Sum of alpha and beta occupations
     if (btest(det0(k,1),i)) n=n+1
     if (btest(det0(k,2),i)) n=n+1
     
  enddo

  if (n/=2*mras1) then
     errmsg='Error in check_ras_input: the RAS1 orbitals have to'&
          //' be doubly-occupied'
     call error_control
  endif

!----------------------------------------------------------------------
! Check the RAS3 orbitals: they have to be unoccupied in the current
! incarnation of the code
!----------------------------------------------------------------------
  n=0
  ! Loop over RAS1 orbitals
  do i1=1,mras3

     ! Current RAS1 orbital
     imo=iras3(i1)

     ! Block index
     k=(imo-1)/64+1

     ! Orbital index in the bit string
     i=imo-(k-1)*64-1

     ! Sum of alpha and beta occupations
     if (btest(det0(k,1),i)) n=n+1
     if (btest(det0(k,2),i)) n=n+1
     
  enddo

  if (n/=0) then
     errmsg='Error in check_ras_input: the RAS3 orbitals have to'&
          //' be unoccupied'
     call error_control
  endif
  
  return
  
end subroutine check_ras_input

!######################################################################
! get_common_ras_data: constructs all the permutation arrays needed
!                      to generate the RAS determinants
!######################################################################
subroutine get_common_ras_data(iras1,iras2,iras3,nras1,mras1,nras2,&
     mras2,nras3,mras3)

  use constants
  use bitglobal
  use bitutils
  use rasperm
  use utils
  use math
  use iomod

  implicit none

  integer(is), intent(inout) :: iras1(nmo),iras2(nmo),iras3(nmo)
  integer(is), intent(in)    :: nras1,mras1,nras2,mras2,nras3,mras3
  integer(is)                :: indx(nmo),aux(nmo)
  integer(is)                :: i,i1,imo,k,ic
  
!----------------------------------------------------------------------
! Sort the iras1, iras2 and iras3 arrays
!----------------------------------------------------------------------
  ! RAS1 orbitals
  if (okras1) then
     call i4sortindx('A',mras1,iras1(1:mras1),indx(1:mras1))
     do i=1,mras1
        aux(i)=iras1(indx(i))
     enddo
     iras1(1:mras1)=aux(1:mras1)
  endif
  
  ! RAS2 orbitals
  if (okras2) then
     call i4sortindx('A',mras2,iras2(1:mras2),indx(1:mras2))
     do i=1,mras2
        aux(i)=iras2(indx(i))
     enddo
     iras2(1:mras2)=aux(1:mras2)
  endif
     
  ! RAS3 orbitals
  if (okras3) then
     call i4sortindx('A',mras3,iras3(1:mras3),indx(1:mras3))
     do i=1,mras3
        aux(i)=iras3(indx(i))
     enddo
     iras3(1:mras3)=aux(1:mras3)
  endif

!----------------------------------------------------------------------
! Get the number of alpha and beta electrons in the RAS2 orbitals
!----------------------------------------------------------------------
  na_ras2=0
  nb_ras2=0
  
  ! Loop over RAS2 orbitals
  do i1=1,mras2
     
     ! Current active space orbital
     imo=iras2(i1)

     ! Block index
     k=(imo-1)/64+1

     ! Orbital index in the bit string
     i=imo-(k-1)*64-1

     ! Cumulative numbers of alpha and beta electrons
     if (btest(det0(k,1),i)) na_ras2=na_ras2+1
     if (btest(det0(k,2),i)) nb_ras2=nb_ras2+1
     
  enddo
  
!----------------------------------------------------------------------
! Determine the number of permumations of electrons in the RAS2
! orbitals
!----------------------------------------------------------------------
  select case(okras2)
     
  case(.true.) ! The number of RAS2 orbitals is non-zero

     ! Permutations of the alpha electrons in RAS2
     nperma_ras2(0)=bico(mras2,na_ras2)

     ! Permutations of the beta electrons in RAS2
     npermb_ras2(0)=bico(mras2,nb_ras2)

     ! Permutations of the alpha electrons in RAS2 with +1 alpha electrons
     if (na_ras2+1<=mras2) then
        nperma_ras2(1)=bico(mras2,na_ras2+1)
     else
        nperma_ras2(1)=0
     endif
        
     ! Permutations of the beta electrons in RAS2 with +1 beta electrons
     if (nb_ras2+1<=mras2) then
        npermb_ras2(1)=bico(mras2,nb_ras2+1)
     else
        npermb_ras2(1)=0
     endif
        
     ! Permutations of the alpha electrons in RAS2 with +2 alpha electrons
     if (na_ras2+2<=mras2) then
        nperma_ras2(2)=bico(mras2,na_ras2+2)
     else
        nperma_ras2(2)=0
     endif
        
     ! Permutations of the beta electrons in RAS2 with +2 beta electrons
     if (nb_ras2+2<=mras2) then
        npermb_ras2(2)=bico(mras2,nb_ras2+2)
     else
        npermb_ras2(2)=0
     endif
        
     ! Permutations of the alpha electrons in RAS2 with -1 alpha electrons
     if (na_ras2>0) then
        nperma_ras2(-1)=bico(mras2,na_ras2-1)
     else
        nperma_ras2(-1)=0
     endif

     ! Permutations of the beta electrons in RAS2 with -1 beta electrons
     if (nb_ras2>0) then
        npermb_ras2(-1)=bico(mras2,nb_ras2-1)
     else
        npermb_ras2(-1)=0
     endif

     ! Permutations of the alpha electrons in RAS2 with -2 alpha electrons
     if (na_ras2>1) then
        nperma_ras2(-2)=bico(mras2,na_ras2-2)
     else
        nperma_ras2(-2)=0
     endif

     ! Permutations of the beta electrons in RAS2 with -2 beta electrons
     if (nb_ras2>1) then
        npermb_ras2(-2)=bico(mras2,nb_ras2-2)
     else
        npermb_ras2(-2)=0
     endif

     ! Allocate the RAS2 electron permutation array
     v_ras2_dim=max(maxval(nperma_ras2),maxval(npermb_ras2))
     allocate(v_ras2(v_ras2_dim,-2:2,2))
     v_ras2=0_ib
          
  case(.false.) ! The number of RAS2 orbitals is zero

     ! We will only have RAS1->RAS3 excitations in the wavefunction.
     ! We need to set things up so that the base determinant
     ! orbital occupations are used in the RAS2 space
     
     ! Set nperma_ras2(0)=npermb_ras2(0)=1
     ! Set all other elements to 0
     nperma_ras2=0
     npermb_ras2=0
     nperma_ras2(0)=1
     npermb_ras2(0)=1

  end select

!----------------------------------------------------------------------
! Compute the permutations of the alpha and beta electrons in the RAS2
! orbitals
!----------------------------------------------------------------------
  if (okras2) then

     ! Permutations of the alpha electrons in RAS2
     call get_permutations(na_ras2,mras2-na_ras2,&
          v_ras2(1:nperma_ras2(0),0,1),nperma_ras2(0))
     
     ! Permutations of the beta electrons in RAS2
     call get_permutations(nb_ras2,mras2-nb_ras2,&
          v_ras2(1:npermb_ras2(0),0,2),npermb_ras2(0))

     ! Permutations of the alpha electrons in RAS2 with +1 alpha electrons
     call get_permutations(na_ras2+1,mras2-(na_ras2+1),&
          v_ras2(1:nperma_ras2(1),1,1),nperma_ras2(1))

     ! Permutations of the beta electrons in RAS2 with +1 beta electrons
     call get_permutations(nb_ras2+1,mras2-(nb_ras2+1),&
          v_ras2(1:npermb_ras2(1),1,2),npermb_ras2(1))

     ! Permutations of the alpha electrons in RAS2 with +2 alpha electrons
     call get_permutations(na_ras2+2,mras2-(na_ras2+2),&
          v_ras2(1:nperma_ras2(2),2,1),nperma_ras2(2))

     ! Permutations of the beta electrons in RAS2 with +2 beta electrons
     call get_permutations(nb_ras2+2,mras2-(nb_ras2+2),&
          v_ras2(1:npermb_ras2(2),2,2),nperma_ras2(2))

     ! Permutations of the alpha electrons in RAS2 with -1 alpha electrons
     if (na_ras2>0) then
        call get_permutations(na_ras2-1,mras2-(na_ras2-1),&
             v_ras2(1:nperma_ras2(-1),-1,1),nperma_ras2(-1))
     endif

     ! Permutations of the beta electrons in RAS2 with -1 beta electrons
     if (nb_ras2>0) then
        call get_permutations(nb_ras2-1,mras2-(nb_ras2-1),&
             v_ras2(1:npermb_ras2(-1),-1,2),npermb_ras2(-1))
     endif

     ! Permutations of the alpha electrons in RAS2 with -2 alpha electrons
     if (na_ras2>1) then
        call get_permutations(na_ras2-2,mras2-(na_ras2-2),&
             v_ras2(1:nperma_ras2(-2),-2,1),nperma_ras2(-2))
     endif

     ! Permutations of the beta electrons in RAS2 with -2 beta electrons
     if (nb_ras2>1) then
        call get_permutations(nb_ras2-2,mras2-(nb_ras2-2),&
             v_ras2(1:npermb_ras2(-2),-2,2),npermb_ras2(-2))
     endif

  endif
  
!----------------------------------------------------------------------
! Determine the number of permumations of electrons in the RAS1
! orbitals
!----------------------------------------------------------------------
  select case(okras1)

  case(.true.) ! The number of RAS1 orbitals is non-zero
     ! Number of alpha and beta electrons in RAS1
     ! Note that the RAS1 orbitals are constrained to be doubly-occupied
     na_ras1=mras1
     nb_ras1=mras1

     ! Number of permutations of the alpha electrons RAS1 with zero holes
     nperma_ras1(0)=1
     
     ! Number of permutations of the beta electrons RAS1 with zero holes
     npermb_ras1(0)=1

     ! Number of permutations of the alpha electrons RAS1 with one hole
     nperma_ras1(-1)=bico(mras1,1)
     
     ! Number of permutations of the beta electrons in RAS1 with one hole
     npermb_ras1(-1)=bico(mras1,1)
     
     ! Number of permutations of the alpha electrons RAS1 with two holes
     if (na_ras1>1) then
        nperma_ras1(-2)=bico(mras1,2)
     else
        nperma_ras1(-2)=0
     endif

     ! Number of permutations of the beta electrons RAS1 with two holes
     if (nb_ras1>1) then
        npermb_ras1(-2)=bico(mras1,2)
     else
        npermb_ras1(-2)=0
     endif

     ! Allocate the RAS1 electron permutation array
     v_ras1_dim=max(maxval(nperma_ras1),maxval(npermb_ras1))
     allocate(v_ras1(v_ras1_dim,-nras1:0,2))
     v_ras1=0_ib
     
  case(.false.) ! The number of RAS1 orbitals is zero

     ! We will only have RAS2->RAS3 excitations in the wavefunction.
     ! We need to set things up so that the base determinant
     ! orbital occupations are used in the RAS1 space

     ! Set nperma_ras1(0)=npermb_ras1(0)=1
     ! Set all other elements to 0
     nperma_ras1=0
     npermb_ras1=0
     nperma_ras1(0)=1
     npermb_ras1(0)=1

  end select
       
!----------------------------------------------------------------------
! Compute the permutations of the alpha and beta electrons in the RAS1
! orbitals
!----------------------------------------------------------------------
  if (okras1) then

     ! Permutations of the alpha electrons in RAS1 with zero holes
     v_ras1(1,0,1)=ishft(1_ib,mras1)-1_ib
    
     ! Permutations of the beta electrons in RAS1 with zero holes
     v_ras1(1,0,2)=ishft(1_ib,mras1)-1_ib
     
     ! Permutations of the alpha electrons in RAS1 with one hole
     call get_permutations(mras1-1,1,v_ras1(1:nperma_ras1(-1),-1,1),&
          nperma_ras1(-1))
     
     ! Permutations of the beta electrons in RAS1 with one hole
     call get_permutations(mras1-1,1,v_ras1(1:npermb_ras1(-1),-1,2),&
          npermb_ras1(-1))
     
     ! Permutations of the alpha electrons in RAS1 with two holes
     if (mras1>1.and.nras1>1) then
        call get_permutations(mras1-2,2,v_ras1(1:nperma_ras1(-2),-2,1),&
             nperma_ras1(-2))
     endif
     
     ! Permutations of the beta electrons in RAS1 with two holes
     !if (mras1-2>0) then
     if (mras1>1.and.nras1>1) then
        call get_permutations(mras1-2,2,v_ras1(1:npermb_ras1(-2),-2,2),&
             npermb_ras1(-2))
     endif
     
  endif
  
!----------------------------------------------------------------------
! Determine the number of permumations of electrons in the RAS3
! orbitals
!----------------------------------------------------------------------
  select case(okras3)
      
  case(.true.) ! The number of RAS3 orbitals is non-zero
  
     ! Number of permutations of the alpha holes RAS3 with zero particles
     nperma_ras3(0)=1
     
     ! Number of permutations of the beta holes RAS3 with zero particles
     npermb_ras3(0)=1
     
     ! Number of permutations of the alpha holes RAS3 with one particle
     nperma_ras3(1)=bico(mras3,1)
     
     ! Number of permutations of the beta electrons in RAS3 with one particle
     npermb_ras3(1)=bico(mras3,1)
     
     ! Number of permutations of the alpha electrons RAS3 with two particles
     if (mras3>1.and.nras3>1) then
        nperma_ras3(2)=bico(mras3,2)
     else
        nperma_ras3(2)=0
     endif
     
     ! Number of permutations of the beta electrons RAS3 with two particles
     if (mras3>1.and.nras3>1) then
        npermb_ras3(2)=bico(mras3,2)
     else
        npermb_ras3(2)=0
     endif

     ! Allocate the RAS3 electron permutation array
     v_ras3_dim=max(maxval(nperma_ras3),maxval(npermb_ras3))
     allocate(v_ras3(v_ras3_dim,0:nras3,2))
     v_ras3=0_ib
     
  case(.false.) ! The number of RAS3 orbitals is non-zero

     ! We will only have RAS1->RAS2 excitations in the wavefunction.
     ! We need to set things up so that the base determinant
     ! orbital occupations are used in the RAS3 space
     
     ! Set nperma_ras3(0)=npermb_ras3(0)=1
     ! Set all other elements to 0
     nperma_ras3=0
     npermb_ras3=0
     nperma_ras3(0)=1
     npermb_ras3(0)=1
     
  end select

!----------------------------------------------------------------------
! Compute the permutations of the alpha and beta electrons in the RAS1
! orbitals
!----------------------------------------------------------------------
  if (okras3) then

     ! Permutations of the alpha electrons in RAS3 with zero particles
     v_ras3(1,0,1)=0_ib

     ! Permutations of the beta electrons in  RAS3 with zero particles
     v_ras3(1,0,2)=0_ib

     ! Permutations of the alpha electrons in RAS3 with one particle
     call get_permutations(1,mras3-1,&
          v_ras3(1:nperma_ras3(1),1,1),nperma_ras3(1))

     ! Permutations of the beta electrons in RAS3 with one particle
     call get_permutations(1,mras3-1,&
          v_ras3(1:npermb_ras3(1),1,2),npermb_ras3(1))

     ! Permutations of the alpha electrons in RAS3 with two particles
     if (mras3>1.and.nras3>1) then
        call get_permutations(2,mras3-2,&
             v_ras3(1:nperma_ras3(2),2,1),nperma_ras3(2))
     endif

     ! Permutations of the beta electrons in RAS3 with two particles
     if (mras3>1.and.nras3>1) then
        call get_permutations(2,mras3-2,&
             v_ras3(1:npermb_ras3(2),2,2),npermb_ras3(2))
     endif

  endif
  
  return
  
end subroutine get_common_ras_data

!######################################################################
! make_ras_dets: computes all the determinants of the RAS wavefunction
!                and writes them to disk
!######################################################################
subroutine make_ras_dets(detras,ndet,iras1,iras2,iras3,nras1,mras1,&
     nras2,mras2,nras3,mras3)

  use constants
  use bitglobal
  use rasperm
  use iomod

  implicit none

  integer(is), intent(in)    :: ndet
  integer(ib), intent(inout) :: detras(n_int,2,ndet)
  integer(is), intent(in)    :: iras1(nmo),iras2(nmo),iras3(nmo)
  integer(is), intent(in)    :: nras1,mras1,nras2,mras2,nras3,mras3
  integer(is)                :: i,i1,i2,imo,k,ic,n,ispin
  integer(is)                :: ia_ras1,ib_ras1,ja_ras1,jb_ras1
  integer(is)                :: ia_ras2,ib_ras2,ja_ras2,jb_ras2
  integer(is)                :: ia_ras3,ib_ras3,ja_ras3,jb_ras3
  integer(is)                :: icount
  integer(ib)                :: det(n_int,2)

!----------------------------------------------------------------------
! Construct the RAS determinants  
!----------------------------------------------------------------------
  ! Initialise the determinant counter
  icount=0
  
  !
  ! Outer loop over the different RAS excitation classes
  !
  !
  ! Loop over alpha and beta RAS1 hole numbers
  do ia_ras1=-nras1,0
     do ib_ras1=-nras1,0
        ! Loop over alpha and beta RAS2 hole/particle numbers
        do ia_ras2=-nras3,nras2
           do ib_ras2=-nras3,nras2
              ! Loop over alpha and beta RAS3 particle numbers
              do ia_ras3=0,nras3
                 do ib_ras3=0,nras3

                    ! Cycle if the current set of excitations is not
                    ! electron-number conserving (this will be for the
                    ! majority of iterations)
                    if (ia_ras1+ib_ras1+ia_ras2+ib_ras2+ia_ras3+ib_ras3&
                         /=0) cycle

                    ! Cycle if the current set of excitations is not
                    ! spin conserving
                    if (ia_ras1+ia_ras2+ia_ras3/=0) cycle
                    if (ib_ras1+ib_ras2+ib_ras3/=0) cycle
                    
                    ! Cycle if the number of holes in RAS1 is too large
                    if (abs(ia_ras1+ib_ras1)>nras1) cycle

                    ! Cycle if the number of particles in RAS3 is too large
                    if (ia_ras3+ib_ras3>nras3) cycle
                    
                    !
                    ! Inner loop over determinants within each RAS
                    ! excitation class
                    !
                    do ja_ras1=1,nperma_ras1(ia_ras1)
                       do jb_ras1=1,npermb_ras1(ib_ras1)
                          do ja_ras2=1,nperma_ras2(ia_ras2)
                             do jb_ras2=1,npermb_ras2(ib_ras2)
                                do ja_ras3=1,nperma_ras3(ia_ras3)
                                   do jb_ras3=1,npermb_ras3(ib_ras3)

                                      ! Increment the determinant counter
                                      icount=icount+1

                                      ! Initialise the determinant to
                                      ! the base determinant
                                      det=det0
                                      
                                      ! Construct the RAS1 alpha-string
                                      if (okras1) &
                                           call fill_orbitals(det,1,&
                                           mras1,&
                                           v_ras1(ja_ras1,ia_ras1,1),&
                                           iras1)
                                      
                                      ! Construct the RAS1 beta-string
                                      if (okras1) &
                                           call fill_orbitals(det,2,&
                                           mras1,&
                                           v_ras1(jb_ras1,ib_ras1,2),&
                                           iras1)
                                      
                                      ! Construct the RAS2 alpha-string
                                      if (okras2) &
                                           call fill_orbitals(det,1,&
                                           mras2,&
                                           v_ras2(ja_ras2,ia_ras2,1),&
                                           iras2)
                                      
                                      ! Construct the RAS2 beta-string
                                      if (okras2) &
                                           call fill_orbitals(det,2,&
                                           mras2,&
                                           v_ras2(jb_ras2,ib_ras2,2),&
                                           iras2)
                                      
                                      ! Construct the RAS3 alpha-string
                                      if (okras3) &
                                           call fill_orbitals(det,1,&
                                           mras3,&
                                           v_ras3(ja_ras3,ia_ras3,1),&
                                           iras3)
                                      
                                      ! Construct the RAS3 beta-string
                                      if (okras3) &
                                           call fill_orbitals(det,2,&
                                           mras3,&
                                           v_ras3(jb_ras3,ib_ras3,2),&
                                           iras3)

                                      ! Save the determinant
                                      detras(:,:,icount)=det
                                      
                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo

                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
    
  return

end subroutine make_ras_dets

!######################################################################
! get_ndet_ras_dumb: determines the total number of RAS determiants
!######################################################################
subroutine get_ndet_ras_dumb(ndet,nras1,nras2,nras3)

  use constants
  use bitglobal
  use rasperm

  integer(is), intent(out)   :: ndet
  integer(is), intent(in)    :: nras1,nras2,nras3
  integer(is)                :: ia_ras1,ib_ras1,ja_ras1,jb_ras1
  integer(is)                :: ia_ras2,ib_ras2,ja_ras2,jb_ras2
  integer(is)                :: ia_ras3,ib_ras3,ja_ras3,jb_ras3

  ! Initialise the determinant counter
  ndet=0

  !
  ! Outer loop over the different RAS excitation classes
  !
  !
  ! Loop over alpha and beta RAS1 hole numbers
  do ia_ras1=-nras1,0
     do ib_ras1=-nras1,0
        ! Loop over alpha and beta RAS2 hole/particle numbers
        do ia_ras2=-nras3,nras2
           do ib_ras2=-nras3,nras2
              ! Loop over alpha and beta RAS3 particle numbers
              do ia_ras3=0,nras3
                 do ib_ras3=0,nras3

                    ! Cycle if the current set of excitations is not
                    ! electron-number conserving (this will be for the
                    ! majority of iterations)
                    if (ia_ras1+ib_ras1+ia_ras2+ib_ras2+ia_ras3+ib_ras3&
                         /=0) cycle

                    ! Cycle if the current set of excitations is not
                    ! spin conserving
                    if (ia_ras1+ia_ras2+ia_ras3/=0) cycle
                    if (ib_ras1+ib_ras2+ib_ras3/=0) cycle
                    
                    ! Cycle if the number of holes in RAS1 is too large
                    if (abs(ia_ras1+ib_ras1)>nras1) cycle

                    ! Cycle if the number of particles in RAS3 is too large
                    if (ia_ras3+ib_ras3>nras3) cycle
                    
                    !
                    ! Inner loop over determinants within each RAS
                    ! excitation class
                    !
                    do ja_ras1=1,nperma_ras1(ia_ras1)
                       do jb_ras1=1,npermb_ras1(ib_ras1)
                          do ja_ras2=1,nperma_ras2(ia_ras2)
                             do jb_ras2=1,npermb_ras2(ib_ras2)
                                do ja_ras3=1,nperma_ras3(ia_ras3)
                                   do jb_ras3=1,npermb_ras3(ib_ras3)

                                      ! Increment the determinant counter
                                      ndet=ndet+1

                                   enddo
                                enddo
                             enddo
                          enddo
                       enddo
                    enddo

                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  
  return
  
end subroutine get_ndet_ras_dumb
  
!######################################################################
! fill_orbitals: sets the ispin (1<->alpha, 2<->beta) bits in the
!                determinant det corresponding to the morb orbitals
!                in the iorb array using the permutation pattern v
!######################################################################
subroutine fill_orbitals(det,ispin,morb,v,iorb)

  use constants
  use bitglobal
  
  implicit none

  integer(ib), intent(out) :: det(n_int,2)
  integer(is), intent(in)  :: ispin
  integer(is), intent(in)  :: morb
  integer(ib), intent(in)  :: v
  integer(is), intent(in)  :: iorb(nmo)
  integer(is)              :: i1,i,k
  
  ! Loop over orbitals
  do i1=1,morb

     ! Block index
     k=(iorb(i1)-1)/64+1

     ! Index of the orbital in the kth block
     i=iorb(i1)-(k-1)*64-1

     ! Set the bit
     if (btest(v,i1-1)) then
        ! Set the bit to 1
        det(k,ispin)=ibset(det(k,ispin),i)
     else
        ! Set the bit to 0
        det(k,ispin)=ibclr(det(k,ispin),i)
     endif
     
  enddo
     
  return

end subroutine fill_orbitals
  
!######################################################################
! write_ras_dets_sorted: Writes the double-sorted RAS determinants and
!                        offset information to disk
!######################################################################
subroutine write_ras_dets_sorted(scrnum,detras_a,detras_b,ndet,nsym,&
     mapab,offset_a,offset_b,offdim_a,offdim_b,nunique_a,nunique_b)

  use constants
  use bitglobal
  use iomod
  
  implicit none

  integer(is), intent(out) :: scrnum
  integer(is), intent(in)  :: ndet,offdim_a,offdim_b
  integer(ib), intent(in)  :: detras_a(n_int,2,ndet),&
                              detras_b(n_int,2,ndet)
  integer(is), intent(in)  :: nsym(0:nirrep-1)
  integer(is), intent(in)  :: mapab(ndet)
  integer(is), intent(in)  :: offset_a(offdim_a,0:nirrep-1),&
                              offset_b(offdim_b,0:nirrep-1)
  integer(is), intent(in)  :: nunique_a(0:nirrep-1),&
                              nunique_b(0:nirrep-1)
  integer(is)              :: iscratch
  character(len=60)        :: rasfile
  character(len=2)         :: amult
  
  !
  ! Register a new scratch file
  !
  write(amult,'(i0)') imult
  call scratch_name('detras_sorted.mult'//trim(amult),rasfile)
  call register_scratch_file(scrnum,rasfile)

  !
  ! Open the scratch file
  !
  iscratch=scrunit(scrnum)
  open(iscratch,file=scrname(scrnum),form='unformatted',&
       status='unknown')

  !
  ! Number of RAS determinants
  !
  write(iscratch) ndet

  !
  ! Leading dimensions of the offset arrays
  !
  write(iscratch) offdim_a
  write(iscratch) offdim_b
  
  !
  ! Number of determinants per irrep
  !
  write(iscratch) nsym

  !
  ! RAS determinants in alpha- and beta-major order
  !
  write(iscratch) detras_a
  write(iscratch) detras_b
  
  !
  ! Number of unique alpha and beta strings per irrep
  !
  write(iscratch) nunique_a
  write(iscratch) nunique_b
  
  !
  ! Offsets for the alpha- and beta-major ordered determinants
  !
  write(iscratch) offset_a
  write(iscratch) offset_b
  
  !
  ! Mapping between the alpha- and beta-major ordered determinants
  !
  write(iscratch) mapab
  
  !
  ! Close the scratch file
  !
  close(iscratch)
  
  return
  
end subroutine write_ras_dets_sorted

!######################################################################
! write_ras_dets_unsorted: Writes the unsorted RAS determinants to disk
!######################################################################
subroutine write_ras_dets_unsorted(scrnum,detras,ndet)

  use constants
  use bitglobal
  use iomod
  
  implicit none

  integer(is), intent(out) :: scrnum
  integer(is), intent(in)  :: ndet
  integer(ib), intent(in)  :: detras(n_int,2,ndet)
  integer(is)              :: iscratch
  character(len=60)        :: rasfile
  character(len=2)         :: amult
  
  !
  ! Register a new scratch file
  !
  write(amult,'(i0)') imult
  call scratch_name('detras.mult'//trim(amult),rasfile)
  call register_scratch_file(scrnum,rasfile)

  !
  ! Open the scratch file
  !
  iscratch=scrunit(scrnum)
  open(iscratch,file=scrname(scrnum),form='unformatted',&
       status='unknown')

  !
  ! Number of RAS determinants
  !
  write(iscratch) ndet

  !
  ! RAS determinants
  !
  write(iscratch) detras

  !
  ! Close the scratch file
  !
  close(iscratch)
  
  return
  
end subroutine write_ras_dets_unsorted
