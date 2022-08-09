!**********************************************************************
! Routines for the generation of the Slater determinants of an MRCI
! wavefunction
!**********************************************************************

!######################################################################
! generate_mrci_dets: Generates single or single+double excitations
!                     out of some set of a given set of reference
!                     space determinants
!######################################################################
#ifdef CBINDING
subroutine generate_mrci_dets(refscr,mrciscr,ndet,order) &
     bind(c,name="generate_mrci_dets")
#else
subroutine generate_mrci_dets(refscr,mrciscr,ndet,order)
#endif

  use constants
  use bitglobal
  use utils
  use dethash
  use detsort
  use iomod
  use timing
  
  implicit none

  ! Scratch file numbers for the reference space and MRCI determinants
  integer(is), intent(in)  :: refscr
  integer(is), intent(out) :: mrciscr
  
  ! MRCI order
  integer(is), intent(in)  :: order
  
  ! Total number of MRCI determinants
  integer(is), intent(out) :: ndet

  ! Determinant hash table
  type(dhtbl)              :: h
  
  ! Determinant scratch file
  integer(is)              :: iscratch
  character(len=250)       :: detfile
  character(len=2)         :: amult

  ! Determinant sorting
  integer(ib), allocatable :: detmrci_a(:,:,:),detmrci_b(:,:,:)
  integer(is)              :: offdim_a,offdim_b
  integer(is), allocatable :: mapab(:)
  integer(is), allocatable :: nsym(:)
  integer(is), allocatable :: offset_a(:,:),offset_b(:,:)
  integer(is), allocatable :: nunique_a(:),nunique_b(:)
  
!----------------------------------------------------------------------
! Check the excitation order
!----------------------------------------------------------------------
  if (order.lt.1.or.order.gt.2) then
     errmsg='Error in generate_mrci_dets: the excitation order can '&
          //'only be 1 or 2'
     call error_control
  endif

!----------------------------------------------------------------------
! Generate the MRCI determinants and store them in the hash table
!----------------------------------------------------------------------
  call make_mrci_dets(h,refscr,order)

!----------------------------------------------------------------------
! Number of MRCI determinants
!----------------------------------------------------------------------
  ndet=h%n_keys_stored
  if (verbose) &
       write(6,'(/,x,a,x,i0)') 'Total no. MRCI determinants:',ndet  

!----------------------------------------------------------------------
! Allocate the MRCI determinant arrays
!----------------------------------------------------------------------
  allocate(detmrci_a(n_int,2,ndet))
  allocate(detmrci_b(n_int,2,ndet))
  allocate(mapab(ndet))
  allocate(nsym(0:nirrep-1))
  allocate(nunique_a(0:nirrep-1))
  allocate(nunique_b(0:nirrep-1))
  
!----------------------------------------------------------------------
! Retrieve the MRCI determinants from the hash table
!----------------------------------------------------------------------
  call h%retrieve_keys(detmrci_a)
  detmrci_b=detmrci_a
  
!----------------------------------------------------------------------
! Sort the MRCI determinants
!----------------------------------------------------------------------
  call det_sort_all(detmrci_a,detmrci_b,ndet,nsym,mapab,offset_a,&
       offset_b,offdim_a,offdim_b,nunique_a,nunique_b)
  
!----------------------------------------------------------------------
! Write the sorted MRCI determinants to disk
!----------------------------------------------------------------------
  call write_mrci_dets(mrciscr,detmrci_a,detmrci_b,ndet,nsym,mapab,&
       offset_a,offset_b,offdim_a,offdim_b,nunique_a,nunique_b)

!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
  call h%delete_table
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(detmrci_a)
  deallocate(detmrci_b)
  deallocate(mapab)
  deallocate(nsym)
  deallocate(nunique_a)
  deallocate(nunique_b)

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine generate_mrci_dets

!######################################################################

subroutine make_mrci_dets(h,refscr,order)

  use constants
  use bitglobal
  use detutils
  use utils
  use dethash
  use iomod
  
  implicit none

  ! Determinant hash table
  type(dhtbl), intent(out) :: h

  ! Scratch file number for the reference space determinants
  integer(is), intent(in)  :: refscr

  ! Excitation order
  integer(is), intent(in)  :: order

  ! Reference space determinants
  integer(is)              :: nref
  integer(ib), allocatable :: dref(:,:,:)
  
  ! Determinant generation variables
  integer(is)              :: ntot
  integer(is), allocatable :: occ(:,:,:),unocc(:,:,:)
  integer(is)              :: nea,neb,nha,nhb
  integer(is)              :: i,j,k,iref
  integer(is)              :: ia1,iap1,ia,iap,ika,ikap
  integer(is)              :: ic1,icp1,ic,icp,ikc,ikcp
  integer(is)              :: iexc1,iexc2
  integer(ib)              :: key(n_int,2)
  integer(is)              :: initial_size

!----------------------------------------------------------------------
! Read the reference space determinants from disk
!----------------------------------------------------------------------
  call read_dets_alpha_major(refscr,nref,dref)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! Indices of the occupied and unoccupied orbitals in the reference
  ! space alpha and beta bit strings
  allocate(occ(nmo,2,nref))
  allocate(unocc(nmo,2,nref))
  occ=0
  unocc=0

!----------------------------------------------------------------------
! Get the indices of the occupied orbitals in the alpha and beta bit
! strings
!----------------------------------------------------------------------
  ! Loop over reference determinants
  do i=1,nref
     
     ! Get the indices of the occupied orbitals
     call get_occ(dref(:,:,i),occ(:,:,i))
     
     ! Get the indices of the unoccupied orbitals
     call get_unocc(dref(:,:,i),unocc(:,:,i))
     
  enddo

!----------------------------------------------------------------------
! Determine the number of occupied and unoccupied alpha and beta
! orbitals
!----------------------------------------------------------------------
  ! Number of occupied alpha and beta orbitals
  nea=n_alpha_electrons(dref(:,:,1))
  neb=n_beta_electrons(dref(:,:,1))

  ! Number of unoccupied alpha and beta orbitals
  nha=nmo-nea
  nhb=nmo-neb

!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------  
  initial_size=10
  call h%initialise_table(initial_size)

!----------------------------------------------------------------------
! Insert the reference space determinants into the hash table
!----------------------------------------------------------------------
  do iref=1,nref
     call h%insert_key(dref(:,:,iref))
     ntot=ntot+1
  enddo

!----------------------------------------------------------------------
! Generate all the singly-excited determinants
!----------------------------------------------------------------------
  ! Initialise the determinant counter
  iexc1=0

  !
  ! Alpha spin excitations
  !
  ! Loop over reference space determinants
  do iref=1,nref

     ! Loop over annihilators
     do ia1=1,nea

        ! Annihilator index
        ia=occ(ia1,1,iref)

        ! Block that the annihilated electron sits in
        ika=(ia-1)/n_bits+1
        
        ! Loop over creators
        do ic1=1,nha

           ! Creator index
           ic=unocc(ic1,1,iref)

           ! Block that the created electron sits in
           ikc=(ic-1)/n_bits+1
           
           ! Increment the determinant counter
           iexc1=iexc1+1
           
           ! Form the singly-excited determinant
           key=dref(:,:,iref)
           key(ika,1)=ibclr(key(ika,1),ia-1)
           key(ikc,1)=ibset(key(ikc,1),ic-1)

           ! Insert the singly-excited determinant into the hash table
           call h%insert_key(key)
           ntot=ntot+1
           
        enddo
        
     enddo
     
  enddo
  
  !
  ! Beta spin excitations
  !
  ! Loop over reference space determinants
  do iref=1,nref

     ! Loop over annihilators
     do ia1=1,neb

        ! Annihilator index
        ia=occ(ia1,2,iref)

        ! Block that the annihilated electron sits in
        ika=(ia-1)/n_bits+1
        
        ! Loop over creators
        do ic1=1,nhb

           ! Creator index
           ic=unocc(ic1,2,iref)

           ! Block that the created electron sits in
           ikc=(ic-1)/n_bits+1
           
           ! Increment the determinant counter
           iexc1=iexc1+1
           
           ! Form the singly-excited determinant
           key=dref(:,:,iref)
           key(ika,2)=ibclr(key(ika,2),ia-1)
           key(ikc,2)=ibset(key(ikc,2),ic-1)

           ! Insert the singly-excited determinant into the hash table
           call h%insert_key(key)
           ntot=ntot+1
           
        enddo
        
     enddo
     
  enddo

!----------------------------------------------------------------------
! Generate all the doubly-excited determinants
!----------------------------------------------------------------------
  if (order==2) then

     ! Initialise the determinant counter
     iexc2=0

     !
     ! alpha alpha -> alpha alpha excitations
     !
     ! Loop over reference space determinants
     do iref=1,nref
        
        ! Loop over the first annihilator
        do ia1=1,nea-1
           
           ! First annihilator index
           ia=occ(ia1,1,iref)
           
           ! Block index for the first annihilator
           ika=(ia-1)/n_bits+1
           
           ! Loop over the second annihilator
           do iap1=ia1+1,nea

              ! Second annihilator index
              iap=occ(iap1,1,iref)

              ! Block index for the first annihilator
              ikap=(iap-1)/n_bits+1
           
              ! Loop over the first creator
              do ic1=1,nha-1
                 
                 ! First creator index
                 ic=unocc(ic1,1,iref)
                 
                 ! Block index for the first creator
                 ikc=(ic-1)/n_bits+1
                 
                 ! Loop over the second creator
                 do icp1=ic1+1,nha
                    
                    ! Second creator index
                    icp=unocc(icp1,1,iref)
                    
                    ! Block index for the second creator
                    ikcp=(icp-1)/n_bits+1
                    
                    ! Increment the determinant counter
                    iexc2=iexc2+1
                    
                    ! Form the doubly-excited determinant
                    key=dref(:,:,iref)
                    key(ika,1)=ibclr(key(ika,1),ia-1)
                    key(ikap,1)=ibclr(key(ikap,1),iap-1)
                    key(ikc,1)=ibset(key(ikc,1),ic-1)
                    key(ikcp,1)=ibset(key(ikcp,1),icp-1)
                    
                    ! Insert the doubly-excited determinant into the
                    ! hash table
                    call h%insert_key(key)

                 enddo
                 
              enddo
              
           enddo
           
        enddo
        
     enddo
     
     !
     ! beta beta -> beta beta excitations
     !
     ! Loop over reference space determinants
     do iref=1,nref
        
        ! Loop over the first annihilator
        do ia1=1,neb-1
           
           ! First annihilator index
           ia=occ(ia1,2,iref)
           
           ! Block index for the first annihilator
           ika=(ia-1)/n_bits+1
           
           ! Loop over the second annihilator
           do iap1=ia1+1,neb
              
              ! Second annihilator index
              iap=occ(iap1,2,iref)
              
              ! Block index for the first annihilator
              ikap=(iap-1)/n_bits+1
              
              ! Loop over the first creator
              do ic1=1,nhb-1
                 
                 ! First creator index
                 ic=unocc(ic1,2,iref)
                 
                 ! Block index for the first creator
                 ikc=(ic-1)/n_bits+1
                 
                 ! Loop over the second creator
                 do icp1=ic1+1,nhb
                    
                    ! Second creator index
                    icp=unocc(icp1,2,iref)
                    
                    ! Block index for the second creator
                    ikcp=(icp-1)/n_bits+1
                 
                    ! Increment the determinant counter
                    iexc2=iexc2+1

                    ! Form the doubly-excited determinant
                    key=dref(:,:,iref)
                    key(ika,2)=ibclr(key(ika,2),ia-1)
                    key(ikap,2)=ibclr(key(ikap,2),iap-1)
                    key(ikc,2)=ibset(key(ikc,2),ic-1)
                    key(ikcp,2)=ibset(key(ikcp,2),icp-1)

                    ! Insert the doubly-excited determinant into the
                    ! hash table
                    call h%insert_key(key)

                 enddo
                 
              enddo
              
           enddo
           
        enddo
        
     enddo
     
     !
     ! alpha beta -> alpha beta excitations
     !
     ! Loop over reference space determinants
     do iref=1,nref
        
        ! Loop over the first annihilator (alpha)
        do ia1=1,nea
           
           ! First annihilator index
           ia=occ(ia1,1,iref)
           
           ! Block index for the first annihilator
           ika=(ia-1)/n_bits+1
        
           ! Loop over the second annihilator (beta)
           do iap1=1,neb

              ! Second annihilator index
              iap=occ(iap1,2,iref)

              ! Block index for the first annihilator
              ikap=(iap-1)/n_bits+1
              
              ! Loop over the first creator (alpha)
              do ic1=1,nha
                 
                 ! First creator index
                 ic=unocc(ic1,1,iref)

                 ! Block index for the first creator
                 ikc=(ic-1)/n_bits+1
              
                 ! Loop over the second creator (beta)
                 do icp1=1,nhb

                    ! Second creator index
                    icp=unocc(icp1,2,iref)

                    ! Block index for the second creator
                    ikcp=(icp-1)/n_bits+1
                 
                    ! Increment the determinant counter
                    iexc2=iexc2+1
                    
                    ! Form the doubly-excited determinant
                    key=dref(:,:,iref)
                    key(ika,1)=ibclr(key(ika,1),ia-1)    ! alpha annihilator
                    key(ikap,2)=ibclr(key(ikap,2),iap-1) ! beta annihilator
                    key(ikc,1)=ibset(key(ikc,1),ic-1)    ! alpha creator
                    key(ikcp,2)=ibset(key(ikcp,2),icp-1) ! beta creator

                    ! Insert the doubly-excited determinant into the
                    ! hash table
                    call h%insert_key(key)
                    
                 enddo
                 
              enddo
              
           enddo
           
        enddo
        
     enddo

  endif

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(dref)
  deallocate(occ,unocc)
  
  return
  
end subroutine make_mrci_dets

!######################################################################

subroutine write_mrci_dets(mrciscr,detmrci_a,detmrci_b,ndet,nsym,&
     mapab,offset_a,offset_b,offdim_a,offdim_b,nunique_a,nunique_b)

  use constants
  use bitglobal
  use iomod
  
  implicit none

  integer(is), intent(out) :: mrciscr
  integer(is), intent(in)  :: ndet,offdim_a,offdim_b
  integer(ib), intent(in)  :: detmrci_a(n_int,2,ndet),&
                              detmrci_b(n_int,2,ndet)
  integer(is), intent(in)  :: nsym(0:nirrep-1)
  integer(is), intent(in)  :: mapab(ndet)
  integer(is), intent(in)  :: offset_a(offdim_a,0:nirrep-1),&
                              offset_b(offdim_b,0:nirrep-1)
  integer(is), intent(in)  :: nunique_a(0:nirrep-1),&
                              nunique_b(0:nirrep-1)
  integer(is)              :: iscratch
  character(len=250)       :: mrcifile
  character(len=2)         :: amult
  
  !
  ! Register a new scratch file
  !
  write(amult,'(i0)') imult
  call scratch_name('detmrci.mult'//trim(amult),mrcifile)
  call register_scratch_file(mrciscr,mrcifile)

  !
  ! Open the scratch file
  !
  iscratch=scrunit(mrciscr)
  open(iscratch,file=scrname(mrciscr),form='unformatted',&
       status='unknown')

  !
  ! Number of MRCI determinants
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
  write(iscratch) detmrci_a
  write(iscratch) detmrci_b
  
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

end subroutine write_mrci_dets

!######################################################################
