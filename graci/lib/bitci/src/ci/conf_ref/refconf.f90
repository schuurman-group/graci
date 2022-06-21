!**********************************************************************
! Routines for the generation of bit string representations of the
! reference space configurations of a RAS wavefunction
!**********************************************************************

module refconf

contains
  
!######################################################################
! generate_ref_confs: generates the set of configurations corresponding
!                     to a RAS active space and the base determinant
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
!----------------------------------------------------------------------
! Output variables:  -nconf:   numbers of reference configurations per
!                              irrep
!                    -scrnum:  the scratch file numbers for the
!                              reference configurations of each irrep
!######################################################################
#ifdef CBINDING
  subroutine generate_ref_confs(iras1,iras2,iras3,nras1,mras1,nras2,&
       mras2,nras3,mras3,icvs,nconf,scrnum) &
       bind(c,name="generate_ref_confs")
#else
  subroutine generate_ref_confs(iras1,iras2,iras3,nras1,mras1,nras2,&
       mras2,nras3,mras3,icvs,nconf,scrnum)
#endif

    use constants
    use bitglobal
    use dethash
    use mrciutils
    use timing
  
    implicit none
    
    ! RAS wavefunction information
    integer(is), intent(inout) :: iras1(nmo),iras2(nmo),iras3(nmo)
    integer(is), intent(in)    :: nras1,mras1,nras2,mras2,nras3,mras3
    integer(is), intent(out)   :: nconf(0:nirrep-1)
    integer(is), intent(out)   :: scrnum(0:nirrep-1)
    integer(is)                :: ndet,detscr

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    integer(is)                :: ncvs
    
    ! SOP hash table
    type(dhtbl)                :: h
    
    ! SOPs
    integer(ib), allocatable   :: sop(:,:,:)
    
    ! Configurations
    integer(is)                :: ntot
    integer(ib), allocatable   :: conf(:,:,:)
    integer(is)                :: nsym(0:nirrep-1)
    
    ! Internal MO space information
    integer(is)                :: nmoI,nmoE
    integer(is)                :: Ilist(nmo),Elist(nmo)
    
    ! Canonical-to-MRCI MO mapping
    integer(is)                :: m2c(nmo),c2m(nmo)
  
    ! Timing
    real(dp)                   :: tcpu_start,tcpu_end,twall_start,&
                                  twall_end
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! CVS-MRCI: determine the dimension of the core MO space in which
! a single hole will be enforced
!----------------------------------------------------------------------
    ncvs=sum(icvs)
    
!----------------------------------------------------------------------
! Generate the RAS determinants
!----------------------------------------------------------------------
    call generate_ras_dets(iras1,iras2,iras3,nras1,mras1,nras2,mras2,&
         nras3,mras3,ndet,detscr,.false.)
    
!----------------------------------------------------------------------
! Determine the unique reference space SOPs
!----------------------------------------------------------------------
    call get_ref_sops(ndet,detscr,h,ncvs,icvs)

!----------------------------------------------------------------------
! Allocate the reference space configuration and SOP arrays
!----------------------------------------------------------------------
    ! Total number of ref confugurations
    ntot=h%n_keys_stored
    
    ! Configurations
    allocate(conf(n_int,2,ntot))
    conf=0_ib
    
    ! SOPs
    allocate(sop(n_int,2,ntot))
    sop=0_ib
    
!----------------------------------------------------------------------
! Retrieve the SOPs from the hash table
!----------------------------------------------------------------------
    call h%retrieve_keys(sop)
  
!----------------------------------------------------------------------  
! Construct the reference space configuration bit string pairs from
! the SOPs
!----------------------------------------------------------------------
    call get_ref_confs(conf,sop,ntot)
  
!----------------------------------------------------------------------
! Get the list of internal MOs
! Here, Ilist is the list of internal MO indices, Elist is the list of
! external MO indices, nmoI is the number of internal MOs, and nmoE
! is the number of external MOs
!----------------------------------------------------------------------
    call get_internal_external_mos(conf,ntot,nmoI,nmoE,Ilist,Elist)

!----------------------------------------------------------------------
! Construct the canonical-to-MRCI MO index mapping array
!----------------------------------------------------------------------
    call get_mo_mapping(nmoI,nmoE,Ilist,Elist,m2c,c2m)
  
!----------------------------------------------------------------------
! Re-arrange the reference space configurations such that the internal
! MOs come before the external MOs
!----------------------------------------------------------------------
    call rearrange_ref_confs(c2m,conf,sop,ntot)

!----------------------------------------------------------------------
! Sort the reference space configurations by symmetry
!----------------------------------------------------------------------
    ! Sort the confs and SOPs
    call sort_by_irrep(nsym,conf,sop,ntot,m2c)

    ! Save the numbers of confs per irrep
    nconf=nsym

!----------------------------------------------------------------------
! Write the reference space configuration information to disk
!----------------------------------------------------------------------
    call write_ref_confs(conf,sop,ntot,m2c,c2m,nmoI,nmoE,nsym,scrnum)
  
!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
    call h%delete_table

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(conf)
    deallocate(sop)
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'generate_ref_confs')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
    flush(6)
    
    return
  
  end subroutine generate_ref_confs

!######################################################################
! get_ref_sops: Determines the set of unique reference space SOPs
!               from the set of reference space determinants
!######################################################################
  subroutine get_ref_sops(ndet,detscr,h,ncvs,icvs)

    use constants
    use bitglobal
    use detutils
    use dethash
    use slater_condon
    use mrciutils
    use iomod
  
    implicit none

    ! Number of ref determinants and the ref determinant scratch
    ! file number
    integer(is), intent(in)  :: ndet,detscr

    ! CVS-MRCI: core MOs
    integer(is), intent(in)  :: icvs(nmo)
    integer(is), intent(in)  :: ncvs
    
    ! Determinant arrays
    integer(is)              :: ndet1
    integer(ib), allocatable :: det(:,:,:)

    ! SOP hash table
    type(dhtbl), intent(out) :: h
    integer(is)              :: initial_size
    integer(ib)              :: key(n_int,2)
    
    ! Number of open shells
    integer(is)              :: nopen
  
    ! Everything else
    integer(is)              :: i,nexci
  
!----------------------------------------------------------------------
! Read the ref determinants from disk
!----------------------------------------------------------------------
    call read_det_file_unsorted(detscr,ndet1,det)

!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------  
    initial_size=10
    call h%initialise_table(initial_size)
  
!----------------------------------------------------------------------
! Insert the ref SOPs into hash table
!----------------------------------------------------------------------
    ! Loop over ref determinants
    do i=1,ndet
       
       ! Calculate the SOP
       call get_sop(det(:,:,i),key)

       ! Cycle if the number of open shells is greater than nomax
       nopen=sop_nopen(key,n_int)
       if (nopen > nomax) cycle

       ! Cycle if the spin multiplicity cannot be supported by the
       ! number of open shells
       if (nopen < imult-1) cycle
       
       ! Cycle if the excitation degree relative to the base
       ! determinant is greater than nexmax
       nexci=exc_degree_det(det(:,:,i),det0)
       if (nexci > nexmax) cycle

       ! If this is a CVS-MRCI calculation, then cycle if the no. holes
       ! in the selected core MO space is not equal to 1
       if (ncvs > 0 .and. nhole_cvs(key,icvs) /= 1) cycle
       
       ! Insert the SOP into the hash table
       call h%insert_key(key)
       
    enddo
  
!----------------------------------------------------------------------
! Dellocate arrays
!----------------------------------------------------------------------
    deallocate(det)
  
    return
  
  end subroutine get_ref_sops

!######################################################################
! nhole_cvs: Given an SOP and a list of CVS core-level MOs, returns
!            the no. holes in the CVS core level space
!######################################################################
  function nhole_cvs(sop,icvs) result(nhole)

    use constants
    use bitglobal
    
    implicit none

    ! Function result
    integer(is)             :: nhole

    ! SOP
    integer(ib), intent(in) :: sop(n_int,2)

    ! CVS core-level MOs
    integer(is), intent(in) :: icvs(nmo)

    ! Everything else
    integer(is)             :: i,k,imo

    ! Initialisation
    nhole=0

    ! Loop over MOs in the CVS core-level space
    do imo=1,nmo
       if (icvs(imo) == 1) then

          ! Block index
          k=(imo-1)/n_bits+1
          
          ! Orbital index in the bit string
          i=imo-(k-1)*n_bits-1

          if (btest(sop(k,2),i)) then
             ! Doubly-occupied MO
             cycle
          else if (btest(sop(k,1),i)) then
             ! Singly-occupied MO
             nhole=nhole+1
          else
             ! Unoccupied MO
             nhole=nhole+2
          endif
          
       endif
    enddo
    
    return

  end function nhole_cvs
  
!######################################################################
! get_ref_confs: Generates the reference space confs from the
!                corresponding SOPs
!######################################################################
  subroutine get_ref_confs(conf,sop,ntot)

    use constants
    use bitglobal
  
    implicit none

    integer(is), intent(in)  :: ntot
    integer(ib), intent(in)  :: sop(n_int,2,ntot)
    integer(ib), intent(out) :: conf(n_int,2,ntot)
    integer(is)              :: i,k

    ! Loop over SOPs
    do i=1,ntot
       
       ! Construct the bit string pair encoding of the configuration
       do k=1,n_int
          conf(k,1,i)=ieor(sop(k,1,i),sop(k,2,i))
          conf(k,2,i)=sop(k,2,i)
       enddo
       
    enddo
    
    return
    
  end subroutine get_ref_confs

!######################################################################
! rearrange_ref_confs: Rearranges the reference configurations to be
!                      in the MRCI ordering. That is, with the
!                      internal MOs before the external MOs
!######################################################################
  subroutine rearrange_ref_confs(c2m,conf,sop,ntot)

    use constants
    use bitglobal
  
    implicit none
    
    integer(is), intent(in)    :: c2m(nmo)
    integer(is), intent(in)    :: ntot
    integer(ib), intent(inout) :: conf(n_int,2,ntot)
    integer(ib), intent(inout) :: sop(n_int,2,ntot)
    integer(ib)                :: conf_new(n_int,2)
    integer(ib)                :: sop_new(n_int,2)
    integer(is)                :: iconf,k,i,imo,n
    integer(is)                :: imo1,k1,i1

!----------------------------------------------------------------------
! Rearrange the reference space configurations and associated SOPs
! to correspond to the MRCI internal-external MO ordering
!----------------------------------------------------------------------
    ! Loop over reference space configurations
    do iconf=1,ntot

       ! Initialise the reordered conf and SOP
       conf_new=0_ib
       sop_new=0_ib
       
       ! Loop over MOs
       do imo=1,nmo
          
          ! Block index
          k=(imo-1)/n_bits+1
          
          ! Orbital index in the bit string
          i=imo-(k-1)*n_bits-1
          
          ! Fill in the new configuration and SOP bit strings
          do n=1,2
             
             ! New MO index
             imo1=c2m(imo)
             
             ! New block index
             k1=(imo1-1)/n_bits+1
             
             ! New orbital index in the bit string
             i1=imo1-(k1-1)*n_bits-1
             
             ! Configuration
             if (btest(conf(k,n,iconf),i)) &
                  conf_new(k1,n)=ibset(conf_new(k1,n),i1)
             
             ! SOP
             if (btest(sop(k,n,iconf),i)) &
                  sop_new(k1,n)=ibset(sop_new(k1,n),i1)
              
          enddo
           
       enddo

       ! Save the configuration and SOP
       conf(:,:,iconf)=conf_new
       sop(:,:,iconf)=sop_new
       
    enddo
     
    return
  
  end subroutine rearrange_ref_confs

!######################################################################
! sort_by_irrep: Sorts the reference configurations by irrep
!######################################################################
  subroutine sort_by_irrep(nsym,conf,sop,ntot,m2c)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    integer(is), intent(out)   :: nsym(0:nirrep-1)
    integer(is), intent(in)    :: ntot
    integer(ib), intent(inout) :: conf(n_int,2,ntot),sop(n_int,2,ntot)
    integer(is), intent(in)    :: m2c(nmo)
    integer(is), allocatable   :: isym(:)
    integer(ib), allocatable   :: conf_sorted(:,:,:),sop_sorted(:,:,:)
    integer(is)                :: i,k,indx
    integer(is), allocatable   :: ioff(:),counter(:)
    integer(is)                :: sum
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(isym(ntot))
    allocate(conf_sorted(n_int,2,ntot))
    allocate(sop_sorted(n_int,2,ntot))
    allocate(ioff(0:nirrep-1))
    allocate(counter(0:nirrep-1))
  
!----------------------------------------------------------------------
! Determine the irreps generated by each configuration
!----------------------------------------------------------------------
    ! Loop over configurations
    do i=1,ntot
       
       ! Compute the irrep generated by the current configuration
       isym(i)=sop_sym_mrci(sop(:,:,i),m2c)
       
    enddo

!----------------------------------------------------------------------
! Determine the number of configurations per irrep
!----------------------------------------------------------------------
    nsym=0
    do i=1,ntot
       nsym(isym(i))=nsym(isym(i))+1
    enddo

!----------------------------------------------------------------------
! Sort the configurations by symmetry
!----------------------------------------------------------------------
    ! Offsets
    ioff=0
    sum=0
    do i=0,nirrep-1
       ioff(i)=sum
       sum=sum+nsym(i)
    enddo
    
    ! Fill in the irrep-sorted configuration and SOP arrays
    counter=0
    do k=1,ntot
       counter(isym(k))=counter(isym(k))+1
       indx=ioff(isym(k))+counter(isym(k))
       conf_sorted(:,:,indx)=conf(:,:,k)
       sop_sorted(:,:,indx)=sop(:,:,k)
    enddo
    
    ! Re-write the conf and SOP arrays with the sorted ones
    conf=conf_sorted
    sop=sop_sorted
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(isym)
    deallocate(conf_sorted)
    deallocate(sop_sorted)
    deallocate(ioff)
    deallocate(counter)
  
    return
  
  end subroutine sort_by_irrep

!######################################################################
! write_ref_confs: Writes the reference configurations, SOPs and
!                  information about the internal-external
!                  partitioning of the MOs to disk.
!######################################################################
  subroutine write_ref_confs(conf,sop,ntot,m2c,c2m,nmoI,nmoE,nsym,&
       scrnum)

    use constants
    use bitglobal
    use iomod
  
    implicit none

    integer(is), intent(in)  :: ntot,nmoI,nmoE
    integer(ib), intent(in)  :: conf(n_int,2,ntot),sop(n_int,2,ntot)
    integer(is), intent(in)  :: m2c(nmo),c2m(nmo)
    integer(is)              :: nsym(0:nirrep-1)
    integer(is), intent(out) :: scrnum(0:nirrep-1)
    integer(is)              :: n_int_I
    integer(is)              :: iscratch
    integer(is)              :: irrep,istart,iend
    character(len=60)        :: scrfile(0:nirrep-1)
    character(len=2)         :: amult,asym
  
!----------------------------------------------------------------------
! Register the scratch files
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1
       
       ! If there are no configurations for the current irrep, then
       ! flag this with a negative scratch file number and cycle
       if (nsym(irrep) == 0) then
          scrnum(irrep)=-1
          cycle
       endif
       
       ! Scratch file name
       write(amult,'(i0)') imult
       write(asym,'(i0)') irrep
       call scratch_name(&
            'refconf.mult'//trim(amult)//'.sym'//trim(asym),&
            scrfile(irrep))
       
       ! Register the scratch file
       call register_scratch_file(scrnum(irrep),scrfile(irrep))
       
    enddo

!----------------------------------------------------------------------
! Number of (n_bits)-bit integers needed to represent each reference
! space configuration and SOP bit string
!----------------------------------------------------------------------
    n_int_I=(nmoI-1)/n_bits+1
  
!----------------------------------------------------------------------
! Write the configuration information to the scratch files
!----------------------------------------------------------------------
! Important: We only save the first n_int_I integers representing each
!            bit string, as this is all that we need
!----------------------------------------------------------------------
    ! Initialise the configuration offset counters
    istart=0
    iend=0
  
    ! Loop over irreps
    do irrep=0,nirrep-1
       
       ! Cycle if there are no configurations for the current irrep
       if (nsym(irrep) == 0) cycle
       
       ! Open the scratch file
       iscratch=scrunit(scrnum(irrep))
       open(iscratch,file=scrname(scrnum(irrep)),form='unformatted',&
            status='unknown')
       
       ! Number of reference space configurations
       write(iscratch) nsym(irrep)
       
       ! Number of (n_bits)-bit integers needed to represent each
       ! reference space configuration and SOP bit string
       write(iscratch) n_int_I
       
       ! Number of internal and external MOs
       write(iscratch) nmoI
       write(iscratch) nmoE
       
       ! Update the offset counters
       istart=iend+1
       iend=iend+nsym(irrep)
       
       ! Configurations
       write(iscratch) conf(1:n_int_I,:,istart:iend)
       
       ! SOPs
       write(iscratch) sop(1:n_int_I,:,istart:iend)
       
       ! MO mapping arrays
       write(iscratch) m2c
       write(iscratch) c2m
     
       ! Close the scratch file
       close(iscratch)
       
    enddo
    
    return
    
  end subroutine write_ref_confs

!######################################################################

end module refconf
