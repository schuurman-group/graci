!**********************************************************************
! Routines for the construction of a Hamiltonian matrix using the
! simple double loop algorithm
!**********************************************************************
module hbuild_double

  implicit none

contains

!######################################################################
! hii_double: Calculation of the on-diagonal Hamiltonian matrix
!             elements and the their spin-coupling averaged values
!######################################################################
  subroutine hii_double(nconf,hdim,offset,conf,sop,n_int_in,m2c,irrep,&
       hii,averageii)

    use constants
    use bitglobal
    use mrciutils
    use bitutils
    use hparam
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    use iomod
    use timing
    
    implicit none

    ! Irrep number
    integer(is), intent(in)  :: irrep
  
    ! Dimensions
    integer(is), intent(in)  :: nconf,hdim
    integer(is), intent(in)  :: offset(nconf+1)
    
    ! Configurations and SOPs
    integer(is), intent(in)  :: n_int_in
    integer(ib), intent(in)  :: conf(n_int_in,2,nconf)
    integer(ib), intent(in)  :: sop(n_int_in,2,nconf)

    ! MO mapping array
    integer(is), intent(in)  :: m2c(nmo)

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(out)    :: hii(hdim)
    
    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(out)    :: averageii(nconf)

    ! Temporary Hij array
    integer(is)              :: arrdim
    real(dp), allocatable    :: harr(:)

    ! Hamiltonian build variables
    integer(is)              :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)              :: nopen,nsocc,ndocc,nunocc
    integer(ib)              :: conf_full(n_int,2)
    integer(ib)              :: sop_full(n_int,2)

    ! Difference configuration information
    integer(is)              :: ndiff
    integer(is)              :: Dw(nmo,2)

    ! Number of open shells preceding each MO
    integer(is)              :: nbefore(nmo)
    
    ! Everything else
    integer(is)              :: nsp,iconf,icsf,iomega
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(harr(maxval(ncsfs(0:nomax))**2))    
    harr=0.0d0

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    hii=0.0d0
    averageii=0.0d0
    
!----------------------------------------------------------------------
! Compute the on-diagonal elements
!----------------------------------------------------------------------
    ! Loop over configurations
    do iconf=1,nconf
       
       ! Number of open shells
       nopen=sop_nopen(sop(:,:,iconf),n_int_in)
       
       ! Configuration and SOP in the full MO space
       conf_full=0_ib
       sop_full=0_ib
       conf_full(1:n_int_in,:)=conf(:,:,iconf)
       sop_full(1:n_int_in,:)=sop(:,:,iconf)
       
       ! Get all the configuration information needed to evaluate
       ! the on-diagonal Hamiltonian matrix element
       call package_confinfo(sop_full,conf_full,unocc,socc,docc,nunocc,&
            nsocc,ndocc,Dw,ndiff,nbefore)
       
       ! Number of CSFs generated by the configuration
       nsp=ncsfs(nopen)
       
       ! On-diagonal Hamiltonian matrix elements for all CSFs generated
       ! by the current configuration
       call hii_mrci(harr(1:nsp),nsp,sop_full,nopen,socc,nsocc, &
            docc,ndocc,unocc,nunocc,nbefore,Dw,ndiff,m2c)
       
       ! Apply any DFT/MRCI corrections to the on-diagonal
       ! matrix elements
       if (ldftmrci) call hii_dftmrci(harr(1:nsp),nsp,Dw,ndiff,&
            nopen,m2c,sop_full,socc,nsocc,nbefore)
       
       ! Fill in the hii array
       iomega=0
       do icsf=offset(iconf),offset(iconf+1)-1
          iomega=iomega+1
          hii(icsf)=harr(iomega)
       enddo

    enddo

!----------------------------------------------------------------------
! Average values of the on-diagonal matrix elements generated by each
! configuration
!----------------------------------------------------------------------
  ! Loop over configurations
  do iconf=1,nconf

     ! Number of configurations generated by the current configuration
     nsp=offset(iconf+1)-offset(iconf)

     ! Average on-diagonal Hamiltonian matrix element value
     do icsf=offset(iconf),offset(iconf+1)-1
        averageii(iconf)=averageii(iconf)+hii(icsf)
     enddo
     averageii(iconf)=averageii(iconf)/nsp

  enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(harr)
    
    return
    
  end subroutine hii_double
    
!######################################################################
! save_hij_double: Saves the non-zero off-diagonal elements of the
!                  Hamiltonian to disk
!######################################################################
  subroutine save_hij_double(nconf,hdim,offset,averageii,conf,sop,&
       n_int_in,m2c,irrep,hscr,nrec,hstem)
    
    use constants
    use bitglobal
    use mrciutils
    use bitutils
    use hparam
    use hbuild_mrci
    use mrci_integrals
    use dftmrci
    use iomod
    use timing
    
    implicit none

    ! Irrep number
    integer(is), intent(in)      :: irrep
  
    ! Dimensions
    integer(is), intent(in)      :: nconf,hdim
    integer(is), intent(in)      :: offset(nconf+1)

    ! Spin-coupling averaged on-diagonal Hamiltonian matrix elements
    real(dp), intent(in)         :: averageii(nconf)
    
    ! Configurations and SOPs
    integer(is), intent(in)      :: n_int_in
    integer(ib), intent(in)      :: conf(n_int_in,2,nconf)
    integer(ib), intent(in)      :: sop(n_int_in,2,nconf)

    ! MO mapping array
    integer(is), intent(in)      :: m2c(nmo)

    ! Hamiltonian scratch file number
    integer(is), intent(out)     :: hscr

    ! Number of Hamiltonian scratch file records
    integer(is), intent(out)     :: nrec

    ! Hamiltonian scratch file stem
    character(len=*), intent(in) :: hstem
    
    ! Hamiltonian build variables
    integer(is)                  :: nexci
    integer(is)                  :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)                  :: nopen,nsocc,ndocc,nunocc
    integer(ib)                  :: conf_full(n_int,2)
    integer(ib)                  :: sop_full(n_int,2)
    integer(is)                  :: knopen,bnopen
    integer(ib)                  :: bconf_full(n_int,2)
    integer(ib)                  :: kconf_full(n_int,2)
    integer(ib)                  :: bsop_full(n_int,2)
    integer(ib)                  :: ksop_full(n_int,2)
    integer(is), parameter       :: maxexci=2
    integer(is)                  :: hlist(maxexci),plist(maxexci)
        
    ! Difference configuration information
    integer(is)                  :: ndiff
    integer(is)                  :: Dw(nmo,2)
    
    ! Number of open shells preceding each MO
    integer(is)                  :: nbefore(nmo)
  
    ! Temporary Hij array
    integer(is)                  :: arrdim,harr2dim
    real(dp), allocatable        :: harr(:),harr2(:)

    ! I/O variables
    integer(is)                  :: iscratch
    character(len=250)           :: hamfile
    character(len=2)             :: amult,airrep

    ! Buffer
    integer(is)                  :: nbuf
    integer(is), allocatable     :: ibuffer(:,:)
    real(dp), allocatable        :: hbuffer(:)
    
    ! Everything else
    integer(is)                  :: iconf,kconf,bconf
    integer(is)                  :: icsf,kcsf,bcsf
    integer(is)                  :: iomega,komega,bomega
    integer(is)                  :: nsp,count,blim1,blim2,&
                                    klim1,klim2,bnsp,knsp
    
    ! Timing variables
    real(dp)                     :: tcpu_start,tcpu_end,twall_start,&
                                    twall_end

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)
    
!----------------------------------------------------------------------
! Register Hamiltonian scratch file
!----------------------------------------------------------------------
    write(amult,'(i0)') imult
    write(airrep,'(i0)') irrep
    call scratch_name(trim(hstem)//'.mult'//trim(amult)//'.sym'&
         //trim(airrep),hamfile)    
    call register_scratch_file(hscr,hamfile)

!----------------------------------------------------------------------
! Open the Hamiltonian scratch file
!----------------------------------------------------------------------
    iscratch=scrunit(hscr)
    open(iscratch,file=scrname(hscr),form='unformatted',&
         status='unknown')
    
!----------------------------------------------------------------------
! Initialise the buffer and working arrays
!----------------------------------------------------------------------
    ! Buffer variables & arrays
    allocate(ibuffer(2,bufsize))
    ibuffer=0

    allocate(hbuffer(bufsize))
    hbuffer=0.0d0

    nrec=0
    nbuf=0

    ! H_ij working arrays
    harr2dim=maxval(ncsfs(0:nomax))**2
    allocate(harr(harr2dim))
    allocate(harr2(harr2dim))
    harr=0.0d0
    harr2=0.0d0

!----------------------------------------------------------------------
! (1) Off-diagonal elements between CSFs with the same spatial
!     configuration but different spin-couplings
!----------------------------------------------------------------------
    ! Loop over configurations
    do iconf=1,nconf

       ! Number of open shells
       nopen=sop_nopen(sop(:,:,iconf),n_int_in)
       
       ! Configuration and SOP in the full MO space
       conf_full=0_ib
       sop_full=0_ib
       conf_full(1:n_int_in,:)=conf(:,:,iconf)
       sop_full(1:n_int_in,:)=sop(:,:,iconf)
       
       ! Get all the configuration information needed to evaluate
       ! the on-diagonal Hamiltonian matrix element
       call package_confinfo(sop_full,conf_full,unocc,socc,docc,&
            nunocc,nsocc,ndocc,Dw,ndiff,nbefore)
       
       ! Number of CSFs generated by the configuration
       nsp=ncsfs(nopen)
          
       ! Off-diagonal Hamiltonian matrix elements
       ! < w omega' | H | w omega > between CSFs with the same
       ! spatial configuration but different spin-couplings
       if (nsp > 1) then

          ! Compute the Hamiltonian matrix elements
          arrdim=nsp*(nsp-1)/2
          call hij_same_mrci(harr(1:arrdim),arrdim,sop_full,socc,&
               nsocc,nbefore,Dw,ndiff,m2c)
          
          ! Save the above threshold matrix elements
          count=0
          ! Loop over ket CSFs
          do kcsf=offset(iconf),offset(iconf+1)-2
             ! Loop over bra CSFs
             do bcsf=kcsf+1,offset(iconf+1)-1
                count=count+1
                if (abs(harr(count)) > epshij) then
                   nbuf=nbuf+1
                   ibuffer(1,nbuf)=bcsf
                   ibuffer(2,nbuf)=kcsf
                   hbuffer(nbuf)=harr(count)
                   if (nbuf == bufsize) then
                      write(iscratch) hbuffer,ibuffer,nbuf
                      nbuf=0
                      nrec=nrec+1
                   endif
                endif
             enddo
          enddo
                
       endif
       
    enddo

!----------------------------------------------------------------------
! (2) Off-diagonal elements between CSFs with different spatial
!     configurations
!----------------------------------------------------------------------
    ! Loop over ket configurations
    do kconf=1,nconf-1
       
       ! Number of open shells in the ket configuration
       knopen=sop_nopen(sop(:,:,kconf),n_int_in)

       ! Number of ket CSFs
       knsp=ncsfs(knopen)
     
       ! Ket configuration and SOP in the full MO space
       kconf_full=0_ib
       ksop_full=0_ib
       kconf_full(1:n_int_in,:)=conf(:,:,kconf)
       ksop_full(1:n_int_in,:)=sop(:,:,kconf)
     
       ! Package the ket configuration information
       call package_confinfo(ksop_full,kconf_full,unocc,socc,docc,&
            nunocc,nsocc,ndocc,Dw,ndiff,nbefore)

       ! Loop over bra configurations
       do bconf=kconf+1,nconf
        
          ! Compute the excitation degree between the two configurations
          nexci=exc_degree_conf(conf(:,:,kconf),conf(:,:,bconf),n_int_in)

          ! Cycle if the excitation degree is greater than 2
          if (nexci > 2) cycle
        
          ! Number of open shells in the bra configuration
          bnopen=sop_nopen(sop(:,:,bconf),n_int_in)

          ! Number of bra CSFs
          bnsp=ncsfs(bnopen)
        
          ! Bra configuration and SOP in the full MO space
          bconf_full=0_ib
          bsop_full=0_ib
          bconf_full(1:n_int_in,:)=conf(:,:,bconf)
          bsop_full(1:n_int_in,:)=sop(:,:,bconf)
          
          ! Get the indices of the MOs involved in the excitation
          hlist=0
          plist=0
          call get_exci_indices(conf(:,:,kconf),conf(:,:,bconf),&
               n_int_in,hlist(1:nexci),plist(1:nexci),nexci)

          ! Compute the matrix elements between the CSFs generated
          ! by the bra and ket configurations
          call hij_mrci(harr2,harr2dim,nexci,bconf,kconf,&
               bsop_full,ksop_full,bnsp,knsp,bnopen,knopen,&
               hlist,plist,m2c,socc,nsocc,nbefore,Dw,ndiff,&
               offset,offset,nconf+1,nconf+1,averageii(bconf),&
               averageii(kconf))
        
          ! Save the above threshold matrix elements
          count=0
          do kcsf=offset(kconf),offset(kconf+1)-1
             do bcsf=offset(bconf),offset(bconf+1)-1
                count=count+1
                if (abs(harr2(count)) > epshij) then
                   nbuf=nbuf+1
                   ibuffer(1,nbuf)=bcsf
                   ibuffer(2,nbuf)=kcsf
                   hbuffer(nbuf)=harr2(count)
                   if (nbuf == bufsize) then
                      write(iscratch) hbuffer,ibuffer,nbuf
                      nbuf=0
                      nrec=nrec+1
                   endif
                endif
             enddo
          enddo
          
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Number of saved Hamiltonian matrix elements
!----------------------------------------------------------------------
    if (verbose) &
         write(6,'(/,x,a,x,i0)') 'Number of saved matrix elements:', &
         nrec*bufsize+nbuf
    
!----------------------------------------------------------------------
! Write the remaining elements in the buffer to disk
!----------------------------------------------------------------------
    if (nbuf > 0) then
       write(iscratch) hbuffer,ibuffer,nbuf
       nrec=nrec+1
    endif

!----------------------------------------------------------------------
! Close the Hamiltonian scratch file
!----------------------------------------------------------------------
    close(iscratch)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ibuffer)
    deallocate(hbuffer)
    deallocate(harr)
    deallocate(harr2)

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'save_hij_double')
    
    return
    
  end subroutine save_hij_double
    
!######################################################################
  
end module hbuild_double
