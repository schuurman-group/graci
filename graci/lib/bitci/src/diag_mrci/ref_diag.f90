!**********************************************************************
! Routines for the diagonalisation if the MRCI reference space
! Hamiltonian
!**********************************************************************

!######################################################################
! ref_diag_mrci: Diagonalisation of the reference space Hamiltonian.
!                For now, we will just assume that the Hamiltonian
!                matrix is small enough to be amenable to full
!                diagonalisation.
!######################################################################
#ifdef CBINDING
subroutine ref_diag_mrci(irrep,nroots,confscr,nconf,vecscr) &
     bind(c,name="ref_diag_mrci")
#else
subroutine ref_diag_mrci(irrep,nroots,confscr,nconf,vecscr)
#endif

  use constants
  use bitglobal
  use utils
  use iomod
  
  implicit none

  ! Irrep number and the requested number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(inout) :: nroots
  
  ! Array of reference configuration scratch file numbers
  integer(is), intent(in)    :: confscr(0:nirrep-1)

  ! Array of numbers of referecne configurations
  integer(is), intent(in)    :: nconf(0:nirrep-1)

  ! Eigenvector scratch file index
  integer(is), intent(out)   :: vecscr
  
  ! Number of configurations found in the scratch file
  integer(is)                :: nconf1

  ! Number of 64-bit integers required to represent the configurations
  integer(is)                :: n_int_I

  ! Number of internal and external MOs
  integer(is)                :: nmoI,nmoE

  ! Reference configurations and SOPs
  integer(ib), allocatable   :: conf(:,:,:),sop(:,:,:)
  
  ! MO mapping arrays
  integer(is)                :: m2c(nmo),c2m(nmo)

  ! Hamiltonian matrix
  real(dp), allocatable      :: hmat(:,:)

  ! Eigenpairs
  real(dp), allocatable      :: heig(:),hvec(:,:)
  
  ! Numbers of CSFs
  integer(is)                :: hdim
  integer(is), allocatable   :: offset(:)

  ! Everything else
  integer(is)                :: i

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,72a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'Reference space diagonalisation in the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(72a)') ('-',i=1,52)

!----------------------------------------------------------------------
! Return if there are no configurations for the current irrep
!----------------------------------------------------------------------
  if (nconf(irrep) == 0) then
     write(6,'(/,x,a)') 'No reference space configurations of '&
          //trim(irreplbl(irrep,ipg))//' symmetry'
     nroots=0
     return
  endif

!----------------------------------------------------------------------
! Read the configurations from file
!----------------------------------------------------------------------
  call read_ref_confs(confscr(irrep),nconf1,n_int_I,nmoI,nmoE,conf,&
       sop,m2c,c2m)

!----------------------------------------------------------------------
! Sanity check on the number of configurations
!----------------------------------------------------------------------
  if (nconf(irrep) /= nconf1) then
     errmsg='Error in ref_diag_mrci: '&
          //'inconsistent configuration numbers'
     call error_control
  endif

!----------------------------------------------------------------------
! Sanity check on the number of roots
!----------------------------------------------------------------------
  if (nconf(irrep) < nroots) then
     errmsg='Error in ref_diag_mrci: nroots > nconf'
     call error_control
  endif
  
!----------------------------------------------------------------------
! Determine the total number of CSFs and the offsets for each
! configuration
!----------------------------------------------------------------------
  allocate(offset(nconf1+1))
  offset=0

  call basis_dimensions(hdim,offset,sop,n_int_I,nconf1)
  
!----------------------------------------------------------------------
! Construct the reference space Hamiltonian matrix
!----------------------------------------------------------------------
  allocate(hmat(hdim,hdim))
  hmat=0.0d0

  call build_href(nconf1,hdim,hmat,offset,conf,sop,n_int_I,m2c,nmoI,&
       irrep)

!----------------------------------------------------------------------
! Diagonalise the reference space Hamiltonian matrix
!----------------------------------------------------------------------
  ! Allocate arrays
  allocate(heig(hdim))
  allocate(hvec(hdim,hdim))
  heig=0.0d0
  hvec=0.0d0

  ! Diagonalisation
  call diag_matrix_real(hmat,heig,hvec,hdim)

  ! Add on E_SCF + E_nuc
  heig=heig+escf+enuc
  
!----------------------------------------------------------------------
! Save the eigenpairs to disk
!----------------------------------------------------------------------
  call write_ref_eigen(vecscr,irrep,nroots,hdim,heig,hvec)
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(offset)
  deallocate(hmat)
  deallocate(heig)
  deallocate(hvec)

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return
  
end subroutine ref_diag_mrci

!######################################################################
! basis_dimensions: Determines: (1) the total number of CSFs
!                   generated by the reference space configurations,
!                   (hdim) and; (2) the starting points for the CSFs
!                   generated by each configuration (offset)
!######################################################################
subroutine basis_dimensions(hdim,offset,sop,n_int_I,nconf)

  use constants
  use bitglobal
  
  implicit none

  integer(is), intent(out) :: hdim
  integer(is), intent(in)  :: n_int_I,nconf
  integer(is), intent(out) :: offset(nconf+1)
  integer(ib), intent(in)  :: sop(n_int_I,2,nconf)
  integer(is), allocatable :: nopen(:)
  integer(is)              :: ndet
  integer(is)              :: i,k,sum

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(nopen(nconf))
  
!----------------------------------------------------------------------
! Number of open shells for each configuration
!----------------------------------------------------------------------
  ! Initialisation
  nopen=0

  ! Loop over configurations
  do i=1,nconf

     ! Number of open shells
     do k=1,n_int_I
        nopen(i)=nopen(i)+popcnt(sop(k,1,i))
     enddo
     
  enddo
  
!----------------------------------------------------------------------
! Total number of CSFs and determinants
!----------------------------------------------------------------------
  ! Initialisation
  hdim=0
  ndet=0
  
  ! Loop over configurations
  do i=1,nconf

     ! Number of CSFs generated by the current configuration
     hdim=hdim+ncsfs(nopen(i))

     ! Number of determinants generated by the current configuration
     ndet=ndet+ndets(nopen(i))
     
  enddo

!----------------------------------------------------------------------
! Offsets
!----------------------------------------------------------------------
  sum=1
  offset=0

  do i=1,nconf
     offset(i)=sum
     sum=sum+ncsfs(nopen(i))
  enddo

  offset(nconf+1)=hdim+1

!----------------------------------------------------------------------
! Output the reference space dimensions
!----------------------------------------------------------------------
  write(6,'(/,x,a,x,i0)') &
       'Number of configurations in the reference space:',nconf
  write(6,'(x,a,x,i0)') &
       'Number of CSFs in the reference space:',hdim
  write(6,'(x,a,x,i0)') &
       'Number of determinants in the reference space:',ndet

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(nopen)
  
  return
  
end subroutine basis_dimensions
  
!######################################################################
! build_href: Constructs the reference space Hamiltonian
!######################################################################
subroutine build_href(nconf,hdim,hmat,offset,conf,sop,n_int_I,m2c,nmoI,&
     irrep)

  use constants
  use bitglobal
  use mrciutils
  use bitutils
  use hparam
  use hbuild_mrci
  use mrci_integrals
  use dftmrci
  use timing
  
  implicit none

  ! Irrep number
  integer(is), intent(in) :: irrep
  
  ! Dimensions
  integer(is), intent(in) :: nconf,hdim,nmoI
  integer(is), intent(in) :: offset(nconf+1)
  
  ! Hamiltonian matrix
  real(dp), intent(out)   :: hmat(hdim,hdim)

  ! Configurations and SOPs
  integer(is), intent(in) :: n_int_I
  integer(ib), intent(in) :: conf(n_int_I,2,nconf)
  integer(ib), intent(in) :: sop(n_int_I,2,nconf)

  ! MO mapping array
  integer(is), intent(in) :: m2c(nmo)
  
  ! Hamiltonian build variables
  integer(is)             :: nexci
  integer(is)             :: socc(nmo),docc(nmo),unocc(nmo)
  integer(is)             :: nopen,nsocc,ndocc,nunocc
  integer(ib)             :: conf_full(n_int,2)
  integer(ib)             :: sop_full(n_int,2)
  integer(is)             :: knopen,bnopen
  integer(ib)             :: bconf_full(n_int,2),kconf_full(n_int,2)
  integer(ib)             :: bsop_full(n_int,2),ksop_full(n_int,2)
  integer(is), parameter  :: maxexci=2
  integer(is)             :: hlist(maxexci),plist(maxexci)
  real(dp), allocatable   :: averageii(:)

  ! Difference configuration information
  integer(is)             :: ndiff
  integer(is)             :: Dw(nmo,2)

  ! Number of open shells preceding each MO
  integer(is)             :: nbefore(nmo)
  
  ! Temporary Hij array
  integer(is)             :: arrdim
  real(dp), allocatable   :: harr(:),harr2(:,:)

  ! Everything else
  integer(is)             :: iconf,kconf,bconf
  integer(is)             :: icsf,kcsf,bcsf
  integer(is)             :: iomega,komega,bomega
  integer(is)             :: nsp,count,count1,count2,blim1,blim2,&
                             klim1,klim2,bnsp,knsp
  
  ! Timing variables
  real(dp)                :: tcpu_start,tcpu_end,twall_start,&
                             twall_end
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  allocate(harr(ncsfs(nomax)**2))
  harr=0.0d0

  allocate(harr2(ncsfs(nomax),ncsfs(nomax)))
  harr2=0.0d0

  allocate(averageii(nconf))
  averageii=0.0d0

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
  hmat=0.0d0
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! (1) On-diagonal elements and off-diagonal elements between CSFs with
!     the same spatial configuration but different spin-couplings
!----------------------------------------------------------------------
  ! Loop over configurations
  do iconf=1,nconf

     ! Number of open shells
     nopen=sop_nopen(sop(:,:,iconf),n_int_I)
  
     ! Configuration and SOP in the full MO space
     conf_full=0_ib
     sop_full=0_ib
     conf_full(1:n_int_I,:)=conf(:,:,iconf)
     sop_full(1:n_int_I,:)=sop(:,:,iconf)
  
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
          nopen,m2c)
     
     ! Fill in the hmat array with the on-diagonal elements
     iomega=0
     do icsf=offset(iconf),offset(iconf+1)-1
        iomega=iomega+1
        hmat(icsf,icsf)=harr(iomega)
     enddo
     
     ! Off-diagonal Hamiltonian matrix elements
     ! < w omega' | H | w omega > between CSFs with the same spatial
     ! configuration but different spin-couplings
     if (nsp > 1) then

        ! Compute the Hamiltonian matrix elements
        arrdim=nsp*(nsp-1)/2
        call hij_same_mrci(harr(1:arrdim),arrdim,sop_full,socc,nsocc,&
             nbefore,m2c)
        
        ! Fill in the hmat array with the off-diagonal elements
        count=0
        ! Loop over ket CSFs
        do kcsf=offset(iconf),offset(iconf+1)-2
           ! Loop over bra CSFs
           do bcsf=kcsf+1,offset(iconf+1)-1
              count=count+1
              hmat(bcsf,kcsf)=harr(count)
              hmat(kcsf,bcsf)=harr(count)
           enddo
        enddo
        
     endif
        
  enddo
    
!----------------------------------------------------------------------
! Pre-compute the average values of the on-diagonal matrix elements
! generated by each configuration
!----------------------------------------------------------------------
  ! Loop over configurations
  do iconf=1,nconf

     ! Number of configurations generated by the current configuration
     nsp=offset(iconf+1)-offset(iconf)

     ! Average on-diagonal Hamiltonian matrix element value
     do icsf=offset(iconf),offset(iconf+1)-1
        averageii(iconf)=averageii(iconf)+hmat(icsf,icsf)
     enddo
     averageii(iconf)=averageii(iconf)/nsp

  enddo
  
!----------------------------------------------------------------------
! (2) Off-diagonal elements between CSFs with different spatial
!     configurations
!----------------------------------------------------------------------
  ! Loop over ket configurations
  do kconf=1,nconf-1

     ! Number of open shells in the ket configuration
     knopen=sop_nopen(sop(:,:,kconf),n_int_I)

     ! Number of ket CSFs
     knsp=ncsfs(knopen)
     
     ! Ket configuration and SOP in the full MO space
     kconf_full=0_ib
     ksop_full=0_ib
     kconf_full(1:n_int_I,:)=conf(:,:,kconf)
     ksop_full(1:n_int_I,:)=sop(:,:,kconf)
     
     ! Package the ket configuration information
     call package_confinfo(ksop_full,kconf_full,unocc,socc,docc,&
          nunocc,nsocc,ndocc,Dw,ndiff,nbefore)

     ! Loop over bra configurations
     do bconf=kconf+1,nconf
        
        ! Compute the excitation degree between the two configurations
        nexci=exc_degree_conf(conf(:,:,kconf),conf(:,:,bconf),n_int_I)

        ! Cycle if the excitation degree is greater than 2
        if (nexci > 2) cycle
        
        ! Number of open shells in the bra configuration
        bnopen=sop_nopen(sop(:,:,bconf),n_int_I)

        ! Number of bra CSFs
        bnsp=ncsfs(bnopen)
        
        ! Bra configuration and SOP in the full MO space
        bconf_full=0_ib
        bsop_full=0_ib
        bconf_full(1:n_int_I,:)=conf(:,:,bconf)
        bsop_full(1:n_int_I,:)=sop(:,:,bconf)

        ! Get the indices of the MOs involved in the excitation
        hlist=0
        plist=0
        call get_exci_indices(conf(:,:,kconf),conf(:,:,bconf),&
             n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

        ! Compute the matrix elements between the CSFs generated
        ! by the bra and ket configurations
        call hij_mrci(harr2,ncsfs(nomax),nexci,bconf,kconf,&
             bsop_full,ksop_full,bnsp,knsp,bnopen,knopen,&
             hlist,plist,m2c,socc,nsocc,nbefore,Dw,ndiff,&
             offset,offset,nconf+1,nconf+1,averageii(bconf),&
             averageii(kconf))
        
        ! Fill in the hmat array
        count1=0
        do kcsf=offset(kconf),offset(kconf+1)-1
           count1=count1+1
           count2=0
           do bcsf=offset(bconf),offset(bconf+1)-1
              count2=count2+1
              hmat(bcsf,kcsf)=harr2(count2,count1)
              hmat(kcsf,bcsf)=hmat(bcsf,kcsf)
           enddo
        enddo

     enddo

  enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(harr)
  deallocate(harr2)
  deallocate(averageii)
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'build_href')
  
  return

end subroutine build_href
  
!######################################################################
! write_ref_eigen: writes the reference space eigenpairs to disk
!######################################################################
subroutine write_ref_eigen(vecscr,irrep,nroots,hdim,heig,hvec)

  use constants
  use bitglobal
  use iomod
  
  implicit none

  ! Eigenpair scratch file index
  integer(is), intent(out)   :: vecscr

  ! Irrep number and the requested number of roots
  integer(is), intent(in)    :: irrep
  integer(is), intent(inout) :: nroots

  ! Hamiltonian dimension (no. CSFs)
  integer(is), intent(in)    :: hdim

  ! Hamiltonian eigenpairs
  real(dp), intent(in)       :: heig(hdim),hvec(hdim,hdim)
  
  ! I/O variables
  integer(is)                :: iscratch
  character(len=60)          :: vecfile
  character(len=2)           :: amult,airrep

  ! Everything else
  integer(is)                :: i
  
!----------------------------------------------------------------------
! Register the scratch file
!----------------------------------------------------------------------
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('refvec.mult'//trim(amult)//&
       '.sym'//trim(airrep),vecfile)
  call register_scratch_file(vecscr,vecfile)

!----------------------------------------------------------------------
! Open the scratch file
!----------------------------------------------------------------------
  iscratch=scrunit(vecscr)
  open(iscratch,file=scrname(vecscr),form='unformatted',&
       status='unknown')

!----------------------------------------------------------------------
! Write the eigenpairs to disk
!----------------------------------------------------------------------
  ! No. CSFs
  write(iscratch) hdim

  ! No. roots
  write(iscratch) nroots

  ! Eigenvalues
  write(iscratch) heig(1:nroots)

  ! Eigenvectors
  do i=1,nroots
     write(iscratch) hvec(:,i)
  enddo
     
!----------------------------------------------------------------------
! Close the scratch file
!----------------------------------------------------------------------
  close(iscratch)
  
  return
    
end subroutine write_ref_eigen

