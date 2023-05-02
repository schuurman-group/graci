!**********************************************************************
! Generation of the MRCI configurations
!**********************************************************************

!######################################################################
! generate_mrci_confs: generates the MRCI configurations for all irreps
!######################################################################
#ifdef CBINDING
subroutine generate_mrci_confs(nroots,conf0scr,confscr,nconf,E0max1,&
     icvs,ddci) bind(c,name="generate_mrci_confs")
#else
subroutine generate_mrci_confs(nroots,conf0scr,confscr,nconf,E0max1,&
       icvs,ddci)
#endif

  use constants
  use bitglobal
  use iomod
  use conftype
  use holeconfs
  use confbuilder
  use filter_confs
  use confinfo
  use timing

  implicit none
  
  ! Number of roots requested
  integer(is), intent(in)    :: nroots(0:nirrep-1)

  ! Reference space configuration scratch file numbers
  integer(is), intent(in)    :: conf0scr(0:nirrep-1)

  ! MRCI space configuration scratch file numbers
  integer(is), intent(inout) :: confscr(0:nirrep-1)

  ! Number of MRCI configurations
  integer(is), intent(inout) :: nconf(0:nirrep-1)
  
  ! Energy of the highest-lying reference space state of interest
  real(dp), intent(in)       :: E0max1
  real(dp)                   :: E0max

  ! CVS-MRCI: core MOs
  integer(is), intent(in)    :: icvs(nmo)

  ! DDCI configuration reduction
  logical, intent(in)        :: ddci
  
  ! Number of reference space configurations for each irrep
  integer(is)                :: nconf0(0:nirrep-1)
  integer(is)                :: maxconf0

  ! Reference space configurations
  integer(ib), allocatable   :: conf0h(:,:,:,:)
  integer(ib), allocatable   :: sop0h(:,:,:,:)
  integer(is)                :: n_int_I,nmoI,nmoE

  ! Configuration derived data types: one per irrep
  type(mrcfg), allocatable   :: cfgM(:)
  
  ! MO mapping arrays
  integer(is)                :: m2c(nmo),c2m(nmo)

  ! Active MO information
  integer(is)                :: nactive
  integer(is)                :: active(nmo),iact(nmo)
  
  ! Scratch file variables
  integer(is)                :: iscratch
  character(len=250)         :: conffile
  character(len=2)           :: amult,airrep

  ! Timing variables
  real(dp)                   :: tcpu_start,tcpu_end,twall_start,&
                                twall_end
  
  ! Everything else
  integer(is)                :: i,n,counter,irrep
  integer(is)                :: ntotal(0:nirrep-1)

!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  if (verbose) then
     write(6,'(/,72a)') ('-',i=1,52)
     write(6,'(3(x,a))') 'MRCI configuration generation for all irreps'
     write(6,'(72a)') ('-',i=1,52)
  endif
     
!----------------------------------------------------------------------
! Allocate the MRCI configuraytion derived types
!----------------------------------------------------------------------
  allocate(cfgM(0:nirrep-1))
  
!----------------------------------------------------------------------
! Subtract E_SCF from the highest reference space energy to obtain the
! true Hamiltonian eigenvalue
!----------------------------------------------------------------------
  E0max=E0max1-Escf

  if (verbose) write(6,'(/,x,a,x,F12.6,x,a)') &
       'Max ref space eigenvalue:',E0max*eh2ev,'eV'
  
!----------------------------------------------------------------------
! Read the reference configurations for all irreps from disk
!----------------------------------------------------------------------
  call read_ref_confs_all(conf0scr,nconf0,maxconf0,conf0h,sop0h,&
       n_int_I,nmo,nmoI,nmoE,m2c,c2m,nroots)

!----------------------------------------------------------------------
! Add the reference space information to the MRCI configuration derived
! data types
!----------------------------------------------------------------------
  ! Total number of ref confs across all irreps
  cfgM%nR=sum(nconf0(0:nirrep-1))
  
  ! Loop over irreps
  do irrep=0,nirrep-1

     ! Cycle if there are no roots for this irrep
     if (nroots(irrep) == 0) cycle
     
     ! Allocate arrays
     allocate(cfgM(irrep)%conf0h(n_int_I,2,nconf0(irrep)))
     allocate(cfgM(irrep)%sop0h(n_int_I,2,nconf0(irrep)))
     allocate(cfgM(irrep)%confR(n_int_I,2,cfgM(irrep)%nR))
     allocate(cfgM(irrep)%sopR(n_int_I,2,cfgM(irrep)%nR))
     allocate(cfgM(irrep)%m2c(nmo))
     allocate(cfgM(irrep)%c2m(nmo))

     ! Irrep number
     cfgM(irrep)%irrep=irrep
  
     ! No. ref. confs of the current irrep
     cfgM(irrep)%n0h=nconf0(irrep)

     ! MO subspace dimensions
     cfgM(irrep)%n_int_I=n_int_I
     cfgM(irrep)%nmoI=nmoI
     cfgM(irrep)%nmoE=nmoE
  
     ! Ref confs and SOPs of the current irrep
     cfgM(irrep)%conf0h=conf0h(:,:,1:nconf0(irrep),irrep)
     cfgM(irrep)%sop0h=sop0h(:,:,1:nconf0(irrep),irrep)

     ! Ref confs and SOPs across all irreps
     counter=0
     do i=0,nirrep-1
        do n=1,nconf0(i)
           counter=counter+1
           cfgM(irrep)%confR(:,:,counter)=conf0h(:,:,n,i)
           cfgM(irrep)%sopR(:,:,counter)=sop0h(:,:,n,i)
        enddo
     enddo

     ! MO mapping arrays
     cfgM(irrep)%m2c=m2c
     cfgM(irrep)%c2m=c2m

  enddo

!----------------------------------------------------------------------
! Generate the 1-hole, 2-hole configurations
!----------------------------------------------------------------------
  call generate_hole_confs(cfgM,icvs,nroots)
  
!----------------------------------------------------------------------
! Generate the configurations with one internal hole and one external
! electron
!----------------------------------------------------------------------
  call generate_1E_confs(E0max,cfgM,nroots)

  STOP
  
!----------------------------------------------------------------------
! Generate the configurations with two internal holes and two external
! electrons
!----------------------------------------------------------------------
  call generate_2E_confs(E0max,cfgM,ddci)
  
!----------------------------------------------------------------------
! Generate the configurations with one internal hole and one internal
! electron
!----------------------------------------------------------------------
  call generate_1I_confs(E0max,cfgM,icvs)

!----------------------------------------------------------------------
! Generate the configurations with:
! (a) Two internal holes and two internal electrons
! (b) Two internal holes, one internal electron and one external
!     electron
!----------------------------------------------------------------------
! Important: this routine must be called *after* the generation of
!            the 1I and 1E confs have been generated
!----------------------------------------------------------------------
  call generate_2I_1I1E_confs(E0max,cfgM,icvs,ddci)
  
!----------------------------------------------------------------------
! Filter out any hole configurations which do not generate any
! full configurations
!----------------------------------------------------------------------
  do irrep=0,nirrep-1
     call filter_hole_confs(cfgM(irrep))
  enddo

!----------------------------------------------------------------------
! Set the total number of configurations per irrep
!----------------------------------------------------------------------
  do irrep=0,nirrep-1
     ntotal(irrep)=cfgM(irrep)%n0h+cfgM(irrep)%n1I+cfgM(irrep)%n1E &
          +cfgM(irrep)%n2E+cfgM(irrep)%n1I1E+cfgM(irrep)%n2I
     nconf(irrep)=ntotal(irrep)
  enddo
     
!----------------------------------------------------------------------
! Write the configurations to disk
!----------------------------------------------------------------------
  ! Loop over irreps
  do irrep=0,nirrep-1

     ! Register the scratch file
     write(amult,'(i0)') imult
     write(airrep,'(i0)') irrep
     call scratch_name('mrciconf.mult'//trim(amult)//&
          '.sym'//trim(airrep),conffile)
     call register_scratch_file(confscr(irrep),conffile)
     
     ! Open the scratch file
     iscratch=scrunit(confscr(irrep))
     open(iscratch,file=scrname(confscr(irrep)),form='unformatted',&
          status='unknown')

     ! Subspace dimensions
     write(iscratch) n_int_I
     write(iscratch) nmoI
     write(iscratch) nmoE
     write(iscratch) ntotal(irrep)
     write(iscratch) cfgM(irrep)%nR
     write(iscratch) cfgM(irrep)%n1h
     write(iscratch) cfgM(irrep)%n2h
     write(iscratch) cfgM(irrep)%n0h
     write(iscratch) cfgM(irrep)%n1I
     write(iscratch) cfgM(irrep)%n2I
     write(iscratch) cfgM(irrep)%n1E
     write(iscratch) cfgM(irrep)%n2E
     write(iscratch) cfgM(irrep)%n1I1E
  
     ! Configuration information
     write(iscratch) cfgM(irrep)%confR
     write(iscratch) cfgM(irrep)%conf1h
     write(iscratch) cfgM(irrep)%a1h
     write(iscratch) cfgM(irrep)%off1h
     write(iscratch) cfgM(irrep)%conf2h
     write(iscratch) cfgM(irrep)%a2h
     write(iscratch) cfgM(irrep)%off2h
     write(iscratch) cfgM(irrep)%conf0h
     if (cfgM(irrep)%n1I > 0) then
        write(iscratch) cfgM(irrep)%a1I
        write(iscratch) cfgM(irrep)%off1I
     endif
     if (cfgM(irrep)%n2I > 0) then
        write(iscratch) cfgM(irrep)%a2I
        write(iscratch) cfgM(irrep)%off2I
     endif
     if (cfgM(irrep)%n1E > 0) then
        write(iscratch) cfgM(irrep)%a1E
        write(iscratch) cfgM(irrep)%off1E
     endif
     if (cfgM(irrep)%n2E > 0) then
        write(iscratch) cfgM(irrep)%a2E
        write(iscratch) cfgM(irrep)%off2E
     endif
     if (cfgM(irrep)%n1I1E > 0 ) then
        write(iscratch) cfgM(irrep)%a1I1E
        write(iscratch) cfgM(irrep)%off1I1E
     endif
     
     ! MO mapping arrays
     write(iscratch) cfgM(irrep)%m2c
     write(iscratch) cfgM(irrep)%c2m

     ! Number of CSFs as a function of the the number of open shells
     write(iscratch) ncsfs
     
     ! Close the scratch file
     close(iscratch)

  enddo
     
!----------------------------------------------------------------------
! Output the various configuration subspace dimensions
!----------------------------------------------------------------------
  if (verbose) then
     
     ! Loop over irreps
     do irrep=0,nirrep-1

        ! Table of subspace dimensions
        write(6,'(/,x,20a)') ('-',i=1,20)
        write(6,'(2x,a)') 'Irrep: '//trim(irreplbl(irrep,ipg))
        write(6,'(x,20a)') ('-',i=1,20)
        write(6,'(2x,a)') 'Class | Nconf'
        write(6,'(x,20a)') ('-',i=1,20)
        write(6,'(2x,a,x,i0)') '1H    |',cfgM(irrep)%n1h
        write(6,'(2x,a,x,i0)') '2H    |',cfgM(irrep)%n2h
        write(6,'(2x,a,x,i0)') ' R    |',cfgM(irrep)%n0h
        write(6,'(2x,a,x,i0)') ' I    |',cfgM(irrep)%n1I
        write(6,'(2x,a,x,i0)') ' E    |',cfgM(irrep)%n1E
        write(6,'(2x,a,x,i0)') 'II    |',cfgM(irrep)%n2I
        write(6,'(2x,a,x,i0)') 'IE    |',cfgM(irrep)%n1I1E
        write(6,'(2x,a,x,i0)') 'EE    |',cfgM(irrep)%n2E
        write(6,'(x,20a)') ('-',i=1,20)
        write(6,'(2x,a,x,i0)') 'Total |',ntotal(irrep)
        write(6,'(x,20a)') ('-',i=1,20)
        
     enddo

  endif
     
!----------------------------------------------------------------------
! Debugging: check for duplicate configurations
!----------------------------------------------------------------------
  !do irrep=0,nirrep-1
  !   call check_confs(ntotal(irrep),n_int_I,nmoI,cfgM(irrep))
  !enddo

!----------------------------------------------------------------------
! Check on the number of 'active' MOs. That is, the number of variably
! occupied MOs across all configurations
!----------------------------------------------------------------------
  ! Initialisation
  nactive=0
  iact=0

  ! Loop over irreps
  do irrep=0,nirrep-1

     ! Construct the full configuration derived type including the
     ! SOPs, etc.
     call cfgM(irrep)%finalise
     call cfgM(irrep)%initialise(irrep,confscr(irrep))
     
     ! Determine the active MO indices for this irrep
     call get_active_mos(cfgM(irrep),n,active)
     
     ! Update the total number of active MOs
     do i=1,n
        if (iact(active(i)) == 0) then
           nactive=nactive+1
           iact(active(i))=1
        endif
     enddo
          
     ! Deallocate the configuration derived type
     call cfgM(irrep)%finalise

  enddo
     
  ! Output the no. active MOs
  if (verbose) &
       write(6,'(/,x,a,x,i0)') 'Total number of active MOs:',nactive
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  if (verbose) &
       call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'generate_mrci_confs')

!----------------------------------------------------------------------
! Flush stdout
!----------------------------------------------------------------------
  flush(6)
  
  return

end subroutine generate_mrci_confs
  
!######################################################################
! check_confs: temporary debugging routine to check the generated
!              configurations for duplicates
!######################################################################
subroutine check_confs(ntotal,n_int_I,nmoI,cfgM)
  
  use constants
  use bitglobal
  use conftype
  use mrciutils
  use dethash
  use iomod
  
  implicit none

  ! Subspace dimensions
  integer(is), intent(in) :: ntotal,n_int_I,nmoI

  ! MRCI configurations
  type(mrcfg), intent(in) :: cfgM
    
  ! Full configurations
  integer(ib)             :: conftmp(n_int,2)
    
  ! Hash table
  type(dhtbl)             :: h
  integer(is)             :: initial_size
  integer(ib)             :: key(n_int,2)

  ! Everything else
  integer(is)             :: iconf,ioff,imo,imo1,imo2
  integer(is)             :: n,counter,k,i
  
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------
  ! Initialisation
  initial_size=ntotal
  call h%initialise_table(initial_size)

!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
  ! Loop over reference configurations
  do iconf=1,cfgM%n0h

     ! Insert the conf into the hash table
     key=0_ib
     key(1:n_int_I,:)=cfgM%conf0h(:,:,iconf)
     call h%insert_key(key)
     
  enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
  if (cfgM%n1I > 0) then
  
     ! Initialise the configuration counter
     counter=0

     ! Loop over 1-hole configurations
     do n=1,cfgM%n1h
        
        ! Loop over the 1I configurations generated by the current
        ! 1-hole configuration
        do ioff=cfgM%off1I(n),cfgM%off1I(n+1)-1

           ! Increment the configuration counter
           counter=counter+1

           ! Construct the configuration
           imo=cfgM%a1I(counter)
           conftmp=0_ib
           conftmp(1:n_int_I,:)=cfgM%conf1h(:,:,n)
           key=create_electron(conftmp,n_int,imo)
           
           ! Insert the configuration into the hash table
           call h%insert_key(key)
        
        enddo
     
     enddo

  endif
     
!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
  if (cfgM%n2I > 0) then

     ! Initialise the configuration counter
     counter=0

     ! Loop over 2-hole configurations
     do n=1,cfgM%n2h
     
        ! Loop over the 2I configurations generated by the current
        ! 2-hole configuration
        do ioff=cfgM%off2I(n),cfgM%off2I(n+1)-1
           
           ! Increment the configuration counter
           counter=counter+1

           ! Construct the configuration
           imo1=cfgM%a2I(1,counter)
           imo2=cfgM%a2I(2,counter)
           conftmp=0_ib
           conftmp(1:n_int_I,:)=cfgM%conf2h(:,:,n)
           key=create_electron(conftmp,n_int,imo1)
           conftmp=key
           key=create_electron(conftmp,n_int,imo2)
           
           ! Insert the configuration into the hash table
           call h%insert_key(key)
           
        enddo
        
     enddo

  endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
  if (cfgM%n1E > 0) then

     ! Initialise the configuration counter
     counter=0
  
     ! Loop over 1-hole configurations
     do n=1,cfgM%n1h

        ! Loop over the 1E configurations generated by the current
        ! 1-hole configuration
        do ioff=cfgM%off1E(n),cfgM%off1E(n+1)-1

           ! Increment the configuration counter
           counter=counter+1
           
           ! Construct the configuration
           imo=cfgM%a1E(counter)
           conftmp=0_ib
           conftmp(1:n_int_I,:)=cfgM%conf1h(:,:,n)
           key=create_electron(conftmp,n_int,imo)

           ! Insert the configuration into the hash table
           call h%insert_key(key)
        
        enddo
     
     enddo

  endif
     
!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
  if (cfgM%n2E > 0) then

     ! Initialise the configuration counter
     counter=0

     ! Loop over 2-hole configurations
     do n=1,cfgM%n2h
     
        ! Loop over the 2E configurations generated by the current
        ! 2-hole configuration
        do ioff=cfgM%off2E(n),cfgM%off2E(n+1)-1
        
           ! Increment the configuration counter
           counter=counter+1

           ! Construct the configuration
           imo1=cfgM%a2E(1,counter)
           imo2=cfgM%a2E(2,counter)
           conftmp=0_ib
           conftmp(1:n_int_I,:)=cfgM%conf2h(:,:,n)
           key=create_electron(conftmp,n_int,imo1)
           conftmp=key
           key=create_electron(conftmp,n_int,imo2)

           ! Insert the configuration into the hash table
           call h%insert_key(key)

        enddo
        
     enddo

  endif
  
!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
  if (cfgM%n1I1E > 0) then

     ! Initialise the configuration counter
     counter=0

     ! Loop over 2-hole configurations
     do n=1,cfgM%n2h
     
        ! Loop over the 1I1E configurations generated by the
        ! current 2-hole configuration
        do ioff=cfgM%off1I1E(n),cfgM%off1I1E(n+1)-1
        
           ! Increment the configuration counter
           counter=counter+1

           ! Construct the configuration
           imo1=cfgM%a1I1E(1,counter)
           imo2=cfgM%a1I1E(2,counter)
           conftmp=0_ib
           conftmp(1:n_int_I,:)=cfgM%conf2h(:,:,n)
           key=create_electron(conftmp,n_int,imo1)
           conftmp=key
           key=create_electron(conftmp,n_int,imo2)

           ! Insert the configuration into the hash table
           call h%insert_key(key)

        enddo
        
     enddo

  endif
  
!----------------------------------------------------------------------
! Did we generate any duplicate configurations?
!----------------------------------------------------------------------
  if (h%n_keys_stored /= ntotal) then
     errmsg='Duplicate configurations found in check_confs'
     call error_control
  endif
    
!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
  call h%delete_table
    
  return
    
end subroutine check_confs
    
!######################################################################
  

