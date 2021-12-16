!**********************************************************************
! Generation of the MRCI configurations
!**********************************************************************

!######################################################################
! generate_mrci_confs: generates the int-ext representation of the
!                      MRCI configurations for a single irrep
!######################################################################
#ifdef CBINDING
subroutine generate_mrci_confs(irrep,nroots,conf0scr,confscr,nconf,&
     E0max1,icvs) bind(c,name="generate_mrci_confs")
#else
subroutine generate_mrci_confs(irrep,nroots,conf0scr,confscr,nconf,&
     E0max1,icvs)
#endif

  use constants
  use bitglobal
  use iomod
  use conftype
  use holeconfs
  use confbuilder
  use confinfo
  use timing

  implicit none
  
  ! Irrep number
  integer(is), intent(in)    :: irrep

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
  
  ! Number of reference space configurations for each irrep
  integer(is)                :: nconf0(0:nirrep-1)
  integer(is)                :: maxconf0

  ! Reference space configurations
  integer(ib), allocatable   :: conf0h(:,:,:,:)
  integer(ib), allocatable   :: sop0h(:,:,:,:)
  integer(is)                :: n_int_I,nmoI,nmoE

  ! 1H1I configurations
  integer(ib), allocatable   :: conf1h1I(:,:,:)
  integer(is)                :: n1h1I
  integer(is), allocatable   :: indx1h1I(:,:)

  ! Configuration derived data type
  type(mrcfg)                :: cfgM
  
  ! MO mapping arrays
  integer(is)                :: m2c(nmo),c2m(nmo)

  ! Active MO information
  integer(is)                :: nactive
  integer(is)                :: active(nmo)
  
  ! Scratch file variables
  integer(is)                :: iscratch
  character(len=60)          :: vecfile
  character(len=2)           :: amult,airrep

  ! Timing variables
  real(dp)                   :: tcpu_start,tcpu_end,twall_start,&
                                twall_end

  ! Everything else
  integer(is)                :: i,n,ntotal,counter


  ! Timing variables
  real(dp)                   :: tc1,tc2,tw1,tw2
  
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
  write(6,'(/,72a)') ('-',i=1,52)
  write(6,'(3(x,a))') 'MRCI configuration generation for the',&
       trim(irreplbl(irrep,ipg)),'subspace'
  write(6,'(72a)') ('-',i=1,52)

!----------------------------------------------------------------------
! Subtract E_SCF from the highest reference space energy to obtain the
! true MRCI Hamiltonian eigenvalue
!----------------------------------------------------------------------
  E0max=E0max1-Escf

!----------------------------------------------------------------------
! Read the reference configurations for all irreps from disk
!----------------------------------------------------------------------
  call read_ref_confs_all(conf0scr,nconf0,maxconf0,conf0h,sop0h,&
       n_int_I,nmoI,nmoE,m2c,c2m,nroots)

!----------------------------------------------------------------------
! Add the reference space information to the MRCI configuration derived
! data type
!----------------------------------------------------------------------
  ! Total number of ref confs across all irreps
  cfgM%nR=sum(nconf0(0:nirrep-1))

  ! Allocate arrays
  allocate(cfgM%conf0h(n_int_I,2,nconf0(irrep)))
  allocate(cfgM%sop0h(n_int_I,2,nconf0(irrep)))
  allocate(cfgM%confR(n_int_I,2,cfgM%nR))
  allocate(cfgM%sopR(n_int_I,2,cfgM%nR))
  allocate(cfgM%m2c(nmo))
  allocate(cfgM%c2m(nmo))

  ! Irrep number
  cfgM%irrep=irrep
  
  ! No. ref. confs of the current irrep
  cfgM%n0h=nconf0(irrep)

  ! MO subspace dimensions
  cfgM%n_int_I=n_int_I
  cfgM%nmoI=nmoI
  cfgM%nmoE=nmoE
  
  ! Ref confs and SOPs of the current irrep
  cfgM%conf0h=conf0h(:,:,1:nconf0(irrep),irrep)
  cfgM%sop0h=sop0h(:,:,1:nconf0(irrep),irrep)

  ! Ref confs and SOPs across all irreps
  counter=0
  do i=0,nirrep-1
     do n=1,nconf0(i)
        counter=counter+1
        cfgM%confR(:,:,counter)=conf0h(:,:,n,i)
        cfgM%sopR(:,:,counter)=sop0h(:,:,n,i)
     enddo
  enddo

  ! MO mapping arrays
  cfgM%m2c=m2c
  cfgM%c2m=c2m
  
!----------------------------------------------------------------------
! Generate the 1-hole, 2-hole configurations
!----------------------------------------------------------------------
  call generate_hole_confs(cfgM,icvs)
  
!----------------------------------------------------------------------
! Apply the internal creation operators to the 2-hole configurations,
! filtering out any 1-hole configurations that have already been
! generated from the reference configurations. This will yield the
! 1H1I configurations, from which the 2I and 1I1E configurations
! will be generated
!----------------------------------------------------------------------
  call generate_1hole_1I_confs(conf1h1I,n1h1I,indx1h1I,cfgM,icvs)

  print*,''
  print*,'n1h1I old:',n1h1I
  

  call generate_2I_1I1E_confs(conf1h1I,n1h1I,indx1h1I,cfgM,icvs,E0max)
  
  print*,''
  print*,'n1h1I new:',n1h1I
  
!----------------------------------------------------------------------
! Generate the configurations with one internal hole and one external
! electron
!----------------------------------------------------------------------

  call get_times(tw1,tc1)

  call generate_1E_confs(irrep,E0max,cfgM)

  call get_times(tw2,tc2)

  print*,'1E:',tw2-tw1
  
!----------------------------------------------------------------------
! Generate the configurations with two internal holes and two external
! electrons
!----------------------------------------------------------------------

  call get_times(tw1,tc1)
  
  call generate_2E_confs(irrep,E0max,cfgM)

  call get_times(tw2,tc2)

  print*,'2E:',tw2-tw1
  
!----------------------------------------------------------------------
! Generate the configurations with one internal hole and one internal
! electron
!----------------------------------------------------------------------

  call get_times(tw1,tc1)

  call generate_1I_confs(irrep,E0max,cfgM,icvs)

  call get_times(tw2,tc2)

  print*,'1I:',tw2-tw1
  
!----------------------------------------------------------------------
! Generate the configurations with two internal holes, one internal
! electron and one external electron
!----------------------------------------------------------------------

  call get_times(tw1,tc1)

  print*,''
  print*,'n1I1E new:',cfgM%n1I1E
  
  call generate_1I_1E_confs(irrep,E0max,conf1h1I,indx1h1I,&
       n_int_I,n1h1I,cfgM)

  print*,'n1I1E old:',cfgM%n1I1E
  
  call get_times(tw2,tc2)

  print*,'1I1E:',tw2-tw1

  
!----------------------------------------------------------------------
! Generate the configurations with two internal holes, two internal
! electrons
!----------------------------------------------------------------------

  call get_times(tw1,tc1)

  print*,''
  print*,'n2I new:',cfgM%n2I
  
  call generate_2I_confs(irrep,E0max,conf1h1I,indx1h1I,n_int_I,&
       n1h1I,cfgM,icvs)

  print*,'n2I old:',cfgM%n2I
  
  call get_times(tw2,tc2)

  print*,'2I:',tw2-tw1

  STOP
  
!----------------------------------------------------------------------
! Filter out any hole configurations which do not generate any
! full configurations
!----------------------------------------------------------------------
  call filter_hole_confs(cfgM)
  
!----------------------------------------------------------------------
! Set the total number of configurations
!----------------------------------------------------------------------
  ntotal=cfgM%n0h+cfgM%n1I+cfgM%n1E+cfgM%n2E+cfgM%n1I1E+cfgM%n2I
  nconf(irrep)=ntotal

!----------------------------------------------------------------------
! Write the configurations to disk
!----------------------------------------------------------------------
  ! Register the scratch file
  write(amult,'(i0)') imult
  write(airrep,'(i0)') irrep
  call scratch_name('mrciconf.mult'//trim(amult)//&
       '.sym'//trim(airrep),vecfile)
  call register_scratch_file(confscr(irrep),vecfile)

  ! Open the scratch file
  iscratch=scrunit(confscr(irrep))
  open(iscratch,file=scrname(confscr(irrep)),form='unformatted',&
       status='unknown')

  ! Subspace dimensions
  write(iscratch) n_int_I
  write(iscratch) nmoI
  write(iscratch) nmoE
  write(iscratch) ntotal
  write(iscratch) cfgM%nR
  write(iscratch) cfgM%n1h
  write(iscratch) cfgM%n2h
  write(iscratch) cfgM%n0h
  write(iscratch) cfgM%n1I
  write(iscratch) cfgM%n2I
  write(iscratch) cfgM%n1E
  write(iscratch) cfgM%n2E
  write(iscratch) cfgM%n1I1E
  
  ! Configuration information
  write(iscratch) cfgM%confR
  write(iscratch) cfgM%conf1h
  write(iscratch) cfgM%a1h
  write(iscratch) cfgM%off1h
  write(iscratch) cfgM%conf2h
  write(iscratch) cfgM%a2h
  write(iscratch) cfgM%off2h
  write(iscratch) cfgM%conf0h
  if (cfgM%n1I > 0) then
     write(iscratch) cfgM%a1I
     write(iscratch) cfgM%off1I
  endif
  if (cfgM%n2I > 0) then
     write(iscratch) cfgM%a2I
     write(iscratch) cfgM%off2I
  endif
  if (cfgM%n1E > 0) then
     write(iscratch) cfgM%a1E
     write(iscratch) cfgM%off1E
  endif
  if (cfgM%n2E > 0) then
     write(iscratch) cfgM%a2E
     write(iscratch) cfgM%off2E
  endif
  if (cfgM%n1I1E > 0 ) then
     write(iscratch) cfgM%a1I1E
     write(iscratch) cfgM%off1I1E
  endif

  ! MO mapping arrays
  write(iscratch) cfgM%m2c
  write(iscratch) cfgM%c2m
  
  ! Close the scratch file
  close(iscratch)

!----------------------------------------------------------------------
! Output the various configuration subspace dimensions
!----------------------------------------------------------------------
  ! Table of subspace dimensions
  write(6,'(/,x,20a)') ('-',i=1,20)
  write(6,'(2x,a)') 'Class | Nconf'
  write(6,'(x,20a)') ('-',i=1,20)
  write(6,'(2x,a,x,i0)') '1H    |',cfgM%n1h
  write(6,'(2x,a,x,i0)') '2H    |',cfgM%n2h
  write(6,'(2x,a,x,i0)') ' R    |',cfgM%n0h
  write(6,'(2x,a,x,i0)') ' I    |',cfgM%n1I
  write(6,'(2x,a,x,i0)') ' E    |',cfgM%n1E
  write(6,'(2x,a,x,i0)') 'II    |',cfgM%n2I
  write(6,'(2x,a,x,i0)') 'IE    |',cfgM%n1I1E
  write(6,'(2x,a,x,i0)') 'EE    |',cfgM%n2E
  write(6,'(x,20a)') ('-',i=1,20)

  ! Total number of configurations
  write(6,'(/,x,a,x,i0)') &
       'Total number of MRCI configurations:',ntotal
  
!----------------------------------------------------------------------
! Debugging: check for duplicate configurations
!----------------------------------------------------------------------
  !call check_confs(ntotal,n_int_I,nmoI,cfgM)
  
!----------------------------------------------------------------------
! Check on the number of 'active' MOs. That is, the number of variably
! occupied MOs across all configurations
!----------------------------------------------------------------------
  ! Construct the full configuration derived type including the SOPs,
  ! etc.
  call cfgM%finalise
  call cfgM%initialise(irrep,confscr(irrep))
  
  ! Determine the active MO indices
  call get_active_mos(cfgM,nactive,active)

  ! Deallocate the configuration derived type
  call cfgM%finalise

  ! Output the no. active MOs
  write(6,'(/,x,a,x,i0)') 'Number of active MOs:',nactive
  
!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
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
  

