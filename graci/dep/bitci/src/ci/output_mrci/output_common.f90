!**********************************************************************
! Useful common MRCI output routines
!**********************************************************************
module output_common

  implicit none

contains

!######################################################################
! sort_energies: returns irrep/state number pairs in order of
!                increasing energy
!######################################################################
  subroutine sort_energies(ener,nroot_max,nroot_tot,nroots,sindx)

    use constants
    use bitglobal
    
    implicit none

    ! Dimensions
    integer(is), intent(in)  :: nroot_max,nroot_tot
    integer(is), intent(in)  :: nroots(0:nirrep-1)
    
    ! Energies
    real(dp), intent(in)     :: ener(nroot_max,0:nirrep-1)

    ! Mapping array
    integer(is), intent(out) :: sindx(nroot_tot,2)

    ! Everything else
    integer(is)              :: irrep,n,k
    real(dp), allocatable    :: etmp(:,:)
    real(dp)                 :: minval
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(etmp(nroot_max,0:nirrep-1))
    etmp=0.0d0
    
!----------------------------------------------------------------------
! Fill in the mapping array
! This is a pretty dumb approach, but speed doesn't matter here
!----------------------------------------------------------------------
    ! Initialisation
    etmp=ener

    ! Loop over the total number of states
    do n=1,nroot_tot

       minval=1e+90_dp
       
       ! Loop over irreps
       do irrep=0,nirrep-1
          
          ! Loop over roots for this irrep
          do k=1,nroots(irrep)

             if (etmp(k,irrep) < minval) then
                minval=etmp(k,irrep)
                sindx(n,2)=irrep
                sindx(n,1)=k
             endif
          
          enddo
          
       enddo

       etmp(sindx(n,1),sindx(n,2))=1e+90_dp
       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(etmp)
    
    return
    
  end subroutine sort_energies

!######################################################################
! dominant_csfs: for a given eigenvector, returns the configuration
!                and spin-coupling information corresponding to the
!                dominant CSFs
!######################################################################
  subroutine dominant_csfs(modus,csfthrsh,cfg,vecscr,k,ndom,refnorm,&
       dconf,domega,dcoe)

    use constants
    use bitglobal
    use conftype
    use iomod
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of dominant
    !                                CSFs
    !                    modus=1 <-> package the configuration and
    !                                spin-coupling information
    integer(is), intent(in)    :: modus
    
    ! CSF coefficient threshold
    real(dp), intent(in)       :: csfthrsh

    ! MRCI configuration derived type
    type(mrcfg), intent(in)    :: cfg

    ! Eigenpair scratch file number
    integer(is), intent(in)    :: vecscr

    ! Root number
    integer(is), intent(in)    :: k

    ! No. dominant CSFs
    integer(is), intent(inout) :: ndom
    
    ! Norm of the wavefunction projected onto the ref space
    real(dp), intent(out)      :: refnorm

    ! Dominant configurations and spin-couplings
    integer(ib), allocatable   :: dconf(:,:,:)
    integer(is), allocatable   :: domega(:)
    real(dp), allocatable      :: dcoe(:)
    
    ! Eigenvector
    real(dp), allocatable      :: vec(:)
    
    ! Everything else
    integer(is)                :: csfdim,refdim,n,ioff,csf,counter
    integer(is)                :: omega
    integer(is)                :: n_int_I
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    csfdim=cfg%csfdim
    allocate(vec(csfdim))
    
!----------------------------------------------------------------------
! Read in the eigenvector
!----------------------------------------------------------------------
    call read_single_vector(vecscr,vec,csfdim,k)

!----------------------------------------------------------------------
! Norm of the wavefunction projected onto the reference space
!----------------------------------------------------------------------
    refdim=cfg%csfs0h(cfg%n0h+1)-1
    refnorm=sqrt(dot_product(vec(1:refdim),vec(1:refdim)))

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    counter=0
    n_int_I=cfg%n_int_I
    
    if (modus == 1) then
       allocate(dconf(n_int,2,ndom))
       allocate(domega(ndom))
       allocate(dcoe(ndom))
       dconf=0_ib
       domega=0
       dcoe=0.0d0
    endif
       
!----------------------------------------------------------------------
! Reference space configurations
!----------------------------------------------------------------------
    ! Loop over ref space configurations
    do n=1,cfg%n0h

       ! Loop over the CSFs generated by this configuration
       omega=0
       do csf=cfg%csfs0h(n),cfg%csfs0h(n+1)-1

          ! Spin-coupling index
          omega=omega+1
          
          ! Cycle if this is not a dominant CSF
          if (vec(csf)**2 < csfthrsh) cycle
          
          ! Increment the dominant CSF counter
          counter=counter+1

          ! Cycle if we are just determinig the no. dominant CSFs
          if (modus == 0) cycle

          ! Save the CSF information
          dconf(1:n_int_I,:,counter)=cfg%conf0h(:,:,n)
          domega(counter)=omega
          dcoe(counter)=vec(csf)
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! 1I configurations
!----------------------------------------------------------------------
    if (cfg%n1I > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1I configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1I(n),cfg%off1I(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs1I(ioff),cfg%csfs1I(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf1I(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo

          enddo
          
       enddo

    endif

!----------------------------------------------------------------------
! 2I configurations
!----------------------------------------------------------------------
    if (cfg%n2I > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2I configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2I(n),cfg%off2I(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs2I(ioff),cfg%csfs2I(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf2I(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1E configurations
!----------------------------------------------------------------------
    if (cfg%n1E > 0) then

       ! Loop over 1-hole configurations
       do n=1,cfg%n1h
          
          ! Loop over the 1E configurations generated by the 1-hole
          ! configuration
          do ioff=cfg%off1E(n),cfg%off1E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs1E(ioff),cfg%csfs1E(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf1E(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 2E configurations
!----------------------------------------------------------------------
    if (cfg%n2E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 2E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off2E(n),cfg%off2E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs2E(ioff),cfg%csfs2E(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf2E(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! 1I1E configurations
!----------------------------------------------------------------------
    if (cfg%n1I1E > 0) then

       ! Loop over 2-hole configurations
       do n=1,cfg%n2h
          
          ! Loop over the 1I1E configurations generated by the 2-hole
          ! configuration
          do ioff=cfg%off1I1E(n),cfg%off1I1E(n+1)-1
             
             ! Loop over the CSFs generated by this configuration
             omega=0
             do csf=cfg%csfs1I1E(ioff),cfg%csfs1I1E(ioff+1)-1

                ! Spin-coupling index
                omega=omega+1
                
                ! Cycle if this is not a dominant CSF
                if (vec(csf)**2 < csfthrsh) cycle
                
                ! Increment the dominant CSF counter
                counter=counter+1
                
                ! Cycle if we are just determinig the no. dominant CSFs
                if (modus == 0) cycle

                ! Save the CSF information
                dconf(:,:,counter)=cfg%conf1I1E(:,:,ioff)
                domega(counter)=omega
                dcoe(counter)=vec(csf)
                
             enddo
             
          enddo

       enddo

    endif

!----------------------------------------------------------------------
! Set the number of dominant CSFs
!----------------------------------------------------------------------
    ndom=counter
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(vec)
    
    return
    
  end subroutine dominant_csfs

!######################################################################
! exci_to_string: given lists of hole and particle indices, returns
!                 a character string of the form i, j,... -> a, b,...
!######################################################################
  subroutine exci_to_string(string,hlist,plist,nexci,m2c)

    use constants
    use bitglobal
    
    implicit none

    ! Creation and annihilation operator indices
    integer(is), intent(in) :: nexci
    integer(is), intent(in) :: hlist(nexci),plist(nexci)

    ! MO mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Ouput character string
    character(len=*)        :: string

    ! Everything else
    integer(is)             :: i,k

    !
    ! Initialisation
    !
    string=''
    k=1

    !
    ! Base configuration case
    !
    if (nexci == 0) then
       string='base conf'
       return
    endif
    
    !
    ! Annihilation operators
    !
    do i=1,nexci
       write(string(k:),'(x,i0)') m2c(hlist(i))
       k=len_trim(string)+1
    enddo

    !
    ! Creation operators
    !
    write(string(k:),'(x,a)') '->'
    k=len_trim(string)+1
    do i=1,nexci
       write(string(k:),'(x,i0)') m2c(plist(i))
       k=len_trim(string)+1
    enddo

    return
    
  end subroutine exci_to_string
  
!######################################################################
! deadwood: deadwood analysis
!######################################################################
  subroutine deadwood(vecscr,nroots,cfg)

    use constants
    use bitglobal
    use conftype
    use iomod
    use utils
    
    implicit none

    ! MRCI eigenpair scratch file numbers
    integer(is), intent(in)  :: vecscr(0:nirrep-1)
    
    ! Number of roots per irrep
    integer(is), intent(in)  :: nroots(0:nirrep-1)

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg(0:nirrep-1)
    
    ! Eigenvectors
    real(dp), allocatable    :: vec(:)

    ! Everything else
    integer(is)              :: irrep,k,i,csfdim,unit
    integer(is), allocatable :: indx(:)
    real(dp)                 :: normsq
    character(len=255)       :: filename
    
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! No. CSFs
       csfdim=cfg(irrep)%csfdim
       
       ! Allocate arrays
       allocate(vec(csfdim),indx(csfdim))
       
       ! Loop over roots
       do k=1,nroots(irrep)

          ! Open the output file
          write(filename,'(a,i0,a)') &
               'normsq_',k,trim(irreplbl(irrep,ipg))//'.dat'
          call freeunit(unit)
          open(unit,file=filename,form='formatted',status='unknown')
          
          ! Read in the eigenvector
          call read_single_vector(vecscr(irrep),vec,csfdim,k)
       
          ! Sort the coefficients
          vec=abs(vec)
          call dsortindxa1('D',csfdim,vec,indx)

          ! Squared norm of the truncated wavefunction
          normsq=0.0d0
          do i=1,csfdim
             normsq=normsq+vec(indx(i))**2
             write(unit,'(i0,x,ES10.4)') i,normsq
             if (normsq > 0.9999d0) exit
          enddo

          ! Close the output file
          close(unit)
          
       enddo
          
       ! Deallocate arrays
       deallocate(vec,indx)
       
    enddo
    
    return
    
  end subroutine deadwood
  
!######################################################################
! print_base_conf: prints the base configuration
!######################################################################
  subroutine print_base_conf

    use constants
    use bitglobal
    use detutils
    use iomod
    
    implicit none

    integer(is)                                 :: ihomo,imo,k,i,occ
    integer(is)                                 :: ind
    character(len=60)                           :: string1,string2
    character(len=1), parameter, dimension(0:2) :: aocc= ['0', '1', '2']
    
    !
    ! Get the idex of the HOMO
    !
    ihomo=homo_index(conf0)

    !
    ! Write the full base conf character string
    !
    string1=''
    ind=0
    do imo=1,ihomo

       ! Block index
       k=(imo-1)/n_bits+1
       
       ! Orbital position with the block
       i=imo-1-(k-1)*n_bits

       ! Orbital occupancy
       occ=0
       if (btest(conf0(k,1),i)) occ=occ+1
       if (btest(conf0(k,2),i)) occ=occ+1

       ! Add the orbital occupancy to the character string
       write(string1(imo:imo), '(a1)') aocc(occ)

       ! Save the position of the first non-doubly-occupied MO
       if (occ /= 2 .and. ind == 0) ind=imo
       
    enddo
    
    if (ind == 0) ind=ihomo+1
    
    !
    ! Write the truncated base conf character string
    !
    string2=''

    ! String of doubly-occupied MOs
    write(string2(1:16),'(a)') 'base conf: |2...'
    if (ind /= 2) then
       write(string2(17:17),'(a)') '2'
       k=18
    else
       k=17
    endif

    ! Remaining MOs
    do imo=ind,ihomo
       write(string2(k:k),'(a)') string1(imo:imo)
       k=k+1
    enddo

    ! Ket end
    write(string2(k:k),'(a)') '>'

    !
    ! Print the base conf
    !
    write(6,'(2/,2x,a)') trim(string2)

    !
    ! Print the last doubly-occupied MO index
    !
    string1=''
    write(string1(1:17),'(16a,a)') (' ', k=1,16), achar(094)
    write(6,'(2x,a)') trim(string1)
    string1=''
    write(string1(1:16),'(16a)') (' ', k=1,16)
    write(string1(17:),'(i0)') ind-1
    write(6,'(2x,a)') trim(string1)
    
    return
    
  end subroutine print_base_conf
  
!######################################################################
 
end module output_common
