!**********************************************************************
! Routines for the generation of the various classes of configurations
! entering into an MRCI calculation from the 1-hole and 2-hole
! configurations previously generated
!**********************************************************************
module confbuilder

  implicit none
  
contains

!######################################################################
! generate_2I_1I1E_confs: for all irreps, generates all the allowable
!                         2I and 1I1E confs
!######################################################################
  subroutine generate_2I_1I1E_confs(E0max,cfgM,icvs)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use iomod

    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max
    
    ! MRCI configurations for all irreps
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs

    ! No. 1H1I, 2I and 1I1E configurations
    integer(is)                :: n1h1I,n2I,n1I1E
    
    ! Configuration scratch files
    integer(is)                :: nrec2I,nrec1I1E
    character(len=60)          :: file2I,file1I1E
    
!----------------------------------------------------------------------
! (1) Generate the 2I and 1I1E configurations for all irreps
!----------------------------------------------------------------------
    call builder_2I_1I1E(n1h1I,n2I,n1I1E,cfgM(0),icvs,E0max,file2I,&
         file1I1E,nrec2I,nrec1I1E)

!----------------------------------------------------------------------
! (2) Sort the 2I and 1I1E configurations by irrep
!----------------------------------------------------------------------
    call sort_2I_1I1E(cfgM,file2I,file1I1E,nrec2I,nrec1I1E,n2I,n1I1E)
    
    return
    
  end subroutine generate_2I_1I1E_confs

!######################################################################
! builder_2I_1I1E: peforms all the heavy lifting involved in the
!                  generation of the 2I and 1I1E configurations
!                  across all irreps
!######################################################################
  subroutine builder_2I_1I1E(n1h1I,n2I,n1I1E,cfgM,icvs,E0max,file2I,&
       file1I1E,nrec2I,nrec1I1E)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hparam
    use utils
    use iomod
    
    implicit none

    ! Number of 1H1I, 2I and 1I1E configurations
    integer(is), intent(out)       :: n1h1I,n2I,n1I1E

    ! MRCI configurations
    type(mrcfg), intent(inout)     :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)        :: icvs(nmo)
    logical                        :: lcvs

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)           :: E0max

    ! Configuration bit string buffers
    integer(is), intent(out)       :: nrec2I,nrec1I1E
    character(len=60), intent(out) :: file2I,file1I1E
    integer(is)                    :: nbuf2I,nbuf1I1E
    integer(is), allocatable       :: ibuf2I(:,:),ibuf1I1E(:,:)
    integer(is)                    :: iscratch2I,iscratch1I1E
    
    ! Orbital classes
    integer(is)                    :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)                    :: nopen,nsocc,ndocc,nunocc

    ! 2-hole SOPs
    integer(ib), allocatable       :: sop2h(:,:,:)

    ! Irreps generated the various different configurations
    integer(is)                    :: irrepR,irrep2H,irrep1H1I,&
                                      irrep1I1E,irrep2I
    integer(ib)                    :: ib1,ib2
    
    ! Lists of hole/particle indices linking ref confs
    integer(is), parameter         :: maxexci=4
    integer(is)                    :: hlist(maxexci),plist(maxexci)
    integer(is)                    :: indx(maxexci)
    
    ! Allowable internal creation operator indices
    integer(is)                    :: ncreate1,ncreate2
    integer(is), allocatable       :: icreate1(:),icreate2(:)

    ! Difference configuration information
    integer(is), allocatable       :: rhp(:,:)

    ! Work conf bit strings
    integer(ib), allocatable       :: conf(:,:),confI(:,:),sop(:,:),&
                                      conf_int(:,:)

    ! Sums of KS orbital energies multipled by occupation
    ! differences relative to the base conf
    real(dp)                       :: esumR,esum2H,esum1H1I,esum1I1E,&
                                      esum2I
    
    ! Everything else
    integer(is)                    :: i,j,k,n,np,i1,i2,i3
    integer(is)                    :: i2h,iext,iint1,iint2
    integer(is)                    :: n_int_I,nmoI
    integer(is)                    :: ia2h,ja2h,itmp,jtmp,nexci,nmatch
    integer(is)                    :: ic,counter
    integer(is)                    :: nacR,nac2H,nac1H1I,nac1I1E,nac2I
    
!----------------------------------------------------------------------
! Is this a CVS-MRCI calculation
!----------------------------------------------------------------------
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    nmoI=cfgM%nmoI
    n_int_I=cfgM%n_int_I

    allocate(icreate1(nmoI))
    icreate1=0

    allocate(icreate2(nmoI))
    icreate2=0
    
    allocate(conf(n_int,2))
    conf=0_ib
    
    allocate(confI(n_int_I,2))
    confI=0_ib

    allocate(conf_int(n_int_I,2))
    conf_int=0_ib
    
    allocate(sop(n_int,2))
    sop=0_ib
    
    allocate(sop2h(n_int_I,2,cfgM%n2h))
    sop2h=0_ib

    allocate(rhp(9,cfgM%nR))
    rhp=0

    allocate(ibuf2I(4,bufsize))
    ibuf2I=0

    allocate(ibuf1I1E(4,bufsize))
    ibuf1I1E=0

!----------------------------------------------------------------------
! Open the configuration scratch files
!----------------------------------------------------------------------
    ! 2I confs and the indices of their parent 2-hole confs
    call freeunit(iscratch2I)
    call scratch_name('2I',file2I)
    open(iscratch2I,file=file2I,form='unformatted',status='unknown')

    ! 1I1E confs and the indices of their parent 2-hole confs
    call freeunit(iscratch1I1E)
    call scratch_name('1I1E',file1I1E)
    open(iscratch1I1E,file=file1I1E,form='unformatted',&
         status='unknown')

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    n1h1I=0
    n2I=0
    n1I1E=0
    counter=0
    nbuf2I=0
    nbuf1I1E=0
    nrec2I=0
    nrec1I1E=0

!----------------------------------------------------------------------
! Generate the 2-hole SOPs
!----------------------------------------------------------------------    
    ! Loop over 2-hole configurations
    do i=1,cfgM%n2h

       ! Generate the next 2-hole SOP
       confI=cfgM%conf2h(:,:,i)
       sop2h(:,:,i)=conf_to_sop(confI,n_int_I)

    enddo

!----------------------------------------------------------------------
! Generate the 1I1E and 2I configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do n=1,cfgM%nR

       ! Sum_p F_pp^(KS) Delta W_p for the reference configuration
       esumR=esum(cfgM%confR(:,:,n),n_int_I,cfgM%m2c)

       ! Irrep generated by the reference configuration
       sop=0_ib
       sop(1:n_int_I,:)=cfgM%sopR(:,:,n)
       irrepR=sop_sym_mrci(sop,cfgM%m2c)

       ! Number of creation and annihilation operators linking the
       ! reference and base configurations
       conf_int=conf0(1:n_int_I,:)
       nacR=n_create_annihilate(cfgM%confR(:,:,n),conf_int,n_int_I)
       
       ! Initialise the array of inter-ref-conf annihilation/creation
       ! operator indices
       rhp=0

       ! Loop over all other reference configurations
       do np=1,cfgM%nR
          if (n == np) cycle
       
          ! Determine the excitation degree between the to reference
          ! configurations
          nexci=exc_degree_conf(cfgM%confR(:,:,n),cfgM%confR(:,:,np),&
               n_int_I)

          ! Cycle if the excitation degree is greater than 3
          if (nexci > 4) cycle
          
          ! Determine the indices of the creation/annihilation
          ! operators linking the two ref confs
          call get_exci_indices(cfgM%confR(:,:,n),cfgM%confR(:,:,np),&
               n_int_I,hlist(1:nexci),plist(1:nexci),nexci)

          ! Excitation degree
          rhp(1,np)=nexci

          ! Creation/annihilation operator indices in ascending order
          call i4sortindx('A',nexci,hlist(1:nexci),indx(1:nexci))
          do i=1,nexci
             rhp(1+i,np)=hlist(indx(i))
          enddo

          call i4sortindx('A',nexci,plist(1:nexci),indx(1:nexci))
          do i=1,nexci
             rhp(5+i,np)=plist(indx(i))
          enddo

       enddo
       
       ! Loop over the 2-hole configurations generated by the current
       ! reference configuration
       do i2h=cfgM%off2h(n),cfgM%off2h(n+1)-1

          ! 2-hole annihilation operator indices
          ia2h=cfgM%a2h(1,i2h)
          ja2h=cfgM%a2h(2,i2h)

          ! Sum_p F_pp^(KS) Delta W_p for the 2-hole configuration
          esum2H=esumR-moen(cfgM%m2c(ia2h))-moen(cfgM%m2c(ja2h))

          ! Irrep generated by the 2-hole configuration
          ib1=irrepR
          ib2=ieor(ib1,mosym(cfgM%m2c(ia2h)))
          ib1=ieor(ib2,mosym(cfgM%m2c(ja2h)))
          irrep2H=ib1

          ! Number of creation and annihilation operators linking the
          ! 2-hole and base configurations
          nac2H=nacR
          ! First annihilation operator
          k=(ia2h-1)/64+1
          i=ia2h-(k-1)*64-1
          if (iocc0(ia2h) == 0) then
             ! Unoccupied MO in the base conf
             nac2H=nac2H-1
          else if (iocc0(ia2h) == 1) then
             ! Doubly-occupied MO in the base conf
             if (btest(cfgM%confR(k,2,n),i)) then
                ! Doubly-occupied MO in the ref conf
                nac2H=nac2H-1
             else
                ! Singly-occupied MO in the ref conf
                nac2H=nac2H+1
             endif
          else if (iocc0(ia2h) == 2) then
             ! Doubly-occupied MO in the base conf
             nac2H=nac2H+1
          endif
          !
          ! Second annihilation operator
          k=(ja2h-1)/64+1
          i=ja2h-(k-1)*64-1
          if (iocc0(ja2h) == 0) then
             ! Unoccupied MO in the base conf
             nac2H=nac2H-1
          else if (iocc0(ja2h) == 1) then
             ! Singly-occupied MO in the base conf
             if (btest(cfgM%confR(k,2,n),i)) then
                ! Doubly-occupied MO in the ref conf
                if (ia2h == ja2h) then
                   ! Both electrons annihilated
                   nac2H=nac2H+1
                else
                   ! One electron annihilated
                   nac2H=nac2H-1
                endif
             else
                ! Singly-occupied MO in the ref conf
                nac2H=nac2H+1
             endif
          else if (iocc0(ja2h) == 2) then
             ! Doubly-occupied MO in the base conf
             nac2H=nac2H+1
          endif

          ! Initialise the creation operator array
          ncreate1=nmoI
          icreate1=1
          
          ! Remove the 2-hole annihilation operator indices from the
          ! list of allowed creation operators: these would lead to
          ! duplicate confs
          icreate1(ia2h)=0
          ncreate1=ncreate1-1
          if (ia2h /= ja2h) then
             icreate1(ja2h)=0
             ncreate1=ncreate1-1
          endif

          ! Remove the doubly-occupied MOs from the list
          ! of allowed creation operators
          call sop_docc_list(sop2h(:,:,i2h),n_int_I,docc,nmo,ndocc)
          do j=1,ndocc
             if (icreate1(docc(j)) /= 0) then
                icreate1(docc(j))=0
                ncreate1=ncreate1-1
             endif
          enddo

          ! Flag creation operators as forbidden if they will
          ! produce duplicate confs based on the excitations linking
          ! the current ref conf to all other ref confs
          call remove_creators_1H1I(n,cfgM%nR,rhp,ia2h,ja2h,icreate1,&
               ncreate1)

          ! Creation of 1H1I confs by the addition of electrons to
          ! the current 2-hole conf
          do iint1=1,nmoI

             ! Cycle if this this creation operator will yield
             ! a duplicate conf
             if (icreate1(iint1) == 0) cycle
             
             ! Cycle if this is a CVS-MRCI calculation and we are creating
             ! an electron in a flagged core MO
             if (lcvs .and. icvs(cfgM%m2c(iint1)) == 1) cycle

             ! Update the no. 1H1I confs
             n1h1I=n1h1I+1

             ! Initialise the list of 2nd creation operators
             icreate2=icreate1
             ncreate2=ncreate1
             
             ! Sum_p F_pp^(KS) Delta W_p for the 1H1I configuration
             esum1H1I=esum2H+moen(cfgM%m2c(iint1))

             ! Irrep generated by the 1H1I configuration
             ib1=irrep2H
             ib2=ieor(ib1,mosym(cfgM%m2c(iint1)))
             irrep1H1I=ib2

             ! Block index
             k=(iint1-1)/64+1

             ! Postion of the external MO within the kth block
             i=iint1-(k-1)*64-1
             
             ! Number of creation and annihilation operators linking the
             ! 1H1I and base configurations
             nac1H1I=nac2H
             if (iocc0(iint1) == 0) then
                ! Unoccupied MO in the base conf
                nac1H1I=nac1H1I+1
             else if (iocc0(iint1) == 1) then
                ! Singly-occupied MO in the base conf
                if (btest(cfgM%conf2h(k,1,i2h),i)) then
                   ! Singly-occupied MO in the 2-hole conf
                   nac1H1I=nac1H1I+11
                else
                   ! Unoccupied MO in the 2-hole conf
                   nac1H1I=nac1H1I-1
                endif
             else if (iocc0(iint1) == 2) then
                ! Doubly-occupied MO in the base conf
                nac1H1I=nac1H1I-1
             endif
             
             ! Generate the 1H1I conf bit string
             confI=create_electron(cfgM%conf2h(:,:,i2h),n_int_I,iint1)

             ! If we have created an electron in a singly-occupied MO
             ! of the 2-hole conf, then remove it from the list of
             ! allowed 2nd creation operators
             if (btest(confI(k,2),i) .and. icreate2(iint1) /= 0) then
                icreate2(iint1)=0
                ncreate2=ncreate2-1
             endif

             ! Before adding the second electron to create the 2I
             ! confs, we need to flag those creation operators that
             ! would lead to duplicate confs based on the excitations
             ! linking the ref confs
             if (ncreate2 /= 0) call remove_creators_2I(n,cfgM%nR,rhp,&
                  ia2h,ja2h,iint1,icreate2,ncreate2)

             !
             ! Generate all the 1I1E confs resulting from the addition
             ! of an external electron to the current 1H1I conf
             !
             ! Loop over external MOs
             do iext=nmoI+1,lastvirt

                ! Sum_p F_pp^(KS) Delta W_p for the 1I1E configuration
                esum1I1E=esum1H1I+moen(cfgM%m2c(iext))
                
                ! If this is a DFT/MRCI calculation, then skip this
                ! configuration if it doesn't satisfy the energy-based
                ! selection criterion
                if (ldftmrci .and. esum1I1E > E0max + desel) cycle

                ! Block index
                k=(iext-1)/64+1
                
                ! Postion of the external MO within the kth block
                i=iext-(k-1)*64-1

                ! Full configuration
                conf=0_ib
                conf(1:n_int_I,:)=confI
                conf(k,1)=ibset(conf(k,1),i)

                ! Number of creation and annihilation operators
                ! linking the 1I1E and base configurations
                nac1I1E=nac1H1I+1
                
                ! Skip this configuration if the excitation degree
                ! relative to the base configuration is too high
                nexci=nac1I1E/2
                if (nexci > nexmax) cycle
                
                ! Irrep generated by the 1I1E configuration
                ib1=irrep1H1I
                ib2=ieor(ib1,mosym(cfgM%m2c(iext)))
                irrep1I1E=ib2
                
                ! Update the number of 1I1E confs
                n1I1E=n1I1E+1

                ! Save the 1I1E conf
                nbuf1I1E=nbuf1I1E+1
                ibuf1I1E(1,nbuf1I1E)=irrep1I1E
                ibuf1I1E(2,nbuf1I1E)=i2h
                ibuf1I1E(3,nbuf1I1E)=iint1
                ibuf1I1E(4,nbuf1I1E)=iext
                if (nbuf1I1E == bufsize) then
                   nrec1I1E=nrec1I1E+1
                   write(iscratch1I1E) bufsize,ibuf1I1E
                   nbuf1I1E=0
                endif
                
             enddo
                
             !
             ! Generate all the 2I confs resulting from the addition
             ! of an internal electron to the current 1H1I conf
             !
             ! Cycle if there are no allowable internal creation
             ! operators
             if (ncreate2 == 0) cycle
             
             ! Loop over the second internal MO
             do iint2=iint1,nmoI

                ! Cycle if this is a non-allowed creation operator
                if (icreate2(iint2) == 0) cycle

                ! Cycle if this is a CVS-MRCI calculation and we are
                ! creating an electron in a flagged core MO
                if (lcvs .and. icvs(cfgM%m2c(iint2)) == 1) cycle

                ! Sum_p F_pp^(KS) Delta W_p for the 2I configuration
                esum2I=esum1H1I+moen(cfgM%m2c(iint2))
                
                ! If this is a DFT/MRCI calculation, then skip this
                ! configuration if it doesn't satisfy the energy-based
                ! selection criterion
                if (ldftmrci .and. esum2I > E0max + desel) cycle
                
                ! Block index
                k=(iint2-1)/64+1
                
                ! Postion of the internal MO within the kth block
                i=iint2-(k-1)*64-1

                ! Full configuration
                conf=0_ib
                conf(1:n_int_I,:)=confI
                if (btest(conf(k,1),i)) then
                   ! Creation of a doubly-occupied internal MO
                   conf(k,2)=ibset(conf(k,2),i)
                else
                   ! Creation of a singly-occupied internal MO
                   conf(k,1)=ibset(conf(k,1),i)
                endif
                                
                ! Number of creation and annihilation operators linking the
                ! 2I and base configurations
                nac2I=nac1H1I
                if (iocc0(iint2) == 0) then
                   ! Unoccupied MO in the base conf
                   nac2I=nac2I+1
                else if (iocc0(iint2) == 1) then
                   ! Singly-occupied MO in the base conf
                   if (btest(confI(k,1),i)) then
                      ! Singly-occupied MO in the 1H1I conf
                      nac2I=nac2I+1
                   else
                      ! Unoccupied MO in the 1H1I conf
                      nac2I=nac2I-1
                   endif
                else if (iocc0(iint2) == 2) then
                   ! Doubly-occupied MO in the base conf
                   nac2I=nac2I-1
                endif
                
                ! Skip this configuration if the excitation degree
                ! relative to the base configuration is too high
                nexci=nac2I/2
                if (nexci > nexmax) cycle
                
                ! Irrep generated by the 2I configuration
                ib1=irrep1H1I
                ib2=ieor(ib1,mosym(cfgM%m2c(iint2)))
                irrep2I=ib2
                
                ! Update the number of 2I confs
                n2I=n2I+1
                
                ! Save the 2I conf
                nbuf2I=nbuf2I+1
                ibuf2I(1,nbuf2I)=irrep2I
                ibuf2I(2,nbuf2I)=i2h
                ibuf2I(3,nbuf2I)=iint1
                ibuf2I(4,nbuf2I)=iint2
                if (nbuf2I == bufsize) then
                   nrec2I=nrec2I+1
                   write(iscratch2I) bufsize,ibuf2I
                   nbuf2I=0
                endif
                
             enddo
             
          enddo
          
          
       enddo
          
    enddo

!----------------------------------------------------------------------
! Write any non-empty conf buffers to disk    
!----------------------------------------------------------------------
    ! 1I1E
    if (nbuf1I1E /= 0) then
       nrec1I1E=nrec1I1E+1
       write(iscratch1I1E) nbuf1I1E,ibuf1I1E
    endif

    ! 2I
    if (nbuf2I /= 0) then
       nrec2I=nrec2I+1
       write(iscratch2I) nbuf2I,ibuf2I
    endif

!----------------------------------------------------------------------
! Close the configuration scratch files
!----------------------------------------------------------------------
    close(iscratch2I)
    close(iscratch1I1E)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(icreate1)
    deallocate(icreate2)
    deallocate(conf)
    deallocate(confI)
    deallocate(conf_int)
    deallocate(sop)
    deallocate(sop2h)
    deallocate(rhp)
    deallocate(ibuf2I)
    deallocate(ibuf1I1E)
    
    return
    
  end subroutine builder_2I_1I1E
    
!######################################################################
! remove_creators_1H1I: given a 2-hole conf (defined by a ref conf n
!                       and annihilation operators ia and ja), flags
!                       those creation operators that, when operating
!                       on the 2-hole conf, will yield duplicate confs
!######################################################################
  subroutine remove_creators_1H1I(n,nR,rhp,ia,ja,icreate,ncreate)

    use constants
    use bitglobal
    
    implicit none
    
    ! Index of the current ref conf
    integer(is), intent(in)    :: n

    ! Excitations linking the current ref conf to all others
    integer(is), intent(in)    :: nR
    integer(is), intent(in)    :: rhp(9,nR)

    ! Indices of the 2-hole annihilation operators
    integer(is), intent(in)    :: ia,ja

    ! List of allowed/disallowed creation operators
    integer(is), intent(inout) :: icreate(lastvirt)

    ! Number of allowed creation operators
    integer(is), intent(inout) :: ncreate
    
    ! Everything else
    integer(is)                :: j,np,itmp,jtmp,nmatch
    
!----------------------------------------------------------------------
! Excitations linking the current ref conf to the preceding ones
!----------------------------------------------------------------------
    ! Loop over preceding ref confs
    do np=1,n-1

       ! Greater than triple excitation between ref confs n and np:
       ! no creation operators need to be removed
       if (rhp(1,np) > 3) cycle
       
       ! 2-hole annihilation operators
       itmp=ia
       jtmp=ja
             
       if (rhp(1,np) == 1) then
                
          ! Single excitation bewtween ref confs n and np:
          ! Creating an electron in the orbital corresponding
          ! to the creation operator of the excitation linking
          ! the two ref confs will lead to duplicate confs
          
          ! Flag the creation operator linking the ref confs
          ! as forbidden
          if (icreate(rhp(6,np)) /= 0) then
             icreate(rhp(6,np))=0
             ncreate=ncreate-1
          endif
             
       else if (rhp(1,np) == 2) then
                
          ! Double excitation bewtween ref confs n and np:
          ! If one of the two annihilation operators linking
          ! the ref confs matches one of the annihilation
          ! operators of the 2-hole conf, then duplicate
          ! confs will be formed from excitations into the
          ! orbitals corresponding to the creation operators
          ! linking the ref confs

          ! How many matched annihilation operators do we have?
          nmatch=0
          do j=1,2
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
             endif
          enddo
          
          ! If we have a single matched annihilation operator,
          ! then flag the creation operators linking the ref
          ! confs as forbidden
          if (nmatch == 1) then
             do j=1,2
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             enddo
          endif
          
       else if (rhp(1,np) == 3) then
          
          ! Triple excitation bewtween ref confs n and np:
          ! If two of the two annihilation operators linking
          ! the ref confs match the annihilation operators
          ! of the 2-hole conf, then duplicate confs will be
          ! formed from excitations into the orbitals
          ! corresponding to the creation operators linking
          ! the ref confs
          
          ! How many matched annihilation operators do we have?
          nmatch=0
          do j=1,3
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
             endif
          enddo
          
          ! If we have two matched annihilation operators,
          ! then flag the creation operators linking the ref
          ! confs as forbidden
          if (nmatch == 2) then
             do j=1,3
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             enddo
          endif
                
       endif
             
    enddo

!----------------------------------------------------------------------
! Excitations linking the current ref conf to the proceding ones
!----------------------------------------------------------------------
    ! Loop over proceding ref confs
    do np=n+1,nR

       ! Greater than double excitation between ref confs n and np:
       ! no creation operators need to be removed
       if (rhp(1,np) > 2) cycle
       
       ! 2-hole annihilation operators
       itmp=ia
       jtmp=ja
             
       if (rhp(1,np) == 1) then

          ! Single excitation between ref confs n and np:
          ! If the annihilation operator linking the two
          ! ref confs matchs either of the 2-hole conf
          ! annihilation operators, then duplicate confs
          ! will be formed by creating an electron in the
          ! orbital corresponding to the creation operators
          ! linking the ref confs
                
          if (ia == rhp(2,np) .or. ja == rhp(2,np)) then
             if (icreate(rhp(6,np)) /= 0) then
                icreate(rhp(6,np))=0
                ncreate=ncreate-1
             endif
          endif
          
       else if (rhp(1,np) == 2) then
          ! Double excitation bewtween ref confs n and np:
          ! If both of the two annihilation operators linking
          ! the ref confs match of the annihilation
          ! operators of the 2-hole conf, then duplicate
          ! confs will be formed from excitations into the
          ! orbitals corresponding to the creation operators
          ! linking the ref confs
          
          ! How many matched annihilation operators do we have?
          nmatch=0
          do j=1,2
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
             endif
          enddo

          ! If we have two matched annihilation operators,
          ! then flag the creation operators linking the ref
          ! confs as forbidden
          if (nmatch == 2) then
             do j=1,2
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             enddo
          endif
          
       endif
       
    enddo
    
    return
    
  end subroutine remove_creators_1H1I

!######################################################################
! remove_creators_2I: given a 1H1I conf (defined by a ref conf n,
!                     annihilation operators ia and ja and a
!                     creation operator ic), flags those creation
!                     operators that, when operating on the 1H1I
!                     conf, will yield duplicate confs
!######################################################################
  subroutine remove_creators_2I(n,nR,rhp,ia,ja,ic,icreate,ncreate)

    use constants
    use bitglobal

    implicit none

    ! Index of the current ref conf
    integer(is), intent(in)    :: n

    ! Excitations linking the current ref conf to all others
    integer(is), intent(in)    :: nR
    integer(is), intent(in)    :: rhp(9,nR)

    ! Indices of the 2-hole annihilation operators
    integer(is), intent(in)    :: ia,ja

    ! Index of the 1H1I creation operator
    integer(is), intent(in)    :: ic
    
    ! List of allowed/disallowed creation operators
    integer(is), intent(inout) :: icreate(lastvirt)

    ! Number of allowed/disallowed creation operators
    integer(is), intent(inout) :: ncreate
    
    ! Everything else
    integer(is)                :: j,np,itmp,jtmp,nmatch,icmatch
    logical                    :: match1
    
!----------------------------------------------------------------------
! Excitations linking the current ref conf to the preceding ones
!----------------------------------------------------------------------
    ! Loop over preceding ref confs
    do np=1,n-1

       ! Greater than quadruple excitation between ref confs n and np:
       ! no creation operators need to be removed
       if (rhp(1,np) > 4) cycle
              
       ! 2-hole annihilation operators
       itmp=ia
       jtmp=ja

       if (rhp(1,np) == 2) then
          ! Double excitation between ref confs n an np:
          !
          ! (i)  If neither of the annihilation operators
          !      linking the ref confs match the annihilation
          !      operators of the 1H1I conf *and* there is
          !      a matched creation operator, then excitations
          !      into the orbitals above the matched creation
          !      operator will give rise to duplicate confs
          !
          ! (ii) If one of the annihilation operators
          !      linking the ref confs matches an annihilation
          !      operator of the 1H1I conf *and* there are no
          !      matched creators, then excitations into the
          !      orbitals of the creation operators linking
          !      the ref confs will give rise to duplicate
          !      configurations
          
          ! Number of matched annihilation operators
          nmatch=0
          
          ! Comparison of the annihilation operators
          ! of the 2-hole conf and those linking the
          ! ref confs
          do j=1,2
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
             endif
          enddo
          
          ! Comparison of the 1H1I creation operator
          ! with those of the excitation linking the
          ! ref confs
          icmatch=0
          do j=1,2
             if (ic == rhp(5+j,np)) then
                icmatch=j
                exit
             endif
          enddo
          
          ! Zero matched annihilation operators and a
          ! matched creation operator: remove all the
          ! creation operator past the matched one
          if (nmatch == 0 .and. icmatch == 1) then
             if (icreate(rhp(7,np)) /= 0) then
                icreate(rhp(7,np))=0
                ncreate=ncreate-1
             endif
          endif
               
          ! One matched annihilation operator and zero
          ! matched creation operators: remove all the
          ! creation operators
          if (nmatch == 1 .and. icmatch == 0) then
             do j=1,2
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             enddo
          endif
             
       else if (rhp(1,np) == 3) then
          ! Triple excitation between ref confs n an np:
          !
          ! (i)  If one of the annihilation operators
          !      linking the ref confs matches an annihilation
          !      operator of the 1H1I conf *and* there is
          !      a matched creation operator, then excitations
          !      into the orbitals above the matched creation
          !      operator will give rise to duplicate confs
          !
          ! (ii) If two of the annihilation operators
          !      linking the ref confs matches an annihilation
          !      operator of the 1H1I conf *and* there are no
          !      matched creators, then excitations into the
          !      orbitals of the creation operators linking
          !      the ref confs will give rise to duplicate
          !      configurations
          
          ! Number of matched annihilation operators
          nmatch=0
          
          ! Comparison of the annihilation operators
          ! of the 2-hole conf and those linking the
          ! ref confs
          do j=1,3
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
                if (nmatch == 2) exit
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
                if (nmatch == 2) exit
             endif
          enddo
          
          ! Comparison of the 1H1I creation operator
          ! with those of the excitation linking the
          ! ref confs
          icmatch=0
          do j=1,3
             if (ic == rhp(5+j,np)) then
                icmatch=j
                exit
             endif
          enddo
          
          ! One matched annihilation operator and a
          ! matched creation operator: remove all the
          ! creation operator past the matched one
          if (nmatch == 1 .and. icmatch /= 0) then
             do j=icmatch+1,3
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             enddo
          endif
          
          ! Two matched annihilation operators and zero
          ! matched creation operators: remove all the
          ! creation operators
          if (nmatch == 2 .and. icmatch == 0) then
             do j=1,3
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             enddo
          endif

       else if (rhp(1,np) == 4) then
          ! Quadruple excitation between ref confs n an np:
          !
          ! If two of the annihilation operators
          ! linking the ref confs matches an annihilation
          ! operator of the 1H1I conf *and* there is a
          ! matched creator, then excitations into the orbitals
          ! above the matched creation operator will give rise
          ! to duplicate confs
          
          ! Number of matched annihilation operators
          nmatch=0
          
          ! Comparison of the annihilation operators
          ! of the 2-hole conf and those linking the
          ! ref confs
          do j=1,4
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
                if (nmatch == 2) exit
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
                if (nmatch == 2) exit
             endif
          enddo

          ! Comparison of the 1H1I creation operator
          ! with those of the excitation linking the
          ! ref confs
          icmatch=0
          do j=1,4
             if (ic == rhp(5+j,np)) then
                icmatch=j
                exit
             endif
          enddo
          
          ! Two matched annihilation operators and a
          ! matched creation operator: remove all the
          ! creation operator past the matched one
          if (nmatch == 2 .and. icmatch /= 0) then
             do j=icmatch+1,4
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             enddo
          endif

       endif
                
    enddo

!----------------------------------------------------------------------
! Excitations linking the current ref conf to the proceding ones
!----------------------------------------------------------------------
    ! Loop over proceding ref confs
    do np=n+1,nR

       ! Greater than triple excitation between ref confs n and np:
       ! no creation operators need to be removed
       if (rhp(1,np) > 3) cycle

       ! 2-hole annihilation operators
       itmp=ia
       jtmp=ja

       if (rhp(1,np) == 1) then
          ! Single excitation linking ref confs n and np:
          ! If any one of the annihilation operators of the 2-hole conf
          ! matches that of the excitation linking the ref confs, then
          ! creating an electron in the orbital corresponding to the
          ! creation operator of the excitation linking the ref confs
          ! will yield duplicate confs
          ! Futhermore, if the creation operator equals that of the
          ! 1H1I conf, then all creation operations must be removed
          
          if (ia == rhp(2,np) .or. ja == rhp(2,np)) then
             if (ic == rhp(6,np)) then
                icreate=0
                ncreate=0
                return
             else
                if (icreate(rhp(6,np)) /= 0) then
                   icreate(rhp(6,np))=0
                   ncreate=ncreate-1
                endif
             endif
          endif

       else if (rhp(1,np) == 2) then
          ! Double excitation linking ref confs n and np:
          ! If there are matched annihilators *and* the 1H1I creator
          ! matches one the the R_n -> R-np excitation, each creator
          ! of the R_n -> R-np excitation has to be removed if:
          ! (1) the 1H1I creator doesn't match the R_n -> R-np creator
          ! (2) the 1H1I creator does match the R_n -> R-np creator
          !     *and* the R_n -> R-np creators are equal
          
          ! How many matched annihilation operators do we have
          nmatch=0
          do j=1,2
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
             endif
          enddo

          ! Cycle if there are no matched annihilation operators
          if (nmatch == 0) cycle

          ! Cycle if the creation operator of the 1H1I conf doesn't
          ! match one of the R_n -> R-np excitation
          if (ic /= rhp(6,np) .and. ic /= rhp(7,np)) cycle

          ! Loop over R_n -> R-np creation operators
          do j=1,2
             if (ic /= rhp(5+j,np)) then
                ! Unmatched creation operators: remove the
                ! R_n -> R-np creation operator from the list
                if (icreate(rhp(5+j,np)) /= 0) then
                   icreate(rhp(5+j,np))=0
                   ncreate=ncreate-1
                endif
             else
                ! Matched creation operators: remove the
                ! R_n -> R-np creation operator from the list if
                ! the two R_n -> R-np creation operator are not
                ! equal
                if (rhp(6,np) == rhp(7,np)) then
                   if (icreate(rhp(5+j,np)) /= 0) then
                      icreate(rhp(5+j,np))=0
                      ncreate=ncreate-1
                   endif
                endif
             endif
          enddo

       else if (rhp(1,np) == 3) then
          ! Triple excitation linking ref confs n and np:
          ! If there are two matched annihilators, then the
          ! creators of the R_n -> R_np excitation
          ! have to be removed, up to the first one that matches
          ! with the creator of the 1H1I conf
          
          ! How many matched annihilation operators do we have
          nmatch=0
          do j=1,3
             if (itmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                itmp=-1
             else if (jtmp == rhp(1+j,np)) then
                nmatch=nmatch+1
                jtmp=-1
             endif
          enddo

          ! Cycle we do not have 2 matched annihilation operators
          if (nmatch /= 2) cycle

          ! Cycle if we don't have a matched creation operator
          if (ic /= rhp(6,np) .and. ic /= rhp(7,np) &
               .and. ic /= rhp(8,np)) cycle
          
          ! Flag the creation operators of the R_n -> R_np
          ! excitation as forbidden, *not* including the first
          ! match with the creation operator of the 1H1I conf
          match1=.false.
          do j=1,3

             ! Skip the first matched creation operator
             if (ic == rhp(5+j,np) .and. .not. match1) then
                match1=.true.
                cycle
             endif

             ! Flag the remaining ones as forbidden
             if (icreate(rhp(5+j,np)) /= 0) then
                icreate(rhp(5+j,np))=0
                ncreate=ncreate-1
             endif
             
          enddo
          
       endif
       
    enddo
    
    return
    
  end subroutine remove_creators_2I

!######################################################################
! sort_2I_1I1E: sorts and saves the 2I and 1I1E confg=igurations by
!               irrep
!######################################################################
  subroutine sort_2I_1I1E(cfgM,file2I,file1I1E,nrec2I,nrec1I1E,&
       n2I_tot,n1I1E_tot)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use iomod
    
    implicit none

    ! MRCI configurations for all irreps
    type(mrcfg), intent(inout)    :: cfgM(0:nirrep-1)

    ! Total number of 2I and 1I1E configurations
    integer(is), intent(in)       :: n2I_tot,n1I1E_tot
    
    ! Configuration scratch files
    integer(is), intent(in)       :: nrec2I,nrec1I1E
    character(len=60), intent(in) :: file2I,file1I1E
    integer(is), allocatable      :: ibuf(:,:)
    integer(is)                   :: iscratch2I,iscratch1I1E

    ! Working arrays
    integer(ib), allocatable      :: sop(:,:)

    ! No. confs generated by each 2-hole conf
    integer(is), allocatable      :: ngen(:,:)
    
    ! Everything else
    integer(is)                   :: nbuf,i,irec,irrep
    integer(is)                   :: counter(0:nirrep-1)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(ibuf(4,bufsize))
    ibuf=0

    allocate(sop(n_int,2))
    sop=0_ib

    allocate(ngen(cfgM(0)%n2h,0:nirrep-1))
    ngen=0
    
!----------------------------------------------------------------------
! Open the configuration scratch files
!----------------------------------------------------------------------
    ! 2I confs and the indices of their parent 2-hole confs
    call freeunit(iscratch2I)
    call scratch_name('2I',file2I)
    open(iscratch2I,file=file2I,form='unformatted',status='unknown')

    ! 1I1E confs and the indices of their parent 2-hole confs
    call freeunit(iscratch1I1E)
    call scratch_name('1I1E',file1I1E)
    open(iscratch1I1E,file=file1I1E,form='unformatted',&
         status='unknown')

!----------------------------------------------------------------------
! Determine the no. confs per irrep
!----------------------------------------------------------------------
    if (nirrep == 1) then
       ! C1 symmetry

       ! Initialisation
       cfgM(0)%n2I=n2I_tot
       cfgM(0)%n1I1E=n1I1E_tot

    else
       ! Non-C1 symmetry

       ! Initialisation
       do irrep=0,nirrep-1
          cfgM(irrep)%n2I=0
          cfgM(irrep)%n1I1E=0
       enddo
       
       !
       ! 2I
       !
       ! Loop over records
       do irec=1,nrec2I
          
          ! Read in the next batch of confs
          read(iscratch2I) nbuf,ibuf

          ! Determine the irreps generated by the confs
          do i=1,nbuf
             irrep=ibuf(1,i)
             cfgM(irrep)%n2I=cfgM(irrep)%n2I+1
          enddo
          
       enddo

       !
       ! 1I1E
       !
       ! Loop over records
       do irec=1,nrec1I1E
          
          ! Read in the next batch of confs
          read(iscratch1I1E) nbuf,ibuf
          
          ! Determine the irreps generated by the confs
          do i=1,nbuf
             irrep=ibuf(1,i)
             cfgM(irrep)%n1I1E=cfgM(irrep)%n1I1E+1
          enddo
          
       enddo

    endif
       
!----------------------------------------------------------------------
! Allocate the offset and creation operator index arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Creation operator indices
       allocate(cfgM(irrep)%a2I(2,cfgM(irrep)%n2I))
       allocate(cfgM(irrep)%a1I1E(2,cfgM(irrep)%n1I1E))
       cfgM(irrep)%a2I=0
       cfgM(irrep)%a1I1E=0
       
       ! Offsets
       allocate(cfgM(irrep)%off2I(cfgM(irrep)%n2h+1))
       allocate(cfgM(irrep)%off1I1E(cfgM(irrep)%n2h+1))
       cfgM(irrep)%off1I1E=0
       cfgM(irrep)%off1I1E=0
       
    enddo

!----------------------------------------------------------------------
! 2I operator indices and offsets
!----------------------------------------------------------------------
    ! Initialise counters
    counter=0
    ngen=0

    ! Rewind the 2I scratch file
    rewind(iscratch2I)
    
    ! Loop over records
    do irec=1,nrec2I
       
       ! Read in the next batch of confs
       read(iscratch2I) nbuf,ibuf
       
       ! Loop over confs in the current batch
       do i=1,nbuf

          ! Symmetry of the conf
          irrep=ibuf(1,i)

          ! Increment the conf counter for this irrep
          counter(irrep)=counter(irrep)+1

          ! Update the no. confs generated by the parent
          ! 2-hole conf
          ngen(ibuf(2,i),irrep)=ngen(ibuf(2,i),irrep)+1

          ! Save the creation operator indices
          cfgM(irrep)%a2I(:,counter(irrep))=ibuf(3:4,i)
          
       enddo

    enddo

    ! Fill in the offset arrays
    do irrep=0,nirrep-1
       cfgM(irrep)%off2I(1)=1
       do i=2,cfgM(irrep)%n2h+1
          cfgM(irrep)%off2I(i)=cfgM(irrep)%off2I(i-1)+ngen((i-1),irrep)
       enddo
    enddo

!----------------------------------------------------------------------
! I1IE operator indices and offsets
!----------------------------------------------------------------------
    ! Initialise counters
    counter=0
    ngen=0

    ! Rewind the 1I1E scratch file
    rewind(iscratch1I1E)
    
    ! Loop over records
    do irec=1,nrec1I1E
       
       ! Read in the next batch of confs
       read(iscratch1I1E) nbuf,ibuf
       
       ! Loop over confs in the current batch
       do i=1,nbuf

          ! Symmetry of the conf
          irrep=ibuf(1,i)

          ! Increment the conf counter for this irrep
          counter(irrep)=counter(irrep)+1

          ! Update the no. confs generated by the parent
          ! 2-hole conf
          ngen(ibuf(2,i),irrep)=ngen(ibuf(2,i),irrep)+1

          ! Save the creation operator indices
          cfgM(irrep)%a1I1E(:,counter(irrep))=ibuf(3:4,i)
          
       enddo

    enddo

    ! Fill in the offset arrays
    do irrep=0,nirrep-1
       cfgM(irrep)%off1I1E(1)=1
       do i=2,cfgM(irrep)%n2h+1
          cfgM(irrep)%off1I1E(i)=&
               cfgM(irrep)%off1I1E(i-1)+ngen((i-1),irrep)
       enddo
    enddo
    
!----------------------------------------------------------------------
! Close the configuration scratch files
!----------------------------------------------------------------------
    close(iscratch2I)
    close(iscratch1I1E)

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ibuf)
    deallocate(sop)
    deallocate(ngen)
    
    return
    
  end subroutine sort_2I_1I1E
  
!######################################################################
! generate_1E_hash: for a given irrep, generates all the allowable
!                   configurations with one internal hole and one
!                   external electron
!######################################################################
  subroutine generate_1E_hash(irrep,E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep
    
    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM
    
    ! Everything else
    integer(is)                :: modus

!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! correct symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_1E_hash(modus,irrep,E0max,cfgM)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the external creation operators
    allocate(cfgM%a1E(cfgM%n1E))
    cfgM%a1E=0

    ! Offsets
    allocate(cfgM%off1E(cfgM%n1h+1))
    cfgM%off1E=0

!----------------------------------------------------------------------
! Fill in the external creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off1E array has the dimension (n1h+1), even though
! not all 1-hole configurations may generate allowable
! 1E configurations. However, this is OK due to how the offset array
! is used when looping over the 1E configurations generated by each
! 1-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_1E_hash(modus,irrep,E0max,cfgM)

    return
    
  end subroutine generate_1E_hash
  
!######################################################################
! 1E_builder_hash: performs all the heavy lifting involved in the
!                  generation of the configurations with one internal
!                  hole and one external electron
!######################################################################
  subroutine builder_1E_hash(modus,irrep,E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM
    
    ! Full configurations and SOPs
    integer(ib)                :: conf_full(n_int,2)
    integer(ib)                :: sop_full(n_int,2)

    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:)
    
    ! Everything else
    integer(is)                :: n_int_I,nmoI,n1h
    integer(is)                :: n,iext,k,i,counter
    logical                    :: ok,checksym

!----------------------------------------------------------------------    
! Allocate arrays
!----------------------------------------------------------------------
    n1h=cfgM%n1h
    n_int_I=cfgM%n_int_I
    nmoI=cfgM%nmoI

    allocate(ngen(n1h))
    ngen=0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n1E=0
    checksym=.true.
    counter=0

!----------------------------------------------------------------------    
! Generate the 1E configurations
!----------------------------------------------------------------------
    ! Loop over the 1-hole configurations
    do n=1,n1h

       ! Loop over external MOs
       do iext=nmoI+1,lastvirt
       
          ! Block index
          k=(iext-1)/64+1
          
          ! Postion of the external MO within the kth block
          i=iext-(k-1)*64-1
          
          ! Full configuration
          conf_full=0_ib
          conf_full(1:n_int_I,:)=cfgM%conf1h(:,:,n)
          conf_full(k,1)=ibset(conf_full(k,1),i)
          
          ! Full SOP
          sop_full=conf_to_sop(conf_full,n_int)

          ! Is this an allowable configuration?
          ok=allowable_conf(conf_full,sop_full,irrep,cfgM%m2c,E0max,&
               nomax,checksym)

          if (ok) then
             ! Update the total no. allowable configurations generated
             counter=counter+1
             ! Update the no. confs generated by the current hole conf
             ngen(n)=ngen(n)+1
          endif
             
          ! Cycle if we are just determining the no. allowable
          ! configurations...
          if (modus == 0) cycle

          ! ...else fill in the creation operator array
          if (ok) cfgM%a1E(counter)=iext

       enddo
       
    enddo
    
!----------------------------------------------------------------------
! If we are just determining the number of allowable configurations,
! then save this number
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n1E=counter

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    if (modus == 1) then
       cfgM%off1E(1)=1
       do n=2,n1h+1
          cfgM%off1E(n)=cfgM%off1E(n-1)+ngen(n-1)
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)
    
    return
    
  end subroutine builder_1E_hash

!######################################################################
! generate_2E_confs_hash: for a given irrep, generates all the
!                         allowable configurations with two internal
!                         holes and two external electrons
!######################################################################
  subroutine generate_2E_hash(irrep,E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM
    
    ! Everything else
    integer(is)                :: modus

!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! correct symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_2E_hash(modus,irrep,E0max,cfgM)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the external creation operators
    allocate(cfgM%a2E(2,cfgM%n2E))
    cfgM%a2E=0
    
    ! Offsets
    allocate(cfgM%off2E(cfgM%n2h+1))
    cfgM%off2E=0

!----------------------------------------------------------------------
! Fill in the external creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off2E array has the dimension (n2h+1), even though
! not all 2-hole configurations may generate allowable
! ext configurations. However, this is OK due to how the offset array
! is used when looping over the 2E configurations generated by each
! 2-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_2E_hash(modus,irrep,E0max,cfgM)
    
    return
    
  end subroutine generate_2E_hash

!######################################################################
! builder_2E_hash: performs all the heavy lifting involved in the
!                  generation of the configurations with two internal
!                  holes and two external electrons
!######################################################################
  subroutine builder_2E_hash(modus,irrep,E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hparam
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM
        
    ! Configurations and SOPs
    integer(ib)                :: conf_full1(n_int,2),&
                                  conf_full2(n_int,2)
    integer(ib)                :: sop_full(n_int,2)

    ! Difference configuration information
    integer(is)                :: Dw(nmo,2)
    integer(is)                :: ndiff

    ! DFT/MRCI energy selection criterion
    real(dp), allocatable      :: sum2h(:)
    real(dp)                   :: elim,sum,sum1E,sum2E

    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:)
    
    ! Everything else
    integer(is)                :: n_int_I,nmoI,n2h
    integer(is)                :: n,i,iext1,iext2,k1,k2,i1,i2,Dwi,&
                                  counter
    logical                    :: ok,checksym

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgM%n_int_I
    nmoI=cfgM%nmoI
    n2h=cfgM%n2h

    allocate(sum2h(n2h))
    sum2h=0.0d0

    allocate(ngen(n2h))
    ngen=0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n2E=0
    counter=0

!----------------------------------------------------------------------
! If this is a DFT/MRCI calculation, then first determine the
! contribution from the internal holes to the energy selection
! sum
!----------------------------------------------------------------------
    if (ldftmrci) then

       ! Loop over 2-hole configurations
       do n=1,n2h

          ! Compute the difference configuration
          call diffconf(cfgM%conf2h(:,:,n),n_int_I,Dw,nmo,ndiff)

          ! Sum_i (F_ii^KS) Dw_i
          do i=1,ndiff

             ! MO index
             i1=cfgM%m2c(Dw(i,1))
             
             ! Delta w_i value
             Dwi=Dw(i,2)

             ! Sum the contribution
             sum2h(n)=sum2h(n)+moen(i1)*Dwi
          
          enddo
          
       enddo

    endif

!----------------------------------------------------------------------    
! Generate the 2E configurations
!----------------------------------------------------------------------
    ! DFT/MRCI energy selection limit
    if (ldftmrci) elim=E0max+desel

    ! Loop over the 2-hole configurations
    do n=1,n2h

       ! If this is a DFT/MRCI calculation, cycle if the current
       ! 2-hole configuration cannot possibly generate a 2E
       ! configuration satisfying the energy selection criterion
       if (ldftmrci .and. sum2h(n) > elim) cycle
       
       ! Loop over the first creation operator
       do iext1=nmoI+1,lastvirt
          
          ! If this is a DFT/MRCI calculation, cycle if the current
          ! 1H1E configuration doesn't satisfy the energy selection
          ! criterion
          sum1E=sum2h(n)+moen(cfgM%m2c(iext1))
          if (ldftmrci .and. sum1E > elim) cycle
          
          ! First block index
          k1=(iext1-1)/64+1

          ! Postion of the first external MO within the k1th block
          i1=iext1-(k1-1)*64-1

          ! Add the first electron
          conf_full1=0_ib
          conf_full1(1:n_int_I,:)=cfgM%conf2h(:,:,n)
          conf_full1(k1,1)=ibset(conf_full1(k1,1),i1)
          sop_full=conf_to_sop(conf_full1,n_int)
          
          ! Check to see if the partial configuration can
          ! generate an allowable 2E configuration
          checksym=.false.
          ok=allowable_conf(conf_full1,sop_full,irrep,cfgM%m2c,E0max,&
               nomax+1,checksym)
          if (.not. ok) cycle

          ! Loop over the second creation operator
          do iext2=iext1,lastvirt
          
             ! If this is a DFT/MRCI calculation, cycle if the current
             ! 2E configuration satisfying the energy selection
             ! criterion
             sum2E=sum1E+moen(cfgM%m2c(iext2))
             if (ldftmrci .and. sum2E > elim) cycle
             
             ! Second block index
             k2=(iext2-1)/64+1

             ! Postion of the seconf external MO within the k2th block
             i2=iext2-(k2-1)*64-1

             ! Add the second electron
             conf_full2=conf_full1
             if (k1 /= k2) then
                ! Two singly-occupied external MOs
                conf_full2(k2,1)=ibset(conf_full2(k2,1),i2)   
             else
                ! One doubly-occupied external MO
                conf_full2(k2,2)=ibset(conf_full2(k2,2),i2)   
             endif

             ! Full SOP
             sop_full=conf_to_sop(conf_full2,n_int)

             ! Is this an allowable configuration?
             checksym=.true.
             ok=allowable_conf(conf_full2,sop_full,irrep,cfgM%m2c,&
                  E0max,nomax,checksym)

             if (ok) then
                ! Update the total no. allowable configurations
                ! generated
                counter=counter+1
                ! Update the no. confs generated by the current hole
                ! conf
                ngen(n)=ngen(n)+1
             endif
             
             ! Fill in the creation operator array
             if (ok .and. modus == 1) then
                cfgM%a2E(1,counter)=iext1
                cfgM%a2E(2,counter)=iext2
             endif
                
          enddo

       enddo

    enddo

!----------------------------------------------------------------------
! If we are just determining the number of allowable configurations,
! then save this number
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n2E=counter

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    if (modus == 1) then
       cfgM%off2E(1)=1
       do n=2,n2h+1
          cfgM%off2E(n)=cfgM%off2E(n-1)+ngen(n-1)
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(sum2h)
    deallocate(ngen)
    
    return
    
  end subroutine builder_2E_hash

!######################################################################
! generate_1I_hash: for a given irrep, generates all the allowable
!                   configurations with one internal hole and one
!                   internal electron
!######################################################################
  subroutine generate_1I_hash(irrep,E0max,cfgM,icvs)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs
    
    ! Everything else
    integer(is)                :: modus

!----------------------------------------------------------------------
! Is this a CVS-MRCI calculation
!----------------------------------------------------------------------
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif
    
!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! correct symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_1I_hash(modus,irrep,E0max,cfgM,icvs,lcvs)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the internal creation operators
    allocate(cfgM%a1I(cfgM%n1I))
    cfgM%a1I=0

    ! Offsets
    allocate(cfgM%off1I(cfgM%n1h+1))
    cfgM%off1I=0

!----------------------------------------------------------------------
! Fill in the internal creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off1I array has the dimension (n1h+1), even though
! not all 1-hole configurations may generate allowable
! int configurations. However, this is OK due to how the offset array
! is used when looping over the 1I configurations generated by each
! 1-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_1I_hash(modus,irrep,E0max,cfgM,icvs,lcvs)
    
    return
    
  end subroutine generate_1I_hash

!######################################################################
! builder_1I_hash: performs all the heavy lifting involved in the
!                  generation of the configurations with one internal
!                  hole and one internal electron
!######################################################################
  subroutine builder_1I_hash(modus,irrep,E0max,cfgM,icvs,lcvs)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use dethash
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical, intent(in)        :: lcvs
    
    ! Full configurations and SOPs
    integer(ib)                :: conf_full(n_int,2)
    integer(ib)                :: sop_full(n_int,2)

    ! Hash table
    type(dhtbl)                :: h
    integer(is)                :: initial_size
    integer(is)                :: nold
    integer(ib)                :: key(n_int,2)

    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:)
    
    ! Everything else
    integer(is)                :: n_int_I,nmoI,n1h
    integer(is)                :: n,iint,k,i,counter,iref
    logical                    :: ok,unique,checksym

!----------------------------------------------------------------------    
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgM%n_int_I
    nmoI=cfgM%nmoI
    n1h=cfgM%n1h

    allocate(ngen(n1h))
    ngen=0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n1I=0
    counter=0
    checksym=.true.

!----------------------------------------------------------------------
! Initialise the hash table and insert the reference configurations
!----------------------------------------------------------------------
    ! Initialisation
    initial_size=2**15
    call h%initialise_table(initial_size)

    ! Insertion of the ref confs
    do iref=1,cfgM%n0h
       key=0_ib
       key(1:n_int_I,:)=cfgM%conf0h(:,:,iref)
       call h%insert_key(key)
    enddo

!----------------------------------------------------------------------
! Generate the 1I configurations
!----------------------------------------------------------------------
    ! Number of confs currently stored in the hash table
    nold=h%n_keys_stored
    
    ! Loop over the 1-hole configurations
    do n=1,n1h

       ! Loop over internal MOs
       do iint=1,nmoI
          
          ! Cycle if this is a CVS-MRCI calculation and we are creating
          ! an electron in a flagged core MO
          if (lcvs .and. icvs(cfgM%m2c(iint)) == 1) cycle

          ! Block index
          k=(iint-1)/64+1

          ! Postion of the internal MO within the kth block
          i=iint-(k-1)*64-1
          
          ! Cycle if this MO is doubly-occupied in the
          ! 1-hole configuration
          if (btest(cfgM%conf1h(k,2,n),i)) cycle

          ! Full configuration
          conf_full=0_ib
          conf_full(1:n_int_I,:)=cfgM%conf1h(:,:,n)
          if (btest(conf_full(k,1),i)) then
             ! Creation of a doubly-occupied internal MO
             conf_full(k,2)=ibset(conf_full(k,2),i)
          else
             ! Creation of a singly-occupied internal MO
             conf_full(k,1)=ibset(conf_full(k,1),i)
          endif
             
          ! Full SOP
          sop_full=conf_to_sop(conf_full,n_int)
          
          ! Is this an allowable configuration?
          ok=allowable_conf(conf_full,sop_full,irrep,cfgM%m2c,E0max,&
               nomax,checksym)

          ! Is this a unique configuration?
          if (ok) then

             ! Attempt an insertion into the hash table
             call h%insert_key(conf_full)

             ! If the number of stored keys has increased, then
             ! this is a new configuration
             if (h%n_keys_stored > nold) then
                unique=.true.
                nold=h%n_keys_stored
             else
                unique=.false.
             endif
             
          endif

          if (ok .and. unique) then
             ! Update the total no. allowable configurations generated
             counter=counter+1
             ! Update the no. confs generated by the current hole conf
             ngen(n)=ngen(n)+1
          endif
             
          ! Cycle if we are just determining the no. allowable
          ! configurations...
          if (modus == 0) cycle

          ! ...else fill in the creation operator array
          if (ok .and. unique) cfgM%a1I(counter)=iint
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
    call h%delete_table

!----------------------------------------------------------------------
! If we are just determining the number of allowable configurations,
! then save this number
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n1I=counter

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    if (modus == 1) then
       cfgM%off1I(1)=1
       do n=2,n1h+1
          cfgM%off1I(n)=cfgM%off1I(n-1)+ngen(n-1)
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)
    
    return
    
  end subroutine builder_1I_hash

!######################################################################
! generate_1I1E_hash: for a given irrep, generates all the allowable
!                      configurations with two internal holes, one
!                      internal electron and one external electron
!######################################################################
  subroutine generate_1I1E_hash(irrep,E0max,conf1h1I,indx1h1I,&
       n_int_I,n1h1I,cfgM)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max
    
    ! 1-hole and 2-hole configurations
    integer(is), intent(in)    :: n_int_I,n1h1I
    integer(ib), intent(in)    :: conf1h1I(n_int_I,2,n1h1I)
    integer(is), intent(in)    :: indx1h1I(2,n1h1I)
    
    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM
    
    ! Everything else
    integer(is)                :: modus
    
!----------------------------------------------------------------------
! Determine the number of allowable configurations of the correct
! symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_1I1E_hash(modus,irrep,E0max,conf1h1I,indx1h1I,&
         n_int_I,n1h1I,cfgM)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the internal and external creation operators
    allocate(cfgM%a1I1E(2,cfgM%n1I1E))
    cfgM%a1I1E=0

    ! Offsets
    allocate(cfgM%off1I1E(cfgM%n2h+1))
    cfgM%off1I1E=0

!----------------------------------------------------------------------
! Fill in the internal and external creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off1I1E array has the dimension (n2h+1), even though
! not all 2-hole configurations may generate allowable
! 1I1E configurations. However, this is OK due to how the
! offset arrayis used when looping over the 1I1E configurations
! generated by each 2-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_1I1E_hash(modus,irrep,E0max,conf1h1I,indx1h1I,&
         n_int_I,n1h1I,cfgM)
    
    return
    
  end subroutine generate_1I1E_hash

!######################################################################
! builder_1I1E_hash: performs all the heavy lifting involved in the
!                    generation of the configurations with with two
!                    internal holes, one internal electron and one
!                    external electron
!######################################################################  
  subroutine builder_1I1E_hash(modus,irrep,E0max,conf1h1I,indx1h1I,&
       n_int_I,n1h1I,cfgM)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hparam
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep
    
    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! 1-hole and 2-hole configurations
    integer(is), intent(in)    :: n_int_I,n1h1I
    integer(ib), intent(in)    :: conf1h1I(n_int_I,2,n1h1I)
    integer(is), intent(in)    :: indx1h1I(2,n1h1I)
    
    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM
        
    ! Configurations and SOPs
    integer(ib)                :: conf_full(n_int,2)
    integer(ib)                :: sop_full(n_int,2)

    ! Difference configuration information
    integer(is)                :: Dw(nmo,2)
    integer(is)                :: ndiff

    ! DFT/MRCI energy selection criterion
    real(dp), allocatable      :: sum1h1I(:)
    real(dp)                   :: elim,sum

    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:)
    
    ! Everything else
    integer(is)                :: nmoI,n2h
    integer(is)                :: n,i,i1,k,iext,counter,Dwi,i2h
    logical                    :: ok,checksym

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    nmoI=cfgM%nmoI
    n2h=cfgM%n2h

    allocate(sum1h1I(n1h1I))
    sum1h1I=0.0d0

    allocate(ngen(n2h))
    ngen=0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n1I1E=0
    counter=0
    checksym=.true.
    
!----------------------------------------------------------------------
! If this is a DFT/MRCI calculation, then first determine the
! contribution from the 1H1I parts of the configurations to the
! energy selection sum
!----------------------------------------------------------------------
    if (ldftmrci) then

       ! Loop over the 1H1I configurations
       do n=1,n1h1I

          ! Compute the difference configuration
          call diffconf(conf1h1I(:,:,n),n_int_I,Dw,nmo,ndiff)

          ! Sum_i (F_ii^KS) Dw_i
          do i=1,ndiff

             ! MO index
             i1=cfgM%m2c(Dw(i,1))
             
             ! Delta w_i value
             Dwi=Dw(i,2)

             ! Sum the contribution
             sum1h1I(n)=sum1h1I(n)+moen(i1)*Dwi
          
          enddo
          
       enddo
       
    endif

!----------------------------------------------------------------------
! Generate the 1H1E configurations
!----------------------------------------------------------------------
    ! DFT/MRCI energy selection limit
    if (ldftmrci) elim=E0max+desel

    ! Loop over the 1H1I configurations
    do n=1,n1h1I

       ! Index of the 2-hole configuration from which the
       ! 1H1I configuration was derived
       i2h=indx1h1I(1,n)
       
       ! If this is a DFT/MRCI calculation, cycle if the current
       ! 1H1I configuration cannot possibly generate a 1I1E
       ! configuration satisfying the energy selection criterion
       if (ldftmrci .and. sum1h1I(n) > elim) cycle

       ! Loop over external MOs
       do iext=nmoI+1,lastvirt
       
          ! Block index
          k=(iext-1)/64+1

          ! Postion of the external MO within the kth block
          i=iext-(k-1)*64-1
          
          ! Full configuration
          conf_full=0_ib
          conf_full(1:n_int_I,:)=conf1h1I(:,:,n)
          conf_full(k,1)=ibset(conf_full(k,1),i)

          ! Full SOP
          sop_full=conf_to_sop(conf_full,n_int)

          ! Is this an allowable configuration?
          ok=allowable_conf(conf_full,sop_full,irrep,cfgM%m2c,E0max,&
               nomax,checksym)

          if (ok) then
             ! Update the total no. allowable configurations generated
             counter=counter+1
             ! Update the no. confs generated by the current hole conf
             ngen(i2h)=ngen(i2h)+1
          endif
          
          ! Fill in the creation operator array
          if (ok .and. modus == 1) then
             cfgM%a1I1E(1,counter)=indx1h1I(2,n)
             cfgM%a1I1E(2,counter)=iext
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! If we are just determining the number of allowable configurations,
! then save this number
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n1I1E=counter

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    if (modus == 1) then
       cfgM%off1I1E(1)=1
       do n=2,n2h+1
          cfgM%off1I1E(n)=cfgM%off1I1E(n-1)+ngen(n-1)
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(sum1h1I)
    deallocate(ngen)
    
    return
    
  end subroutine builder_1I1E_hash

!######################################################################
! generate_2I_hash: for a given irrep, generates all the allowable
!                   configurations with two internal holes, two
!                   internal electrons
!######################################################################
  subroutine generate_2I_hash(irrep,E0max,conf1h1I,indx1h1I,&
       n_int_I,n1h1I,cfgM,icvs)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max
    
    ! 1H1I configurations
    integer(is), intent(in)    :: n_int_I,n1h1I
    integer(ib), intent(in)    :: conf1h1I(n_int_I,2,n1h1I)
    integer(is), intent(in)    :: indx1h1I(2,n1h1I)
    
    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs
    
    ! Everything else
    integer(is)                :: modus

!----------------------------------------------------------------------
! Is this a CVS-MRCI calculation
!----------------------------------------------------------------------
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif
    
!----------------------------------------------------------------------
! Determine the number of allowable configurations of the correct
! symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_2I_hash(modus,irrep,E0max,conf1h1I,indx1h1I,n_int_I,&
         n1h1I,cfgM,icvs,lcvs)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Indices of the internal creation operators
    allocate(cfgM%a2I(2,cfgM%n2I))
    cfgM%a2I=0

    ! Offsets
    allocate(cfgM%off2I(cfgM%n2h+1))
    cfgM%off2I=0

!----------------------------------------------------------------------
! Fill in the internal creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off2I array has the dimension (n2h+1), even though
! not all 2-hole configurations may generate allowable
! 2I configurations. However, this is OK due to how the offset array
! is used when looping over the 2I configurations generated by each
! 2-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_2I_hash(modus,irrep,E0max,conf1h1I,indx1h1I,n_int_I,&
         n1h1I,cfgM,icvs,lcvs)
    
    return
    
  end subroutine generate_2I_hash

!######################################################################
! builder_2I_hash: performs all the heavy lifting involved in the
!                  generation of the configurations with with two
!                  internal holes, two internal electrons
!######################################################################  
  subroutine builder_2I_hash(modus,irrep,E0max,conf1h1I,indx1h1I,&
       n_int_I,n1h1I,cfgM,icvs,lcvs)
    
    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hparam
    use dethash

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Symmetry of the subspace
    integer(is), intent(in)    :: irrep

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max
        
    ! 1H1I configurations
    integer(is), intent(in)    :: n_int_I,n1h1I
    integer(ib), intent(in)    :: conf1h1I(n_int_I,2,n1h1I)
    integer(is), intent(in)    :: indx1h1I(2,n1h1I)

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical, intent(in)        :: lcvs
    
    ! Configurations and SOPs
    integer(ib)                :: conf_full(n_int,2)
    integer(ib)                :: sop_full(n_int,2)
    integer(ib)                :: conftmp(n_int,2)
    
    ! Difference configuration information
    integer(is)                :: Dw(nmo,2)
    integer(is)                :: ndiff

    ! DFT/MRCI energy selection criterion
    real(dp), allocatable      :: sum1h1I(:)
    real(dp)                   :: elim,sum
    
    ! Hash table
    type(dhtbl)                :: h
    integer(is)                :: initial_size
    integer(is)                :: nold
    integer(ib)                :: key(n_int,2)

    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:)
    
    ! Everything else
    integer(is)                :: nmoI,n1h,n2h
    integer(is)                :: n,i,i1,k,iint,counter,Dwi,i2h,iref
    integer(is)                :: ioff,i1I,old2h,imo,ilow
    logical                    :: ok,unique,checksym

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    nmoI=cfgM%nmoI
    n1h=cfgM%n1h
    n2h=cfgM%n2h

    allocate(sum1h1I(n1h1I))
    sum1h1I=0.0d0

    allocate(ngen(n2h))
    ngen=0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n2I=0
    counter=0
    checksym=.true.

!----------------------------------------------------------------------
! Initialise the hash table and insert the reference and 1I
! configurations
!----------------------------------------------------------------------
    ! Initialisation
    initial_size=2**15
    call h%initialise_table(initial_size)

    ! Insertion of the ref confs
    do iref=1,cfgM%n0h
       key=0_ib
       key(1:n_int_I,:)=cfgM%conf0h(:,:,iref)
       call h%insert_key(key)
    enddo

    ! Insertion of the 1I confs
    i1I=0
    ! Loop over 1-hole configurations
    do n=1,n1h
       
       ! Loop over the 1I configurations generated by the current
       ! 1-hole configuration
       do ioff=cfgM%off1I(n),cfgM%off1I(n+1)-1

          ! Increment the configuration counter
          i1I=i1I+1

          ! Construct the configuration
          imo=cfgM%a1I(i1I)
          conftmp=0_ib
          conftmp(1:n_int_I,:)=cfgM%conf1h(:,:,n)
          key=create_electron(conftmp,n_int,imo)

          ! Insert the configuration into the hash table
          call h%insert_key(key)
        
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! If this is a DFT/MRCI calculation, then first determine the
! contribution from the 1H1I parts of the configurations to the
! energy selection sum
!----------------------------------------------------------------------
    if (ldftmrci) then

       ! Loop over the 1H1I configurations
       do n=1,n1h1I

          ! Compute the difference configuration
          call diffconf(conf1h1I(:,:,n),n_int_I,Dw,nmo,ndiff)

          ! Sum_i (F_ii^KS) Dw_i
          do i=1,ndiff

             ! MO index
             i1=cfgM%m2c(Dw(i,1))
             
             ! Delta w_i value
             Dwi=Dw(i,2)

             ! Sum the contribution
             sum1h1I(n)=sum1h1I(n)+moen(i1)*Dwi
          
          enddo
          
       enddo
       
    endif

!----------------------------------------------------------------------
! Generate the 2I configurations
!----------------------------------------------------------------------
    ! Number of confs currently stored in the hash table
    nold=h%n_keys_stored

    ! DFT/MRCI energy selection limit
    if (ldftmrci) elim=E0max+desel

    ! Loop over the 1H1I configurations
    do n=1,n1h1I
       
       ! Index of the 2-hole configuration from which the
       ! 1H1I configuration was derived
       i2h=indx1h1I(1,n)

       ! Index of the lowest-lying unoccupied or singly-occupied MO
       ilow=lowest_particle_index(conf1h1I(:,:,n),n_int_I)

       ! If this is a DFT/MRCI calculation, cycle if the current
       ! 1H1I configuration cannot possibly generate an
       ! a 2I configuration satisfying the energy selection criterion
       if (ldftmrci) then
          sum=sum1h1I(n)+moen(cfgM%m2c(ilow))
          if (sum > elim) cycle
       endif
       
       ! Loop over internal MOs
       do iint=1,nmoI

          ! Cycle if this is a CVS-MRCI calculation and we are creating
          ! an electron in a flagged core MO
          if (lcvs .and. icvs(cfgM%m2c(iint)) == 1) cycle
          
          ! Block index
          k=(iint-1)/64+1

          ! Postion of the external MO within the kth block
          i=iint-(k-1)*64-1

          ! Cycle if the internal MO is doubly-occupied in the
          ! 1H1I configuration
          if (btest(conf1h1I(k,2,n),i)) cycle
          
          ! Full configuration
          conf_full=0_ib
          conf_full(1:n_int_I,:)=conf1h1I(:,:,n)
          if (btest(conf_full(k,1),i)) then
             ! Creation of a doubly-occupied internal MO
             conf_full(k,2)=ibset(conf_full(k,2),i)
          else
             ! Creation of a singly-occupied internal MO
             conf_full(k,1)=ibset(conf_full(k,1),i)
          endif
          
          ! Full SOP
          sop_full=conf_to_sop(conf_full,n_int)

          ! Is this an allowable configuration?
          ok=allowable_conf(conf_full,sop_full,irrep,cfgM%m2c,E0max,&
               nomax,checksym)

          ! Is this a unique configuration?
          if (ok) then

             ! Attempt an insertion into the hash table
             call h%insert_key(conf_full)

             ! If the number of stored keys has increased, then
             ! this is a new configuration
             if (h%n_keys_stored > nold) then
                unique=.true.
                nold=h%n_keys_stored
             else
                unique=.false.
             endif
             
          endif

          if (ok .and. unique) then
             ! Update the total no. allowable configurations generated
             counter=counter+1
             ! Update the no. confs generated by the current hole conf
             ngen(i2h)=ngen(i2h)+1
          endif
             
          ! Cycle if we are just determining the no. allowable
          ! configurations...
          if (modus == 0) cycle

          ! ...else fill in the creation operator array
          if (ok .and. unique) then
             cfgM%a2I(1,counter)=indx1h1I(2,n)
             cfgM%a2I(2,counter)=iint
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Delete the hash table
!----------------------------------------------------------------------
    call h%delete_table
    
!----------------------------------------------------------------------
! If we are just determining the number of allowable configurations,
! then save this number
!----------------------------------------------------------------------
    if (modus == 0) cfgM%n2I=counter

!----------------------------------------------------------------------
! Fill in the offset array
!----------------------------------------------------------------------
    if (modus == 1) then
       cfgM%off2I(1)=1
       do n=2,n2h+1
          cfgM%off2I(n)=cfgM%off2I(n-1)+ngen(n-1)
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(sum1h1I)
    deallocate(ngen)
    
    return
    
  end subroutine builder_2I_hash
    
!######################################################################
! allowable_conf: Determines whether a given configuration is
!                 is allowable. That is, if it satisfies the following
!                 criteria:
!
!                 (1) Generates the correct irrep
!
!                 (2) Has an excitation degree less than the maximum
!                     allowable value relative to the base
!                     configuration
!
!                 (3) Doesn't have more than the maximum allowed
!                     number of open shells
!
!                 (4) For a DFT/MRCI run, satisfies the energy
!                     selection criterion
!######################################################################
  function allowable_conf(conf,sop,irrep,m2c,E0max,nolim,checksym) &
       result(allowable)

    use constants
    use bitglobal
    use mrciutils
    use hparam
    
    implicit none

    ! Function result
    logical                 :: allowable

    ! Configuration and SOP bitstrings
    integer(ib), intent(in) :: conf(n_int,2),sop(n_int,2)

    ! Symmetry of the subspace
    integer(is), intent(in) :: irrep

    ! MO mapping array
    integer(is), intent(in) :: m2c(nmo)

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)    :: E0max

    ! Maximum number of open shells
    integer(is), intent(in) :: nolim

    ! Perform the symmetry check?
    logical, intent(in)     :: checksym
    
    ! Difference configuration information
    integer(is)             :: ndiff
    integer(is)             :: Dw(nmo,2)
    
    ! Everything else
    integer(is)             :: isym,nopen,nexci,i,i1,Dwi
    real(dp)                :: sum
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    allowable=.true.

!----------------------------------------------------------------------
! Symmetry
!----------------------------------------------------------------------
    if (checksym) then
       
       isym=sop_sym_mrci(sop,m2c)
       
       if (isym /= irrep) then
          allowable=.false.
          return
       endif

    endif
    
!----------------------------------------------------------------------
! Number of open shells
!----------------------------------------------------------------------
    nopen=sop_nopen(sop,n_int)

    if (nopen > nolim) then
       allowable=.false.
       return
    endif

!----------------------------------------------------------------------
! Excitation degree relative to the base configuration
!----------------------------------------------------------------------
    nexci=exc_degree_conf(conf,conf0,n_int)
            
    if (nexci > nexmax) then
       allowable=.false.
       return
    endif
    
!----------------------------------------------------------------------
! DFT/MRCI energy selection criterion
!----------------------------------------------------------------------
    if (ldftmrci) then

       ! Generate the difference configuration information relative
       ! to the base configuration
       call diffconf(conf,n_int,Dw,nmo,ndiff)

       ! Sum_i (F_ii^KS) Dw_i
       sum=0.0d0
       do i=1,ndiff

          ! MO index
          i1=m2c(Dw(i,1))
    
          ! Delta w_i value
          Dwi=Dw(i,2)

          ! Sum the contribution
          sum=sum+moen(i1)*Dwi
          
       enddo
       
       ! Does this configuration satisfy the energy selection
       ! criterion?
       if (sum > E0max + desel) then
          allowable=.false.
          return
       endif
       
    endif
    
    return
    
  end function allowable_conf

!######################################################################
! lowest_particle_index: for a given configuration, returns the index
!                        of the lowest-lying unoccupied or singly-
!                        occupied MO
!######################################################################
  function lowest_particle_index(conf,ldc) result(indx)

    use constants
    use bitglobal
    
    implicit none

    ! Function result
    integer(is)             :: indx

    ! Input configuration
    integer(is), intent(in) :: ldc
    integer(ib), intent(in) :: conf(ldc,2)

    ! Everything else
    integer(is)             :: k,i

    ! Initialisation
    indx=0
    
    ! Loop over blocks
    do k=1,ldc

       ! Cycle if this block doesn't contain a singly-occupied
       ! or unoccupied MO
       if (popcnt(conf(k,2)) == 64) then
          indx=indx+64
          cycle
       endif
          
       ! Index of the lowest lying singly-occupied
       ! or unoccupied MO
       indx=indx+trailz(not(conf(k,2)))+1
       exit
       
    enddo
    
    return
    
  end function lowest_particle_index
  
!######################################################################
! esum: given a configuration conf, returns Sum_p F_pp^(KS) Delta w_p
!       relative to the base configuration
!######################################################################
  function esum(conf,ldc,m2c)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    ! Function result
    real(dp) :: esum

    ! Configuration
    integer(is), intent(in) :: ldc
    integer(ib), intent(in) :: conf(ldc,2)

    ! Difference configuration information
    integer(is)             :: Dw(nmo,2)
    integer(is)             :: ndiff

    ! MO mapping array
    integer(is), intent(in) :: m2c(nmo)
    
    ! Everything else
    integer(is)             :: i,i1,Dwi
    
    !
    ! Compute the difference configuration relative to the base
    ! configuration
    !
    call diffconf(conf,ldc,Dw,nmo,ndiff)

    !
    ! Sum_i (F_ii^KS) Dw_i
    !
    esum=0.0d0
    do i=1,ndiff

       ! MO index
       i1=m2c(Dw(i,1))
       
       ! Delta w_i value
       Dwi=Dw(i,2)
       
       ! Sum the contribution
       esum=esum+moen(i1)*Dwi
       
    enddo
    
    return
    
  end function esum
  
!######################################################################
! generate_1E_confs: for all irreps, generates all the allowable
!                    configurations with one internal hole and one
!                    external electron
!######################################################################
  subroutine generate_1E_confs(E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)
    
    ! Everything else
    integer(is)                :: modus,irrep

!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! each symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_1E(modus,E0max,cfgM)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Indices of the external creation operators
       allocate(cfgM(irrep)%a1E(cfgM(irrep)%n1E))
       cfgM(irrep)%a1E=0
       
       ! Offsets
       allocate(cfgM(irrep)%off1E(cfgM(irrep)%n1h+1))
       cfgM(irrep)%off1E=0

    enddo

!----------------------------------------------------------------------
! Fill in the external creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off1E array has the dimension (n1h+1), even though
! not all 1-hole configurations may generate allowable
! 1E configurations. However, this is OK due to how the offset array
! is used when looping over the 1E configurations generated by each
! 1-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_1E(modus,E0max,cfgM)
    
    return
    
  end subroutine generate_1E_confs

!######################################################################
! 1E_builder: performs all the heavy lifting involved in the
!             generation of the configurations with one internal hole
!             and one external electron
!######################################################################
  subroutine builder_1E(modus,E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hparam
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)
    
    ! Full SOP
    integer(ib)                :: sop(n_int,2)

    ! Base configuration, internal MOs only
    integer(ib), allocatable   :: baseconf(:,:)
    
    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:,:)

    ! Sums of KS orbital energies multipled by occupation
    ! differences relative to the base conf
    real(dp)                   :: esumR,esum1H,esum1E

    ! Irreps generated the various different configurations
    integer(is)                :: irrepR,irrep1H,irrep1E
    integer(ib)                :: ib1,ib2
    
    ! Everything else
    integer(is)                :: n_int_I,nmoI,nR,n1h
    integer(is)                :: n,i1h,ia1h,iext,irrep
    integer(is)                :: counter(0:nirrep-1)
    integer(is)                :: nexci1H

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n1h=cfgM(0)%n1h
    nR=cfgM(0)%nR
    n_int_I=cfgM(0)%n_int_I
    nmoI=cfgM(0)%nmoI

    allocate(ngen(n1h,0:nirrep-1))
    ngen=0

    allocate(baseconf(n_int_I,2))
    baseconf=0_ib
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) then
       do irrep=0,nirrep-1
          cfgM(irrep)%n1E=0
       enddo
    endif

    counter=0

    baseconf=conf0(1:n_int_I,:)
    
!----------------------------------------------------------------------    
! Generate the 1E configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do n=1,nR

       ! Sum_p F_pp^(KS) Delta W_p for the reference configuration
       esumR=esum(cfgM(0)%confR(:,:,n),n_int_I,cfgM(0)%m2c)

       ! Irrep generated by the reference configuration
       sop=0_ib
       sop(1:n_int_I,:)=cfgM(0)%sopR(:,:,n)
       irrepR=sop_sym_mrci(sop,cfgM(0)%m2c)

       ! Loop over the 1-hole configurations generated by the current
       ! reference configuration
       do i1h=cfgM(0)%off1h(n),cfgM(0)%off1h(n+1)-1

          ! 1-hole annihilation operator index
          ia1h=cfgM(0)%a1h(i1h)

          ! Sum_p F_pp^(KS) Delta W_p for the 1-hole configuration
          esum1H=esumR-moen(cfgM(0)%m2c(ia1h))

          ! If this is a DFT/MRCI calculation, then skip this
          ! 1-hole configuration cannot possibly generate 1E
          ! configurations satisfying the energy-based
          ! selection criterion
          if (ldftmrci .and. esum1H > E0max + desel) cycle

          ! Cycle if the excitation degree of the 1E configurations
          ! generated by this 1-hole configuration will be above
          ! threshold
          nexci1H=exc_degree_conf(cfgM(0)%conf1h(:,:,i1h),baseconf,&
               n_int_I)
          if (nexci1H+1 > nexmax) cycle
          
          ! Irrep generated by the 1-hole configuration
          ib1=irrepR
          ib2=ieor(ib1,mosym(cfgM(0)%m2c(ia1h)))
          irrep1H=ib2

          ! Loop over external MOs
          do iext=nmoI+1,lastvirt

             ! Sum_p F_pp^(KS) Delta W_p for the 1E configuration
             esum1E=esum1H+moen(cfgM(0)%m2c(iext))

             ! If this is a DFT/MRCI calculation, then skip this
             ! configuration if it doesn't satisfy the energy-based
             ! selection criterion
             if (ldftmrci .and. esum1E > E0max + desel) cycle

             ! Irrep generated by the 1I1E configuration
             ib1=irrep1H
             ib2=ieor(ib1,mosym(cfgM(0)%m2c(iext)))
             irrep1E=ib2

             ! Update the no. 1E confs
             if (modus == 0) cfgM(irrep1E)%n1E=cfgM(irrep1E)%n1E+1
             
             if (modus == 1) then

                ! Fill in the creation operator array
                counter(irrep1E)=counter(irrep1E)+1
                cfgM(irrep1E)%a1E(counter(irrep1E))=iext

                ! Update the no. confs generated by the current hole
                ! conf
                ngen(i1h,irrep1E)=ngen(i1h,irrep1E)+1
                
             endif
             
          enddo
          
       enddo
          
    enddo

!----------------------------------------------------------------------
! Fill in the offset arrays
!----------------------------------------------------------------------
    if (modus == 1) then
       do irrep=0,nirrep-1
          cfgM(irrep)%off1E(1)=1
          do i1h=2,n1h+1
             cfgM(irrep)%off1E(i1h)=cfgM(irrep)%off1E(i1h-1) &
                  +ngen(i1h-1,irrep)
          enddo
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)
    deallocate(baseconf)
    
    return
    
  end subroutine builder_1E

!######################################################################
! generate_2E_confs: for all irreps, generates all the allowable
!                    configurations with two internal holes and two
!                    external electrons
!######################################################################
  subroutine generate_2E_confs(E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)
    
    ! Everything else
    integer(is)                :: modus,irrep

!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! each symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_2E(modus,E0max,cfgM)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Indices of the external creation operators
       allocate(cfgM(irrep)%a2E(2,cfgM(irrep)%n2E))
       cfgM(irrep)%a2E=0
       
       ! Offsets
       allocate(cfgM(irrep)%off2E(cfgM(irrep)%n2h+1))
       cfgM(irrep)%off2E=0

    enddo

!----------------------------------------------------------------------
! Fill in the external creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off2E array has the dimension (n2h+1), even though
! not all 2-hole configurations may generate allowable
! ext configurations. However, this is OK due to how the offset array
! is used when looping over the 2E configurations generated by each
! 2-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_2E(modus,E0max,cfgM)
    
    return
    
  end subroutine generate_2E_confs

!######################################################################
! 2E_builder: performs all the heavy lifting involved in the
!             generation of the configurations with two internal holes
!             and two external electrons
!######################################################################
  subroutine builder_2E(modus,E0max,cfgM)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hparam
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)
    
    ! Full SOP
    integer(ib)                :: sop(n_int,2)

    ! Base configuration, internal MOs only
    integer(ib), allocatable   :: baseconf(:,:)
    
    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:,:)

    ! Sums of KS orbital energies multipled by occupation
    ! differences relative to the base conf
    real(dp)                   :: esumR,esum2H,esum1H1E,esum2E

    ! Irreps generated the various different configurations
    integer(is)                :: irrepR,irrep2H,irrep1H1E,irrep2E
    integer(ib)                :: ib1,ib2
    
    ! Everything else
    integer(is)                :: n_int_I,nmoI,nR,n2h
    integer(is)                :: n,i2h,ia2h,ja2h,iext1,iext2,k,i,irrep
    integer(is)                :: counter(0:nirrep-1)
    integer(is)                :: nexci2H

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n2h=cfgM(0)%n2h
    nR=cfgM(0)%nR
    n_int_I=cfgM(0)%n_int_I
    nmoI=cfgM(0)%nmoI

    allocate(ngen(n2h,0:nirrep-1))
    ngen=0

    allocate(baseconf(n_int_I,2))
    baseconf=0_ib

!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) then
       do irrep=0,nirrep-1
          cfgM(irrep)%n2E=0
       enddo
    endif

    counter=0

    baseconf=conf0(1:n_int_I,:)

!----------------------------------------------------------------------    
! Generate the 2E configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do n=1,nR

       ! Sum_p F_pp^(KS) Delta W_p for the reference configuration
       esumR=esum(cfgM(0)%confR(:,:,n),n_int_I,cfgM(0)%m2c)

       ! Irrep generated by the reference configuration
       sop=0_ib
       sop(1:n_int_I,:)=cfgM(0)%sopR(:,:,n)
       irrepR=sop_sym_mrci(sop,cfgM(0)%m2c)

       ! Loop over the 2-hole configurations generated by the current
       ! reference configuration
       do i2h=cfgM(0)%off2h(n),cfgM(0)%off2h(n+1)-1

          ! 2-hole annihilation operator indices
          ia2h=cfgM(0)%a2h(1,i2h)
          ja2h=cfgM(0)%a2h(2,i2h)

          ! Sum_p F_pp^(KS) Delta W_p for the 2-hole configuration
          esum2H=esumR-moen(cfgM(0)%m2c(ia2h))-moen(cfgM(0)%m2c(ja2h))

          ! If this is a DFT/MRCI calculation, then skip this
          ! 2-hole configuration cannot possibly generate 2E
          ! configurations satisfying the energy-based
          ! selection criterion
          if (ldftmrci .and. esum2H > E0max + desel) cycle

          ! Cycle if the excitation degree of the 2E configurations
          ! generated by this 2-hole configuration will be above
          ! threshold
          nexci2H=exc_degree_conf(cfgM(0)%conf2h(:,:,i2h),baseconf,&
               n_int_I)
          if (nexci2H+2 > nexmax) cycle

          ! Irrep generated by the 2-hole configuration
          ib1=irrepR
          ib2=ieor(ib1,mosym(cfgM(0)%m2c(ia2h)))
          ib1=ieor(ib2,mosym(cfgM(0)%m2c(ja2h)))
          irrep2H=ib1
          
          ! Loop over the first creation operator
          do iext1=nmoI+1,lastvirt

             ! Sum_p F_pp^(KS) Delta W_p for the 1H1E configuration
             esum1H1E=esum2H+moen(cfgM(0)%m2c(iext1))

             ! If this is a DFT/MRCI calculation, then skip this
             ! 1H1E configuration cannot possibly generate 2E
             ! configurations satisfying the energy-based
             ! selection criterion
             if (ldftmrci .and. esum1H1E > E0max + desel) cycle

             ! Irrep generated by the 1H1E configuration
             ib1=irrep2H
             ib2=ieor(ib1,mosym(cfgM(0)%m2c(iext1)))
             irrep1H1E=ib2

             ! Loop over the second creation operator
             do iext2=iext1,lastvirt

                ! Sum_p F_pp^(KS) Delta W_p for the 2E configuration
                esum2E=esum1H1E+moen(cfgM(0)%m2c(iext2))
                
                ! If this is a DFT/MRCI calculation, then skip this
                ! 2E configuration if it does not satisfy the
                ! energy-basedselection criterion
                if (ldftmrci .and. esum2E > E0max + desel) cycle
                
                ! Irrep generated by the 1H1E configuration
                ib1=irrep1H1E
                ib2=ieor(ib1,mosym(cfgM(0)%m2c(iext2)))
                irrep2E=ib2

                ! Update the no. 2E confs
                if (modus == 0) cfgM(irrep2E)%n2E=cfgM(irrep2E)%n2E+1

                if (modus == 1) then

                   ! Fill in the creation operator array
                   counter(irrep2E)=counter(irrep2E)+1
                   cfgM(irrep2E)%a2E(1,counter(irrep2E))=iext1
                   cfgM(irrep2E)%a2E(2,counter(irrep2E))=iext2

                   ! Update the no. confs generated by the current hole
                   ! conf
                   ngen(i2h,irrep2E)=ngen(i2h,irrep2E)+1
                   
                endif
                
             enddo
                
          enddo
             
       enddo
          
    enddo

!----------------------------------------------------------------------
! Fill in the offset arrays
!----------------------------------------------------------------------
    if (modus == 1) then
       do irrep=0,nirrep-1
          cfgM(irrep)%off2E(1)=1
          do i2h=2,n2h+1
             cfgM(irrep)%off2E(i2h)=cfgM(irrep)%off2E(i2h-1) &
                  +ngen(i2h-1,irrep)
          enddo
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ngen)
    deallocate(baseconf)
    
    return

  end subroutine builder_2E

!######################################################################
! generate_1I_confs: for all irreps, generates all the allowable
!                    configurations with one internal hole and one
!                    external electron
!######################################################################
  subroutine generate_1I_confs(E0max,cfgM,icvs)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs
    
    ! Everything else
    integer(is)                :: modus,irrep

!----------------------------------------------------------------------
! Is this a CVS-MRCI calculation
!----------------------------------------------------------------------
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif
    
!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! each symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_1I(modus,E0max,cfgM,icvs,lcvs)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Indices of the external creation operators
       allocate(cfgM(irrep)%a1I(cfgM(irrep)%n1I))
       cfgM(irrep)%a1I=0
       
       ! Offsets
       allocate(cfgM(irrep)%off1I(cfgM(irrep)%n1h+1))
       cfgM(irrep)%off1I=0

    enddo

!----------------------------------------------------------------------
! Fill in the external creation operator and offset arrays
!----------------------------------------------------------------------
! Note that the off1I array has the dimension (n1h+1), even though
! not all 1-hole configurations may generate allowable
! 1I configurations. However, this is OK due to how the offset array
! is used when looping over the 1I configurations generated by each
! 1-hole configuration.
!----------------------------------------------------------------------
    modus=1
    call builder_1I(modus,E0max,cfgM,icvs,lcvs)
    
    return
    
  end subroutine generate_1I_confs

!######################################################################
! 1I_builder: performs all the heavy lifting involved in the
!             generation of the configurations with one internal hole
!             and one internal electron
!######################################################################
  subroutine builder_1I(modus,E0max,cfgM,icvs,lcvs)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use hparam
    use iomod
    use dethash
    
    implicit none

    ! Mode of operation: modus=0 <-> determine the number of allowable
    !                                configurations
    !                    modus=1 <-> build all of the allowable
    !                                configurations
    integer(is), intent(in)    :: modus

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical, intent(in)        :: lcvs
    
    ! Base configuration, internal MOs only
    integer(ib), allocatable   :: baseconf(:,:)
    
    ! No. confs generated by each hole conf
    integer(is), allocatable   :: ngen(:,:)

    ! Allowable internal creation operator indices
    integer(is)                :: ncreate1
    integer(is), allocatable   :: icreate1(:)
    
    ! Difference configuration information
    integer(is), allocatable   :: rhp(:,:)
    integer(is), parameter     :: maxexci=2
    integer(is)                :: hlist(maxexci),plist(maxexci)
    
    ! Work bit strings
    integer(ib), allocatable   :: sop(:,:),conf(:,:),confI(:,:)

    ! Orbital classes
    integer(is)                :: socc(nmo),docc(nmo),unocc(nmo)
    integer(is)                :: nopen,nsocc,ndocc,nunocc
    
    ! 1-hole SOPs
    integer(ib), allocatable   :: sop1h(:,:,:)
        
    ! Sums of KS orbital energies multipled by occupation
    ! differences relative to the base conf
    real(dp)                   :: esumR,esum1H,esum1I

    ! Irreps generated the various different configurations
    integer(is)                :: irrepR,irrep1H,irrep1I
    integer(ib)                :: ib1,ib2

    ! Hash table
    type(dhtbl)                :: h
    integer(is)                :: initial_size
    integer(is)                :: nold
    
    ! Everything else
    integer(is)                :: n,np,i1h,irrep,n1h,nR,n_int_I,nmoI
    integer(is)                :: ia1h,iint,k,i,j
    integer(is)                :: counter(0:nirrep-1)
    integer(is)                :: nacR,nac1H,nac1I
    integer(is)                :: nexci,iaR,jaR,icR,jcR
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n1h=cfgM(0)%n1h
    nR=cfgM(0)%nR
    n_int_I=cfgM(0)%n_int_I
    nmoI=cfgM(0)%nmoI

    allocate(icreate1(nmoI))
    icreate1=0
    
    allocate(ngen(n1h,0:nirrep-1))
    ngen=0

    allocate(baseconf(n_int_I,2))
    baseconf=0_ib

    allocate(sop(n_int,2))
    sop=0_ib

    allocate(sop1h(n_int_I,2,n1h))
    sop1h=0_ib

    allocate(confI(n_int_I,2))
    confI=0_ib

    allocate(conf(n_int,2))
    conf=0_ib
    
    allocate(rhp(5,nR))
    rhp=0
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    if (modus == 0) then
       do irrep=0,nirrep-1
          cfgM(irrep)%n1I=0
       enddo
    endif

    counter=0

    baseconf=conf0(1:n_int_I,:)

!----------------------------------------------------------------------
! For now, we will check for duplicate reference configurations by
! using a hash table-based approach 
!----------------------------------------------------------------------
! Initialise the hash table and insert the reference configurations
!----------------------------------------------------------------------
    ! Initialisation
    initial_size=2**15
    call h%initialise_table(initial_size)

    ! Insertion of the ref confs
    do n=1,cfgM(0)%nR
       conf=0_ib
       conf(1:n_int_I,:)=cfgM(0)%confR(:,:,n)
       call h%insert_key(conf)
    enddo

    ! No. keys inserted
    nold=h%n_keys_stored

!----------------------------------------------------------------------
! Generate the 1-hole SOPs
!----------------------------------------------------------------------    
    ! Loop over 1-hole configurations
    do i1h=1,n1h

       ! Generate the next 1-hole SOP
       confI=cfgM(0)%conf1h(:,:,i1h)
       sop1h(:,:,i1h)=conf_to_sop(confI,n_int_I)

    enddo
    
!----------------------------------------------------------------------
! Generate the 1I configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do n=1,nR

       ! Sum_p F_pp^(KS) Delta W_p for the reference configuration
       esumR=esum(cfgM(0)%confR(:,:,n),n_int_I,cfgM(0)%m2c)
       
       ! Irrep generated by the reference configuration
       sop=0_ib
       sop(1:n_int_I,:)=cfgM(0)%sopR(:,:,n)
       irrepR=sop_sym_mrci(sop,cfgM(0)%m2c)

       ! Number of creation and annihilation operators linking the
       ! reference and base configurations
       nacR=n_create_annihilate(cfgM(0)%confR(:,:,n),baseconf,n_int_I)

       ! Initialise the array of inter-ref-conf annihilation/creation
       ! operator indices
       rhp=0

       ! Loop over preceding reference configurations
       do np=1,n-1

          ! Determine the excitation degree between the two reference
          ! configurations
          nexci=exc_degree_conf(cfgM(0)%confR(:,:,n),&
               cfgM(0)%confR(:,:,np),n_int_I)

          ! Cycle if the excitation degree is greater than 2
          if (nexci > 2) cycle

          ! Determine the indices of the creation/annihilation
          ! operators linking the two ref confs
          call get_exci_indices(cfgM(0)%confR(:,:,n),&
               cfgM(0)%confR(:,:,np),n_int_I,hlist(1:nexci),&
               plist(1:nexci),nexci)
          
          ! Excitation degree
          rhp(1,np)=nexci

          ! Creation/annihilation operator indices in ascending order
          if (nexci == 1) then
             rhp(2,np)=hlist(1)
             rhp(4,np)=plist(1)
          else if (nexci == 2) then
             rhp(2,np)=min(hlist(1),hlist(2))
             rhp(3,np)=max(hlist(1),hlist(2))
             rhp(4,np)=min(plist(1),plist(2))
             rhp(5,np)=max(plist(1),plist(2))
          endif
             
       enddo

       ! Loop over the 1-hole configurations generated by the current
       ! reference configuration
       do i1h=cfgM(0)%off1h(n),cfgM(0)%off1h(n+1)-1

          ! 1-hole annihilation operator index
          ia1h=cfgM(0)%a1h(i1h)

          ! Sum_p F_pp^(KS) Delta W_p for the 1-hole configuration
          esum1H=esumR-moen(cfgM(0)%m2c(ia1h))

          ! Irrep generated by the 1-hole configuration
          ib1=irrepR
          ib2=ieor(ib1,mosym(cfgM(0)%m2c(ia1h)))
          irrep1H=ib2

          ! Number of creation and annihilation operators linking the
          ! 1-hole and base configurations
          nac1H=nacR
          k=(ia1h-1)/64+1
          i=ia1h-(k-1)*64-1
          if (iocc0(ia1h) == 0) then
             ! Unoccupied MO in the base conf
             nac1H=nac1H-1
          else if (iocc0(ia1h) == 1) then
             ! Singly-occupied MO in the base conf
             if (btest(cfgM(0)%confR(k,2,n),i)) then
                ! Doubly-occupied MO in the ref conf
                nac1H=nac1H-1
             else
                ! Singly-occupied MO in the ref conf
                nac1H=nac1H+1
             endif
          else if (iocc0(ia1h) == 2) then
             ! Doubly-occupied MO in the base conf
             nac1H=nac1H+1
          endif

          ! Initialise the creation operator array
          ncreate1=nmoI
          icreate1=1
          
          ! Remove the 1-hole annihilation operator index from the
          ! list of allowed creation operators
          icreate1(ia1h)=0
          ncreate1=ncreate1-1
          
          ! Remove the doubly-occupied MOs from the list
          ! of allowed creation operators
          call sop_docc_list(sop1h(:,:,i1h),n_int_I,docc,nmo,ndocc)
          do j=1,ndocc
             if (icreate1(docc(j)) /= 0) then
                icreate1(docc(j))=0
                ncreate1=ncreate1-1
             endif
          enddo
          
          !
          ! Flag creation operators as forbidden if they will
          ! produce duplicate confs based on the excitations linking
          ! the current ref conf to the preceding ref confs
          !
          ! For ref confs linked by a single excitation:
          ! The creation operator of the excitation will yield
          ! a duplicate conf when acting on the 1-hole conf
          !
          ! For ref confs linked by a double excitation:
          ! If there exists a matched annihilation operator between
          ! the 1-hole conf and the excitation, then the creation
          ! operators of the excitation will yield duplicate confs
          ! when acting on the 1-hole conf
          !
          ! Loop over preceding ref confs
          do np=1,n-1
             
             ! Cycle if the excitation degree is greater than 2
             nexci=rhp(1,np)
             if (nexci > 2) cycle

             ! Creation and annhilation operators linking the two
             ! ref confs
             iaR=rhp(2,np)
             jaR=rhp(3,np)
             icR=rhp(4,np)
             jcR=rhp(5,np)
             
             ! Single excitation linking the ref confs
             if (nexci == 1) then
                if (icreate1(icR) /= 0) then
                   icreate1(icR)=0
                   ncreate1=ncreate1-1
                endif
             endif

             ! Double excitation linking the ref confs
             if (nexci == 2) then
                if (ia1h == iaR .or. ia1h == jaR) then
                   ! Matched annihilation operators
                   if (icreate1(icR) /= 0) then
                      icreate1(icR)=0
                      ncreate1=ncreate1-1
                   endif
                   if (icreate1(jcR) /= 0) then
                      icreate1(jcR)=0
                      ncreate1=ncreate1-1
                   endif
                endif
             endif
             
          enddo

          ! Skip this 1-hole conf if all creation operators are
          ! flagged as disallowed
          if (ncreate1 == 0) cycle

          ! Loop over internal creation operators
          do iint=1,nmoI
             
             ! Cycle if this this creation operator will yield
             ! a duplicate conf
             if (icreate1(iint) == 0) cycle
             
             ! Cycle if this is a CVS-MRCI calculation and we are creating
             ! an electron in a flagged core MO
             if (lcvs .and. icvs(cfgM(0)%m2c(iint)) == 1) cycle

             ! Sum_p F_pp^(KS) Delta W_p for the 1I configuration
             esum1I=esum1H+moen(cfgM(0)%m2c(iint))
             
             ! If this is a DFT/MRCI calculation, then skip this
             ! 1I configuration if it does not satisfy the
             ! energy-basedselection criterion
             if (ldftmrci .and. esum1I > E0max + desel) cycle
             
             ! Irrep generated by the 1I configuration
             ib1=irrep1H
             ib2=ieor(ib1,mosym(cfgM(0)%m2c(iint)))
             irrep1I=ib2

             ! Block index
             k=(iint-1)/64+1

             ! Postion of the external MO within the kth block
             i=iint-(k-1)*64-1

             ! Number of creation and annihilation operators linking the
             ! 1I and base configurations
             nac1I=nac1H
             if (iocc0(iint) == 0) then
                ! Unoccupied MO in the base conf
                nac1I=nac1I+1
             else if (iocc0(iint) == 1) then
                ! Singly-occupied MO in the base conf
                if (btest(cfgM(0)%conf1h(k,1,i1h),i)) then
                   ! Singly-occupied MO in the 1-hole conf
                   nac1I=nac1I+11
                else
                   ! Unoccupied MO in the 1-hole conf
                   nac1I=nac1I-1
                endif
             else if (iocc0(iint) == 2) then
                ! Doubly-occupied MO in the base conf
                nac1I=nac1I-1
             endif

             ! Cycle if the excitation degree relative to the base
             ! conf is too high
             if (nac1I/2 > nexmax) cycle

             ! Full configuration
             conf=0_ib
             confI=create_electron(cfgM(0)%conf1h(:,:,i1h),&
                  n_int_I,iint)
             conf(1:n_int_I,:)=confI

             ! Cycle if we have generated a ref conf
             if (h%key_exists(conf)) cycle
             
             ! Update the no. 1I confs
             if (modus == 0) cfgM(irrep1I)%n1I=cfgM(irrep1I)%n1I+1

             if (modus == 1) then

                ! Fill in the creation operator array
                counter(irrep1I)=counter(irrep1I)+1
                cfgM(irrep1I)%a1I(counter(irrep1I))=iint

                ! Update the no. confs generated by the current hole
                ! conf
                ngen(i1h,irrep1I)=ngen(i1h,irrep1I)+1

             endif
             
          enddo
          
       enddo
          
    enddo

!----------------------------------------------------------------------
! Fill in the offset arrays
!----------------------------------------------------------------------
    if (modus == 1) then
       do irrep=0,nirrep-1
          cfgM(irrep)%off1I(1)=1
          do i1h=2,n1h+1
             cfgM(irrep)%off1I(i1h)=cfgM(irrep)%off1I(i1h-1) &
                  +ngen(i1h-1,irrep)
          enddo
       enddo
    endif
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(icreate1)
    deallocate(ngen)
    deallocate(baseconf)
    deallocate(sop)
    deallocate(sop1h)
    deallocate(conf)
    deallocate(confI)
    deallocate(rhp)
    
    return
    
  end subroutine builder_1I
    
!######################################################################
  
end module confbuilder
