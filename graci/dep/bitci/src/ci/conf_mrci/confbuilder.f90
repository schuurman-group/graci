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
  subroutine generate_2I_1I1E_confs(E0max,cfgM,ddci,nroots)

    use constants
    use bitglobal
    use conftype
    use confinfo
    use mrciutils
    use iomod
    
    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max
    
    ! MRCI configurations for all irreps
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)

    ! DDCI2 configuration reduction
    logical, intent(in)        :: ddci

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)

    ! Doubly-occupied MOs in the reference space
    integer(is)                :: ndoccR(0:nirrep-1)
    integer(is)                :: doccR(nmo,0:nirrep-1)
    integer(is)                :: idoccR(nmo,0:nirrep-1)
    
    ! No. 1H1I, 2I and 1I1E configurations
    integer(is)                :: n1h1I,n2I,n1I1E
    
    ! Configuration scratch files
    integer(is)                :: nrec2I,nrec1I1E
    character(len=250)         :: file2I,file1I1E

    ! Everything else
    integer(is)                :: irrep,i,istart

!----------------------------------------------------------------------
! First index of an irrep with a non-zero number of roots
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       if (nroots(irrep) > 0) then
          istart=irrep
          exit
       endif
    enddo
    
!----------------------------------------------------------------------
! Get the list of MOs doubly-occupied in all reference configurations
! (this will be used for DDCI2 configuration reduction)
!----------------------------------------------------------------------
    if (ddci) then
       ! Get the lists of docc space MOs: one for each irrep
       do irrep=0,nirrep-1
          if (nroots(irrep) == 0) cycle
          call get_ref_docc_mos(cfgM(irrep),ndoccR(irrep),doccR(:,irrep))
       enddo
       
       ! Put the docc lists into a more useful form
       idoccR=0
       do irrep=0,nirrep-1
          if (nroots(irrep) == 0) cycle
          do i=1,ndoccR(irrep)
             idoccR(doccR(i,irrep),irrep)=1
          enddo
       enddo
    endif
    
!----------------------------------------------------------------------
! (1) Generate the 2I and 1I1E configurations for all irreps
!     *** including duplicates ***
!----------------------------------------------------------------------
    call builder_2I_1I1E(n1h1I,n2I,n1I1E,cfgM(istart),ddci,&
         idoccR,E0max,file2I,file1I1E,nrec2I,nrec1I1E,nroots)

!----------------------------------------------------------------------
! (2) Remove the duplicate 2I and 1I1E configurations
!----------------------------------------------------------------------
    call remove_duplicates_2I(n2I,cfgM,file2I,nrec2I,nroots,istart)
    call remove_duplicates_1I1E(n1I1E,cfgM,file1I1E,nrec1I1E,nroots,&
         istart)

!----------------------------------------------------------------------
! (3) Sort the 2I and 1I1E configurations by irrep
!----------------------------------------------------------------------
    call sort_2I_1I1E(cfgM,file2I,file1I1E,nrec2I,nrec1I1E,n2I,n1I1E,&
         nroots,istart)

    return
    
  end subroutine generate_2I_1I1E_confs

!######################################################################
! builder_2I_1I1E: peforms all the heavy lifting involved in the
!                  generation of the 2I and 1I1E configurations
!                  across all irreps
!######################################################################
  subroutine builder_2I_1I1E(n1h1I,n2I,n1I1E,cfgM,ddci,idoccR,&
       E0max,file2I,file1I1E,nrec2I,nrec1I1E,nroots)

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

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)           :: E0max

    ! DDCI2 configuration reduction
    logical, intent(in)            :: ddci
    integer(is), intent(in)        :: idoccR(nmo,0:nirrep-1)
    
    ! Configuration scratch files
    integer(is), intent(out)       :: nrec2I,nrec1I1E
    character(len=250), intent(out):: file2I,file1I1E

    ! Number of roots per irrep
    integer(is), intent(in)        :: nroots(0:nirrep-1)

    ! Configuration bit string buffers
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
    
    ! Allowable internal creation operator indices
    integer(is)                    :: ncreate1,ncreate2
    integer(is), allocatable       :: icreate1(:),icreate2(:)

    ! Work conf bit strings
    integer(ib), allocatable       :: conf(:,:),confI(:,:),sop(:,:),&
                                      conf_int(:,:)

    ! Sums of KS orbital energies multipled by occupation
    ! differences relative to the base conf
    real(dp)                       :: esumR,esum2H,esum1H1I,esum1I1E,&
                                      esum2I

    ! Number of annihilation/creation operators linking the various
    ! confs to the base conf
    integer(is)                    :: nacR,nac2H,nac1H1I,nac1I1E,nac2I

    ! Numbers of open shells
    integer(is)                    :: noR,no2H,no1H1I,no1I1E,no2I
    
    ! Everything else
    integer(is)                    :: i,j,k,n,np,i1,i2,i3
    integer(is)                    :: i2h,iext,iint1,iint2
    integer(is)                    :: n_int_I,nmoI
    integer(is)                    :: ia2h,ja2h,itmp,jtmp,nexci,nmatch
    integer(is)                    :: ic,counter

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

       ! Number of open shells in the reference configuration
       noR=sop_nopen(cfgM%sopR(:,:,n),n_int_I)
       
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
          k=(ia2h-1)/n_bits+1
          i=ia2h-(k-1)*n_bits-1
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
          k=(ja2h-1)/n_bits+1
          i=ja2h-(k-1)*n_bits-1
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

          ! Number of open shells in the 2-hole configuration
          !
          ! 1st annihilation operator
          !
          k=(ia2h-1)/n_bits+1
          i=ia2h-(k-1)*n_bits-1
          if (btest(cfgM%confR(k,2,n),i)) then
             ! Doubly-occupied MO in the ref conf: nopen -> nopen+1
             no2H=noR+1
          else
             ! Singly-occupied MO in the ref conf: nopen -> nopen-1
             no2H=noR-1
          endif
          !
          ! 2nd annihilation operator
          !
          k=(ja2h-1)/n_bits+1
          i=ja2h-(k-1)*n_bits-1
          if (btest(cfgM%confR(k,2,n),i)) then
             ! Doubly-occupied MO in the ref conf:
             if (ia2h == ja2h) then
                ! Equal annihilation operator indices: nopen -> nopen-1
                no2h=no2H-1
             else
                ! Unequal annihilation operator indices: nopen -> nopen+1
                no2h=no2H+1
             endif
          else
             ! Singly-occupied MO in the ref conf: nopen -> nopen-1
             no2H=no2H-1
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
             k=(iint1-1)/n_bits+1

             ! Postion of the external MO within the kth block
             i=iint1-(k-1)*n_bits-1
             
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

             ! Number of open shells in the 1H1I configuration
             if (btest(cfgM%conf2h(k,1,i2h),i)) then
                ! Singly-occupied MO in the 2-hole conf:
                ! nopen -> nopen-1
                no1H1I=no2H-1
             else
                ! Unoccupied MO in the 2-hole conf:
                ! nopen -> nopen+1
                no1H1I=no2H+1
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
                k=(iext-1)/n_bits+1
                
                ! Postion of the external MO within the kth block
                i=iext-(k-1)*n_bits-1

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

                ! Number of open shells in the 1I1E configuration
                no1I1E=no1H1I+1

                ! Skip this configuration if the number of open
                ! shells is above the maximum value
                if (no1I1E > nomax) cycle

                ! Skip this configuration if the number of open
                ! shells doesn't support the spin multiplicity
                if (no1I1E < imult-1) cycle
                
                ! Irrep generated by the 1I1E configuration
                ib1=irrep1H1I
                ib2=ieor(ib1,mosym(cfgM%m2c(iext)))
                irrep1I1E=ib2

                ! Cycle if there are no roots for this irrep
                if (nroots(irrep1I1E) == 0) cycle
                
                ! DDCI2 configuration reduction: skip this
                ! configuration if both annihilation operators
                ! correspond to the docc space
                if (ddci .and. idoccR(ia2h,irrep1I1E) == 1 &
                     .and. idoccR(ja2h,irrep1I1E) == 1) cycle
                
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
                k=(iint2-1)/n_bits+1
                
                ! Postion of the internal MO within the kth block
                i=iint2-(k-1)*n_bits-1

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

                ! Number of open shells in the 2I configuration
                if (btest(confI(k,1),i)) then
                   ! Singly-occupied MO in the 1H1I conf:
                   ! nopen -> nopen-1
                   no2I=no1H1I-1
                else
                   ! Unoccupied MO in the 1H1I conf:
                   ! nopen -> nopen+1
                   no2I=no1H1I+1
                endif

                ! Skip this configuration if the number of open shells
                ! is above the maximum value
                if (no2I > nomax) cycle

                ! Skip this configuration if the number of open
                ! shells doesn't support the spin multiplicity
                if (no2I < imult-1) cycle
                
                ! Irrep generated by the 2I configuration
                ib1=irrep1H1I
                ib2=ieor(ib1,mosym(cfgM%m2c(iint2)))
                irrep2I=ib2

                ! Cycle if there are no roots for this irrep
                if (nroots(irrep2I) == 0) cycle
                
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
    deallocate(ibuf2I)
    deallocate(ibuf1I1E)

    return
    
  end subroutine builder_2I_1I1E
    
!######################################################################
! remove_duplicates_2I: removes the duplicate 2I configurations using
!                       a hash table-based approach
!######################################################################
  subroutine remove_duplicates_2I(n2I,cfgM,file2I,nrec2I,nroots,istart)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use dethash
    use iomod    

    implicit none

    ! Number of 2I configurations
    integer(is), intent(inout)     :: n2I

    ! MRCI configurations
    type(mrcfg), intent(inout)     :: cfgM(0:nirrep-1)

    ! Original 2I configuration scratch file
    integer(is), intent(inout)     :: nrec2I
    character(len=250), intent(in) :: file2I

    ! Number of roots per irrep
    integer(is), intent(in)        :: nroots(0:nirrep-1)    

    ! First index of an irrep with a non-zero number of roots
    integer(is), intent(in)        :: istart
    
    ! Scratch file I/O
    integer(is), allocatable       :: ibuf(:,:)
    integer(is)                    :: iscratch

    ! New 2I configuration scratch file
    integer(is)                    :: nrec2I_new
    integer(is), allocatable       :: ibuf_new(:,:)
    integer(is)                    :: iscratch_new
    character(len=250)             :: file2I_new
    
    ! Hash table
    type(dhtbl)                    :: h
    integer(is)                    :: initial_size
    integer(is)                    :: nold
    integer(ib), allocatable       :: key(:,:)

    ! Configuration bit strings
    integer(ib), allocatable       :: conf(:,:),confI(:,:)
    
    ! Everything else
    integer(is)                    :: irrep,i,n,ioff,imo,ic,jc,i2h
    integer(is)                    :: irec,nbuf,nbuf_new
    integer(is)                    :: n_int_I,n2I_new
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgM(istart)%n_int_I

    allocate(ibuf(4,bufsize))
    ibuf=0

    allocate(ibuf_new(4,bufsize))
    ibuf_new=0
    
    allocate(key(n_int,2))
    key=0_ib
    
    allocate(conf(n_int,2))
    conf=0_ib

    allocate(confI(n_int_I,2))
    confI=0_ib
    
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------    
    ! Rough stab at an appropriate hash table size
    initial_size=(nrec2I*bufsize)/3
    do irrep=0,nirrep-1
       if (nroots(irrep) /= 0) then
          initial_size=initial_size+cfgM(irrep)%n0h
          initial_size=initial_size+cfgM(irrep)%n1I
       endif
    enddo

    ! Initialise the hash table
    call h%initialise_table(initial_size)

!----------------------------------------------------------------------
! Insert the reference and 1I configurations into the hash table
!----------------------------------------------------------------------
    !
    ! Ref confs
    !
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(irrep) == 0) cycle

       ! Loop over ref confs for this irrsp
       do i=1,cfgM(irrep)%n0h
          key=0_ib
          key(1:n_int_I,:)=cfgM(irrep)%conf0h(:,:,i)
          call h%insert_key(key)
       enddo

    enddo
       
    !
    ! 1I confs
    !
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(irrep) == 0) cycle
       
       ! Loop over 1-hole confs
       do n=1,cfgM(irrep)%n1h

          ! Loop over the 1I confs generated by this 1-hole conf
          do ioff=cfgM(irrep)%off1I(n),cfgM(irrep)%off1I(n+1)-1
             
             ! Construct the configuration
             imo=cfgM(irrep)%a1I(ioff)
             conf=0_ib
             conf(1:n_int_I,:)=cfgM(irrep)%conf1h(:,:,n)
             key=create_electron(conf,n_int,imo)
             
             ! Insert the configuration into the hash table
             call h%insert_key(key)
             
          enddo
          
       enddo

    enddo
       
    !
    ! Number of keys stored
    !
    nold=h%n_keys_stored
    
!----------------------------------------------------------------------
! Open the 2I configuration scratch files
!----------------------------------------------------------------------
    ! Old scratch file including duplicates
    call freeunit(iscratch)
    open(iscratch,file=file2I,form='unformatted',status='old')

    ! New scratch file with duplicates removed
    call freeunit(iscratch_new)
    call scratch_name('2I_new',file2I_new)
    open(iscratch_new,file=file2I_new,form='unformatted',&
         status='unknown')
    
!----------------------------------------------------------------------    
! Find the unique 2I configurations
!----------------------------------------------------------------------
    ! Number of unique 2I confs
    n2I_new=0

    ! Initialise the buffer for the unique confs
    nbuf_new=0
    nrec2I_new=0
    
    ! Loop over records
    do irec=1,nrec2I

       ! Read in the next batch of confs
       read(iscratch) nbuf,ibuf

       ! Loop over the confs in the current batch
       do i=1,nbuf

          ! Irrep
          irrep=ibuf(1,i)

          ! Cycle if there are no roots for this irrep
          if (nroots(irrep) == 0) cycle
          
          ! 2-hole conf index
          i2h=ibuf(2,i)
          
          ! Creation operator indices
          ic=ibuf(3,i)
          jc=ibuf(4,i)

          ! Full 2I conf
          confI=create_electron(cfgM(irrep)%conf2h(:,:,i2h),n_int_I,ic)
          conf=0_ib
          conf(1:n_int_I,:)=confI
          key=create_electron(conf,n_int,jc)

          ! Attempt an insertion into the hash table
          call h%insert_key(key)

          ! Save the conf if it is unique
          if (h%n_keys_stored > nold) then

             nold=h%n_keys_stored

             n2I_new=n2I_new+1

             nbuf_new=nbuf_new+1

             ibuf_new(:,nbuf_new)=ibuf(:,i)
             
          endif

          ! Dump the buffer to disk
          if (nbuf_new == bufsize) then
             nrec2I_new=nrec2I_new+1
             write(iscratch_new) bufsize,ibuf_new
             nbuf_new=0
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Write the remaining buffer to disk
!----------------------------------------------------------------------
    if (nbuf_new /= 0) then
       nrec2I_new=nrec2I_new+1
       write(iscratch_new) nbuf_new,ibuf_new
    endif

!----------------------------------------------------------------------
! New number of 2I confs
!----------------------------------------------------------------------
    n2I=n2I_new
    
!----------------------------------------------------------------------
! New number of records
!----------------------------------------------------------------------
    nrec2I=nrec2I_new
    
!----------------------------------------------------------------------
! Close the 2I configuration scratch files
!----------------------------------------------------------------------
    close(iscratch)
    close(iscratch_new)

!----------------------------------------------------------------------
! Rename the new 2I scratch file
!----------------------------------------------------------------------
    call system('mv '//trim(file2I_new)//' '//trim(file2I))
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ibuf)
    deallocate(key)
    deallocate(conf)
    deallocate(confI)
    
    return
    
  end subroutine remove_duplicates_2I

!######################################################################
! remove_duplicates_1I1E: removes the duplicate 1I1E configurations
!                         using a hash table-based approach
!######################################################################
  subroutine remove_duplicates_1I1E(n1I1E,cfgM,file1I1E,nrec1I1E,&
       nroots,istart)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use dethash
    use iomod    

    implicit none

    ! Number of 1I1E configurations
    integer(is), intent(inout)     :: n1I1E

    ! MRCI configurations
    type(mrcfg), intent(inout)     :: cfgM(0:nirrep-1)

    ! Original 1I1E configuration scratch file
    integer(is), intent(inout)     :: nrec1I1E
    character(len=250), intent(in) :: file1I1E

    ! Number of roots per irrep
    integer(is), intent(in)        :: nroots(0:nirrep-1)    

    ! First index of an irrep with a non-zero number of roots
    integer(is), intent(in)        :: istart

    ! Scratch file I/O
    integer(is), allocatable       :: ibuf(:,:)
    integer(is)                    :: iscratch

    ! New 1I1E configuration scratch file
    integer(is)                    :: nrec1I1E_new
    integer(is), allocatable       :: ibuf_new(:,:)
    integer(is)                    :: iscratch_new
    character(len=250)             :: file1I1E_new
    
    ! Hash table
    type(dhtbl)                    :: h
    integer(is)                    :: initial_size
    integer(is)                    :: nold
    integer(ib), allocatable       :: key(:,:)

    ! Configuration bit strings
    integer(ib), allocatable       :: conf(:,:),conf1(:,:)
    
    ! Everything else
    integer(is)                    :: irrep,i,n,ioff,imo,ic,jc,i2h
    integer(is)                    :: irec,nbuf,nbuf_new
    integer(is)                    :: n_int_I,n1I1E_new

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n_int_I=cfgM(istart)%n_int_I

    allocate(ibuf(4,bufsize))
    ibuf=0

    allocate(ibuf_new(4,bufsize))
    ibuf_new=0
    
    allocate(key(n_int,2))
    key=0_ib
    
    allocate(conf(n_int,2))
    conf=0_ib

    allocate(conf1(n_int,2))
    conf1=0_ib
    
!----------------------------------------------------------------------
! Initialise the hash table
!----------------------------------------------------------------------    
    ! Rough stab at an appropriate hash table size
    initial_size=(nrec1I1E*bufsize)/3
    do irrep=0,nirrep-1
       if (nroots(irrep) /= 0) &
            initial_size=initial_size+cfgM(irrep)%n1E
    enddo
    
    ! Initialise the hash table
    call h%initialise_table(initial_size)

!----------------------------------------------------------------------
! Insert the 1E configurations into the hash table
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(irrep) == 0) cycle
       
       ! Loop over 1-hole confs
       do n=1,cfgM(irrep)%n1h
          
          ! Loop over the 1E confs generated by this 1-hole conf
          do ioff=cfgM(irrep)%off1E(n),cfgM(irrep)%off1E(n+1)-1
             
             ! Construct the configuration
             imo=cfgM(irrep)%a1E(ioff)
             conf=0_ib
             conf(1:n_int_I,:)=cfgM(irrep)%conf1h(:,:,n)
             key=create_electron(conf,n_int,imo)
             
             ! Insert the configuration into the hash table
             call h%insert_key(key)
             
          enddo
          
       enddo

    enddo
       
    ! Number of keys stored
    nold=h%n_keys_stored

!----------------------------------------------------------------------
! Open the 1I1E configuration scratch files
!----------------------------------------------------------------------
    ! Old scratch file including duplicates
    call freeunit(iscratch)
    open(iscratch,file=file1I1E,form='unformatted',status='old')

    ! New scratch file with duplicates removed
    call freeunit(iscratch_new)
    call scratch_name('1I1E_new',file1I1E_new)
    open(iscratch_new,file=file1I1E_new,form='unformatted',&
         status='unknown')

!----------------------------------------------------------------------    
! Find the unique 1I1E configurations
!----------------------------------------------------------------------
    ! Number of unique 1I1E confs
    n1I1E_new=0

    ! Initialise the buffer for the unique confs
    nbuf_new=0
    nrec1I1E_new=0
    
    ! Loop over records
    do irec=1,nrec1I1E

       ! Read in the next batch of confs
       read(iscratch) nbuf,ibuf

       ! Loop over the confs in the current batch
       do i=1,nbuf

          ! Irrep
          irrep=ibuf(1,i)

          ! Cycle if there are no roots for this irrep
          if (nroots(irrep) == 0) cycle
          
          ! 2-hole conf index
          i2h=ibuf(2,i)
          
          ! Creation operator indices
          ic=ibuf(3,i)
          jc=ibuf(4,i)

          ! Full 1I1E conf
          conf=0_ib
          conf(1:n_int_I,:)=cfgM(irrep)%conf2h(:,:,i2h)
          conf1=create_electron(conf,n_int,ic)
          key=create_electron(conf1,n_int,jc)

          ! Attempt an insertion into the hash table
          call h%insert_key(key)

          ! Save the conf if it is unique
          if (h%n_keys_stored > nold) then

             nold=h%n_keys_stored

             n1I1E_new=n1I1E_new+1

             nbuf_new=nbuf_new+1

             ibuf_new(:,nbuf_new)=ibuf(:,i)
             
          endif

          ! Dump the buffer to disk
          if (nbuf_new == bufsize) then
             nrec1I1E_new=nrec1I1E_new+1
             write(iscratch_new) bufsize,ibuf_new
             nbuf_new=0
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Write the remaining buffer to disk
!----------------------------------------------------------------------
    if (nbuf_new /= 0) then
       nrec1I1E_new=nrec1I1E_new+1
       write(iscratch_new) nbuf_new,ibuf_new
    endif

!----------------------------------------------------------------------
! New number of 1I1E confs
!----------------------------------------------------------------------
    n1I1E=n1I1E_new
    
!----------------------------------------------------------------------
! New number of records
!----------------------------------------------------------------------
    nrec1I1E=nrec1I1E_new
    
!----------------------------------------------------------------------
! Close the 1I1E configuration scratch files
!----------------------------------------------------------------------
    close(iscratch)
    close(iscratch_new)

!----------------------------------------------------------------------
! Rename the new 1I1E scratch file
!----------------------------------------------------------------------
    call system('mv '//trim(file1I1E_new)//' '//trim(file1I1E))

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(ibuf)
    deallocate(key)
    deallocate(conf)
    deallocate(conf1)
    
    return

  end subroutine remove_duplicates_1I1E
    
!######################################################################
! sort_2I_1I1E: sorts and saves the 2I and 1I1E configurations by
!               irrep
!######################################################################
  subroutine sort_2I_1I1E(cfgM,file2I,file1I1E,nrec2I,nrec1I1E,&
       n2I_tot,n1I1E_tot,nroots,istart)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    use iomod
    
    implicit none

    ! MRCI configurations for all irreps
    type(mrcfg), intent(inout)     :: cfgM(0:nirrep-1)

    ! Total number of 2I and 1I1E configurations
    integer(is), intent(in)        :: n2I_tot,n1I1E_tot
    
    ! Configuration scratch files
    integer(is), intent(in)        :: nrec2I,nrec1I1E
    character(len=250), intent(in) :: file2I,file1I1E

    ! Number of roots per irrep
    integer(is), intent(in)        :: nroots(0:nirrep-1)    

    ! First index of an irrep with a non-zero number of roots
    integer(is), intent(in)        :: istart
    
    ! Scratch file I/O
    integer(is), allocatable       :: ibuf(:,:)
    integer(is)                    :: iscratch2I,iscratch1I1E

    ! Working arrays
    integer(ib), allocatable       :: sop(:,:)

    ! No. confs generated by each 2-hole conf
    integer(is), allocatable       :: ngen(:,:)
    
    ! Everything else
    integer(is)                    :: nbuf,i,irec,irrep
    integer(is)                    :: counter(0:nirrep-1)
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(ibuf(4,bufsize))
    ibuf=0

    allocate(sop(n_int,2))
    sop=0_ib

    allocate(ngen(cfgM(istart)%n2h,0:nirrep-1))
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
             if (nroots(irrep) /= 0) &
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
             if (nroots(irrep) /= 0) &
                  cfgM(irrep)%n1I1E=cfgM(irrep)%n1I1E+1
          enddo
          
       enddo

    endif
       
!----------------------------------------------------------------------
! Allocate the offset and creation operator index arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(irrep) == 0) cycle
       
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

          ! Cycle if there are no roots for this irrep
          if (nroots(irrep) == 0) cycle
          
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
       if (nroots(irrep) == 0) cycle
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

          ! Cycle if there are no roots for this irrep
          if (nroots(irrep) == 0) cycle

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
       if (nroots(irrep) == 0) cycle
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
! generate_1E_confs: for all irreps, generates all the allowable
!                    configurations with one internal hole and one
!                    external electron
!######################################################################
  subroutine generate_1E_confs(E0max,cfgM,nroots)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)
    
    ! Everything else
    integer(is)                :: modus,irrep

!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! each symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_1E(modus,E0max,cfgM,nroots)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(irrep) == 0) cycle

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
    call builder_1E(modus,E0max,cfgM,nroots)

    return
    
  end subroutine generate_1E_confs

!######################################################################
! 1E_builder: performs all the heavy lifting involved in the
!             generation of the configurations with one internal hole
!             and one external electron
!######################################################################
  subroutine builder_1E(modus,E0max,cfgM,nroots)

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

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)
    
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
    integer(is)                :: k,i,istart
    integer(is)                :: counter(0:nirrep-1)
    integer(is)                :: nexci1H,noR,no1H

!----------------------------------------------------------------------
! First index of an irrep with a non-zero number of roots
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       if (nroots(irrep) > 0) then
          istart=irrep
          exit
       endif
    enddo
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n1h=cfgM(istart)%n1h
    nR=cfgM(istart)%nR
    n_int_I=cfgM(istart)%n_int_I
    nmoI=cfgM(istart)%nmoI

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
       esumR=esum(cfgM(istart)%confR(:,:,n),n_int_I,cfgM(istart)%m2c)

       ! Irrep generated by the reference configuration
       sop=0_ib
       sop(1:n_int_I,:)=cfgM(istart)%sopR(:,:,n)
       irrepR=sop_sym_mrci(sop,cfgM(istart)%m2c)

       ! Number of open shells in the reference configuration
       noR=sop_nopen(cfgM(istart)%sopR(:,:,n),n_int_I)
       
       ! Loop over the 1-hole configurations generated by the current
       ! reference configuration
       do i1h=cfgM(istart)%off1h(n),cfgM(istart)%off1h(n+1)-1

          ! 1-hole annihilation operator index
          ia1h=cfgM(istart)%a1h(i1h)

          ! Sum_p F_pp^(KS) Delta W_p for the 1-hole configuration
          esum1H=esumR-moen(cfgM(istart)%m2c(ia1h))

          ! If this is a DFT/MRCI calculation, then skip this
          ! 1-hole configuration cannot possibly generate 1E
          ! configurations satisfying the energy-based
          ! selection criterion
          if (ldftmrci .and. esum1H > E0max + desel) cycle

          ! Cycle if the excitation degree of the 1E configurations
          ! generated by this 1-hole configuration will be above
          ! threshold
          nexci1H=exc_degree_conf(cfgM(istart)%conf1h(:,:,i1h),&
               baseconf,n_int_I)
          if (nexci1H+1 > nexmax) cycle

          ! Number of open shells in the 1-hole configuration
          k=(ia1h-1)/n_bits+1
          i=ia1h-(k-1)*n_bits-1
          if (btest(cfgM(istart)%confR(k,2,n),i)) then
             ! Doubly-occupied MO in the ref conf: nopen -> nopen+1
             no1H=noR+1
          else
             ! Singly-occupied MO in the ref conf: nopen -> nopen-1
             no1H=noR-1
          endif

          ! Cycle if the no. open shells in the 1E configurations
          ! will be above the maximum value
          if (no1H+1 > nomax) cycle
          
          ! Irrep generated by the 1-hole configuration
          ib1=irrepR
          ib2=ieor(ib1,mosym(cfgM(istart)%m2c(ia1h)))
          irrep1H=ib2

          ! Loop over external MOs
          do iext=nmoI+1,lastvirt

             ! Sum_p F_pp^(KS) Delta W_p for the 1E configuration
             esum1E=esum1H+moen(cfgM(istart)%m2c(iext))

             ! If this is a DFT/MRCI calculation, then skip this
             ! configuration if it doesn't satisfy the energy-based
             ! selection criterion
             if (ldftmrci .and. esum1E > E0max + desel) cycle

             ! Irrep generated by the 1E configuration
             ib1=irrep1H
             ib2=ieor(ib1,mosym(cfgM(istart)%m2c(iext)))
             irrep1E=ib2

             ! Cycle if there are no roots for this irrep
             if (nroots(irrep1E) == 0) cycle
             
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
          if (nroots(irrep) == 0) cycle
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
  subroutine generate_2E_confs(E0max,cfgM,ddci,nroots)

    use constants
    use bitglobal
    use conftype
    use confinfo
    
    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)
    
    ! DDCI2 configuration reduction
    logical, intent(in)        :: ddci

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)

    ! Doubly-occupied MOs in the reference space
    integer(is)                :: ndoccR(0:nirrep-1)
    integer(is)                :: doccR(nmo,0:nirrep-1)
    integer(is)                :: idoccR(nmo,0:nirrep-1)
    
    ! Everything else
    integer(is)                :: modus,irrep,i

!----------------------------------------------------------------------
! Get the list of MOs doubly-occupied in all reference configurations
! (this will be used for DDCI2 configuration reduction)
!----------------------------------------------------------------------
    if (ddci) then
       ! Get the lists of docc space MOs: one for each irrep
       do irrep=0,nirrep-1
          if (nroots(irrep) == 0) cycle
          call get_ref_docc_mos(cfgM(irrep),ndoccR(irrep),doccR(:,irrep))
       enddo
       
       ! Put the docc lists into a more useful form
       idoccR=0
       do irrep=0,nirrep-1
          if (nroots(irrep) == 0) cycle
          do i=1,ndoccR(irrep)
             idoccR(doccR(i,irrep),irrep)=1
          enddo
       enddo
    endif
       
!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! each symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_2E(modus,E0max,cfgM,ddci,idoccR,nroots)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(irrep) == 0) cycle
       
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
    call builder_2E(modus,E0max,cfgM,ddci,idoccR,nroots)
    
    return
    
  end subroutine generate_2E_confs

!######################################################################
! 2E_builder: performs all the heavy lifting involved in the
!             generation of the configurations with two internal holes
!             and two external electrons
!######################################################################
  subroutine builder_2E(modus,E0max,cfgM,ddci,idoccR,nroots)

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

    ! DDCI2 configuration reduction
    logical, intent(in)        :: ddci
    integer(is), intent(in)    :: idoccR(nmo,0:nirrep-1)

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)
    
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
    integer(is)                :: istart
    integer(is)                :: counter(0:nirrep-1)
    integer(is)                :: nexci2H
    integer(is)                :: noR,no2H,no2E

!----------------------------------------------------------------------
! First index of an irrep with a non-zero number of roots
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       if (nroots(irrep) > 0) then
          istart=irrep
          exit
       endif
    enddo
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n2h=cfgM(istart)%n2h
    nR=cfgM(istart)%nR
    n_int_I=cfgM(istart)%n_int_I
    nmoI=cfgM(istart)%nmoI

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
       esumR=esum(cfgM(istart)%confR(:,:,n),n_int_I,cfgM(istart)%m2c)

       ! Irrep generated by the reference configuration
       sop=0_ib
       sop(1:n_int_I,:)=cfgM(istart)%sopR(:,:,n)
       irrepR=sop_sym_mrci(sop,cfgM(istart)%m2c)

       ! Number of open shells in the reference configuration
       noR=sop_nopen(cfgM(istart)%sopR(:,:,n),n_int_I)
       
       ! Loop over the 2-hole configurations generated by the current
       ! reference configuration
       do i2h=cfgM(istart)%off2h(n),cfgM(istart)%off2h(n+1)-1

          ! 2-hole annihilation operator indices
          ia2h=cfgM(istart)%a2h(1,i2h)
          ja2h=cfgM(istart)%a2h(2,i2h)
          
          ! Sum_p F_pp^(KS) Delta W_p for the 2-hole configuration
          esum2H=esumR-moen(cfgM(istart)%m2c(ia2h))&
               -moen(cfgM(istart)%m2c(ja2h))

          ! If this is a DFT/MRCI calculation, then skip this
          ! 2-hole configuration cannot possibly generate 2E
          ! configurations satisfying the energy-based
          ! selection criterion
          if (ldftmrci .and. esum2H > E0max + desel) cycle

          ! Cycle if the excitation degree of the 2E configurations
          ! generated by this 2-hole configuration will be above
          ! threshold
          nexci2H=exc_degree_conf(cfgM(istart)%conf2h(:,:,i2h),&
               baseconf,n_int_I)
          if (nexci2H+2 > nexmax) cycle

          ! Number of open shells in the 2-hole configuration
          !
          ! 1st annihilation operator
          !
          k=(ia2h-1)/n_bits+1
          i=ia2h-(k-1)*n_bits-1
          if (btest(cfgM(istart)%confR(k,2,n),i)) then
             ! Doubly-occupied MO in the ref conf:
             ! nopen -> nopen+1
             no2H=noR+1
          else
             ! Singly-occupied MO in the ref conf:
             ! nopen -> nopen-1
             no2H=noR-1
          endif
          !
          ! 2nd annihilation operator
          !
          k=(ja2h-1)/n_bits+1
          i=ja2h-(k-1)*n_bits-1
          if (btest(cfgM(istart)%confR(k,2,n),i)) then
             ! Doubly-occupied MO in the ref conf:
             if (ia2h == ja2h) then
                ! Equal annihilation operator indices:
                ! nopen -> nopen-1
                no2H=no2H-1
             else
                ! Unequal annihilation operator indices:
                ! nopen -> nopen+1
                no2H=no2H+1
             endif
          else
             ! Singly-occupied MO in the ref conf:
             ! nopen -> nopen-1
             no2H=no2H-1
          endif

          ! Cycle if the number open shells of the 2E configurations
          ! generated by this 2-hole configuration cannot be below
          ! the maximum value
          if (no2H > nomax) cycle
          
          ! Irrep generated by the 2-hole configuration
          ib1=irrepR
          ib2=ieor(ib1,mosym(cfgM(istart)%m2c(ia2h)))
          ib1=ieor(ib2,mosym(cfgM(istart)%m2c(ja2h)))
          irrep2H=ib1
          
          ! Loop over the first creation operator
          do iext1=nmoI+1,lastvirt

             ! Sum_p F_pp^(KS) Delta W_p for the 1H1E configuration
             esum1H1E=esum2H+moen(cfgM(istart)%m2c(iext1))

             ! If this is a DFT/MRCI calculation, then skip this
             ! 1H1E configuration cannot possibly generate 2E
             ! configurations satisfying the energy-based
             ! selection criterion
             if (ldftmrci .and. esum1H1E > E0max + desel) cycle

             ! Irrep generated by the 1H1E configuration
             ib1=irrep2H
             ib2=ieor(ib1,mosym(cfgM(istart)%m2c(iext1)))
             irrep1H1E=ib2

             ! Loop over the second creation operator
             do iext2=iext1,lastvirt

                ! Sum_p F_pp^(KS) Delta W_p for the 2E configuration
                esum2E=esum1H1E+moen(cfgM(istart)%m2c(iext2))
                
                ! If this is a DFT/MRCI calculation, then skip this
                ! 2E configuration if it does not satisfy the
                ! energy-based selection criterion
                if (ldftmrci .and. esum2E > E0max + desel) cycle

                ! Number of open shells
                if (iext1 == iext2) then
                   no2E=no2H
                else
                   no2E=no2H+2
                endif
                
                ! Skip this 2E configuration if the no. open shells
                ! is too high or cannot support the spin multiplicity
                if (no2E > nomax) cycle
                if (no2E < imult-1) cycle
                
                ! Irrep generated by the 2E configuration
                ib1=irrep1H1E
                ib2=ieor(ib1,mosym(cfgM(istart)%m2c(iext2)))
                irrep2E=ib2

                ! Cycle if there are no roots for this irrep
                if (nroots(irrep2E) == 0) cycle
                
                ! DDCI2 configuration reduction: skip this 2E conf if
                ! one of the annihilation operators correspond to the
                ! docc space
                if (ddci) then
                   if (idoccR(ia2h,irrep2E) == 1 &
                        .or. idoccR(ja2h,irrep2E) == 1) cycle
                endif
                   
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
          if (nroots(irrep) == 0) cycle
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
  subroutine generate_1I_confs(E0max,cfgM,nroots)

    use constants
    use bitglobal
    use conftype
    
    implicit none

    ! Energy of the highest-lying reference space state of interest
    real(dp), intent(in)       :: E0max

    ! MRCI configurations
    type(mrcfg), intent(inout) :: cfgM(0:nirrep-1)

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)
    
    ! Everything else
    integer(is)                :: modus,irrep

!----------------------------------------------------------------------
! First, determine the number of allowable configurations of the
! each symmetry
!----------------------------------------------------------------------
    modus=0
    call builder_1I(modus,E0max,cfgM,nroots)

!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    ! Loop over irreps
    do irrep=0,nirrep-1

       ! Cycle if there are no roots for this irrep
       if (nroots(irrep) == 0) cycle
       
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
    call builder_1I(modus,E0max,cfgM,nroots)
    
    return
    
  end subroutine generate_1I_confs

!######################################################################
! 1I_builder: performs all the heavy lifting involved in the
!             generation of the configurations with one internal hole
!             and one internal electron
!######################################################################
  subroutine builder_1I(modus,E0max,cfgM,nroots)

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

    ! Number of roots per irrep
    integer(is), intent(in)    :: nroots(0:nirrep-1)
    
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
    integer(is)                :: ia1h,iint,k,i,j,istart
    integer(is)                :: counter(0:nirrep-1)
    integer(is)                :: nacR,nac1H,nac1I
    integer(is)                :: noR,no1H,no1I
    integer(is)                :: nexci,iaR,jaR,icR,jcR

!----------------------------------------------------------------------
! First index of an irrep with a non-zero number of roots
!----------------------------------------------------------------------
    do irrep=0,nirrep-1
       if (nroots(irrep) > 0) then
          istart=irrep
          exit
       endif
    enddo
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    n1h=cfgM(istart)%n1h
    nR=cfgM(istart)%nR
    n_int_I=cfgM(istart)%n_int_I
    nmoI=cfgM(istart)%nmoI

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
    do n=1,cfgM(istart)%nR
       conf=0_ib
       conf(1:n_int_I,:)=cfgM(istart)%confR(:,:,n)
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
       confI=cfgM(istart)%conf1h(:,:,i1h)
       sop1h(:,:,i1h)=conf_to_sop(confI,n_int_I)

    enddo
    
!----------------------------------------------------------------------
! Generate the 1I configurations
!----------------------------------------------------------------------
    ! Loop over reference configurations
    do n=1,nR

       ! Sum_p F_pp^(KS) Delta W_p for the reference configuration
       esumR=esum(cfgM(istart)%confR(:,:,n),n_int_I,cfgM(istart)%m2c)
       
       ! Irrep generated by the reference configuration
       sop=0_ib
       sop(1:n_int_I,:)=cfgM(istart)%sopR(:,:,n)
       irrepR=sop_sym_mrci(sop,cfgM(istart)%m2c)

       ! Number of creation and annihilation operators linking the
       ! reference and base configurations
       nacR=n_create_annihilate(cfgM(istart)%confR(:,:,n),baseconf,&
            n_int_I)

       ! Number of open shells in the reference configuration
       noR=sop_nopen(cfgM(istart)%sopR(:,:,n),n_int_I)
       
       ! Initialise the array of inter-ref-conf annihilation/creation
       ! operator indices
       rhp=0

       ! Loop over preceding reference configurations
       do np=1,n-1

          ! Determine the excitation degree between the two reference
          ! configurations
          nexci=exc_degree_conf(cfgM(istart)%confR(:,:,n),&
               cfgM(istart)%confR(:,:,np),n_int_I)

          ! Cycle if the excitation degree is greater than 2
          if (nexci > 2) cycle

          ! Determine the indices of the creation/annihilation
          ! operators linking the two ref confs
          call get_exci_indices(cfgM(istart)%confR(:,:,n),&
               cfgM(istart)%confR(:,:,np),n_int_I,hlist(1:nexci),&
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
       do i1h=cfgM(istart)%off1h(n),cfgM(istart)%off1h(n+1)-1

          ! 1-hole annihilation operator index
          ia1h=cfgM(istart)%a1h(i1h)

          ! Sum_p F_pp^(KS) Delta W_p for the 1-hole configuration
          esum1H=esumR-moen(cfgM(istart)%m2c(ia1h))

          ! Irrep generated by the 1-hole configuration
          ib1=irrepR
          ib2=ieor(ib1,mosym(cfgM(istart)%m2c(ia1h)))
          irrep1H=ib2

          ! Number of creation and annihilation operators linking the
          ! 1-hole and base configurations
          nac1H=nacR
          k=(ia1h-1)/n_bits+1
          i=ia1h-(k-1)*n_bits-1
          if (iocc0(ia1h) == 0) then
             ! Unoccupied MO in the base conf
             nac1H=nac1H-1
          else if (iocc0(ia1h) == 1) then
             ! Singly-occupied MO in the base conf
             if (btest(cfgM(istart)%confR(k,2,n),i)) then
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

          ! Number of open shells in the 1-hole configuration
          k=(ia1h-1)/n_bits+1
          i=ia1h-(k-1)*n_bits-1
          if (btest(cfgM(istart)%confR(k,2,n),i)) then
             ! Doubly-occupied MO in the ref conf: nopen -> nopen+1
             no1H=noR+1
          else
             ! Singly-occupied MO in the ref conf: nopen -> nopen-1
             no1H=noR-1
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
             
             ! Cycle if this is a CVS-MRCI calculation and we are
             ! creating an electron in a flagged core MO
             if (lcvs .and. icvs(cfgM(istart)%m2c(iint)) == 1) cycle

             ! Sum_p F_pp^(KS) Delta W_p for the 1I configuration
             esum1I=esum1H+moen(cfgM(istart)%m2c(iint))
             
             ! If this is a DFT/MRCI calculation, then skip this
             ! 1I configuration if it does not satisfy the
             ! energy-basedselection criterion
             if (ldftmrci .and. esum1I > E0max + desel) cycle
             
             ! Irrep generated by the 1I configuration
             ib1=irrep1H
             ib2=ieor(ib1,mosym(cfgM(istart)%m2c(iint)))
             irrep1I=ib2

             ! Cycle if there are no roots for this irrep
             if (nroots(irrep1I) == 0) cycle
             
             ! Block index
             k=(iint-1)/n_bits+1

             ! Postion of the external MO within the kth block
             i=iint-(k-1)*n_bits-1

             ! Number of creation and annihilation operators linking the
             ! 1I and base configurations
             nac1I=nac1H
             if (iocc0(iint) == 0) then
                ! Unoccupied MO in the base conf
                nac1I=nac1I+1
             else if (iocc0(iint) == 1) then
                ! Singly-occupied MO in the base conf
                if (btest(cfgM(istart)%conf1h(k,1,i1h),i)) then
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

             ! Number of open shells in the 1I configuration
             if (btest(cfgM(istart)%conf1h(k,1,i1h),i)) then
                ! Singly-occupied MO in the 1-hole conf:
                ! nopen -> nopen-1
                no1I=no1H-1
             else
                ! Unoccupied MO in the 1-hole conf:
                ! nopen -> nopen+1
                no1I=no1H+1
             endif

             ! Cycle if the number of open shells is too high
             if (no1I > nomax) cycle

             ! Skip this configuration if the number of open
             ! shells doesn't support the spin multiplicity
             if (no1I < imult-1) cycle
             
             ! Full configuration
             conf=0_ib
             confI=create_electron(cfgM(istart)%conf1h(:,:,i1h),&
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
          if (nroots(irrep) == 0) cycle
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
  
end module confbuilder
