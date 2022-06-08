!**********************************************************************
! Routines for the calculation of the determinant expansions of MRCI
! wave functions
!**********************************************************************
module csf2det

  implicit none

contains

!######################################################################
! get_detdim: determines the number of determinants generated by a
!             given configuration basis
!######################################################################
  subroutine get_detdim(cfg,detdim)

    use constants
    use bitglobal
    use conftype
    use mrciutils

    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Size of the determinant basis
    integer(is), intent(out) :: detdim

    ! Everyting else
    integer(is)              :: nopen,iconf

    !
    ! Initialisation
    !
    detdim=0

    !
    ! Reference space confs
    !
    do iconf=1,cfg%n0h
       nopen=sop_nopen(cfg%sop0h(:,:,iconf),cfg%n_int_I)
       detdim=detdim+ndets(nopen)
    enddo
    
    !
    ! 1I confs
    !
    do iconf=1,cfg%n1I
       nopen=sop_nopen(cfg%sop1I(:,:,iconf),n_int)
       detdim=detdim+ndets(nopen)       
    enddo

    !
    ! 2I confs
    !
    do iconf=1,cfg%n2I
       nopen=sop_nopen(cfg%sop2I(:,:,iconf),n_int)
       detdim=detdim+ndets(nopen)       
    enddo

    !
    ! 1E confs
    !
    do iconf=1,cfg%n1E
       nopen=sop_nopen(cfg%sop1E(:,:,iconf),n_int)
       detdim=detdim+ndets(nopen)       
    enddo

    !
    ! 2E confs
    !
    do iconf=1,cfg%n2E
       nopen=sop_nopen(cfg%sop2E(:,:,iconf),n_int)
       detdim=detdim+ndets(nopen)       
    enddo

    !
    ! 1I1E confs
    !
    do iconf=1,cfg%n1I1E
       nopen=sop_nopen(cfg%sop1I1E(:,:,iconf),n_int)
       detdim=detdim+ndets(nopen)       
    enddo

    return
    
  end subroutine get_detdim
    
!######################################################################
! eigenvectors_detbas: computes a set of eigenvectors in the
!                      determinant basis
!######################################################################
  subroutine eigenvectors_detbas(cfg,nroots,csfdim,detdim,vec_csf,&
       vec_det)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! Dimensions
    integer(is), intent(in)  :: nroots,csfdim,detdim

    ! Eigenvectors in the CSF and determinant bases
    real(dp), intent(in)     :: vec_csf(csfdim,nroots)
    real(dp), intent(out)    :: vec_det(detdim,nroots)

    ! CSF -> det transformation matrix
    real(dp), allocatable    :: transmat(:,:,:)
    
    ! Everything else
    integer(is)              :: n,nopen,ncsf,ndet,iconf
    integer(is)              :: icsf,jcsf,idet,jdet

!----------------------------------------------------------------------
! CSF -> det transformation matrix
!----------------------------------------------------------------------
    ! Allocate array
    allocate(transmat(maxdet,maxcsf,0:nomax))
    transmat=0.0d0

    ! Nopen = 0 case
    transmat(1,1,0)=1.0d0

    ! Nopen > 0 cases
    do n=1,nomax
       transmat(:,:,n)=transpose(csfcoe(:,:,n))
    enddo
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Eigenvectors in the determinant basis
    vec_det=0.0d0

    ! CSF and determinant counters
    icsf=1
    jcsf=0
    idet=1
    jdet=0
    
!----------------------------------------------------------------------
! Ref space confs
!----------------------------------------------------------------------
    do iconf=1,cfg%n0h

       ! No. open shells
       nopen=sop_nopen(cfg%sop0h(:,:,iconf),cfg%n_int_I)
       
       ! No. CSFs & determinants
       ndet=ndets(nopen)
       ncsf=ncsfs(nopen)

       ! End points in the CSF and determinant basis eigenvectors
       jcsf=icsf+ncsf-1
       jdet=idet+ndet-1

       ! Determinant coefficicients
       do n=1,nroots
          vec_det(idet:jdet,n)=matmul(transmat(1:ndet,1:ncsf,nopen),&
               vec_csf(icsf:jcsf,n))
       enddo
          
       ! Start points in CSF and determinant basis eigenvectors for
       ! the next iteration
       icsf=jcsf+1
       idet=jdet+1
       
    enddo

!----------------------------------------------------------------------
! 1I confs
!----------------------------------------------------------------------
    do iconf=1,cfg%n1I

       ! No. open shells
       nopen=sop_nopen(cfg%sop1I(:,:,iconf),n_int)
       
       ! No. CSFs & determinants
       ndet=ndets(nopen)
       ncsf=ncsfs(nopen)

       ! End points in the CSF and determinant basis eigenvectors
       jcsf=icsf+ncsf-1
       jdet=idet+ndet-1

       ! Determinant coefficicients
       do n=1,nroots
          vec_det(idet:jdet,n)=matmul(transmat(1:ndet,1:ncsf,nopen),&
               vec_csf(icsf:jcsf,n))
       enddo
       
       ! Start points in CSF and determinant basis eigenvectors for
       ! the next iteration
       icsf=jcsf+1
       idet=jdet+1
       
    enddo

!----------------------------------------------------------------------
! 2I confs
!----------------------------------------------------------------------
    do iconf=1,cfg%n2I

       ! No. open shells
       nopen=sop_nopen(cfg%sop2I(:,:,iconf),n_int)
       
       ! No. CSFs & determinants
       ndet=ndets(nopen)
       ncsf=ncsfs(nopen)

       ! End points in the CSF and determinant basis eigenvectors
       jcsf=icsf+ncsf-1
       jdet=idet+ndet-1

       ! Determinant coefficicients
       do n=1,nroots
          vec_det(idet:jdet,n)=matmul(transmat(1:ndet,1:ncsf,nopen),&
               vec_csf(icsf:jcsf,n))
       enddo
       
       ! Start points in CSF and determinant basis eigenvectors for
       ! the next iteration
       icsf=jcsf+1
       idet=jdet+1
       
    enddo

!----------------------------------------------------------------------
! 1E confs
!----------------------------------------------------------------------
    do iconf=1,cfg%n1E

       ! No. open shells
       nopen=sop_nopen(cfg%sop1E(:,:,iconf),n_int)
       
       ! No. CSFs & determinants
       ndet=ndets(nopen)
       ncsf=ncsfs(nopen)

       ! End points in the CSF and determinant basis eigenvectors
       jcsf=icsf+ncsf-1
       jdet=idet+ndet-1

       ! Determinant coefficicients
       do n=1,nroots
          vec_det(idet:jdet,n)=matmul(transmat(1:ndet,1:ncsf,nopen),&
               vec_csf(icsf:jcsf,n))
       enddo
       
       ! Start points in CSF and determinant basis eigenvectors for
       ! the next iteration
       icsf=jcsf+1
       idet=jdet+1
       
    enddo

!----------------------------------------------------------------------
! 2E confs
!----------------------------------------------------------------------
    do iconf=1,cfg%n2E

       ! No. open shells
       nopen=sop_nopen(cfg%sop2E(:,:,iconf),n_int)
       
       ! No. CSFs & determinants
       ndet=ndets(nopen)
       ncsf=ncsfs(nopen)

       ! End points in the CSF and determinant basis eigenvectors
       jcsf=icsf+ncsf-1
       jdet=idet+ndet-1

       ! Determinant coefficicients
       do n=1,nroots
          vec_det(idet:jdet,n)=matmul(transmat(1:ndet,1:ncsf,nopen),&
               vec_csf(icsf:jcsf,n))
       enddo
       
       ! Start points in CSF and determinant basis eigenvectors for
       ! the next iteration
       icsf=jcsf+1
       idet=jdet+1
       
    enddo

!----------------------------------------------------------------------
! 1I1E confs
!----------------------------------------------------------------------
    do iconf=1,cfg%n1I1E

       ! No. open shells
       nopen=sop_nopen(cfg%sop1I1E(:,:,iconf),n_int)
       
       ! No. CSFs & determinants
       ndet=ndets(nopen)
       ncsf=ncsfs(nopen)

       ! End points in the CSF and determinant basis eigenvectors
       jcsf=icsf+ncsf-1
       jdet=idet+ndet-1

       ! Determinant coefficicients
       do n=1,nroots
          vec_det(idet:jdet,n)=matmul(transmat(1:ndet,1:ncsf,nopen),&
               vec_csf(icsf:jcsf,n))
       enddo
       
       ! Start points in CSF and determinant basis eigenvectors for
       ! the next iteration
       icsf=jcsf+1
       idet=jdet+1
       
    enddo

    return
    
  end subroutine eigenvectors_detbas
  
!######################################################################
! bitstrings_detbas: returns the determinant basis bit strings given
!                    a set of MRCI configurations
!######################################################################
  subroutine bitstrings_detbas(cfg,detdim,det)

    use constants
    use bitglobal
    use conftype
    use mrciutils
    
    implicit none

    ! MRCI configuration derived type
    type(mrcfg), intent(in)  :: cfg

    ! No. determinants
    integer(is), intent(in)  :: detdim

    ! Determinant bit strings
    integer(ib), intent(out) :: det(n_int,2,detdim)

    ! Everything else
    integer(is)              :: nopen,ndet,iconf,idet,counter
    integer(is)              :: socc(nmo)


    ! TEST
    integer(is) :: i,j,k,n
    ! TEST
    
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! Determinant counter
    counter=0

    ! Determinant bit strings
    det=0_ib
    
!----------------------------------------------------------------------
! Ref space determinants
!----------------------------------------------------------------------
    ! Loop over ref space confs
    do iconf=1,cfg%n0h

       ! Determine the singly-occupied MOs
       call sop_socc_list(cfg%sop0h(:,:,iconf),cfg%n_int_I,socc,nmo,&
            nopen)
       
       ! No. determinants
       ndet=ndets(nopen)

       ! Loop over the determinants generated by this conf
       do idet=1,ndet

          ! Increment the determinant counter
          counter=counter+1

          ! Fill in the doubly-occupied MOs
          det(1:cfg%n_int_I,1,counter)=cfg%sop0h(:,2,iconf)
          det(1:cfg%n_int_I,2,counter)=cfg%sop0h(:,2,iconf)
          
          ! Fill in the singly-occupied MOs
          call fill_socc_mos(nopen,socc,idet,det(:,:,counter))

       enddo
       
    enddo

!----------------------------------------------------------------------
! 1I determinants
!----------------------------------------------------------------------
    ! Loop over 1I confs
    do iconf=1,cfg%n1I

       ! Determine the singly-occupied MOs
       call sop_socc_list(cfg%sop1I(:,:,iconf),n_int,socc,nmo,nopen)
       
       ! No. determinants
       ndet=ndets(nopen)

       ! Loop over the determinants generated by this conf
       do idet=1,ndet

          ! Increment the determinant counter
          counter=counter+1

          ! Fill in the doubly-occupied MOs
          det(:,1,counter)=cfg%sop1I(:,2,iconf)
          det(:,2,counter)=cfg%sop1I(:,2,iconf)
          
          ! Fill in the singly-occupied MOs
          call fill_socc_mos(nopen,socc,idet,det(:,:,counter))

       enddo
       
    enddo

!----------------------------------------------------------------------
! 2I determinants
!----------------------------------------------------------------------
    ! Loop over 2I confs
    do iconf=1,cfg%n2I

       ! Determine the singly-occupied MOs
       call sop_socc_list(cfg%sop2I(:,:,iconf),n_int,socc,nmo,nopen)
       
       ! No. determinants
       ndet=ndets(nopen)

       ! Loop over the determinants generated by this conf
       do idet=1,ndet

          ! Increment the determinant counter
          counter=counter+1

          ! Fill in the doubly-occupied MOs
          det(:,1,counter)=cfg%sop2I(:,2,iconf)
          det(:,2,counter)=cfg%sop2I(:,2,iconf)
          
          ! Fill in the singly-occupied MOs
          call fill_socc_mos(nopen,socc,idet,det(:,:,counter))

       enddo
       
    enddo
    
!----------------------------------------------------------------------
! 1E determinants
!----------------------------------------------------------------------
    ! Loop over 1E confs
    do iconf=1,cfg%n1E

       ! Determine the singly-occupied MOs
       call sop_socc_list(cfg%sop1E(:,:,iconf),n_int,socc,nmo,nopen)
       
       ! No. determinants
       ndet=ndets(nopen)

       ! Loop over the determinants generated by this conf
       do idet=1,ndet

          ! Increment the determinant counter
          counter=counter+1

          ! Fill in the doubly-occupied MOs
          det(:,1,counter)=cfg%sop1E(:,2,iconf)
          det(:,2,counter)=cfg%sop1E(:,2,iconf)
          
          ! Fill in the singly-occupied MOs
          call fill_socc_mos(nopen,socc,idet,det(:,:,counter))

       enddo
       
    enddo

!----------------------------------------------------------------------
! 2E determinants
!----------------------------------------------------------------------
    ! Loop over 2E confs
    do iconf=1,cfg%n2E

       ! Determine the singly-occupied MOs
       call sop_socc_list(cfg%sop2E(:,:,iconf),n_int,socc,nmo,nopen)
       
       ! No. determinants
       ndet=ndets(nopen)

       ! Loop over the determinants generated by this conf
       do idet=1,ndet

          ! Increment the determinant counter
          counter=counter+1

          ! Fill in the doubly-occupied MOs
          det(:,1,counter)=cfg%sop2E(:,2,iconf)
          det(:,2,counter)=cfg%sop2E(:,2,iconf)
          
          ! Fill in the singly-occupied MOs
          call fill_socc_mos(nopen,socc,idet,det(:,:,counter))

       enddo
       
    enddo

!----------------------------------------------------------------------
! 1I1E determinants
!----------------------------------------------------------------------
    ! Loop over 1I1E confs
    do iconf=1,cfg%n1I1E

       ! Determine the singly-occupied MOs
       call sop_socc_list(cfg%sop1I1E(:,:,iconf),n_int,socc,nmo,nopen)
       
       ! No. determinants
       ndet=ndets(nopen)

       ! Loop over the determinants generated by this conf
       do idet=1,ndet

          ! Increment the determinant counter
          counter=counter+1

          ! Fill in the doubly-occupied MOs
          det(:,1,counter)=cfg%sop1I1E(:,2,iconf)
          det(:,2,counter)=cfg%sop1I1E(:,2,iconf)
          
          ! Fill in the singly-occupied MOs
          call fill_socc_mos(nopen,socc,idet,det(:,:,counter))

       enddo
       
    enddo

    return
    
  end subroutine bitstrings_detbas

!######################################################################
! fill_socc_mos: fills in the singly-occupied MOs in a determinant
!                bit string given the no. open shells, the list of
!                open shell MOs and the determinant index
!                idet \in {1,...,ndets(nopen)}
!######################################################################
  subroutine fill_socc_mos(nopen,socc,idet,det)

    use constants
    use bitglobal
    
    implicit none

    ! Open shells
    integer(is), intent(in)    :: nopen
    integer(is), intent(in)    :: socc(nmo)

    ! Determinant bit string
    integer(ib), intent(inout) :: det(n_int,2)

    ! Determinant index
    integer(is), intent(in)    :: idet

    ! Everything else
    integer(is)                :: i,i1,k

    ! Loop over the singly-occupied MOs
    do i1=1,nopen
       
       ! Block index
       k=(socc(i1)-1)/n_bits+1

       ! Position of the MO within the block
       i=socc(i1)-(k-1)*n_bits-1

       ! Add the electron
       if (btest(detvec(idet,nopen),i1-1)) then
          ! Alpha-spin electron
          det(k,1)=ibset(det(k,1),i)
       else
          ! Beta-spin electron
          det(k,2)=ibset(det(k,2),i)
       endif
       
    enddo
    
    return
    
  end subroutine fill_socc_mos

!######################################################################
  
end module csf2det
