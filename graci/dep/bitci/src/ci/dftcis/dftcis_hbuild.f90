!**********************************************************************
! DFT/CIS Hamiltonian build routines
!**********************************************************************
module dftcis_hbuild

  use constants
  
  implicit none

  real(dp), allocatable :: Goo(:,:),Gvv(:,:),Gov(:,:)
  integer(is)           :: nocc,nvirt
    
contains
  
!######################################################################
! save_hij_dftcis: Saves the non-zero DFT/CIS Hamiltonian matrix
!                  elements to disk
!######################################################################
  subroutine save_hij_dftcis(irrep,hscr,nrec,ncsf,iph,loose)

    use constants
    use bitglobal
    use timing
    use iomod

    implicit none

    ! Irrep number
    integer(is), intent(in)  :: irrep
    
    ! No. CIS CSFs
    integer(is), intent(in)  :: ncsf
    
    ! CIS particle/hole indices
    integer(is), intent(in)  :: iph(2,ncsf)
    
    ! Hamiltonian scratch file number
    integer(is), intent(out) :: hscr

    ! Number of Hamiltonian scratch file records
    integer(is), intent(out) :: nrec

    ! Aggressively loose integral screening: only to be used for
    ! the generation of RAS spaces
    logical, intent(in)      :: loose
    
    ! I/O variables
    integer(is)              :: iscratch
    character(len=250)       :: hamfile
    character(len=2)         :: amult,airrep

    ! Buffer
    integer(is)              :: nbuf
    integer(is), allocatable :: ibuffer(:,:)
    real(dp), allocatable    :: hbuffer(:)

    ! Everthing else
    integer(is)              :: nbra,nket,ibra,iket,abra,aket
    integer(is)              :: i,j,a,b
    real(dp)                 :: hij
        
    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
    call get_times(twall_start,tcpu_start)

!----------------------------------------------------------------------
! Precompute G_pq = (p,q|p,q)^1/2
!----------------------------------------------------------------------
    if (loose) call precompute_gmat
    
!----------------------------------------------------------------------
! Register Hamiltonian scratch file
!----------------------------------------------------------------------
    write(amult,'(i0)') imult
    write(airrep,'(i0)') irrep
    call scratch_name('hij_dftcis.mult'//trim(amult)//'.sym'&
         //trim(airrep),hamfile)
    call register_scratch_file(hscr,hamfile)

!----------------------------------------------------------------------
! Open the Hamiltonian scratch file
!----------------------------------------------------------------------
    iscratch=scrunit(hscr)
    open(iscratch,file=scrname(hscr),form='unformatted',&
         status='unknown')

!----------------------------------------------------------------------
! Initialise the buffer
!----------------------------------------------------------------------
    allocate(ibuffer(2,bufsize))
    ibuffer=0

    allocate(hbuffer(bufsize))
    hbuffer=0.0d0

    nrec=0
    nbuf=0

!----------------------------------------------------------------------
! Compute the off-diagonal DFT/CIS Hamiltonian matrix elements
!----------------------------------------------------------------------
    ! Loop over bra CSFs
    do nbra=1,ncsf-1

       ! Bra particle index
       abra=iph(1,nbra)

       ! Bra hole index
       ibra=iph(2,nbra)

       ! Loop over ket CSFs
       do nket=nbra+1,ncsf

          ! Ket particle index
          aket=iph(1,nket)
          
          ! Key hole index
          iket=iph(2,nket)

          ! Loose integral screening
          if (loose)  then
             i=ibra
             a=abra-nocc
             j=iket
             b=aket-nocc
             if (0.5*Goo(i,j)*Gvv(a,b) < 1.2e-2_dp &
                  .and. 2*Gov(i,a)*Gov(j,b) < 1.2e-2_dp) cycle
          endif
                    
          ! Matrix element
          hij=hij_dftcis(ibra,abra,iket,aket)

          ! Save the matrix element if it's above threshold
          if (abs(hij) > epshij) then
             nbuf=nbuf+1
             ibuffer(1,nbuf)=nbra
             ibuffer(2,nbuf)=nket
             hbuffer(nbuf)=hij
             if (nbuf == bufsize) then
                write(iscratch) hbuffer,ibuffer,nbuf
                nbuf=0
                nrec=nrec+1
             endif
          endif
          
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
! Stop timing and print report
!----------------------------------------------------------------------
    call get_times(twall_end,tcpu_end)
    if (verbose) &
         call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
         'save_hij_dftcis')

    return
    
  end subroutine save_hij_dftcis
  
!######################################################################
! hij_dftcis: Computes a single off-diagonal element of the DFT/CIS
!             Hamiltonian matrix
!######################################################################
  function hij_dftcis(i,a,j,b) result(hij)

    use constants
    use bitglobal
    use dftcis_param
    
    implicit none

    ! Function result
    real(dp)                :: hij

    ! Hole indices
    integer(is), intent(in) :: a,b

    ! Particle indices
    integer(is), intent(in) :: i,j

    select case(ihamiltonian)

    case(1)

       ! Grimme's original parameterisation
       hij =-hpar(1)*bitci_ints%mo_int(i,j,a,b)&
            +2.0d0*bitci_ints%mo_int(i,a,j,b)

    case(2)

       ! Our singlet/BHLYP parameterisation
       hij =-hpar(1)*bitci_ints%mo_int(i,j,a,b)&
            +2.0d0*bitci_ints%mo_int(i,a,j,b)
            
    end select
    
    return
    
  end function hij_dftcis
  
!######################################################################
! hii_dftcis: Calculation of the on-diagonal elements of the DFT/CIS
!             Hamiltonian matrix
!######################################################################
  subroutine hii_dftcis(ncsf,iph,hii)

    use constants
    use bitglobal
    use dftcis_param
        
    implicit none

    ! No. CSFs
    integer(is), intent(in) :: ncsf

    ! Particle/hole indices
    integer(is), intent(in) :: iph(2,ncsf)

    ! On-diagonal Hamiltonian matrix elements
    real(dp), intent(out)   :: hii(ncsf)

    ! Everything else
    integer(is)             :: n,i,a

    ! Loop over CSFs
    do n=1,ncsf

       ! Particle index
       a=iph(1,n)

       ! Hole index
       i=iph(2,n)

       ! Compute the Hamiltonian matrix element
       select case(ihamiltonian)

       case(1)

          ! Grimme's original parameterisation
          hii(n)=moen(a)-moen(i) &
               -hpar(1)*Vc(i,a) &
               +2.0d0*Vx(i,a) &
               -0.025d0*moen(i)+hpar(2)*exp(-hpar(3)*Vx(i,a)**4)

       case(2)

          ! Our singlet/BHLYP parameterisation
          hii(n)=moen(a)-moen(i) &
               -hpar(1)*Vc(i,a) &
               +2.0d0*Vx(i,a) &
               -0.025d0*moen(i)+hpar(2)*exp(-hpar(3)*Vx(i,a)**4)

       end select
       
    enddo

    return
    
  end subroutine hii_dftcis
    
!######################################################################

  subroutine precompute_gmat

    use constants
    use bitglobal
        
    implicit none

    integer(is) :: i,j,a,b,a1,b1

    !
    ! MO subspace dimensions
    !
    nocc = nel/2
    nvirt = nmo - nocc

    !
    ! Allocate arrays
    !
    if (allocated(Goo)) deallocate(Goo)
    if (allocated(Gvv)) deallocate(Gvv)
    if (allocated(Gov)) deallocate(Gov)
    allocate(Goo(nocc,nocc), Gvv(nvirt,nvirt), Gov(nocc,nvirt))
    Goo=0.0d0; Gvv=0.0d0; Gov=0.0d0

    !
    ! Occ - occ elements
    !
    do i=1,nocc
       do j=i,nocc
          Goo(i,j)=sqrt(bitci_ints%mo_int(i,j,i,j))
          Goo(j,i)=Goo(i,j)
       enddo
    enddo

    !
    ! Virt - virt elements
    !
    do a1=1,nvirt
       a=nocc+a1
       do b1=a1,nvirt
          b=nocc+b1
          Gvv(a1,b1)=sqrt(bitci_ints%mo_int(a,b,a,b))
          Gvv(b1,a1)=Gvv(a1,b1)
       enddo
    enddo

    !
    ! Occ - virt elements
    !
    do i=1,nocc
       do a1=1,nvirt
          a=nocc+a1
          Gov(i,a1)=sqrt(bitci_ints%mo_int(i,a,i,a))
       enddo
    enddo
    
    return
    
  end subroutine precompute_gmat
    
!######################################################################
  
end module dftcis_hbuild
