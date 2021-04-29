!**********************************************************************
! DFT/CIS Hamiltonian build routines
!**********************************************************************
module dftcis_hbuild

  implicit none

contains
  
!######################################################################
! save_hij_dftcis: Saves the non-zero DFT/CIS Hamiltonian matrix
!                  elements to disk
!######################################################################
  subroutine save_hij_dftcis(irrep,hscr,nrec,ncsf,iph)

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
    
    ! I/O variables
    integer(is)              :: iscratch
    character(len=60)        :: hamfile
    character(len=2)         :: amult,airrep

    ! Buffer
    integer(is)              :: nbuf
    integer(is), allocatable :: ibuffer(:,:)
    real(dp), allocatable    :: hbuffer(:)

    ! Everthing else
    integer(is)              :: nbra,nket,ibra,iket,abra,aket
    real(dp)                 :: hij
    
    ! Timing variables
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
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
    use int_pyscf
    
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
       hij =-hpar(1)*mo_integral(i,j,a,b)&
            +2.0d0*mo_integral(i,a,j,b)

    case(2)
       ! Our singlet/BHLYP parameterisation
       hij =-hpar(1)*mo_integral(i,j,a,b)&
            +2.0d0*mo_integral(i,a,j,b)
       
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
  
end module dftcis_hbuild