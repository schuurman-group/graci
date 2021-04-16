!**********************************************************************
! Routines for the generation of the CIS configuration information
!**********************************************************************
module dftcis_conf

  implicit none

contains

!######################################################################
! generate_cis_confs: Fills in the array of particle/hole indices
!######################################################################
  subroutine generate_cis_confs(irrep,ncsf,icvs,iph)

    use constants
    use bitglobal
    
    implicit none

    ! Irrep number
    integer(is), intent(in)  :: irrep
    
    ! No. CIS CSFs
    integer(is), intent(out) :: ncsf

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    
    ! CIS particle/hole indices
    integer(is), allocatable :: iph(:,:)

    ! Everything else
    integer(is)              :: modus
    
!----------------------------------------------------------------------
! Determine the no. CIS CSFs and allocate arrays
!----------------------------------------------------------------------
    modus=0
    call particle_hole_indices(modus,irrep,ncsf,icvs,iph)

    allocate(iph(2,ncsf))
    iph=0

!----------------------------------------------------------------------
! Fill in the particle/hole indices
!----------------------------------------------------------------------
    modus=1
    call particle_hole_indices(modus,irrep,ncsf,icvs,iph)
    
    return
    
  end subroutine generate_cis_confs
    
!######################################################################
! particle_hole_indices: determines the particle/hole indices of the
!                        allowable CIS CSFs
!######################################################################
  subroutine particle_hole_indices(modus,irrep,ncsf,icvs,iph)

    use constants
    use bitglobal
    use mrciutils
    
    implicit none

    ! Mode of operation:
    !
    ! modus = 0 <-> Determine the no. allowable CSFs
    ! modus = 1 <-> Fill in the particle/hole indices
    integer(is), intent(in)    :: modus

    ! Irrep number
    integer(is), intent(in)    :: irrep
    
    ! No. CSFs
    integer(is), intent(inout) :: ncsf

    ! CVS-MRCI: core MOs
    integer(is), intent(in)    :: icvs(nmo)
    logical                    :: lcvs
    
    ! Particle/hole indices
    integer(is), allocatable   :: iph(:,:)

    ! Working array
    integer(ib)                :: conf(n_int,2),sop(n_int,2)
    integer(is)                :: map(nmo)
    
    ! Everything else
    integer(is)                :: i,a,counter
    integer(is)                :: nocc,isym
    
!----------------------------------------------------------------------
! Initialisation
!----------------------------------------------------------------------
    ! No. CSFs
    if (modus == 0) ncsf=0

    ! Allowable CSF counter
    counter=0

    ! MO mapping array
    do i=1,nmo
       map(i)=i
    enddo

    ! Is this a CVS-DFT/CIS calculation?
    if (sum(icvs) > 0) then
       lcvs=.true.
    else
       lcvs=.false.
    endif
    
!----------------------------------------------------------------------
! Determine the allowable CSFs
!----------------------------------------------------------------------
    ! No. occupied MOs in the base determinant
    nocc=nel/2+mod(nel,2)

    ! Loop over occupied MOs
    do i=1,nocc

       ! If this is a CVS-DFT/CIS calculation, then cycle if we
       ! are at a non-flagged core MO
       if (lcvs .and. icvs(i) == 0) cycle
       
       ! Loop over unoccupied MOs
       do a=nocc+1,nmo

          ! Conf/SOP
          conf=annihilate_electron(conf0,n_int,i)
          conf=create_electron(conf,n_int,a)
          sop=conf_to_sop(conf,n_int)
          
          ! Irrep generated by this CSF
          isym=sop_sym_mrci(sop,map)

          ! Cycle if this is not the correct irrep
          if (isym /= irrep) cycle

          ! Increment the allowable CSF counter
          counter=counter+1

          ! Fill in the particle/hole indices
          if (modus == 1) then
             iph(1,counter)=a
             iph(2,counter)=i
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Set the no. allowable CSFs
!----------------------------------------------------------------------
    if (modus == 0) ncsf=counter
    
    return
    
  end subroutine particle_hole_indices
    
!######################################################################
  
end module dftcis_conf
