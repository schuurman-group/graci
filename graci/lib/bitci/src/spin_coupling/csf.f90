!**********************************************************************
! Routines for the generation of CSFs using the genalogical spin
! coupling scheme.
!**********************************************************************
! The approach taken here is based directly on Section 2.6 of
! T. Helgaker, P. Jorgensen, and J. Olsen, Molecular Electronic 
! Structure Theory (Wiley, Chichester, 2000).
!**********************************************************************
! Important: The determinants entering into the CSF expansions are
!            implicitly in alpha-string, beta-string order.
!            The determinant encodings are stored in the detvec array
!            for later use. These bit strings correspond to the
!            lowercase p-vectors in Helgaker's book.
!**********************************************************************

module csf

  implicit none
  
contains

!######################################################################
! generate_csfs: Generates all the CSFs for the current spin
!                multiplicity (imult) up to the maximum number of open
!                shells + 2 (nocase2=nomax+2).
!######################################################################
  subroutine generate_csfs

    use constants
    use bitglobal
    use timing
    
    implicit none

    integer(ib), allocatable :: csfvec(:,:)
    integer(is)              :: i
    real(dp)                 :: tcpu_start,tcpu_end,twall_start,&
                                twall_end

!----------------------------------------------------------------------
! Output what we are doing
!----------------------------------------------------------------------
    write(6,'(/,52a)') ('-',i=1,52)
    write(6,'(x,a)') 'CSF generation'
    write(6,'(52a)') ('-',i=1,52)
    
!----------------------------------------------------------------------
! Start timing
!----------------------------------------------------------------------
  call get_times(twall_start,tcpu_start)
  
!----------------------------------------------------------------------
! Initialise and allocate arrays
!----------------------------------------------------------------------
    call init_csfs

!----------------------------------------------------------------------
! Output the dimensions of the CSF and deteminant bases
!----------------------------------------------------------------------
    call write_csf_info
    
!----------------------------------------------------------------------
! Generate the determinant bit string encodings for all numbers of open
! shells
!----------------------------------------------------------------------
    call get_det_bit_strings

!----------------------------------------------------------------------
! Generate the CSF bit string encodings for all numbers of open
! shells
!----------------------------------------------------------------------
    ! Initialise the array holding the CSF encodings
    allocate(csfvec(maxcsf,nocase2))
    csfvec=0_ib

    ! Generate the CSF encodings
    call get_csf_bit_strings(csfvec)

!----------------------------------------------------------------------
! Compute the CSF expansion coefficients in the determinant basis
!----------------------------------------------------------------------
    call get_csf_coeffs(csfvec)

!----------------------------------------------------------------------
! At this point, the CSF expansion coefficients correspond to
! determinants with spin-orbitals ordered according to the elements of
! the detvec array (the lowercase p-vectors in Helgaker's book).
! However, we are going to need the determinants in alpha-string,
! beta-string order. So, we will now multiply the CSF expansion
! coefficients by the phase factors associated with this reordering.
!----------------------------------------------------------------------
    call rephase_csf_coeffs
    
!----------------------------------------------------------------------
! Check that the CSFs are eigenfunctions of the total spin operator,
! S^2
!----------------------------------------------------------------------
    call check_csfs(csfvec)
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(csfvec)

!----------------------------------------------------------------------
! Stop timing and print report
!----------------------------------------------------------------------
  call get_times(twall_end,tcpu_end)
  call report_times(twall_end-twall_start,tcpu_end-tcpu_start,&
       'generate_csfs')
    
    return
    
  end subroutine generate_csfs
    
!######################################################################
! init_csfs: Allocation and initialisation of the arrays that will
!            hold the CSF information.
!######################################################################
  subroutine init_csfs

    use constants
    use bitglobal
    use math
    
    implicit none

    integer  :: nopen,na,nb
    real(dp) :: S,Shigh
    
!----------------------------------------------------------------------
! Number of CSFs and determinants as a function of the number of
! open shells
!----------------------------------------------------------------------
    ! Initialise the ncsfs and ndets arrays to zero
    ncsfs=0
    ndets=0

    ! Closed shell case
    ncsfs(0)=1
    ndets(0)=1
    
    ! Total spin
    S=dble(imult-1)/2.0d0

    ! Loop over open shells
    do nopen=1,nocase2

       ! Highest allowed total spin
       Shigh=dble(nopen)/2.0d0

       ! Cycle if the highest allowed total spin is less than
       ! the total spin value under consideration
       if (Shigh < S) cycle
       
       ! Cycle if there are no CSFs for the current (nopen,S) pair
       if (modulo(imult,2) == 0 .and. modulo(nopen,2) == 0) cycle
       if (modulo(imult,2) /= 0 .and. modulo(nopen,2) /= 0) cycle

       ! Number of alpha and beta electrons for S=M
       na=int(0.5d0*nopen+S)
       nb=int(0.5d0*nopen-S)
       
       ! Number of CSFs
       ncsfs(nopen)=int(bico(nopen,int(Shigh-S))*(2.0d0*S+1)&
            /(Shigh+S+1.0d0))

       ! Number of determinants
       ndets(nopen)=int(bico(nopen,int(Shigh+S)))

    enddo

!----------------------------------------------------------------------
! Allocate the CSF and determinant arrays
!----------------------------------------------------------------------
    ! Maximum number of CSFs across all numbers of open shells
    maxcsf=maxval(ncsfs)

    ! Maximum number of determinants across all numbers of open shells
    maxdet=maxval(ndets)
    
    ! CSF expansion coefficients
    allocate(csfcoe(maxcsf,maxdet,nocase2))
    csfcoe=0.0d0

    ! Determinant bit string vectors
    allocate(detvec(maxdet,nocase2))
    detvec=0_ib
    
    return
    
  end subroutine init_csfs

!######################################################################
! write_csf_info: outputs the dimensions of the CSF and determinant
!                 bases  
!######################################################################
  subroutine write_csf_info

    use constants
    use bitglobal
    use math
    
    implicit none

    integer(is) :: n,i
    
!----------------------------------------------------------------------
! Spin multiplicity
!----------------------------------------------------------------------
    write(6,'(/,x,a,x,i0)') 'Spin multiplicity:',imult
    
!----------------------------------------------------------------------
! Numbers of CSFs and determinants
!----------------------------------------------------------------------
    write(6,'(/,x,21a)') ('-',i=1,21)
    write(6,'(2x,a)') 'nopen | ncsf | ndet'
    write(6,'(x,21a)') ('-',i=1,21)
    do n=1,nocase2
       if (ncsfs(n) /= 0) then
          write(6,'(2x,i2,4x,a,x,i4,x,a,x,i4)') n,'|',ncsfs(n),'|',ndets(n)
       endif
    enddo
    write(6,'(x,21a)') ('-',i=1,21)
    
    return
    
  end subroutine write_csf_info
    
!######################################################################
! get_det_bit_strings: Generates the bit strings encoding the
!                      determinants in the expansions of the CSFs
!######################################################################
  subroutine get_det_bit_strings

    use constants
    use bitglobal
    use bitutils
    
    implicit none

    integer(is) :: nopen,na,nb
    real(dp)    :: S

!----------------------------------------------------------------------
! Generate the bit strings encoding the determinants entering into
! the expansion of the CSFs. Here, a set bit (1) is used to encode
! an alpha spin-orbital and an unset bit (0) is used to encode a beta
! spin-orbital.
!
! The determinants will be generates and stored in lexicographical
! order, starting with 11...100...0
!----------------------------------------------------------------------
    ! Total spin
    S=dble(imult-1)/2.0d0

    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Cycle if there are no CSFs for the current number of open
       ! shells
       if (ncsfs(nopen) == 0) cycle

       ! Number of alpha and beta electrons for S=M
       na=int(0.5d0*nopen+S)
       nb=int(0.5d0*nopen-S)
       
       ! Generate the determinant encodings for nopen open shells
       call get_permutations(na,nb,detvec(1:ndets(nopen),nopen),&
            ndets(nopen))

    enddo

    return
    
  end subroutine get_det_bit_strings
  
!######################################################################
! get_csf_bit_strings: Generates the bit strings encoding the CSFs
!######################################################################
  subroutine get_csf_bit_strings(csfvec)

    use constants
    use bitglobal
    use bitutils
    use math
    use iomod
    
    implicit none

    integer(ib), intent(out) :: csfvec(maxcsf,nocase2)
    integer(is)              :: nopen,na,nb,i
    integer(is)              :: stringdim
    integer(ib), allocatable :: bitstring(:)
    integer(is)              :: nvalid
    real(dp)                 :: S
    
!----------------------------------------------------------------------
! Generate the bit strings encoding the CSFs.
!
! Here, a set bit (1) is used to encode an increase in the total spin
! when coupling the associated electron, and an unset bit (0) is used
! to encode a decrease in total spin.
!
! Note that the first bit of each bit string encoding is always set.
!
! The CSFs will be stored in lexicographical order (minus those that
! correspond to forbidden paths through the branching diagram).
!----------------------------------------------------------------------
    ! Total spin
    S=dble(imult-1)/2.0d0

    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Cycle if there are no CSFs for the current number of open
       ! shells
       if (ncsfs(nopen) == 0) cycle

       ! Number of alpha and beta electrons for S=M
       na=int(0.5d0*nopen+S)
       nb=int(0.5d0*nopen-S)

       ! Number of permutations of na-1 alpha electrons in na-1+nb
       ! orbitals
       stringdim=int(bico(na-1+nb,na-1))

       ! Allocate the array of permutations
       allocate(bitstring(stringdim))
       bitstring=0_ib
       
       ! Generate the permutations of na-1 set bits and nb unset bits
       call get_permutations(na-1,nb,bitstring,stringdim)

       ! Right bit shift by 1 position and set the first bit to get
       ! the proper CSF encodings including the common first increase
       ! in S
       do i=1,stringdim
          bitstring(i)=ishft(bitstring(i),1)
          bitstring(i)=ibset(bitstring(i),0)
       enddo
       
       ! Filter out the paths that involve intermediate values of S < 0
       ! and save the valid paths in the csfvec array
       !
       ! Initialise the number of valid paths
       nvalid=0
       ! Loop over all permutations
       do i=1,stringdim
          ! Store the current permutation if it corresponds to a
          ! valid path
          if (valid_path(bitstring(i),na+nb)) then
             nvalid=nvalid+1
             csfvec(nvalid,nopen)=bitstring(i)
          endif
       enddo
       
       ! Check on the number of valid paths found
       if (nvalid /= ncsfs(nopen)) then
          write(errmsg,'(a)') 'Error in get_csf_bit_strings: '&
               //'incorrect number of valid paths found'
          call error_control
       endif

       ! Deallocate the array of permutations
       deallocate(bitstring)
       
    enddo

    return
    
  end subroutine get_csf_bit_strings

!######################################################################
! valid_path: Takes as it's argument a bit string corresponding to a
!             path through the genealogical spin coupling branching
!             diagram. Returns .true. if it is a valid path, and
!             .false. otherwise. Note that here a 'valid path' means
!             that all intermediate total spin values are greater than
!             or equal to zero.
!######################################################################
  function valid_path(bitstring,nbits)

    use constants
    
    implicit none

    logical                 :: valid_path
    integer(ib), intent(in) :: bitstring
    integer(is), intent(in) :: nbits
    integer(is)             :: k
    real(dp)                :: Sint

    ! Initialise S_int
    Sint=0.0d0

    ! Initialise valid_path to .true.
    valid_path=.true.
    
    ! Loop over bits
    do k=1,nbits
       ! Update the total spin
       if (btest(bitstring,k-1)) then
          Sint=Sint+0.5d0
       else
          Sint=Sint-0.5d0
       endif
       ! If the total spin has dropped below zero, then this is not a
       ! valid path
       if (Sint < 0.0d0) then
          valid_path=.false.
          return
       endif
    enddo
    
    return
    
  end function valid_path
  
!######################################################################
! get_csf_coeffs: Computes the coefficients of the CSFs expanded in
!                 the determinant basis
!######################################################################
  subroutine get_csf_coeffs(csfvec)

    use constants
    use bitglobal
    
    implicit none

    integer(ib), intent(in) :: csfvec(maxcsf,nocase2)
    integer(is)             :: nopen,icsf,idet,n
    real(dp)                :: S,M,dS,dM
    
!----------------------------------------------------------------------
! Compute the CSF expansion coefficients
!----------------------------------------------------------------------
    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Loop over CSFs for the current value of nopen
       do icsf=1,ncsfs(nopen)

          ! Loop over determinants for the current value of nopen
          do idet=1,ndets(nopen)

             ! Compute the projection of the determinant onto
             ! the CSF
             csfcoe(icsf,idet,nopen)=1.0d0
             S=0.0d0
             M=0.0d0
             do n=1,nopen

                ! Change in total spin
                if (btest(csfvec(icsf,nopen),n-1)) then
                   dS=0.5d0
                else
                   dS=-0.5d0
                endif
                S=S+dS
                
                ! Change in projected spin
                if (btest(detvec(idet,nopen),n-1)) then
                   dM=0.5d0
                else
                   dM=-0.5d0
                endif
                M=M+dM

                ! Cumulative multiplication by the Clebsch-Gordan
                ! coefficients
                csfcoe(icsf,idet,nopen)=csfcoe(icsf,idet,nopen)&
                     *clebsch_gordan(S,M,dS,dM)
                
             enddo

          enddo

       enddo
       
    enddo

    return
    
  end subroutine get_csf_coeffs

!######################################################################
! clebsch_gordan: Computes the value of the Clebsch-Gordan coefficient
!                 C_{dS,dM}^{S,M}
!######################################################################
  function clebsch_gordan(S,M,dS,dM) result(C)

    use constants
    use iomod
    
    implicit none

    real(dp)             :: C
    real(dp), intent(in) :: S,M,dS,dM

!----------------------------------------------------------------------
! Zero coefficient if |M| > S
!----------------------------------------------------------------------
    if (abs(M) > S) then
       C=0.0d0
       return
    endif
    
!----------------------------------------------------------------------
! Calculate the Clebsch-Gordan coefficient
!----------------------------------------------------------------------
    select case(int(dS*2))

    case(1)

       C=(S+2.0d0*dM*M)/(2.0d0*S)
       C=sqrt(C)

    case(-1)

       C=(S+1.0d0-2.0d0*dM*M)/(2.0d0*(S+1.0d0))
       C=-2.0d0*dM*sqrt(C)

    case default

       errmsg='Error in clebsch_gordan: illegal value of dS'
       call error_control

    end select
       
    return
    
  end function clebsch_gordan

!######################################################################
! rephase_csf_coeffs: Multiplication of the CSF expansion coefficients
!                     by the phase factors corresponding to putting
!                     the determinants into alpha-string, beta-string
!                     order.
!######################################################################
  subroutine rephase_csf_coeffs

    use constants
    use bitglobal
    
    implicit none

    integer(is), allocatable :: abfac(:,:)
    integer(is)              :: nopen,icsf,idet
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(abfac(maxdet,nocase2))
    
!----------------------------------------------------------------------
! Calculate the phase factors corresponding to putting the determinant
! bit string encodings into alpha-string, beta-string order. This has
! to be done because the Slater-Condon rule routines assume this
! ordering.
!----------------------------------------------------------------------
    call get_phases(abfac)
    
!----------------------------------------------------------------------
! Multiply the CSF expansion coefficients by the phase factors
!----------------------------------------------------------------------
    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Cycle if there are no CSFs for this number of open shells
       if (ncsfs(nopen) == 0) cycle

       ! Loop over CSFs for the current number of open shells
       do icsf=1,ncsfs(nopen)

          ! Loop over determinants for the current CSF
          do idet=1,ndets(nopen)

             ! Multiplication by the phase factor
             csfcoe(icsf,idet,nopen)=csfcoe(icsf,idet,nopen) &
                  *abfac(idet,nopen)
             
          enddo
          
       enddo
       
    enddo
    
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(abfac)
    
    return
    
  end subroutine rephase_csf_coeffs
    
!######################################################################
! check_csfs: Check on the CSFs. Specifically, on whether they are
!             normalised, and are eigenfunctions of the total and
!             projected spin operators.
!######################################################################
  subroutine check_csfs(csfvec)

    use constants
    use bitglobal
    use utils
    use iomod
    
    implicit none

    integer(ib), intent(in)  :: csfvec(maxcsf,nocase2)
    integer(is)              :: nopen,icsf,idet,jdet,n
    real(dp)                 :: norm
    real(dp), parameter      :: tiny=5e-12_dp
    real(dp), allocatable    :: s2mat(:,:)
    real(dp)                 :: S,S2,expec

!----------------------------------------------------------------------    
! Target S^2 expectation value
!----------------------------------------------------------------------
    ! Total spin
    S=dble(imult-1)/2.0d0

    ! <S^2>
    S2=S*(S+1.0d0)
    
!----------------------------------------------------------------------
! (1) Check on the normalisation of the CSFs
!----------------------------------------------------------------------
    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Loop over CSFs for the current number of open shells
       do icsf=1,ncsfs(nopen)

          ! Compute the CSF norm
          norm=sqrt(dot_product(csfcoe(icsf,:,nopen),&
               csfcoe(icsf,:,nopen)))

          ! Exit if the deviation of the norm from unity is above
          ! threshold
          if (abs(norm-1.0d0) > tiny) then
             errmsg='Error in check_csfs: Non-unit norm found.'
             call error_control
          endif
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! (2) Check on the S^2 expectation values
!----------------------------------------------------------------------
    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Cycle if there are no CSFs for this number of open shells
       if (ncsfs(nopen) == 0) cycle

       ! Allocate the determinant S^2 matrix array
       allocate(s2mat(ndets(nopen),ndets(nopen)))
       
       ! Compute the determinant representation of and S^2
       call determinant_s2mat(s2mat,nopen)

       ! Loop over CSFs for the current number of open shells
       do icsf=1,ncsfs(nopen)

          ! Calculate <S^2> (accounting for the phase factors)
          expec=0.0d0
          do idet=1,ndets(nopen)
             do jdet=1,ndets(nopen)
                expec=expec+csfcoe(icsf,idet,nopen)*s2mat(idet,jdet)&
                     *csfcoe(icsf,jdet,nopen)
             enddo
          enddo

          ! Exit if the deviation of <S^2> from S*(S+1) is above
          ! threshold
          if (abs(expec-S2) > tiny) then
             errmsg='Error in check_csfs: Incorrect <S^2> value found.'
             call error_control
          endif
          
       enddo

       ! Deallocate the determinant S^2 matrix array
       deallocate(s2mat)
       
    enddo

    return
    
  end subroutine check_csfs

!######################################################################
! get_phases: Computes the phase factors corresponding to putting the
!             determinants into alpha-string, beta-string order.
!######################################################################
  subroutine get_phases(abfac)
  
    use constants
    use bitglobal
    use slater_condon
    
    implicit none

    integer(is), intent(out) :: abfac(maxdet,nocase2)
    integer(is)              :: nopen,idet,i
    integer(is)              :: na
    real(dp)                 :: S

    ! Total spin
    S=dble(imult-1)/2.0d0
    
    ! Loop over numbers of open shells
    do nopen=1,nocase2

       ! Number of alpha electrons for S=M
       na=int(0.5d0*nopen+S)
              
       ! Loop over determinants for nopen open shells
       do idet=1,ndets(nopen)
          
          ! Sign change incurred by the spin-orbital permutations
          ! required to put the determinant into alpha-string,
          ! beta-string order
          abfac(idet,nopen)=ab_phase(detvec(idet,nopen))
          
       enddo
       
    enddo    
    
    return
    
  end subroutine get_phases

!######################################################################
! ab_phase: Computes the phase factor associated with putting a given
!           determinant (in ... + ... - ... form) into alpha-string,
!           beta-string order.
!           Here, vec is the bit string encoding of the determinant.
!######################################################################
  function ab_phase(vec)

    use constants
    use bitglobal
    
    implicit none

    integer(is)             :: ab_phase
    integer(ib), intent(in) :: vec
    integer(is)             :: i,nperm,na
    integer(ib)             :: mask,abpos,apos,bpos
    integer(is)             :: ma(nocase2),mb(nocase2)
    integer(is)             :: namv,nbmv
    
!----------------------------------------------------------------------
! Bit mask for the target order of the spin-orbitals
!----------------------------------------------------------------------
    ! Number of alpha electrons
    na=popcnt(vec)
    
    ! Set the bit mask
    mask=0_ib
    do i=1,na
       mask=ibset(mask,i-1)
    enddo

!----------------------------------------------------------------------
! Positions of the alpha and beta spin-orbitals that are in the 'wrong'
! place
!----------------------------------------------------------------------
    ma=0
    mb=0
    abpos=ieor(mask,vec)
    apos=iand(abpos,vec)
    bpos=iand(abpos,not(vec))
    namv=popcnt(apos)
    nbmv=popcnt(bpos)
    call set_bit_indices(apos,ma(1:namv),namv)
    call set_bit_indices(bpos,mb(1:nbmv),nbmv)

!----------------------------------------------------------------------
! Determine the number of permutations required to but the determinant
! into alpha-string, beta-string order
!----------------------------------------------------------------------
    nperm=0
    do i=1,namv
       nperm=nperm+ma(i)-mb(i)
    enddo

!----------------------------------------------------------------------
! Phase factor
!----------------------------------------------------------------------
    ab_phase=(-1)**nperm
    
    return
    
  end function ab_phase
    
!######################################################################
! set_bit_indices: Determines the indices of the set bits in the bit
!                  string I and returns them in the array list.
!######################################################################
  subroutine set_bit_indices(I,list,listdim)

    use constants
    use bitglobal
    
    implicit none

    integer(is), intent(in)    :: listdim
    integer(ib), intent(inout) :: I
    integer(is), intent(out)   :: list(listdim)
    integer(is)                :: k,n
    integer(is)                :: e
    
    ! Initialisation
    list=0
    n=1
    
    ! Determine the indices of any set bits
    do while (I /= 0_ib)
       e=trailz(I)
       I=ibclr(I,e)
       list(n)=e+1
       n=n+1
    enddo
    
    return
    
  end subroutine set_bit_indices
  
!######################################################################
! determinant_s2mat: Computes the determinant representation of S^2
!                    for a given number of open shells.
!                    Note that here it is assumed that the determinants
!                    are in alpha-string, beta-string order.
!###################################################################### 
  subroutine determinant_s2mat(s2mat,nopen)

    use constants
    use bitglobal
    use bitutils
    use slater_condon
        
    implicit none

    integer(is), intent(in) :: nopen
    real(dp), intent(out)   :: s2mat(ndets(nopen),ndets(nopen))

    integer(is)              :: idet,ibra,iket,n
    integer(ib), allocatable :: d(:,:)
    integer(is)              :: n_int_save
    integer(ib)              :: phase_mask(2)
    integer(is)              :: nexci
    integer(ib)              :: p(2),h(2)
    integer(is), parameter   :: maxex=2
    integer(is)              :: plist(maxex,2),hlist(maxex,2)
    integer(is)              :: phase
    
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
    allocate(d(2,ndets(nopen)))

!----------------------------------------------------------------------
! Save the actual value of n_int and then set this to 1 for use in the
! following. This allows us to use the slater_condon module to
! calculate the S^2 matrix.
!----------------------------------------------------------------------
    n_int_save=n_int
    n_int=1
    
!----------------------------------------------------------------------
! Construct the determinant bit strings
!----------------------------------------------------------------------
    ! Initialise to all unset bits
    d=0_ib

    ! Loop over determinants
    do idet=1,ndets(nopen)
       ! Loop over open shells
       do n=1,nopen
          if (btest(detvec(idet,nopen),n-1)) then
             ! Occupied alpha spin-orbital
             d(1,idet)=ibset(d(1,idet),n-1)
          else
             ! Occupied alpha spin-orbital
             d(2,idet)=ibset(d(2,idet),n-1)
          endif
       enddo
    enddo

!----------------------------------------------------------------------
! Initialise the S^2 matrix to zero
!----------------------------------------------------------------------
    s2mat=0.0d0
    
!----------------------------------------------------------------------
! On-diagonal elements
!----------------------------------------------------------------------
    do idet=1,ndets(nopen)
       s2mat(idet,idet)=s2ii(d(:,idet))
    enddo

!----------------------------------------------------------------------
! Off-diagonal elements
!----------------------------------------------------------------------
    ! Loop over ket determinants
    do iket=1,ndets(nopen)-1

       ! Compute the phase mask
       call phasemask(d(:,iket),phase_mask)

       ! Loop over bra determinants
       do ibra=iket+1,ndets(nopen)

          ! Cycle if the no. excitations connecting the bra and
          ! ket determinants is greater than 2
          nexci=exc_degree_det(d(:,ibra),d(:,iket))
          if (nexci > 2) cycle

          ! Get the indices of the spin-orbitals involved in the
          ! excitations linking the bra and ket determinants
          call exc(d(:,iket),d(:,ibra),p,h)
          call list_from_bitstring(p(ialpha),plist(:,ialpha),maxex)
          call list_from_bitstring(h(ialpha),hlist(:,ialpha),maxex)
          call list_from_bitstring(p(ibeta),plist(:,ibeta),maxex)
          call list_from_bitstring(h(ibeta),hlist(:,ibeta),maxex)

          ! Compute the phase factor
          phase=phase_pure_exc(phase_mask(ialpha),hlist(:,ialpha),&
               plist(:,ialpha),maxex,nexci)&
               *phase_pure_exc(phase_mask(ibeta),hlist(:,ibeta),&
               plist(:,ibeta),maxex,nexci)

          ! Compute the S^2 matrix element
          s2mat(ibra,iket)=s2ij(phase,hlist,plist)
          s2mat(iket,ibra)=s2mat(ibra,iket)
          
       enddo
       
    enddo

!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
    deallocate(d)

!----------------------------------------------------------------------
! Reset n_int
!----------------------------------------------------------------------
    n_int=n_int_save
    
    return
    
  end subroutine determinant_s2mat
  
!######################################################################
  
end module csf
