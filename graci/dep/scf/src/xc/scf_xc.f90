!**********************************************************************
! scf_xc: exchange-correlation functionals through libxc's C API (the
! libxc build shipped with PySCF, so the parameters are PySCF's). The
! functional is the list of (libxc id, family, coefficient) that PySCF's
! parse_xc returns; the components are evaluated and summed with their
! coefficients (as pyscf.dft.libxc.eval_xc does), without density
! thresholds or external parameters of our own. Outputs are PySCF's
! potential weights per grid point:
!   polarised:   w_s(1) = 1/2 w vrho_s,
!                w_s(1+c) = w (2 vsigma_ss d_c rho_s + vsigma_ab d_c rho_s')
!   unpolarised: w(1) = 1/2 w vrho, w(1+c) = 2 w vsigma d_c rho
! (the 1/2 on the density weight goes with the V = T + T^T assembly).
! Evaluations with the shared, initialised functionals are made from
! several OpenMP threads, as PySCF's libxc interface does.
!**********************************************************************
module scf_xc

  use scf_constants
  use iso_c_binding

  implicit none

  integer(is)                 :: nfunc = 0, xc_fam = 0, xc_nspin = 0   ! xc_fam: 0 none, 1 LDA, 2 GGA
  integer(c_int), allocatable :: fid(:), ffam(:)
  real(dp), allocatable       :: fcoef(:)
  type(c_ptr), allocatable    :: fptr(:)

  interface
     function xc_func_alloc() bind(c, name='xc_func_alloc') result(p)
       import :: c_ptr
       type(c_ptr) :: p
     end function xc_func_alloc
     function xc_func_init(p, id, nspin) bind(c, name='xc_func_init') result(ierr)
       import :: c_ptr, c_int
       type(c_ptr), value    :: p
       integer(c_int), value :: id, nspin
       integer(c_int)        :: ierr
     end function xc_func_init
     subroutine xc_func_end(p) bind(c, name='xc_func_end')
       import :: c_ptr
       type(c_ptr), value :: p
     end subroutine xc_func_end
     subroutine xc_func_free(p) bind(c, name='xc_func_free')
       import :: c_ptr
       type(c_ptr), value :: p
     end subroutine xc_func_free
     subroutine xc_lda_exc_vxc(p, np, rho, zk, vrho) bind(c, name='xc_lda_exc_vxc')
       import :: c_ptr, c_size_t, c_double
       type(c_ptr), value       :: p
       integer(c_size_t), value :: np
       real(c_double)           :: rho(*), zk(*), vrho(*)
     end subroutine xc_lda_exc_vxc
     subroutine xc_gga_exc_vxc(p, np, rho, sigma, zk, vrho, vsigma) bind(c, name='xc_gga_exc_vxc')
       import :: c_ptr, c_size_t, c_double
       type(c_ptr), value       :: p
       integer(c_size_t), value :: np
       real(c_double)           :: rho(*), sigma(*), zk(*), vrho(*), vsigma(*)
     end subroutine xc_gga_exc_vxc
  end interface

contains

!######################################################################
! xc_setup: (re)initialise the functional components for nspin (1 or 2);
! fam(i): 1 LDA, 2 GGA; ierr 1: unsupported family, 2: libxc init error
!######################################################################
  subroutine xc_setup(nf, ids, fam, coefs, nspin, ierr)

    integer(is), intent(in)  :: nf, ids(nf), fam(nf), nspin
    real(dp), intent(in)     :: coefs(nf)
    integer(is), intent(out) :: ierr
    integer(is)              :: i

    ierr = 0
    call xc_release()
    nfunc = nf
    xc_fam = 0
    if (nf == 0) return
    allocate(fid(nf), ffam(nf), fcoef(nf), fptr(nf))
    do i = 1, nf
       if (fam(i) /= 1 .and. fam(i) /= 2) then
          ierr = 1
          return
       endif
       fid(i) = int(ids(i), c_int)
       ffam(i) = int(fam(i), c_int)
       fcoef(i) = coefs(i)
       xc_fam = max(xc_fam, fam(i))
       fptr(i) = xc_func_alloc()
       if (xc_func_init(fptr(i), fid(i), int(nspin, c_int)) /= 0) then
          ierr = 2
          return
       endif
    enddo
    xc_nspin = nspin

  end subroutine xc_setup

  subroutine xc_release()
    integer(is) :: i
    if (allocated(fptr)) then
       do i = 1, size(fptr)
          if (c_associated(fptr(i))) then
             call xc_func_end(fptr(i))
             call xc_func_free(fptr(i))
          endif
       enddo
       deallocate(fptr, fid, ffam, fcoef)
    endif
    nfunc = 0
    xc_fam = 0
    xc_nspin = 0
  end subroutine xc_release

!######################################################################
! xc_eval_pol: spin-polarised evaluation on np points. ra, rb (np,4):
! density and its x, y, z derivatives per spin (only column 1 for LDA);
! w: weights; wa, wb (np,4) potential weights; exc += sum w (ra+rb) exc
!######################################################################
  subroutine xc_eval_pol(np, ra, rb, w, wa, wb, exc)

    integer(is), intent(in) :: np
    real(dp), intent(in)    :: ra(np,4), rb(np,4), w(np)
    real(dp), intent(out)   :: wa(np,4), wb(np,4)
    real(dp), intent(inout) :: exc

    real(dp), allocatable :: rho(:,:), sig(:,:), zk(:), vr(:,:), vs(:,:), zt(:), vrt(:,:), vst(:,:)
    integer(is) :: i, g, c
    integer(c_size_t) :: n

    n = int(np, c_size_t)
    allocate(rho(2,np), sig(3,np), zk(np), vr(2,np), vs(3,np), zt(np), vrt(2,np), vst(3,np))
    do g = 1, np
       rho(1,g) = ra(g,1)
       rho(2,g) = rb(g,1)
    enddo
    if (xc_fam == 2) then
       do g = 1, np
          sig(1,g) = ra(g,2)*ra(g,2) + ra(g,3)*ra(g,3) + ra(g,4)*ra(g,4)
          sig(2,g) = ra(g,2)*rb(g,2) + ra(g,3)*rb(g,3) + ra(g,4)*rb(g,4)
          sig(3,g) = rb(g,2)*rb(g,2) + rb(g,3)*rb(g,3) + rb(g,4)*rb(g,4)
       enddo
    endif
    zk = 0.0_dp; vr = 0.0_dp; vs = 0.0_dp
    do i = 1, nfunc
       if (ffam(i) == 2) then
          call xc_gga_exc_vxc(fptr(i), n, rho, sig, zt, vrt, vst)
          vs = vs + fcoef(i)*vst
       else
          call xc_lda_exc_vxc(fptr(i), n, rho, zt, vrt)
       endif
       zk = zk + fcoef(i)*zt
       vr = vr + fcoef(i)*vrt
    enddo
    do g = 1, np
       exc = exc + w(g)*(rho(1,g) + rho(2,g))*zk(g)
       wa(g,1) = 0.5_dp*w(g)*vr(1,g)
       wb(g,1) = 0.5_dp*w(g)*vr(2,g)
    enddo
    if (xc_fam == 2) then
       do c = 2, 4
          do g = 1, np
             wa(g,c) = w(g)*(2.0_dp*vs(1,g)*ra(g,c) + vs(2,g)*rb(g,c))
             wb(g,c) = w(g)*(2.0_dp*vs(3,g)*rb(g,c) + vs(2,g)*ra(g,c))
          enddo
       enddo
    else
       wa(:,2:4) = 0.0_dp
       wb(:,2:4) = 0.0_dp
    endif
    deallocate(rho, sig, zk, vr, vs, zt, vrt, vst)

  end subroutine xc_eval_pol

!######################################################################
! xc_eval_unpol: closed-shell evaluation; r (np,4) total density and its
! derivatives; wv (np,4) potential weights; exc += sum w r exc
!######################################################################
  subroutine xc_eval_unpol(np, r, w, wv, exc)

    integer(is), intent(in) :: np
    real(dp), intent(in)    :: r(np,4), w(np)
    real(dp), intent(out)   :: wv(np,4)
    real(dp), intent(inout) :: exc

    real(dp), allocatable :: rho(:), sig(:), zk(:), vr(:), vs(:), zt(:), vrt(:), vst(:)
    integer(is) :: i, g, c
    integer(c_size_t) :: n

    n = int(np, c_size_t)
    allocate(rho(np), sig(np), zk(np), vr(np), vs(np), zt(np), vrt(np), vst(np))
    do g = 1, np
       rho(g) = r(g,1)
    enddo
    if (xc_fam == 2) then
       do g = 1, np
          sig(g) = r(g,2)*r(g,2) + r(g,3)*r(g,3) + r(g,4)*r(g,4)
       enddo
    endif
    zk = 0.0_dp; vr = 0.0_dp; vs = 0.0_dp
    do i = 1, nfunc
       if (ffam(i) == 2) then
          call xc_gga_exc_vxc(fptr(i), n, rho, sig, zt, vrt, vst)
          vs = vs + fcoef(i)*vst
       else
          call xc_lda_exc_vxc(fptr(i), n, rho, zt, vrt)
       endif
       zk = zk + fcoef(i)*zt
       vr = vr + fcoef(i)*vrt
    enddo
    do g = 1, np
       exc = exc + w(g)*rho(g)*zk(g)
       wv(g,1) = 0.5_dp*w(g)*vr(g)
    enddo
    if (xc_fam == 2) then
       do c = 2, 4
          do g = 1, np
             wv(g,c) = 2.0_dp*w(g)*vs(g)*r(g,c)
          enddo
       enddo
    else
       wv(:,2:4) = 0.0_dp
    endif
    deallocate(rho, sig, zk, vr, vs, zt, vrt, vst)

  end subroutine xc_eval_unpol

end module scf_xc
