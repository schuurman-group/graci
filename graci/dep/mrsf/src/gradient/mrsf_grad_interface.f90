!**********************************************************************
! mrsf_grad_interface: C-bound entry points of the gradient kernels.
! All arrays are passed by address in Fortran (column-major) order with
! the dimensions of the loaded MRSF space (nocca, nvirb, naux, nmo).
!**********************************************************************
module mrsf_grad_interface

  use iso_c_binding
  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_gradient

  implicit none

contains

  subroutine mrsf_grad_init_c(erifile1, ierr1) bind(c, name='mrsf_grad_init')
    character(kind=c_char), intent(in) :: erifile1(*)
    integer(is), intent(out)           :: ierr1
    character(len=255)                 :: fname
    call c2fstr(erifile1, fname)
    call grad_init(fname, ierr1)
  end subroutine mrsf_grad_init_c

  subroutine mrsf_grad_free_c() bind(c, name='mrsf_grad_free')
    call grad_free()
  end subroutine mrsf_grad_free_c

  subroutine mrsf_grad_bso_c(bso1) bind(c, name='mrsf_grad_bso')
    real(dp), intent(out) :: bso1(naux,nmo,2)
    call grad_bso(bso1)
  end subroutine mrsf_grad_bso_c

  subroutine mrsf_gfock_c(nvec1, cx1, flags1, Ahh1, Za1, Zb1, jq1, Ghpa1, Ghpb1, &
       Ghha1, Ghhb1) bind(c, name='mrsf_gfock')
    integer(is), intent(in) :: nvec1, flags1
    real(dp), intent(in)    :: cx1
    real(dp), intent(in)    :: Ahh1(nocca,nocca,nvec1), Za1(nvirb,nocca,nvec1), &
                               Zb1(nvirb,nocca,nvec1), jq1(naux,nvec1)
    real(dp), intent(out)   :: Ghpa1(nocca,nvirb,nvec1), Ghpb1(nocca,nvirb,nvec1), &
                               Ghha1(nocca,nocca,nvec1), Ghhb1(nocca,nocca,nvec1)
    call gfock(nvec1, cx1, flags1, Ahh1, Za1, Zb1, jq1, Ghpa1, Ghpb1, Ghha1, Ghhb1)
  end subroutine mrsf_gfock_c

  subroutine mrsf_jblocks_c(nvec1, jq1, JHH1, JHP1, JPP1) bind(c, name='mrsf_jblocks')
    integer(is), intent(in) :: nvec1
    real(dp), intent(in)    :: jq1(naux,nvec1)
    real(dp), intent(out)   :: JHH1(nocca,nocca,nvec1), JHP1(nocca,nvirb,nvec1), &
                               JPP1(nvirb,nvirb,nvec1)
    call jblocks(nvec1, jq1, JHH1, JHP1, JPP1)
  end subroutine mrsf_jblocks_c

  subroutine mrsf_grad_state_c(cx1, Xt1, LaH1, LaP1, LbH1, LbP1, Ghh1, Fhp1, Yf1, &
       Sq1, gpp1, KbTHP1, KbTHH1, jT1, dq1) bind(c, name='mrsf_grad_state')
    real(dp), intent(in)  :: cx1, Xt1(nvirb,nocca)
    real(dp), intent(out) :: dq1(naux)
    real(dp), intent(out) :: LaH1(nocca,nocca), LaP1(nvirb,nocca), LbH1(nocca,nvirb), &
                             LbP1(nvirb,nvirb)
    real(dp), intent(out) :: Ghh1(nocca,nocca,naux), Fhp1(nocca,nvirb,naux), &
                             Yf1(nvirb,nocca,naux), Sq1(nocca,nocca,naux), gpp1(naux,naux)
    real(dp), intent(out) :: KbTHP1(nocca,nvirb), KbTHH1(nocca,nocca), jT1(naux)
    call grad_state(cx1, Xt1, LaH1, LaP1, LbH1, LbP1, Ghh1, Fhp1, Yf1, Sq1, gpp1, &
         KbTHP1, KbTHH1, jT1, dq1)
  end subroutine mrsf_grad_state_c

  subroutine mrsf_grad_finish_c(cx1, dq1, jT1, Ta1, Xt1, Za1, Zb1, Sq1, occH1, Ghh1, Fhp1, &
       gpp1, g1, pq1) bind(c, name='mrsf_grad_finish')
    real(dp), intent(in)    :: cx1, dq1(naux), jT1(naux), Ta1(nocca,nocca), Xt1(nvirb,nocca)
    real(dp), intent(in)    :: Za1(nvirb,nocca), Zb1(nvirb,nocca), Sq1(nocca,nocca,naux)
    real(dp), intent(in)    :: occH1(nocca)
    real(dp), intent(inout) :: Ghh1(nocca,nocca,naux), Fhp1(nocca,nvirb,naux)
    real(dp), intent(in)    :: gpp1(naux,naux)
    real(dp), intent(out)   :: g1(naux,naux), pq1(naux)
    call grad_finish(cx1, dq1, jT1, Ta1, Xt1, Za1, Zb1, Sq1, occH1, Ghh1, Fhp1, gpp1, g1, pq1)
  end subroutine mrsf_grad_finish_c

end module mrsf_grad_interface
