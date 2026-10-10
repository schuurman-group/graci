!**********************************************************************
! mrsf_xcgrid_interface: C-bound entry points of the grid module
!**********************************************************************
module mrsf_xcgrid_interface

  use iso_c_binding
  use mrsf_constants
  use mrsf_global
  use mrsf_io
  use mrsf_xcgrid

  implicit none

contains

  subroutine mrsf_xc_init_c(nao1, ngrid1, ncomp1, nv1, nblocks1, gbeg1, gend1, wgt1, fxc1, vxc1, &
       CH1, CP1, cache1) bind(c, name='mrsf_xc_init')
    integer(is), intent(in) :: nao1, ngrid1, ncomp1, nv1, nblocks1
    integer(is), intent(in) :: gbeg1(nblocks1), gend1(nblocks1)
    real(dp), intent(in)    :: wgt1(ngrid1), fxc1(ngrid1,nv1,2,nv1,2), vxc1(ngrid1,nv1,2)
    real(dp), intent(in)    :: CH1(nao1,nocca), CP1(nao1,nvirb)
    logical(c_bool), intent(in) :: cache1
    call xc_init(nao1, ngrid1, ncomp1, nv1, nblocks1, gbeg1, gend1, wgt1, fxc1, vxc1, CH1, CP1, &
         logical(cache1))
  end subroutine mrsf_xc_init_c

  subroutine mrsf_xc_add_block_c(ib1, ao1) bind(c, name='mrsf_xc_add_block')
    integer(is), intent(in) :: ib1
    real(dp), intent(in)    :: ao1(*)
    call xc_add_block(ib1, ao1)
  end subroutine mrsf_xc_add_block_c

  subroutine mrsf_xc_free_c() bind(c, name='mrsf_xc_free')
    call xc_free()
  end subroutine mrsf_xc_free_c

  subroutine mrsf_xc_begin_c(nvec1, k1, Lfac1, Rfac1, rch1) bind(c, name='mrsf_xc_begin')
    integer(is), intent(in) :: nvec1, k1, rch1(2,nvec1)
    real(dp), intent(in)    :: Lfac1(nao_g,k1,2,nvec1), Rfac1(nao_g,k1,2,nvec1)
    call xc_begin(nvec1, k1, Lfac1, Rfac1, rch1)
  end subroutine mrsf_xc_begin_c

  subroutine mrsf_xc_block_c(ib1, ao1) bind(c, name='mrsf_xc_block')
    integer(is), intent(in) :: ib1
    real(dp), intent(in)    :: ao1(*)
    call xc_block(ib1, ao1)
  end subroutine mrsf_xc_block_c

  subroutine mrsf_xc_cached_c() bind(c, name='mrsf_xc_cached')
    call xc_cached()
  end subroutine mrsf_xc_cached_c

  subroutine mrsf_xc_end_c(wanthp1, wanthh1, VHP1, VHH1) bind(c, name='mrsf_xc_end')
    logical(c_bool), intent(in) :: wanthp1, wanthh1
    real(dp), intent(out)       :: VHP1(nocca,nvirb,2,nvcur), VHH1(nocca,nocca,2,nvcur)
    call xc_end(logical(wanthp1), logical(wanthh1), VHP1, VHH1)
  end subroutine mrsf_xc_end_c

  subroutine mrsf_xc_probe_begin_c(nst1, k1, Lfac1, Rfac1, rch1) bind(c, name='mrsf_xc_probe_begin')
    integer(is), intent(in) :: nst1, k1, rch1(2,nst1)
    real(dp), intent(in)    :: Lfac1(nao_g,k1,2,nst1), Rfac1(nao_g,k1,2,nst1)
    call xc_probe_begin(nst1, k1, Lfac1, Rfac1, rch1)
  end subroutine mrsf_xc_probe_begin_c

  subroutine mrsf_xc_probe_block_c(ib1, nc2, ao1) bind(c, name='mrsf_xc_probe_block')
    integer(is), intent(in) :: ib1, nc2
    real(dp), intent(in)    :: ao1(*)
    call xc_probe_block(ib1, nc2, ao1)
  end subroutine mrsf_xc_probe_block_c

  subroutine mrsf_xc_probe_end_c(t1) bind(c, name='mrsf_xc_probe_end')
    real(dp), intent(out) :: t1(nao_g,3,nvcur+1)
    call xc_probe_end(t1)
  end subroutine mrsf_xc_probe_end_c

  subroutine mrsf_xc_g_set_c(ngrid1, nv1, vxc1, fxc1, kxc1, havek1) bind(c, name='mrsf_xc_g_set')
    integer(is), intent(in)     :: ngrid1, nv1
    real(dp), intent(in)        :: vxc1(ngrid1,nv1), fxc1(ngrid1,nv1,nv1), kxc1(*)
    logical(c_bool), intent(in) :: havek1
    call xc_g_set(ngrid1, nv1, vxc1, fxc1, kxc1, logical(havek1))
  end subroutine mrsf_xc_g_set_c

  subroutine mrsf_xc_g_free_c() bind(c, name='mrsf_xc_g_free')
    call xc_g_free()
  end subroutine mrsf_xc_g_free_c

  subroutine mrsf_xc_g_potential_c(kM1, LM1, RM1, kS1, LS1, MM1, MK1, V11) bind(c, name='mrsf_xc_g_potential')
    integer(is), intent(in) :: kM1, kS1
    real(dp), intent(in)    :: LM1(nao_g,kM1), RM1(nao_g,kM1), LS1(nao_g,kS1)
    real(dp), intent(out)   :: MM1(nocca,nao_g), MK1(nocca,nao_g), V11(nao_g,nao_g)
    call xc_g_potential(kM1, LM1, RM1, kS1, LS1, MM1, MK1, V11)
  end subroutine mrsf_xc_g_potential_c

  subroutine mrsf_xc_gprobe_set_c(nst1, kM1, LM1, RM1, kS1, LS1) bind(c, name='mrsf_xc_gprobe_set')
    integer(is), intent(in) :: nst1, kM1, kS1
    real(dp), intent(in)    :: LM1(nao_g,kM1,nst1), RM1(nao_g,kM1,nst1), LS1(nao_g,kS1,nst1)
    call xc_gprobe_set(nst1, kM1, LM1, RM1, kS1, LS1)
  end subroutine mrsf_xc_gprobe_set_c

  subroutine mrsf_xc_gprobe_get_c(tG1) bind(c, name='mrsf_xc_gprobe_get')
    real(dp), intent(out) :: tG1(nao_g,3,nstg)
    call xc_gprobe_get(tG1)
  end subroutine mrsf_xc_gprobe_get_c

  subroutine mrsf_xc_timings_c() bind(c, name='mrsf_xc_timings')
    call xc_timings()
  end subroutine mrsf_xc_timings_c

end module mrsf_xcgrid_interface
