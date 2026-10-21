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

  subroutine mrsf_xc_gprobe_set_c(nst1, kM1, LM1, RM1, kS1, LS1, usek1) bind(c, name='mrsf_xc_gprobe_set')
    integer(is), intent(in) :: nst1, kM1, kS1
    real(dp), intent(in)    :: LM1(nao_g,kM1,nst1), RM1(nao_g,kM1,nst1), LS1(nao_g,kS1,nst1)
    logical(c_bool), intent(in) :: usek1
    call xc_gprobe_set(nst1, kM1, LM1, RM1, kS1, LS1, logical(usek1))
  end subroutine mrsf_xc_gprobe_set_c

  subroutine mrsf_xc_gprobe_get_c(tG1) bind(c, name='mrsf_xc_gprobe_get')
    real(dp), intent(out) :: tG1(nao_g,3,nstg)
    call xc_gprobe_get(tG1)
  end subroutine mrsf_xc_gprobe_get_c

  subroutine mrsf_xc_init_grid_c(nao1, ngrid1, ncomp1, nv1, nblocks1, gbeg1, gend1, wgt1, CH1, CP1, cache1) &
       bind(c, name='mrsf_xc_init_grid')
    integer(is), intent(in) :: nao1, ngrid1, ncomp1, nv1, nblocks1
    integer(is), intent(in) :: gbeg1(nblocks1), gend1(nblocks1)
    real(dp), intent(in)    :: wgt1(ngrid1)
    real(dp), intent(in)    :: CH1(nao1,nocca), CP1(nao1,nvirb)
    logical(c_bool), intent(in) :: cache1
    call xc_init_grid(nao1, ngrid1, ncomp1, nv1, nblocks1, gbeg1, gend1, wgt1, CH1, CP1, logical(cache1))
  end subroutine mrsf_xc_init_grid_c

  subroutine mrsf_xc_set_ref_kernel_c(ngrid1, nv1, vxc1, fxc1) bind(c, name='mrsf_xc_set_ref_kernel')
    integer(is), intent(in) :: ngrid1, nv1
    real(dp), intent(in)    :: vxc1(ngrid1,nv1,2), fxc1(ngrid1,nv1,2,nv1,2)
    call xc_set_ref_kernel(ngrid1, nv1, vxc1, fxc1)
  end subroutine mrsf_xc_set_ref_kernel_c

  subroutine mrsf_xc_g_density_c(rho1) bind(c, name='mrsf_xc_g_density')
    real(dp), intent(out) :: rho1(ngrid,ndc)
    call xc_g_density(rho1)
  end subroutine mrsf_xc_g_density_c

  subroutine mrsf_xc_g_vxc_ao_c(V1) bind(c, name='mrsf_xc_g_vxc_ao')
    real(dp), intent(out) :: V1(nao_g,nao_g)
    call xc_g_vxc_ao(V1)
  end subroutine mrsf_xc_g_vxc_ao_c

  subroutine mrsf_xc_cs_from_g_c(fac1) bind(c, name='mrsf_xc_cs_from_g')
    real(dp), intent(in) :: fac1
    call xc_cs_from_g(fac1)
  end subroutine mrsf_xc_cs_from_g_c

  subroutine mrsf_xc_probe_begin_q_c(nst1, nexp1, Bexp1, maxnb1, nb1, idx1, S1, withg1, withk1) &
       bind(c, name='mrsf_xc_probe_begin_q')
    integer(is), intent(in) :: nst1, nexp1, maxnb1, nb1(4,nst1), idx1(maxnb1,4,nst1)
    real(dp), intent(in)    :: Bexp1(nao_g,max(nexp1,1)), S1(maxnb1,maxnb1,4,nst1)
    logical(c_bool), intent(in) :: withg1, withk1
    call xc_probe_begin_q(nst1, nexp1, Bexp1, maxnb1, nb1, idx1, S1, logical(withg1), logical(withk1))
  end subroutine mrsf_xc_probe_begin_q_c

  subroutine mrsf_xc_probe_super_c(ib01, nsub1, ng1, nc21, ao1) bind(c, name='mrsf_xc_probe_super')
    integer(is), intent(in) :: ib01, nsub1, ng1, nc21
    real(dp), intent(in)    :: ao1(*)
    call xc_probe_super(ib01, nsub1, ng1, nc21, ao1)
  end subroutine mrsf_xc_probe_super_c

  subroutine mrsf_xc_g_potential_q_c(nexp1, Bexp1, maxnb1, nbM1, idxM1, SM1, nbS1, idxS1, SS1, wantk1, &
       MM1, MK1, VE1, VH1) bind(c, name='mrsf_xc_g_potential_q')
    integer(is), intent(in) :: nexp1, maxnb1, nbM1, nbS1, idxM1(maxnb1), idxS1(maxnb1)
    real(dp), intent(in)    :: Bexp1(nao_g,max(nexp1,1)), SM1(maxnb1,maxnb1), SS1(maxnb1,maxnb1)
    logical(c_bool), intent(in) :: wantk1
    real(dp), intent(out)   :: MM1(nocca,nao_g), MK1(nocca,nao_g), VE1(nao_g,max(nexp1,1)), VH1(nao_g,nocca)
    call xc_g_potential_q(nexp1, Bexp1, maxnb1, nbM1, idxM1, SM1, nbS1, idxS1, SS1, logical(wantk1), &
         MM1, MK1, VE1, VH1)
  end subroutine mrsf_xc_g_potential_q_c

  subroutine mrsf_xc_timings_c() bind(c, name='mrsf_xc_timings')
    call xc_timings()
  end subroutine mrsf_xc_timings_c

  ! ---- Part XX: base-density weights, Dgamma, the extended potential and probe
  subroutine mrsf_xc_set_hole_weights_c(nocca1, lam1) bind(c, name='mrsf_xc_set_hole_weights')
    integer(is), intent(in) :: nocca1
    real(dp), intent(in)    :: lam1(nocca1)
    call xc_set_hole_weights(nocca1, lam1)
  end subroutine mrsf_xc_set_hole_weights_c

  subroutine mrsf_xc_clear_hole_weights_c() bind(c, name='mrsf_xc_clear_hole_weights')
    call xc_clear_hole_weights()
  end subroutine mrsf_xc_clear_hole_weights_c

  subroutine mrsf_xc_set_dgamma_c(dg1) bind(c, name='mrsf_xc_set_dgamma')
    real(dp), intent(in) :: dg1(2,2)
    call xc_set_dgamma(dg1)
  end subroutine mrsf_xc_set_dgamma_c

  subroutine mrsf_xc_clear_dgamma_c() bind(c, name='mrsf_xc_clear_dgamma')
    call xc_clear_dgamma()
  end subroutine mrsf_xc_clear_dgamma_c

  subroutine mrsf_xc_probe_begin_q2_c(nst1, nexp1, Bexp1, maxnb1, nb1, idx1, S1, withg1, withk1, withw1) &
       bind(c, name='mrsf_xc_probe_begin_q2')
    integer(is), intent(in) :: nst1, nexp1, maxnb1, nb1(5,nst1), idx1(maxnb1,5,nst1)
    real(dp), intent(in)    :: Bexp1(nao_g,max(nexp1,1)), S1(maxnb1,maxnb1,5,nst1)
    logical(c_bool), intent(in) :: withg1, withk1, withw1
    call xc_probe_begin_q2(nst1, nexp1, Bexp1, maxnb1, nb1, idx1, S1, logical(withg1), logical(withk1), &
         logical(withw1))
  end subroutine mrsf_xc_probe_begin_q2_c

  subroutine mrsf_xc_g_potential_q2_c(nexp1, Bexp1, maxnb1, nbM1, idxM1, SM1, nbS1, idxS1, SS1, nbT1, idxT1, ST1, &
       wantk1, wantt1, MM1, MK1, MT1, MK21, VE1, VH1) bind(c, name='mrsf_xc_g_potential_q2')
    integer(is), intent(in) :: nexp1, maxnb1, nbM1, nbS1, nbT1, idxM1(maxnb1), idxS1(maxnb1), idxT1(maxnb1)
    real(dp), intent(in)    :: Bexp1(nao_g,max(nexp1,1)), SM1(maxnb1,maxnb1), SS1(maxnb1,maxnb1), ST1(maxnb1,maxnb1)
    logical(c_bool), intent(in) :: wantk1, wantt1
    real(dp), intent(out)   :: MM1(nocca,nao_g), MK1(nocca,nao_g), MT1(nocca,nao_g), MK21(nocca,nao_g)
    real(dp), intent(out)   :: VE1(nao_g,max(nexp1,1)), VH1(nao_g,nocca)
    call xc_g_potential_q2(nexp1, Bexp1, maxnb1, nbM1, idxM1, SM1, nbS1, idxS1, SS1, nbT1, idxT1, ST1, &
         logical(wantk1), logical(wantt1), MM1, MK1, MT1, MK21, VE1, VH1)
  end subroutine mrsf_xc_g_potential_q2_c

end module mrsf_xcgrid_interface
