program test_simple
  !============================================================================
  ! Simple test: 2 vortices far apart
  !============================================================================
  use multipole_2d_log
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  type(multipole_system) :: sys
  real(dp) :: fx_direct, fy_direct, fx_fmm, fy_fmm
  real(dp) :: err_fx, err_fy, err_rel

  print '(A)', '================================================'
  print '(A)', 'Simple Test: 2 vortices at (-3,0) and (+3,0)'
  print '(A)', '================================================'

  ! Initialize with large domain and fine grid
  call init_system(sys, 5.0_dp, 5.0_dp, 10, 10, 20, 10)

  ! Add two vortices far apart
  call add_particle(sys, -3.0_dp, 0.0_dp, 1.0_dp)   ! Vortex +1
  call add_particle(sys,  3.0_dp, 0.0_dp, 1.0_dp)   ! Vortex +1

  print '(A,I0)', 'Particles: ', sys%n_particles
  print '(A,I0,A,I0)', 'Grid: ', sys%nx_cells, ' x ', sys%ny_cells
  print '(A,I0)', 'p_order: ', sys%p_max
  print *

  ! Compute with direct
  call compute_forces_direct(sys)
  fx_direct = sys%particles(1)%fx
  fy_direct = sys%particles(1)%fy

  ! Compute with FMM
  call compute_all_forces(sys)
  fx_fmm = sys%particles(1)%fx
  fy_fmm = sys%particles(1)%fy

  print '(A)', 'Force on particle 1:'
  print '(A,2ES14.6)', '  Direct: ', fx_direct, fy_direct
  print '(A,2ES14.6)', '  FMM:    ', fx_fmm, fy_fmm

  err_fx = abs(fx_fmm - fx_direct)
  err_fy = abs(fy_fmm - fy_direct)
  err_rel = sqrt(err_fx**2 + err_fy**2) / sqrt(fx_direct**2 + fy_direct**2)

  print *
  print '(A,ES14.6)', 'Relative error: ', err_rel
  print *

  if (err_rel < 1.0e-8_dp) then
    print '(A)', '✓ TEST PASSED (excellent)'
  else if (err_rel < 1.0e-6_dp) then
    print '(A)', '✓ TEST PASSED (good)'
  else if (err_rel < 1.0e-3_dp) then
    print '(A)', '? TEST MARGINAL'
  else
    print '(A)', '✗ TEST FAILED'
  end if

  call destroy_system(sys)
end program test_simple
