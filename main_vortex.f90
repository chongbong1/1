program main
  !============================================================================
  ! 2D Vortex Dynamics Simulation using Fast Multipole Method
  !
  ! Simulates superconductor vortex interactions using logarithmic potential
  ! Output: CSV files with particle positions for visualization
  !============================================================================
  use multipole_2d_log
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  ! System
  type(multipole_system) :: sys
  real(dp) :: box_size
  integer :: n_cells, p_order
  integer :: n_particles, max_particles

  ! Simulation parameters
  integer :: n_steps, write_interval
  real(dp) :: dt, eta

  ! Validation and diagnostics
  real(dp) :: max_force_err, avg_force_err, force_err
  real(dp) :: force_magnitude
  real(dp) :: start_time, end_time, time_direct, time_fmm
  real(dp), allocatable :: fx_direct(:), fy_direct(:)
  real(dp), allocatable :: fx_fmm(:), fy_fmm(:)

  ! Loop variables
  integer :: i, step, n_per_side
  character(len=512) :: filename, output_dir
  real(dp) :: x, y, charge

  !----------------------------------------------------------------------------
  ! CONFIGURATION
  !----------------------------------------------------------------------------
  ! Domain: [-box_size, box_size] × [-box_size, box_size]
  box_size = 1.0_dp

  ! Grid parameters
  n_cells = 7               ! 7×7 uniform grid
  p_order = 20              ! Multipole expansion order (increased for accuracy)

  ! Particles - uniform distribution
  n_per_side = 10           ! 10×10 = 100 particles (moderate test)
  n_particles = n_per_side * n_per_side
  max_particles = n_particles

  ! Time integration
  n_steps = 500
  dt = 0.001_dp             ! Timestep
  eta = 0.1_dp              ! Noise amplitude
  write_interval = 50

  !----------------------------------------------------------------------------
  ! INITIALIZE
  !----------------------------------------------------------------------------
  print '(A)', '============================================'
  print '(A)', '2D Vortex Dynamics - Fast Multipole Method'
  print '(A)', '============================================'
  print '(A,F8.3)', ' Box size:        ', box_size
  print '(A,I0,A,I0)', ' Grid:            ', n_cells, ' × ', n_cells
  print '(A,I0)', ' Multipole order: ', p_order
  print '(A,I0,A,I0,A,I0,A)', ' Particles:       ', n_particles, &
    ' (', n_per_side, ' × ', n_per_side, ')'
  print '(A,I0)', ' Time steps:      ', n_steps
  print '(A)', '============================================'
  print *

  call init_system(sys, box_size, box_size, n_cells, n_cells, p_order, max_particles)

  ! Add particles in uniform grid with alternating charges
  do i = 1, n_per_side
    do step = 1, n_per_side
      ! Position (avoid exact boundaries)
      x = -box_size * 0.9_dp + 2.0_dp * box_size * 0.9_dp * (i - 1) / (n_per_side - 1)
      y = -box_size * 0.9_dp + 2.0_dp * box_size * 0.9_dp * (step - 1) / (n_per_side - 1)

      ! Checkerboard pattern: vortices (+1) and antivortices (-1)
      charge = real((-1)**(i + step), dp)

      call add_particle(sys, x, y, charge)
    end do
  end do

  print '(A,I0,A)', 'Added ', sys%n_particles, ' particles in uniform grid'
  print *

  !----------------------------------------------------------------------------
  ! ACCURACY VALIDATION: FMM vs Direct
  !----------------------------------------------------------------------------
  print '(A)', '============================================'
  print '(A)', 'Force Accuracy Validation'
  print '(A)', '============================================'

  allocate(fx_direct(n_particles), fy_direct(n_particles))
  allocate(fx_fmm(n_particles), fy_fmm(n_particles))

  ! Compute with FMM
  call cpu_time(start_time)
  call compute_all_forces(sys)
  call cpu_time(end_time)
  time_fmm = end_time - start_time

  do i = 1, n_particles
    fx_fmm(i) = sys%particles(i)%fx
    fy_fmm(i) = sys%particles(i)%fy
  end do

  print '(A,F10.6,A)', 'FMM time:    ', time_fmm, ' s'

  ! Compute with direct summation
  call cpu_time(start_time)
  call compute_forces_direct(sys)
  call cpu_time(end_time)
  time_direct = end_time - start_time

  do i = 1, n_particles
    fx_direct(i) = sys%particles(i)%fx
    fy_direct(i) = sys%particles(i)%fy
  end do

  print '(A,F10.6,A)', 'Direct time: ', time_direct, ' s'
  if (time_fmm > 1.0e-10_dp) then
    print '(A,F10.2,A)', 'Speedup:     ', time_direct / time_fmm, ' x'
  end if
  print *

  ! Print first 5 forces in detail
  print '(A)', '--------------------------------------------'
  print '(A)', 'First 5 Particle Forces:'
  print '(A)', '--------------------------------------------'
  print '(A)', '  i     FMM_x        FMM_y        Direct_x     Direct_y     Rel.Error'
  do i = 1, min(5, n_particles)
    force_magnitude = sqrt(fx_direct(i)**2 + fy_direct(i)**2)
    if (force_magnitude > 1.0e-14_dp) then
      force_err = sqrt((fx_fmm(i) - fx_direct(i))**2 + &
                      (fy_fmm(i) - fy_direct(i))**2) / force_magnitude
    else
      force_err = 0.0_dp
    end if
    print '(I3,5(ES13.5))', i, fx_fmm(i), fy_fmm(i), fx_direct(i), fy_direct(i), force_err
  end do
  print '(A)', '--------------------------------------------'
  print *

  ! Compute overall accuracy
  max_force_err = 0.0_dp
  avg_force_err = 0.0_dp
  do i = 1, n_particles
    force_magnitude = sqrt(fx_direct(i)**2 + fy_direct(i)**2)
    if (force_magnitude > 1.0e-14_dp) then
      force_err = sqrt((fx_fmm(i) - fx_direct(i))**2 + &
                      (fy_fmm(i) - fy_direct(i))**2)
      force_err = force_err / force_magnitude
      max_force_err = max(max_force_err, force_err)
      avg_force_err = avg_force_err + force_err
    end if
  end do
  avg_force_err = avg_force_err / n_particles

  print '(A)', 'Force Accuracy Summary:'
  print '(A,ES12.5)', '  Max relative error: ', max_force_err
  print '(A,ES12.5)', '  Avg relative error: ', avg_force_err
  print *

  ! Accuracy assessment
  if (max_force_err > 1.0e-3_dp) then
    print '(A)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print '(A)', 'WARNING: Force accuracy is LOW'
    print '(A)', 'Recommendations:'
    print '(A,I0)', '  - Increase p_order (current: ', p_order, ')'
    print '(A,I0,A,I0)', '  - Refine grid (current: ', n_cells, ' × ', n_cells, ')'
    print '(A)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print *
  else if (max_force_err > 1.0e-6_dp) then
    print '(A)', 'NOTE: Force accuracy is moderate (acceptable)'
    print *
  else
    print '(A)', 'EXCELLENT: High force accuracy achieved!'
    print *
  end if

  deallocate(fx_direct, fy_direct, fx_fmm, fy_fmm)

  !----------------------------------------------------------------------------
  ! OUTPUT DIRECTORY
  !----------------------------------------------------------------------------
  write(output_dir, '(A,I0,A,I0,A,I0)') &
    'output_N', n_particles, '_G', n_cells, '_P', p_order
  call execute_command_line('mkdir -p ' // trim(output_dir), wait=.true.)

  print '(A,A)', 'Output directory: ', trim(output_dir)
  print *

  !----------------------------------------------------------------------------
  ! TIME INTEGRATION
  !----------------------------------------------------------------------------
  print '(A)', '============================================'
  print '(A)', 'Starting Time Integration'
  print '(A)', '============================================'

  ! Restore forces for time integration
  call compute_all_forces(sys)

  ! Main simulation loop
  do step = 1, n_steps
    ! Move particles
    call move_particles(sys, dt, eta)

    ! Compute forces for next step
    call compute_all_forces(sys)

    ! Write output
    if (mod(step, write_interval) == 0) then
      write(filename, '(A,A,I6.6,A)') trim(output_dir), '/particles_', step, '.csv'
      call write_frame(sys, filename)
      print '(A,I6,A,I6)', 'Step ', step, ' / ', n_steps
    end if
  end do

  print *
  print '(A)', '============================================'
  print '(A)', 'Simulation Complete!'
  print '(A)', '============================================'
  print '(A,A)', 'Output: ', trim(output_dir)
  print '(A,I0)', 'N_particles = ', n_particles
  print '(A,I0,A,I0)', 'Grid        = ', n_cells, ' × ', n_cells
  print '(A,I0)', 'P_order     = ', p_order
  print '(A)', '============================================'
  print *

  ! Cleanup
  call destroy_system(sys)

  print '(A)', 'Done!'
end program main
