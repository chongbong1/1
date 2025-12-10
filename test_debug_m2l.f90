program test_debug_m2l
  use multipole_2d_log
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 307)

  type(multipole_system) :: sys
  integer :: i

  print '(A)', '========================================'
  print '(A)', 'DEBUG: M2L Translation'
  print '(A)', '========================================'

  ! Simple setup: 2 vortices in different cells
  call init_system(sys, 5.0_dp, 5.0_dp, 10, 10, 5, 10)

  ! Vortex 1 at (-3, 0) - will be in cell (2,5)
  call add_particle(sys, -3.0_dp, 0.0_dp, 1.0_dp)
  ! Vortex 2 at (+3, 0) - will be in cell (9,5)
  call add_particle(sys,  3.0_dp, 0.0_dp, 1.0_dp)

  print '(A)', 'Particles:'
  do i = 1, sys%n_particles
    print '(A,I0,A,2F8.3,A,I0,A,I0,A)', '  Part', i, ': (', &
      sys%particles(i)%x, sys%particles(i)%y, ') in cell (', &
      sys%particles(i)%cell_i, ',', sys%particles(i)%cell_j, ')'
  end do
  print *

  ! Compute multipole moments
  call p2m_all(sys)

  print '(A)', 'Multipole moments:'
  print '(A,I0,A,I0,A)', '  Cell (', sys%particles(1)%cell_i, ',', &
    sys%particles(1)%cell_j, '):'
  print '(A,2ES14.6)', '    M(0) = ', sys%cells(sys%particles(1)%cell_i, sys%particles(1)%cell_j)%M(0)
  print '(A,2ES14.6)', '    M(1) = ', sys%cells(sys%particles(1)%cell_i, sys%particles(1)%cell_j)%M(1)
  print '(A,2ES14.6)', '    M(2) = ', sys%cells(sys%particles(1)%cell_i, sys%particles(1)%cell_j)%M(2)

  print '(A,I0,A,I0,A)', '  Cell (', sys%particles(2)%cell_i, ',', &
    sys%particles(2)%cell_j, '):'
  print '(A,2ES14.6)', '    M(0) = ', sys%cells(sys%particles(2)%cell_i, sys%particles(2)%cell_j)%M(0)
  print '(A,2ES14.6)', '    M(1) = ', sys%cells(sys%particles(2)%cell_i, sys%particles(2)%cell_j)%M(1)
  print '(A,2ES14.6)', '    M(2) = ', sys%cells(sys%particles(2)%cell_i, sys%particles(2)%cell_j)%M(2)
  print *

  ! M2L translation
  call m2l_all(sys)

  print '(A)', 'Local expansion (particle 1 cell):'
  print '(A,2ES14.6)', '  L(0) = ', sys%cells(sys%particles(1)%cell_i, sys%particles(1)%cell_j)%L(0)
  print '(A,2ES14.6)', '  L(1) = ', sys%cells(sys%particles(1)%cell_i, sys%particles(1)%cell_j)%L(1)
  print '(A,2ES14.6)', '  L(2) = ', sys%cells(sys%particles(1)%cell_i, sys%particles(1)%cell_j)%L(2)
  print *

  ! L2P evaluation
  call l2p_all(sys)

  ! Also compute direct
  call compute_forces_direct(sys)

  print '(A)', 'Forces on particle 1:'
  print '(A,2ES14.6)', '  Direct: ', sys%particles(1)%fx, sys%particles(1)%fy

  ! Reset and compute with FMM
  sys%particles(1)%fx = 0.0_dp
  sys%particles(1)%fy = 0.0_dp
  sys%particles(2)%fx = 0.0_dp
  sys%particles(2)%fy = 0.0_dp

  ! Redo FMM
  call p2m_all(sys)
  call m2l_all(sys)
  call l2p_all(sys)

  print '(A,2ES14.6)', '  FMM:    ', sys%particles(1)%fx, sys%particles(1)%fy

  call destroy_system(sys)
end program test_debug_m2l
