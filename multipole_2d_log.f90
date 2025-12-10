module multipole_2d_log
  !============================================================================
  ! Fast Multipole Method for 2D logarithmic potential
  !
  ! Potential: φ = -q·log|r| (vortex interaction energy)
  ! Force:     F = -∇φ = q·r/r²
  !
  ! Reference: Greengard & Rokhlin (1987), Beatson & Greengard FMM course
  !
  ! For 2D logarithmic kernel:
  !   φ(z) = -log|z-z₀| = -log|z| - Re[Σ_{k=1}^∞ (z₀/z)^k / k]  for |z₀| < |z|
  !
  ! Multipole moments (about center z_c):
  !   M_0 = Σ q_i                    (monopole - total charge)
  !   M_k = Σ q_i · z_i^k / k        (higher moments, k ≥ 1)
  !   where z_i = (particle_i - center) as complex number
  !
  ! Local expansion coefficients:
  !   L_k for Taylor series φ(z) = Σ L_k · z^k
  !============================================================================
  implicit none
  private
  public :: multipole_system
  public :: init_system, destroy_system
  public :: add_particle, compute_all_forces
  public :: compute_energy_direct, compute_energy_multipole
  public :: move_particles, write_frame, compute_forces_direct
  public :: p2m_all, m2l_all, l2p_all, p2p_all  ! For debugging

  ! Double precision kind
  integer, parameter :: dp = selected_real_kind(15, 307)

  ! Mathematical constants
  real(dp), parameter :: PI = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: EPS = 1.0e-14_dp

  ! ============================================================================
  ! DATA STRUCTURES
  ! ============================================================================

  type :: particle
    real(dp) :: q           ! Charge (vortex strength)
    real(dp) :: x, y        ! Position
    real(dp) :: vx, vy      ! Velocity
    real(dp) :: fx, fy      ! Force
    real(dp) :: phi         ! Potential
    integer  :: cell_i, cell_j  ! Cell indices
  end type

  type :: cell
    integer :: n_particles, capacity
    integer, allocatable :: particle_ids(:)
    real(dp) :: center_x, center_y
    ! Complex moments and coefficients for log kernel
    ! M(0) = monopole (total charge)
    ! M(k) = Σ q_i · (z_i/rscale)^k / k  for k ≥ 1 (normalized)
    complex(dp), allocatable :: M(:)  ! Multipole moments
    complex(dp), allocatable :: L(:)  ! Local expansion coefficients
  end type

  type :: multipole_system
    real(dp) :: xmax, ymax          ! Domain half-widths
    integer  :: nx_cells, ny_cells  ! Number of cells
    real(dp) :: dx_cell, dy_cell    ! Cell dimensions
    integer  :: p_max               ! Expansion order
    real(dp) :: rscale              ! Global length scale for normalization

    type(cell), allocatable :: cells(:,:)
    type(particle), allocatable :: particles(:)
    integer :: n_particles, max_particles

    ! Binomial coefficients for translations
    real(dp), allocatable :: binomial(:,:)
  end type

contains

  ! ============================================================================
  ! INITIALIZATION
  ! ============================================================================

  subroutine init_system(sys, xmax, ymax, nx, ny, p_max, max_part)
    type(multipole_system), intent(out) :: sys
    real(dp), intent(in) :: xmax, ymax
    integer, intent(in) :: nx, ny, p_max, max_part
    integer :: i, j, n, k
    integer :: cell_capacity

    sys%xmax = xmax
    sys%ymax = ymax
    sys%nx_cells = nx
    sys%ny_cells = ny
    sys%dx_cell = 2.0_dp * xmax / nx
    sys%dy_cell = 2.0_dp * ymax / ny
    sys%p_max = p_max
    sys%max_particles = max_part
    sys%n_particles = 0

    ! Global rscale: use cell size as characteristic length
    sys%rscale = max(sys%dx_cell, sys%dy_cell)

    ! Allocate particle array
    allocate(sys%particles(max_part))

    ! Calculate cell capacity based on particle count and grid size
    cell_capacity = calculate_cell_capacity(max_part, nx*ny)

    print '(A,I0)', ' Calculated cell capacity: ', cell_capacity
    print '(A,I0)', ' Average particles/cell:   ', max_part/(nx*ny)

    ! Initialize cells
    allocate(sys%cells(nx, ny))
    do j = 1, ny
      do i = 1, nx
        sys%cells(i,j)%center_x = -xmax + (i - 0.5_dp) * sys%dx_cell
        sys%cells(i,j)%center_y = -ymax + (j - 0.5_dp) * sys%dy_cell
        sys%cells(i,j)%n_particles = 0
        sys%cells(i,j)%capacity = cell_capacity
        allocate(sys%cells(i,j)%particle_ids(cell_capacity))
        allocate(sys%cells(i,j)%M(0:p_max))
        allocate(sys%cells(i,j)%L(0:p_max))
        sys%cells(i,j)%M = cmplx(0.0_dp, 0.0_dp, dp)
        sys%cells(i,j)%L = cmplx(0.0_dp, 0.0_dp, dp)
      end do
    end do

    ! Precompute binomial coefficients
    allocate(sys%binomial(0:2*p_max+2, 0:2*p_max+2))
    sys%binomial = 0.0_dp
    sys%binomial(0,0) = 1.0_dp
    do n = 1, 2*p_max+2
      sys%binomial(n,0) = 1.0_dp
      sys%binomial(n,n) = 1.0_dp
      do k = 1, n-1
        sys%binomial(n,k) = sys%binomial(n-1,k-1) + sys%binomial(n-1,k)
      end do
    end do
  end subroutine init_system

  ! ============================================================================

  function calculate_cell_capacity(n_particles, n_cells) result(capacity)
    integer, intent(in) :: n_particles, n_cells
    integer :: capacity

    ! Simple solution: for small grids, allocate generously
    ! For large grids, use average * 5 to handle clustering
    if (n_cells <= 100) then
      ! Small grid: be generous, allow up to half of all particles per cell
      capacity = n_particles / 2 + 100
    else
      ! Large grid: assume better distribution, use 5x average
      capacity = (n_particles / n_cells) * 5 + 100
    end if

    ! Ensure minimum
    capacity = max(capacity, 200)
  end function calculate_cell_capacity

  ! ============================================================================

  subroutine destroy_system(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i, j

    if (allocated(sys%particles)) deallocate(sys%particles)
    if (allocated(sys%binomial)) deallocate(sys%binomial)

    if (allocated(sys%cells)) then
      do j = 1, sys%ny_cells
        do i = 1, sys%nx_cells
          if (allocated(sys%cells(i,j)%particle_ids)) &
            deallocate(sys%cells(i,j)%particle_ids)
          if (allocated(sys%cells(i,j)%M)) deallocate(sys%cells(i,j)%M)
          if (allocated(sys%cells(i,j)%L)) deallocate(sys%cells(i,j)%L)
        end do
      end do
      deallocate(sys%cells)
    end if
  end subroutine destroy_system

  ! ============================================================================

  subroutine add_particle(sys, x, y, q)
    type(multipole_system), intent(inout) :: sys
    real(dp), intent(in) :: x, y, q
    integer :: ci, cj, ip

    sys%n_particles = sys%n_particles + 1
    ip = sys%n_particles

    sys%particles(ip)%x = x
    sys%particles(ip)%y = y
    sys%particles(ip)%q = q
    sys%particles(ip)%vx = 0.0_dp
    sys%particles(ip)%vy = 0.0_dp
    sys%particles(ip)%fx = 0.0_dp
    sys%particles(ip)%fy = 0.0_dp
    sys%particles(ip)%phi = 0.0_dp

    ! Find cell (with bounds checking)
    ci = floor((x + sys%xmax) / sys%dx_cell) + 1
    cj = floor((y + sys%ymax) / sys%dy_cell) + 1
    ci = max(1, min(sys%nx_cells, ci))
    cj = max(1, min(sys%ny_cells, cj))

    sys%particles(ip)%cell_i = ci
    sys%particles(ip)%cell_j = cj

    ! Add to cell with capacity check
    if (sys%cells(ci,cj)%n_particles >= sys%cells(ci,cj)%capacity) then
      print '(A)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print '(A)', 'ERROR: Cell capacity exceeded in add_particle!'
      print '(A,I0,A,I0)', '  Cell (', ci, ',', cj, ')'
      print '(A,I0)', '  Current particles in cell: ', sys%cells(ci,cj)%n_particles
      print '(A,I0)', '  Cell capacity:             ', sys%cells(ci,cj)%capacity
      print '(A)', '  Solutions:'
      print '(A)', '    1. Increase grid size (nx, ny)'
      print '(A)', '    2. Reduce dt if particles move too far'
      print '(A)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      stop 1
    end if
    sys%cells(ci,cj)%n_particles = sys%cells(ci,cj)%n_particles + 1
    sys%cells(ci,cj)%particle_ids(sys%cells(ci,cj)%n_particles) = ip
  end subroutine add_particle

  ! ============================================================================
  ! P2M: PARTICLE TO MULTIPOLE
  ! ============================================================================
  ! Compute multipole moments for all cells from particle positions
  !
  ! For log kernel in 2D with rscale normalization:
  !   φ(z) = -log|z-z_c| = -log|z| - Re[Σ_{k=1}^∞ z_c^k / (k·z^k)]
  !
  ! Multipole expansion about z_c with rscale:
  !   M_0 = Σ q_i  (total charge)
  !   M_k = -Σ q_i · (z_i/rscale)^k / k   for k ≥ 1
  !
  ! The rscale normalization ensures numerical stability by keeping
  ! intermediate values O(1).
  ! ============================================================================

  subroutine p2m_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, i, ip, k
    real(dp) :: dx, dy, q
    complex(dp) :: z, z_scaled, z_pow

    ! Reset all moments
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        sys%cells(ci,cj)%M = cmplx(0.0_dp, 0.0_dp, dp)
      end do
    end do

    ! Compute moments for each cell
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        do i = 1, sys%cells(ci,cj)%n_particles
          ip = sys%cells(ci,cj)%particle_ids(i)
          q = sys%particles(ip)%q

          ! Complex position relative to cell center
          dx = sys%particles(ip)%x - sys%cells(ci,cj)%center_x
          dy = sys%particles(ip)%y - sys%cells(ci,cj)%center_y
          z = cmplx(dx, dy, dp)

          ! Normalize by global rscale
          z_scaled = z / sys%rscale

          ! M_0 = Σ q_i (monopole - not scaled)
          sys%cells(ci,cj)%M(0) = sys%cells(ci,cj)%M(0) + q

          ! M_k = -Σ q_i · (z_i/rscale)^k / k  (k ≥ 1)
          ! CRITICAL: Minus sign as in fmm2d library (l2dformmpc)
          ! Formula: mpole_n = -sum charge * (1/n) * (z0/rscale)^n
          z_pow = z_scaled
          do k = 1, sys%p_max
            sys%cells(ci,cj)%M(k) = sys%cells(ci,cj)%M(k) - q * z_pow / real(k, dp)
            z_pow = z_pow * z_scaled
          end do
        end do
      end do
    end do
  end subroutine p2m_all

  ! ============================================================================
  ! M2L: MULTIPOLE TO LOCAL TRANSLATION
  ! ============================================================================
  ! Convert multipole expansion from source cell to local expansion at target
  !
  ! For 2D log kernel, the key formula is:
  !   -log|z-z₀| = -log|z₀| - Re[Σ_{k=1}^∞ z^k / (k·z₀^k)]
  !
  ! This gives M2L translation from source center z_s to target center z_t:
  !   z₀ = z_t - z_s
  !
  ! Monopole (M_0) contribution:
  !   L_0 += M_0 · (-log|z₀|)
  !   L_j += M_0 · (-1/j) / z₀^j    for j ≥ 1
  !
  ! Higher moments (M_k, k ≥ 1) contribution:
  !   L_0 += Σ_k M_k · conj(1/z₀^k)
  !   L_j += Σ_k M_k · conj(1/z₀^k) / z₀^j    for j ≥ 1
  ! ============================================================================

  subroutine m2l_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, si, sj

    ! Reset local expansions
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        sys%cells(ci,cj)%L = cmplx(0.0_dp, 0.0_dp, dp)
      end do
    end do

    ! M2L for all well-separated cell pairs
    ! Use 2-cell separation (skip if within 2 cells)
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        do sj = 1, sys%ny_cells
          do si = 1, sys%nx_cells
            ! Skip near neighbors (P2P region)
            if (abs(si-ci) <= 2 .and. abs(sj-cj) <= 2) cycle
            call m2l_single(sys, si, sj, ci, cj)
          end do
        end do
      end do
    end do
  end subroutine m2l_all

  ! ============================================================================

  subroutine m2l_single(sys, si, sj, ti, tj)
    type(multipole_system), intent(inout) :: sys
    integer, intent(in) :: si, sj, ti, tj
    integer :: k, j, kk
    real(dp) :: dx, dy, rtmp, binom_coef
    complex(dp) :: z0, z0inv, z0pow1, z0pow2
    complex(dp) :: temp_M(0:sys%p_max), temp_L(0:sys%p_max)
    complex(dp) :: ztemp1, ztemp2, ztemp3

    ! CORRECT IMPLEMENTATION from Flatiron Institute fmm2d library
    ! File: laprouts2d.f, subroutine l2dmploc (line 1130)
    !
    ! Key differences from previous version:
    ! 1. z0 = -(center2 - center1) where center2=target, center1=source
    ! 2. Separate z0pow1 and z0pow2 with different normalizations
    ! 3. Pre-multiply M by z0pow1, then use for L computation

    ! Vector from source to target
    dx = sys%cells(ti,tj)%center_x - sys%cells(si,sj)%center_x
    dy = sys%cells(ti,tj)%center_y - sys%cells(si,sj)%center_y
    z0 = -cmplx(dx, dy, dp)  ! NEGATIVE! (as in fmm2d)

    if (abs(z0) < EPS) return

    z0inv = 1.0_dp / z0

    ! Precompute power sequences (as in fmm2d)
    ztemp1 = z0inv
    ztemp2 = z0inv * sys%rscale
    ztemp3 = -z0inv * sys%rscale

    ! Transform multipole: temp_M(k) = M(k) * z0pow1(k)
    ! z0pow1(k) = (-1/z0 * rscale)^k
    temp_M(0) = sys%cells(si,sj)%M(0)
    z0pow1 = ztemp3
    do k = 1, sys%p_max
      temp_M(k) = sys%cells(si,sj)%M(k) * z0pow1
      z0pow1 = -z0pow1 * ztemp1 * sys%rscale
    end do

    ! Compute L_0: M_0*log|z0| + Σ_{k=1}^p temp_M(k)
    rtmp = log(abs(z0))
    temp_L(0) = temp_M(0) * rtmp

    do k = 1, sys%p_max
      temp_L(0) = temp_L(0) + temp_M(k)
    end do

    ! Compute L_j (j ≥ 1)
    do j = 1, sys%p_max
      ! Start with monopole term
      temp_L(j) = -temp_M(0) / real(j, dp)

      ! Add contributions from higher moments
      do k = 1, sys%p_max
        if (j + k - 1 <= ubound(sys%binomial, 1) .and. &
            k - 1 <= ubound(sys%binomial, 2)) then
          binom_coef = sys%binomial(j+k-1, k-1)
          temp_L(j) = temp_L(j) + temp_M(k) * binom_coef
        end if
      end do

      ! Multiply by z0pow2(j) = (1/z0 * rscale)^j
      z0pow2 = ztemp2
      do kk = 2, j
        z0pow2 = z0pow2 * ztemp1 * sys%rscale
      end do
      temp_L(j) = temp_L(j) * z0pow2
    end do

    ! INCREMENT local expansion (INCREMENTS, not replaces!)
    do j = 0, sys%p_max
      sys%cells(ti,tj)%L(j) = sys%cells(ti,tj)%L(j) + temp_L(j)
    end do
  end subroutine m2l_single

  ! ============================================================================
  ! L2P: LOCAL TO PARTICLE
  ! ============================================================================
  ! Evaluate local expansion at particle positions with rscale normalization
  !
  ! Local expansion: φ(z) = Σ_{k=0}^p L_k · (z/rscale)^k
  ! Gradient: ∇φ = (1/rscale) · Σ_{k=1}^p k · L_k · (z/rscale)^(k-1) · d(z)/dz̄
  !
  ! In complex form with rscale:
  !   φ(z) = L_0 + Σ_{k=1}^p L_k · (z/rscale)^k
  !   dφ/dz = (1/rscale) · Σ_{k=1}^p k · L_k · (z/rscale)^(k-1)
  !   Force = -q · (Re(dφ/dz), Im(dφ/dz))
  !
  ! CRITICAL: The gradient picks up a factor of 1/rscale from chain rule!
  ! ============================================================================

  subroutine l2p_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, i, ip, k
    real(dp) :: dx, dy
    complex(dp) :: z, z_scaled, z_pow, phi_c, grad_phi

    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        do i = 1, sys%cells(ci,cj)%n_particles
          ip = sys%cells(ci,cj)%particle_ids(i)

          ! Position relative to cell center
          dx = sys%particles(ip)%x - sys%cells(ci,cj)%center_x
          dy = sys%particles(ip)%y - sys%cells(ci,cj)%center_y
          z = cmplx(dx, dy, dp)

          ! Normalize by GLOBAL rscale
          z_scaled = z / sys%rscale

          ! Evaluate potential: φ = L_0 + Σ L_k · (z/rscale)^k
          phi_c = sys%cells(ci,cj)%L(0)
          z_pow = z_scaled
          do k = 1, sys%p_max
            phi_c = phi_c + sys%cells(ci,cj)%L(k) * z_pow
            z_pow = z_pow * z_scaled
          end do
          sys%particles(ip)%phi = sys%particles(ip)%phi + real(phi_c, dp)

          ! Evaluate gradient: dφ/dz = (1/rscale) · Σ k·L_k·(z/rscale)^(k-1)
          ! Start with L_1 term (k=1)
          grad_phi = sys%cells(ci,cj)%L(1)
          z_pow = z_scaled  ! Start with (z/rscale)^1 for k=2 term
          do k = 2, sys%p_max
            grad_phi = grad_phi + real(k, dp) * sys%cells(ci,cj)%L(k) * z_pow
            z_pow = z_pow * z_scaled
          end do

          ! CRITICAL: Multiply by 1/rscale due to chain rule d/dz[(z/rscale)^k]
          grad_phi = grad_phi / sys%rscale

          ! Force = -q · ∇φ
          ! For complex gradient dφ/dz = ∂φ/∂x - i·∂φ/∂y
          ! We have: Fx = -q·∂φ/∂x = -q·Re(dφ/dz)
          !          Fy = -q·∂φ/∂y = +q·Im(dφ/dz)  <- PLUS sign!
          sys%particles(ip)%fx = sys%particles(ip)%fx - &
            sys%particles(ip)%q * real(grad_phi, dp)
          sys%particles(ip)%fy = sys%particles(ip)%fy + &
            sys%particles(ip)%q * aimag(grad_phi)
        end do
      end do
    end do
  end subroutine l2p_all

  ! ============================================================================
  ! P2P: DIRECT PARTICLE-PARTICLE INTERACTION
  ! ============================================================================
  ! Compute forces directly for nearby cells (within 2-cell stencil)
  !
  ! For log potential: φ = -q·log(r)
  ! Force: F = -∇φ = q·r/r²
  ! ============================================================================

  subroutine p2p_all(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: ci, cj, si, sj

    ! P2P for near neighbors (5×5 stencil: within 2 cells)
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        do sj = max(1,cj-2), min(sys%ny_cells,cj+2)
          do si = max(1,ci-2), min(sys%nx_cells,ci+2)
            call p2p_cells(sys, ci, cj, si, sj)
          end do
        end do
      end do
    end do
  end subroutine p2p_all

  ! ============================================================================

  subroutine p2p_cells(sys, ci, cj, si, sj)
    type(multipole_system), intent(inout) :: sys
    integer, intent(in) :: ci, cj, si, sj
    integer :: i, j, ip, jp
    real(dp) :: dx, dy, r, r2, qi, qj, fx, fy, log_r
    logical :: same_cell

    same_cell = (ci == si .and. cj == sj)

    do i = 1, sys%cells(ci,cj)%n_particles
      ip = sys%cells(ci,cj)%particle_ids(i)
      qi = sys%particles(ip)%q

      do j = 1, sys%cells(si,sj)%n_particles
        jp = sys%cells(si,sj)%particle_ids(j)

        ! Skip self-interaction and avoid double-counting
        if (same_cell .and. ip >= jp) cycle

        qj = sys%particles(jp)%q

        dx = sys%particles(jp)%x - sys%particles(ip)%x
        dy = sys%particles(jp)%y - sys%particles(ip)%y
        r2 = dx*dx + dy*dy
        if (r2 < EPS*EPS) cycle

        r = sqrt(r2)
        log_r = log(r)

        ! Force: F = q_i·q_j · r/r²
        fx = qi * qj * dx / r2
        fy = qi * qj * dy / r2

        sys%particles(ip)%fx = sys%particles(ip)%fx + fx
        sys%particles(ip)%fy = sys%particles(ip)%fy + fy
        sys%particles(ip)%phi = sys%particles(ip)%phi - qj * log_r

        if (same_cell) then
          sys%particles(jp)%fx = sys%particles(jp)%fx - fx
          sys%particles(jp)%fy = sys%particles(jp)%fy - fy
          sys%particles(jp)%phi = sys%particles(jp)%phi - qi * log_r
        end if
      end do
    end do
  end subroutine p2p_cells

  ! ============================================================================
  ! MAIN FMM FORCE COMPUTATION
  ! ============================================================================

  subroutine compute_all_forces(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i

    ! Zero out forces and potentials
    do i = 1, sys%n_particles
      sys%particles(i)%fx = 0.0_dp
      sys%particles(i)%fy = 0.0_dp
      sys%particles(i)%phi = 0.0_dp
    end do

    ! FMM algorithm
    call p2m_all(sys)    ! Particle to Multipole
    call m2l_all(sys)    ! Multipole to Local (far-field)
    call l2p_all(sys)    ! Local to Particle
    call p2p_all(sys)    ! Direct for near-field
  end subroutine compute_all_forces

  ! ============================================================================
  ! DIRECT COMPUTATION (for validation)
  ! ============================================================================

  subroutine compute_forces_direct(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i, j
    real(dp) :: dx, dy, r, r2, qi, qj, fx, fy, log_r

    ! Zero out forces
    do i = 1, sys%n_particles
      sys%particles(i)%fx = 0.0_dp
      sys%particles(i)%fy = 0.0_dp
      sys%particles(i)%phi = 0.0_dp
    end do

    ! O(N²) direct summation
    do i = 1, sys%n_particles
      qi = sys%particles(i)%q
      do j = i+1, sys%n_particles
        qj = sys%particles(j)%q

        dx = sys%particles(j)%x - sys%particles(i)%x
        dy = sys%particles(j)%y - sys%particles(i)%y
        r2 = dx*dx + dy*dy
        if (r2 < EPS*EPS) cycle

        r = sqrt(r2)
        log_r = log(r)

        ! Force for log potential
        fx = qi * qj * dx / r2
        fy = qi * qj * dy / r2

        sys%particles(i)%fx = sys%particles(i)%fx + fx
        sys%particles(i)%fy = sys%particles(i)%fy + fy
        sys%particles(i)%phi = sys%particles(i)%phi - qj * log_r

        sys%particles(j)%fx = sys%particles(j)%fx - fx
        sys%particles(j)%fy = sys%particles(j)%fy - fy
        sys%particles(j)%phi = sys%particles(j)%phi - qi * log_r
      end do
    end do
  end subroutine compute_forces_direct

  ! ============================================================================
  ! ENERGY COMPUTATION
  ! ============================================================================

  function compute_energy_direct(sys) result(energy)
    type(multipole_system), intent(in) :: sys
    real(dp) :: energy
    integer :: i, j
    real(dp) :: dx, dy, r, log_r

    energy = 0.0_dp
    do i = 1, sys%n_particles
      do j = i+1, sys%n_particles
        dx = sys%particles(j)%x - sys%particles(i)%x
        dy = sys%particles(j)%y - sys%particles(i)%y
        r = sqrt(dx*dx + dy*dy)
        if (r > EPS) then
          log_r = log(r)
          energy = energy - sys%particles(i)%q * sys%particles(j)%q * log_r
        end if
      end do
    end do
  end function compute_energy_direct

  ! ============================================================================

  function compute_energy_multipole(sys) result(energy)
    type(multipole_system), intent(in) :: sys
    real(dp) :: energy
    integer :: i

    energy = 0.0_dp
    do i = 1, sys%n_particles
      energy = energy + 0.5_dp * sys%particles(i)%q * sys%particles(i)%phi
    end do
  end function compute_energy_multipole

  ! ============================================================================
  ! TIME INTEGRATION
  ! ============================================================================

  subroutine move_particles(sys, dt, eta)
    type(multipole_system), intent(inout) :: sys
    real(dp), intent(in) :: dt, eta
    integer :: i
    real(dp) :: x_new, y_new, vx_new, vy_new

    ! Simple Euler method: v(t+dt) = v(t) + F*dt, x(t+dt) = x(t) + v*dt
    do i = 1, sys%n_particles
      ! Update velocity: v_new = v + F*dt
      vx_new = sys%particles(i)%vx + sys%particles(i)%fx * dt
      vy_new = sys%particles(i)%vy + sys%particles(i)%fy * dt

      ! Update position: x_new = x + v_new*dt
      x_new = sys%particles(i)%x + vx_new * dt
      y_new = sys%particles(i)%y + vy_new * dt

      ! Reflective boundary conditions with velocity reversal
      if (x_new > sys%xmax) then
        x_new = 2.0_dp * sys%xmax - x_new
        vx_new = -vx_new
      else if (x_new < -sys%xmax) then
        x_new = -2.0_dp * sys%xmax - x_new
        vx_new = -vx_new
      end if

      if (y_new > sys%ymax) then
        y_new = 2.0_dp * sys%ymax - y_new
        vy_new = -vy_new
      else if (y_new < -sys%ymax) then
        y_new = -2.0_dp * sys%ymax - y_new
        vy_new = -vy_new
      end if

      ! Update particle state
      sys%particles(i)%x = x_new
      sys%particles(i)%y = y_new
      sys%particles(i)%vx = vx_new
      sys%particles(i)%vy = vy_new
    end do

    ! Update cell assignments
    call update_cell_assignments(sys)
  end subroutine move_particles

  ! ============================================================================

  subroutine update_cell_assignments(sys)
    type(multipole_system), intent(inout) :: sys
    integer :: i, ci, cj
    integer :: overflow_count
    logical :: capacity_exceeded

    ! Clear all cells
    do cj = 1, sys%ny_cells
      do ci = 1, sys%nx_cells
        sys%cells(ci,cj)%n_particles = 0
      end do
    end do

    overflow_count = 0
    capacity_exceeded = .false.

    ! Reassign particles
    do i = 1, sys%n_particles
      ci = floor((sys%particles(i)%x + sys%xmax) / sys%dx_cell) + 1
      cj = floor((sys%particles(i)%y + sys%ymax) / sys%dy_cell) + 1
      ci = max(1, min(sys%nx_cells, ci))
      cj = max(1, min(sys%ny_cells, cj))

      sys%particles(i)%cell_i = ci
      sys%particles(i)%cell_j = cj

      sys%cells(ci,cj)%n_particles = sys%cells(ci,cj)%n_particles + 1
      if (sys%cells(ci,cj)%n_particles <= sys%cells(ci,cj)%capacity) then
        sys%cells(ci,cj)%particle_ids(sys%cells(ci,cj)%n_particles) = i
      else
        overflow_count = overflow_count + 1
        capacity_exceeded = .true.
      end if
    end do

    ! Report overflow error
    if (capacity_exceeded) then
      print '(A)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print '(A)', 'ERROR: Cell capacity exceeded during update!'
      print '(A,I0)', '  Particles lost:   ', overflow_count
      print '(A,I0)', '  Total particles:  ', sys%n_particles
      print '(A)', ''
      print '(A)', '  This usually happens when:'
      print '(A)', '    - Particles cluster in one region'
      print '(A)', '    - Time step is too large'
      print '(A)', '    - Grid is too coarse'
      print '(A)', ''
      print '(A)', '  Solutions:'
      print '(A)', '    1. Reduce time step (dt)'
      print '(A)', '    2. Use finer grid (increase nx, ny)'
      print '(A)', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      stop 1
    end if
  end subroutine update_cell_assignments

  ! ============================================================================
  ! OUTPUT
  ! ============================================================================

  subroutine write_frame(sys, filename)
    type(multipole_system), intent(in) :: sys
    character(len=*), intent(in) :: filename
    integer :: i, unit, ios

    ! Use newunit for thread-safe file handling
    open(newunit=unit, file=trim(filename), status='replace', iostat=ios)
    if (ios /= 0) then
      print *, 'ERROR: Cannot open file ', trim(filename)
      print *, 'IO status:', ios
      return
    end if

    write(unit, '(A)', iostat=ios) 'x,y,charge'
    if (ios /= 0) then
      print *, 'ERROR: Cannot write header to ', trim(filename)
      close(unit)
      return
    end if

    do i = 1, sys%n_particles
      write(unit, '(ES22.15,",",ES22.15,",",F5.1)', iostat=ios) &
        sys%particles(i)%x, sys%particles(i)%y, sys%particles(i)%q
      if (ios /= 0) then
        print *, 'ERROR: Cannot write particle', i, ' to ', trim(filename)
        close(unit)
        return
      end if
    end do

    close(unit)
  end subroutine write_frame

end module multipole_2d_log
