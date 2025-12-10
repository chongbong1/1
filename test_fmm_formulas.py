#!/usr/bin/env python3
"""
Test script to validate FMM formulas for 2D logarithmic potential
Compares multipole expansion with direct calculation
"""
import numpy as np
from scipy.special import binom
import matplotlib.pyplot as plt

def log_potential(x, y, x0, y0, q):
    """Direct calculation: φ = -q·log(r)"""
    r = np.sqrt((x - x0)**2 + (y - y0)**2)
    if r < 1e-14:
        return 0.0
    return -q * np.log(r)

def log_force(x, y, x0, y0, q):
    """Direct calculation: F = q·r/r²"""
    dx = x - x0
    dy = y - y0
    r2 = dx**2 + dy**2
    if r2 < 1e-28:
        return 0.0, 0.0
    fx = q * dx / r2
    fy = q * dy / r2
    return fx, fy

class FMM2D_Log:
    """2D FMM for logarithmic kernel with CORRECTED formulas"""

    def __init__(self, p_max, rscale=1.0):
        self.p_max = p_max
        self.rscale = rscale

        # Precompute binomial coefficients
        self.binomial = np.zeros((2*p_max+3, 2*p_max+3))
        for n in range(2*p_max+3):
            for k in range(n+1):
                self.binomial[n, k] = binom(n, k)

    def p2m(self, particles, center):
        """Particle to Multipole (P2M)"""
        M = np.zeros(self.p_max + 1, dtype=complex)
        cx, cy = center

        for x, y, q in particles:
            # Complex position relative to center
            z = complex(x - cx, y - cy)
            z_scaled = z / self.rscale

            # M_0 = Σ q_i (monopole)
            M[0] += q

            # M_k = -Σ q_i · (z_i/rscale)^k / k  (k ≥ 1)
            z_pow = z_scaled
            for k in range(1, self.p_max + 1):
                M[k] += -q * z_pow / k
                z_pow *= z_scaled

        return M

    def m2l(self, M_source, z0):
        """Multipole to Local (M2L) - CORRECTED FORMULAS"""
        L = np.zeros(self.p_max + 1, dtype=complex)

        # Normalize translation vector
        z0_scaled = z0 / self.rscale

        if abs(z0_scaled) < 1e-14:
            return L

        z0_inv_scaled = 1.0 / z0_scaled

        # ====================================================================
        # CORRECTED L_0 formula:
        # L_0 = M_0 * (-log|z0|) + Σ_{k=1}^p M_k * z0^k
        # ====================================================================

        # Monopole contribution: NEGATIVE SIGN!
        L[0] = -M_source[0] * np.log(abs(z0_scaled))

        # Higher moments: with powers of z0^k
        z0_pow = z0_scaled
        for k in range(1, self.p_max + 1):
            L[0] += M_source[k] * z0_pow
            z0_pow *= z0_scaled

        # ====================================================================
        # CORRECTED L_j formula (j ≥ 1):
        # L_j = -M_0/(j·z0^j) + Σ_{k=1}^p M_k · C(j+k-1,k-1) / z0^(j+k)
        # ====================================================================

        z0_inv_pow_j = z0_inv_scaled  # (1/z0)^1
        for j in range(1, self.p_max + 1):
            # Monopole contribution: -M_0 / (j * z0^j)
            temp_sum = -M_source[0] * z0_inv_pow_j / j

            # Higher moment contributions: Σ M_k · C(j+k-1,k-1) / z0^(j+k)
            z0_inv_pow_jk = z0_inv_pow_j * z0_inv_scaled  # (1/z0)^(j+1)
            for k in range(1, self.p_max + 1):
                binom_coef = self.binomial[j+k-1, k-1]
                temp_sum += M_source[k] * binom_coef * z0_inv_pow_jk
                z0_inv_pow_jk *= z0_inv_scaled

            L[j] = temp_sum
            z0_inv_pow_j *= z0_inv_scaled

        return L

    def l2p(self, L, z):
        """Local to Particle (L2P)"""
        z_scaled = z / self.rscale

        # Evaluate potential: φ = L_0 + Σ L_k · (z/rscale)^k
        phi = L[0].real
        z_pow = z_scaled
        for k in range(1, self.p_max + 1):
            phi += (L[k] * z_pow).real
            z_pow *= z_scaled

        # Evaluate gradient: dφ/dz = (1/rscale) · Σ k·L_k·(z/rscale)^(k-1)
        grad_phi = L[1]
        z_pow = z_scaled
        for k in range(2, self.p_max + 1):
            grad_phi += k * L[k] * z_pow
            z_pow *= z_scaled

        grad_phi /= self.rscale

        # Force = -∇φ
        # Fx = -Re(dφ/dz), Fy = +Im(dφ/dz)
        fx = -grad_phi.real
        fy = grad_phi.imag

        return phi, fx, fy


def test_single_vortex():
    """Test 1: Single vortex - should give exact result"""
    print("=" * 60)
    print("TEST 1: Single Vortex")
    print("=" * 60)

    # Setup
    q = 1.0
    center = (0.0, 0.0)
    particles = [(0.0, 0.0, q)]

    # Test points
    test_points = [
        (1.0, 0.0),
        (0.0, 1.0),
        (1.0, 1.0),
        (2.0, 0.0),
    ]

    fmm = FMM2D_Log(p_max=20, rscale=1.0)
    M = fmm.p2m(particles, center)

    print(f"Multipole moments: M_0 = {M[0]:.6f}, M_1 = {M[1]:.6f}")
    print()

    max_err_phi = 0.0
    max_err_force = 0.0

    for x, y in test_points:
        # Direct calculation
        phi_direct = log_potential(x, y, 0.0, 0.0, q)
        fx_direct, fy_direct = log_force(x, y, 0.0, 0.0, q)

        # FMM calculation
        z0 = complex(x - center[0], y - center[1])
        L = fmm.m2l(M, z0)
        phi_fmm, fx_fmm, fy_fmm = fmm.l2p(L, complex(0, 0))

        # Errors
        err_phi = abs(phi_fmm - phi_direct)
        err_force = np.sqrt((fx_fmm - fx_direct)**2 + (fy_fmm - fy_direct)**2)
        force_mag = np.sqrt(fx_direct**2 + fy_direct**2)
        if force_mag > 1e-14:
            err_force_rel = err_force / force_mag
        else:
            err_force_rel = 0.0

        max_err_phi = max(max_err_phi, err_phi)
        max_err_force = max(max_err_force, err_force_rel)

        print(f"Point ({x:.1f}, {y:.1f}):")
        print(f"  Potential: Direct={phi_direct:.6f}, FMM={phi_fmm:.6f}, Err={err_phi:.2e}")
        print(f"  Force:     Direct=({fx_direct:.6f}, {fy_direct:.6f}), "
              f"FMM=({fx_fmm:.6f}, {fy_fmm:.6f}), RelErr={err_force_rel:.2e}")

    print()
    print(f"Max potential error: {max_err_phi:.2e}")
    print(f"Max force rel error: {max_err_force:.2e}")

    if max_err_force < 1e-8:
        print("✓ TEST PASSED (excellent accuracy)")
    elif max_err_force < 1e-6:
        print("✓ TEST PASSED (good accuracy)")
    else:
        print("✗ TEST FAILED (poor accuracy)")
    print()


def test_two_vortices():
    """Test 2: Two vortices - dipole configuration"""
    print("=" * 60)
    print("TEST 2: Two Vortices (Dipole)")
    print("=" * 60)

    # Setup: +q at (-0.5, 0), -q at (+0.5, 0)
    q = 1.0
    particles_source = [(-0.5, 0.0, q), (0.5, 0.0, -q)]
    center_source = (0.0, 0.0)

    # Test points far away
    test_points = [
        (2.0, 0.0),
        (0.0, 2.0),
        (2.0, 2.0),
        (3.0, 0.0),
    ]

    fmm = FMM2D_Log(p_max=20, rscale=1.0)
    M = fmm.p2m(particles_source, center_source)

    print(f"Multipole moments:")
    print(f"  M_0 (monopole) = {M[0]:.6e} (should be ≈ 0)")
    print(f"  M_1 (dipole)   = {M[1]:.6e} (should be ≠ 0)")
    print()

    max_err_force = 0.0

    for x, y in test_points:
        # Direct calculation
        phi_direct = sum(log_potential(x, y, px, py, pq)
                        for px, py, pq in particles_source)
        fx_direct = sum(log_force(x, y, px, py, pq)[0]
                       for px, py, pq in particles_source)
        fy_direct = sum(log_force(x, y, px, py, pq)[1]
                       for px, py, pq in particles_source)

        # FMM calculation
        z0 = complex(x - center_source[0], y - center_source[1])
        L = fmm.m2l(M, z0)
        phi_fmm, fx_fmm, fy_fmm = fmm.l2p(L, complex(0, 0))

        # Errors
        err_phi = abs(phi_fmm - phi_direct)
        err_force = np.sqrt((fx_fmm - fx_direct)**2 + (fy_fmm - fy_direct)**2)
        force_mag = np.sqrt(fx_direct**2 + fy_direct**2)
        if force_mag > 1e-14:
            err_force_rel = err_force / force_mag
        else:
            err_force_rel = 0.0

        max_err_force = max(max_err_force, err_force_rel)

        print(f"Point ({x:.1f}, {y:.1f}):")
        print(f"  Potential: Direct={phi_direct:.6f}, FMM={phi_fmm:.6f}, Err={err_phi:.2e}")
        print(f"  Force:     Direct=({fx_direct:.6f}, {fy_direct:.6f}), "
              f"FMM=({fx_fmm:.6f}, {fy_fmm:.6f}), RelErr={err_force_rel:.2e}")

    print()
    print(f"Max force rel error: {max_err_force:.2e}")

    if max_err_force < 1e-8:
        print("✓ TEST PASSED (excellent accuracy)")
    elif max_err_force < 1e-6:
        print("✓ TEST PASSED (good accuracy)")
    else:
        print("✗ TEST FAILED (poor accuracy)")
    print()


def test_convergence():
    """Test 3: Convergence with p_order"""
    print("=" * 60)
    print("TEST 3: Convergence with p_order")
    print("=" * 60)

    # Setup
    q = 1.0
    particles = [(0.0, 0.0, q)]
    center = (0.0, 0.0)
    test_point = (3.0, 0.0)

    # Direct calculation
    phi_direct = log_potential(*test_point, 0.0, 0.0, q)
    fx_direct, fy_direct = log_force(*test_point, 0.0, 0.0, q)
    force_direct = np.sqrt(fx_direct**2 + fy_direct**2)

    print(f"Test point: ({test_point[0]:.1f}, {test_point[1]:.1f})")
    print(f"Direct: φ={phi_direct:.6f}, F=({fx_direct:.6f}, {fy_direct:.6f})")
    print()
    print("p_order | Potential Error | Force Rel Error")
    print("-" * 50)

    p_orders = [5, 10, 15, 20, 25, 30]
    errors = []

    for p in p_orders:
        fmm = FMM2D_Log(p_max=p, rscale=1.0)
        M = fmm.p2m(particles, center)

        z0 = complex(test_point[0] - center[0], test_point[1] - center[1])
        L = fmm.m2l(M, z0)
        phi_fmm, fx_fmm, fy_fmm = fmm.l2p(L, complex(0, 0))

        err_phi = abs(phi_fmm - phi_direct)
        err_force = np.sqrt((fx_fmm - fx_direct)**2 + (fy_fmm - fy_direct)**2) / force_direct

        errors.append((p, err_phi, err_force))
        print(f"  {p:3d}   |   {err_phi:.2e}      |   {err_force:.2e}")

    print()

    # Check convergence
    improving = True
    for i in range(1, len(errors)):
        if errors[i][2] >= errors[i-1][2]:
            improving = False
            break

    if improving and errors[-1][2] < 1e-8:
        print("✓ TEST PASSED (error decreases with p_order)")
    else:
        print("✗ TEST FAILED (error should decrease with p_order)")
    print()


def test_rotational_symmetry():
    """Test 4: Rotational symmetry for ring of vortices"""
    print("=" * 60)
    print("TEST 4: Rotational Symmetry")
    print("=" * 60)

    # Setup: N vortices on a circle
    N = 8
    R = 1.0
    q = 1.0

    particles = []
    for i in range(N):
        theta = 2 * np.pi * i / N
        x = R * np.cos(theta)
        y = R * np.sin(theta)
        particles.append((x, y, q))

    center = (0.0, 0.0)

    fmm = FMM2D_Log(p_max=20, rscale=1.0)

    # For each vortex, compute force and check it's radial
    print(f"Testing {N} vortices on circle of radius {R}")
    print()

    max_tangential = 0.0

    for i in range(N):
        x, y, _ = particles[i]
        theta = np.arctan2(y, x)

        # Direct calculation
        fx_total = 0.0
        fy_total = 0.0
        for j, (xj, yj, qj) in enumerate(particles):
            if i == j:
                continue
            fx, fy = log_force(x, y, xj, yj, qj)
            fx_total += fx
            fy_total += fy

        # Decompose force into radial and tangential
        F_radial = fx_total * np.cos(theta) + fy_total * np.sin(theta)
        F_tangent = -fx_total * np.sin(theta) + fy_total * np.cos(theta)

        F_mag = np.sqrt(fx_total**2 + fy_total**2)
        if F_mag > 1e-14:
            F_tangent_rel = abs(F_tangent) / F_mag
        else:
            F_tangent_rel = 0.0

        max_tangential = max(max_tangential, F_tangent_rel)

        print(f"Vortex {i}: θ={np.degrees(theta):6.1f}°, "
              f"F_r={F_radial:8.5f}, F_t={F_tangent:.2e} "
              f"(rel: {F_tangent_rel:.2e})")

    print()
    print(f"Max tangential component (relative): {max_tangential:.2e}")

    if max_tangential < 1e-10:
        print("✓ TEST PASSED (excellent symmetry)")
    elif max_tangential < 1e-8:
        print("✓ TEST PASSED (good symmetry)")
    else:
        print("✗ TEST FAILED (symmetry broken)")
    print()


def main():
    """Run all tests"""
    print()
    print("=" * 60)
    print("FMM 2D Logarithmic Potential - Formula Validation")
    print("Testing CORRECTED M2L formulas")
    print("=" * 60)
    print()

    test_single_vortex()
    test_two_vortices()
    test_convergence()
    test_rotational_symmetry()

    print("=" * 60)
    print("ALL TESTS COMPLETED")
    print("=" * 60)


if __name__ == "__main__":
    main()
