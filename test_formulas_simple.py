#!/usr/bin/env python3
"""
Simple test of FMM formulas without external dependencies
Tests the mathematical correctness of M2L translation
"""
import math
import cmath

def binomial(n, k):
    """Compute binomial coefficient C(n,k)"""
    if k > n or k < 0:
        return 0
    if k == 0 or k == n:
        return 1
    result = 1
    for i in range(min(k, n - k)):
        result = result * (n - i) // (i + 1)
    return result

class FMM_Test:
    """Simple FMM test for 2D log kernel"""

    def __init__(self, p_max):
        self.p_max = p_max
        self.rscale = 1.0

    def p2m(self, particles, center):
        """Compute multipole moments"""
        M = [0.0 + 0.0j] * (self.p_max + 1)
        cx, cy = center

        for x, y, q in particles:
            z = complex(x - cx, y - cy)
            z_scaled = z / self.rscale

            # M_0 = Σ q_i
            M[0] += q

            # M_k = -Σ q_i · (z_i/rscale)^k / k
            z_pow = z_scaled
            for k in range(1, self.p_max + 1):
                M[k] += -q * z_pow / k
                z_pow *= z_scaled

        return M

    def m2l_corrected(self, M, z0):
        """M2L with CORRECTED formulas"""
        L = [0.0 + 0.0j] * (self.p_max + 1)
        z0_scaled = z0 / self.rscale

        if abs(z0_scaled) < 1e-14:
            return L

        z0_inv = 1.0 / z0_scaled

        # CORRECTED L_0: M_0 * (-log|z0|) + Σ M_k * z0^k
        L[0] = -M[0] * cmath.log(abs(z0_scaled))  # NEGATIVE SIGN

        z0_pow = z0_scaled
        for k in range(1, self.p_max + 1):
            L[0] += M[k] * z0_pow  # z0^k powers
            z0_pow *= z0_scaled

        # CORRECTED L_j: -M_0/(j·z0^j) + Σ M_k·C(j+k-1,k-1)/z0^(j+k)
        z0_inv_j = z0_inv
        for j in range(1, self.p_max + 1):
            temp = -M[0] * z0_inv_j / j

            z0_inv_jk = z0_inv_j * z0_inv
            for k in range(1, self.p_max + 1):
                bc = binomial(j + k - 1, k - 1)
                temp += M[k] * bc * z0_inv_jk  # z0^(j+k) in denominator
                z0_inv_jk *= z0_inv

            L[j] = temp
            z0_inv_j *= z0_inv

        return L

    def m2l_old_buggy(self, M, z0):
        """M2L with OLD BUGGY formulas (for comparison)"""
        L = [0.0 + 0.0j] * (self.p_max + 1)
        z0_scaled = z0 / self.rscale

        if abs(z0_scaled) < 1e-14:
            return L

        z0_inv = 1.0 / z0_scaled

        # BUGGY L_0: M_0 * (+log|z0|) + Σ M_k (no z0^k!)
        L[0] = M[0] * cmath.log(abs(z0_scaled))  # WRONG SIGN

        for k in range(1, self.p_max + 1):
            L[0] += M[k]  # MISSING z0^k!

        # BUGGY L_j: wrong powers
        z0_inv_j = z0_inv
        for j in range(1, self.p_max + 1):
            temp = -M[0] / j  # Missing z0_inv_j factor

            for k in range(1, self.p_max + 1):
                bc = binomial(j + k - 1, k - 1)
                temp += M[k] * bc  # MISSING z0^(j+k) division!

            L[j] = temp * z0_inv_j  # Wrong - should be inside loop
            z0_inv_j *= z0_inv

        return L

    def l2p(self, L, z):
        """Evaluate local expansion"""
        z_scaled = z / self.rscale

        # Potential
        phi = L[0].real
        z_pow = z_scaled
        for k in range(1, self.p_max + 1):
            phi += (L[k] * z_pow).real
            z_pow *= z_scaled

        # Gradient
        grad = L[1]
        z_pow = z_scaled
        for k in range(2, self.p_max + 1):
            grad += k * L[k] * z_pow
            z_pow *= z_scaled
        grad /= self.rscale

        fx = -grad.real
        fy = grad.imag

        return phi, fx, fy

    def direct_potential(self, x, y, particles):
        """Direct calculation"""
        phi = 0.0
        for px, py, q in particles:
            r = math.sqrt((x - px)**2 + (y - py)**2)
            if r > 1e-14:
                phi += -q * math.log(r)
        return phi

    def direct_force(self, x, y, particles):
        """Direct force calculation"""
        fx, fy = 0.0, 0.0
        for px, py, q in particles:
            dx = x - px
            dy = y - py
            r2 = dx*dx + dy*dy
            if r2 > 1e-28:
                fx += q * dx / r2
                fy += q * dy / r2
        return fx, fy


def test_single_vortex():
    """Test 1: Single vortex at origin"""
    print("=" * 70)
    print("TEST 1: Single Vortex")
    print("=" * 70)

    fmm = FMM_Test(p_max=20)
    particles = [(0.0, 0.0, 1.0)]
    center = (0.0, 0.0)

    test_point = (2.0, 0.0)
    print(f"Test point: {test_point}")
    print()

    # Direct
    phi_direct = fmm.direct_potential(*test_point, particles)
    fx_direct, fy_direct = fmm.direct_force(*test_point, particles)
    print(f"Direct: φ={phi_direct:.8f}, F=({fx_direct:.8f}, {fy_direct:.8f})")

    # FMM (corrected)
    M = fmm.p2m(particles, center)
    z0 = complex(test_point[0] - center[0], test_point[1] - center[1])
    L_correct = fmm.m2l_corrected(M, z0)
    phi_fmm, fx_fmm, fy_fmm = fmm.l2p(L_correct, complex(0, 0))
    print(f"FMM:    φ={phi_fmm:.8f}, F=({fx_fmm:.8f}, {fy_fmm:.8f})")

    # Errors
    err_phi = abs(phi_fmm - phi_direct)
    err_f = math.sqrt((fx_fmm - fx_direct)**2 + (fy_fmm - fy_direct)**2)
    f_mag = math.sqrt(fx_direct**2 + fy_direct**2)
    err_f_rel = err_f / f_mag if f_mag > 1e-14 else 0.0

    print()
    print(f"Absolute errors: φ={err_phi:.2e}, F={err_f:.2e}")
    print(f"Relative force error: {err_f_rel:.2e}")

    if err_f_rel < 1e-8:
        print("✓ PASSED (excellent)")
    elif err_f_rel < 1e-6:
        print("✓ PASSED (good)")
    else:
        print("✗ FAILED")
    print()


def test_comparison_buggy_vs_correct():
    """Test 2: Compare buggy vs corrected formulas"""
    print("=" * 70)
    print("TEST 2: Buggy vs Corrected Formulas")
    print("=" * 70)

    fmm = FMM_Test(p_max=20)
    particles = [(0.0, 0.0, 1.0)]
    center = (0.0, 0.0)
    test_point = (2.0, 0.0)

    print(f"Single vortex at origin, test point at {test_point}")
    print()

    # Direct
    phi_direct = fmm.direct_potential(*test_point, particles)
    fx_direct, fy_direct = fmm.direct_force(*test_point, particles)

    # FMM corrected
    M = fmm.p2m(particles, center)
    z0 = complex(test_point[0], test_point[1])
    L_correct = fmm.m2l_corrected(M, z0)
    phi_correct, fx_correct, fy_correct = fmm.l2p(L_correct, complex(0, 0))

    # FMM buggy
    L_buggy = fmm.m2l_old_buggy(M, z0)
    phi_buggy, fx_buggy, fy_buggy = fmm.l2p(L_buggy, complex(0, 0))

    # Print results
    print(f"Direct:    φ={phi_direct:.8f}, F=({fx_direct:.6f}, {fy_direct:.6f})")
    print(f"Corrected: φ={phi_correct:.8f}, F=({fx_correct:.6f}, {fy_correct:.6f})")
    print(f"Buggy:     φ={phi_buggy:.8f}, F=({fx_buggy:.6f}, {fy_buggy:.6f})")
    print()

    # Errors
    f_mag = math.sqrt(fx_direct**2 + fy_direct**2)
    err_correct = math.sqrt((fx_correct - fx_direct)**2 +
                           (fy_correct - fy_direct)**2) / f_mag
    err_buggy = math.sqrt((fx_buggy - fx_direct)**2 +
                         (fy_buggy - fy_direct)**2) / f_mag

    print(f"Corrected formula error: {err_correct:.2e}")
    print(f"Buggy formula error:     {err_buggy:.2e}")
    print()

    if err_correct < 1e-8 and err_buggy > 1e-4:
        print("✓ PASSED: Corrected formulas are MUCH better!")
    else:
        print("? UNCLEAR: Check the formulas")
    print()


def test_dipole():
    """Test 3: Dipole configuration"""
    print("=" * 70)
    print("TEST 3: Dipole (Two Opposite Vortices)")
    print("=" * 70)

    fmm = FMM_Test(p_max=20)
    particles = [(-0.5, 0.0, 1.0), (0.5, 0.0, -1.0)]
    center = (0.0, 0.0)
    test_point = (3.0, 0.0)

    print(f"Vortex +1 at (-0.5, 0), vortex -1 at (+0.5, 0)")
    print(f"Test point: {test_point}")
    print()

    # Direct
    phi_direct = fmm.direct_potential(*test_point, particles)
    fx_direct, fy_direct = fmm.direct_force(*test_point, particles)

    # FMM
    M = fmm.p2m(particles, center)
    print(f"Multipole moments: M_0={M[0].real:.6f} (should be ≈0)")
    print(f"                   M_1={abs(M[1]):.6f} (dipole, should be ≠0)")
    print()

    z0 = complex(test_point[0] - center[0], test_point[1] - center[1])
    L = fmm.m2l_corrected(M, z0)
    phi_fmm, fx_fmm, fy_fmm = fmm.l2p(L, complex(0, 0))

    print(f"Direct: φ={phi_direct:.8f}, F=({fx_direct:.8f}, {fy_direct:.8f})")
    print(f"FMM:    φ={phi_fmm:.8f}, F=({fx_fmm:.8f}, {fy_fmm:.8f})")

    # Error
    f_mag = math.sqrt(fx_direct**2 + fy_direct**2)
    err_f_rel = math.sqrt((fx_fmm - fx_direct)**2 +
                         (fy_fmm - fy_direct)**2) / f_mag if f_mag > 1e-14 else 0.0

    print()
    print(f"Relative force error: {err_f_rel:.2e}")

    if err_f_rel < 1e-8:
        print("✓ PASSED (excellent)")
    elif err_f_rel < 1e-6:
        print("✓ PASSED (good)")
    else:
        print("✗ FAILED")
    print()


def test_convergence():
    """Test 4: Convergence with p_order"""
    print("=" * 70)
    print("TEST 4: Convergence with p_order")
    print("=" * 70)

    particles = [(0.0, 0.0, 1.0)]
    center = (0.0, 0.0)
    test_point = (3.0, 0.0)

    print(f"Single vortex, test point at distance {test_point[0]}")
    print()
    print("p_order | Force Rel Error")
    print("-" * 30)

    prev_error = 1.0
    converging = True

    for p in [5, 10, 15, 20, 25]:
        fmm = FMM_Test(p_max=p)

        # Direct
        fx_direct, fy_direct = fmm.direct_force(*test_point, particles)
        f_mag = math.sqrt(fx_direct**2 + fy_direct**2)

        # FMM
        M = fmm.p2m(particles, center)
        z0 = complex(test_point[0], test_point[1])
        L = fmm.m2l_corrected(M, z0)
        _, fx_fmm, fy_fmm = fmm.l2p(L, complex(0, 0))

        err = math.sqrt((fx_fmm - fx_direct)**2 +
                       (fy_fmm - fy_direct)**2) / f_mag

        print(f"  {p:3d}   |   {err:.2e}")

        if err >= prev_error:
            converging = False
        prev_error = err

    print()
    if converging:
        print("✓ PASSED: Error decreases with p_order")
    else:
        print("✗ FAILED: Error should decrease with p_order")
    print()


def main():
    print()
    print("=" * 70)
    print("  FMM 2D Logarithmic Kernel - Formula Validation")
    print("  Testing CORRECTED M2L formulas")
    print("=" * 70)
    print()

    test_single_vortex()
    test_comparison_buggy_vs_correct()
    test_dipole()
    test_convergence()

    print("=" * 70)
    print("  ALL TESTS COMPLETED")
    print("=" * 70)
    print()


if __name__ == "__main__":
    main()
