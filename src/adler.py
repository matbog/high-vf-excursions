"""Adler-Taylor expected LKC formulas for Gaussian, Student-t, and margin hitting sets.

The Adler-Taylor kinematic formula (Adler & Taylor 2007, eq. 12.4.2) for a smooth
stationary isotropic Gaussian-related random field on a 3D cube of side a is:

    E[L0(A_u)] = rho^3 a^3 M3 + 3 rho^2 a^2 M2 + 3 rho a M1 + M0

where rho = sqrt(lambda2 / 2pi), lambda2 = -C''(0) is the second spectral moment,
and M_j are the Gaussian Minkowski functionals of the hitting set D_u in value space.

What changes between field families (Gaussian, t_nu, margin) is ONLY M_j.
The spatial side (lambda2, domain LKCs) is identical as long as the correlation
structure is unchanged.
"""

from __future__ import annotations

from math import comb

import numpy as np
from scipy import optimize, stats
from scipy.integrate import quad
from scipy.special import erfc, ndtr

from .fields import GridSpec

Array = np.ndarray


# =============================================================================
# Spectral moment
# =============================================================================

def lambda2(lc: float, sigma: float = 1.0) -> float:
    """Second spectral moment for C(r) = sigma^2 * exp(-r^2 / lc^2).

    lambda2 = -C''(0) = 2 sigma^2 / lc^2.
    """
    return 2.0 * sigma ** 2 / lc ** 2


def _rho(lc: float, sigma: float = 1.0) -> float:
    """rho = sqrt(lambda2 / 2pi) = sigma / (lc * sqrt(pi))."""
    return np.sqrt(lambda2(lc, sigma) / (2.0 * np.pi))


# =============================================================================
# Gaussian hitting set  D_u = [u, +inf)
# =============================================================================

def gaussian_m0(u: Array, sigma: float = 1.0) -> Array:
    """M0 = P(Z >= u) = (1/2) erfc(u / (sigma sqrt(2)))."""
    return 0.5 * erfc(np.asarray(u) / (sigma * np.sqrt(2.0)))


def gaussian_m1(u: Array, sigma: float = 1.0) -> Array:
    """M1 = phi(u) = (1 / sigma sqrt(2pi)) exp(-u^2 / 2 sigma^2)."""
    u = np.asarray(u)
    return np.exp(-0.5 * (u / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi))


def gaussian_m2(u: Array, sigma: float = 1.0) -> Array:
    """M2 = (u / sigma^2) phi(u)."""
    u = np.asarray(u)
    return (u / sigma ** 2) * gaussian_m1(u, sigma)


def gaussian_m3(u: Array, sigma: float = 1.0) -> Array:
    """M3 = ((u^2 - sigma^2) / sigma^4) phi(u)."""
    u = np.asarray(u)
    return ((u ** 2 - sigma ** 2) / sigma ** 4) * gaussian_m1(u, sigma)


def expected_lkcs_gaussian(
    grid: GridSpec, lc: float, kappas: Array, sigma: float = 1.0
) -> dict[str, Array]:
    """Expected L0, L1, L2, L3 and phi for a Gaussian RF on a 3D cube.

    Convention (Adler-Taylor / Roubin System 2.17):
      L0 = Euler characteristic
      L1 = integrated mean curvature / pi  (theoretically exact; not measured numerically)
      L2 = half surface area
      L3 = volume

    Returns dict with keys 'L0', 'L1', 'L2', 'L3', 'phi'.
    """
    kappas = np.asarray(kappas, dtype=float)
    a = grid.length
    rho = _rho(lc, sigma)
    M0 = gaussian_m0(kappas, sigma)
    M1 = gaussian_m1(kappas, sigma)
    M2 = gaussian_m2(kappas, sigma)
    M3 = gaussian_m3(kappas, sigma)
    L0 = rho ** 3 * a ** 3 * M3 + 3.0 * rho ** 2 * a ** 2 * M2 + 3.0 * rho * a * M1 + M0
    L1 = rho ** 2 * a ** 3 * M2 + 3.0 * rho * a ** 2 * M1 + 3.0 * a * M0
    L2 = 2.0 * rho * a ** 3 * M1 + 3.0 * a ** 2 * M0
    L3 = a ** 3 * M0
    return {"L0": L0, "L1": L1, "L2": L2, "L3": L3, "phi": M0}


# =============================================================================
# Student-t hitting set  D_u = [u, +inf) under t_nu marginal
# (Worsley 1994; derivation in theory/student_t_derivation.md)
# =============================================================================

def t_m0(u: Array, nu: float) -> Array:
    """M0 = P(T_nu >= u) = scipy.stats.t.sf(u, df=nu)."""
    return stats.t.sf(np.asarray(u), df=nu)


def t_m1(u: Array, nu: float) -> Array:
    """M1 = p_nu(u), the t_nu PDF."""
    return stats.t.pdf(np.asarray(u), df=nu)


def t_m2(u: Array, nu: float) -> Array:
    """M2 = p_nu(u) * (nu+1) u / (nu + u^2)."""
    u = np.asarray(u)
    return stats.t.pdf(u, df=nu) * (nu + 1.0) * u / (nu + u ** 2)


def t_m3(u: Array, nu: float) -> Array:
    """M3 = p_nu(u) * [(nu+1)(nu+3) u^2 / (nu+u^2)^2 - (nu+1) / (nu+u^2)]."""
    u = np.asarray(u)
    denom = nu + u ** 2
    return stats.t.pdf(u, df=nu) * (
        (nu + 1.0) * (nu + 3.0) * u ** 2 / denom ** 2 - (nu + 1.0) / denom
    )


def expected_lkcs_t(
    grid: GridSpec, lc: float, nu: float, kappas: Array, sigma: float = 1.0
) -> dict[str, Array]:
    """Expected L0 and phi for a Worsley t_nu field on a 3D cube.

    Returns dict with keys 'L0', 'phi'.
    (L2 and L3 formulas exist but L2 has not been independently validated
    in this project — use Gaussian L2 formula as an approximation only.)
    """
    kappas = np.asarray(kappas, dtype=float)
    a = grid.length
    rho = _rho(lc, sigma)
    M0 = t_m0(kappas, nu)
    M1 = t_m1(kappas, nu)
    M2 = t_m2(kappas, nu)
    M3 = t_m3(kappas, nu)
    L0 = rho ** 3 * a ** 3 * M3 + 3.0 * rho ** 2 * a ** 2 * M2 + 3.0 * rho * a * M1 + M0
    return {"L0": L0, "phi": M0}


def find_percolation_threshold(
    grid: GridSpec,
    lc: float,
    nu: float | None = None,
    sigma: float = 1.0,
    u_search: tuple[float, float] = (0.5, 5.0),
) -> tuple[float, float]:
    """Find kappa_c+ (upper zero-crossing of E[L0]) and the corresponding phi_c+.

    phi_c+ is the maximum achievable disconnected volume fraction according to
    the Adler EC criterion (Roubin & Colliat 2016).

    Args:
        grid:     Grid specification (for domain LKCs).
        lc:       Correlation length.
        nu:       t_nu degrees of freedom, or None for Gaussian.
        sigma:    Field standard deviation.
        u_search: (lower, upper) bracket for root search.

    Returns:
        (kappa_c+, phi_c+).  Returns (nan, nan) if no zero-crossing found.
    """
    def ec(u_scalar: float) -> float:
        u = np.asarray([u_scalar])
        if nu is None:
            return float(expected_lkcs_gaussian(grid, lc, u, sigma)["L0"][0])
        return float(expected_lkcs_t(grid, lc, nu, u, sigma)["L0"][0])

    u_lo, u_hi = u_search
    if ec(u_lo) * ec(u_hi) > 0:
        return float("nan"), float("nan")
    kc = float(optimize.brentq(ec, u_lo, u_hi, xtol=1e-4))
    if nu is None:
        phi_c = float(gaussian_m0(np.asarray([kc]), sigma)[0])
    else:
        phi_c = float(stats.t.sf(kc, df=nu))
    return kc, phi_c


# =============================================================================
# Margin tessellation hitting set  D_delta = {z: max_i z_i - secondmax_i z_i >= delta}
# (derived in theory/gaussian_minkowski_functionals.md)
# =============================================================================

def _normal_density(x: float) -> float:
    return float(np.exp(-0.5 * x * x) / np.sqrt(2.0 * np.pi))


def margin_m0(q: int, delta: float) -> float:
    """M0 = q integral phi(t) Phi(t-delta)^(q-1) dt = P(A_delta)."""
    integrand = lambda t: _normal_density(t) * ndtr(t - delta) ** (q - 1)
    val, _ = quad(integrand, -10.0, 10.0, epsabs=2e-11, epsrel=2e-11, limit=200)
    return float(q * val)


def margin_m1(q: int, delta: float) -> float:
    """M1: Gaussian boundary measure of the q polyhedral winner chambers."""
    integrand = (
        lambda t: _normal_density(t) * _normal_density(t - delta) * ndtr(t - delta) ** (q - 2)
    )
    val, _ = quad(integrand, -10.0, 10.0, epsabs=2e-11, epsrel=2e-11, limit=200)
    return float(np.sqrt(2.0) * q * (q - 1) * val)


def margin_stratum_measure(q: int, delta: float, codimension: int) -> float:
    """Gaussian measure of codimension-k boundary strata.

    One winning coordinate, k runner-up coordinates exactly delta below it,
    all remaining below those. The sqrt(k+1) factor is the induced surface Jacobian.
    """
    k = codimension
    if not 1 <= k <= min(3, q - 1):
        return 0.0
    integrand = (
        lambda t: _normal_density(t + delta) * _normal_density(t) ** k * ndtr(t) ** (q - 1 - k)
    )
    val, _ = quad(integrand, -10.0, 10.0, epsabs=2e-11, epsrel=2e-11, limit=300)
    multiplicity = q * comb(q - 1, k)
    surface_jacobian = np.sqrt(k + 1.0)
    return float(multiplicity * surface_jacobian * val)


def margin_gmfs(q: int, delta: float) -> np.ndarray:
    """Closed one-dimensional integral formulas for M0..M3 of D_delta.

    Derived from the Gaussian tube expansion of the polyhedral stratification.
    Normal-cone constants: pairwise wedge angle pi/3; trihedral solid angle
    3*arccos(1/3) - pi.

    NOTE: These formulas are numerically validated for q=2 (exact half-space
    identities) and for q=4,8,16 at phi=0.50,0.70,0.80. A formal proof of
    the stratified-tube decomposition has not yet been published.

    Returns:
        Array [M0, M1, M2, M3].
    """
    a = delta / np.sqrt(2.0)
    S1 = margin_stratum_measure(q, delta, 1)
    S2 = margin_stratum_measure(q, delta, 2)
    S3 = margin_stratum_measure(q, delta, 3)
    wedge = np.pi / 3.0
    trihedral = 3.0 * np.arccos(1.0 / 3.0) - np.pi
    M0 = margin_m0(q, delta)
    M1 = S1
    M2 = a * S1 + wedge * S2
    M3 = (a ** 2 - 1.0) * S1 + (4.0 * a / np.sqrt(3.0)) * S2 + 2.0 * trihedral * S3
    return np.array([M0, M1, M2, M3], dtype=float)


def margin_expected_euler(
    grid: GridSpec, lc: float, q: int, delta: float, sigma: float = 1.0
) -> float:
    """Expected Euler characteristic for the margin tessellation on a 3D cube."""
    a = grid.length
    rho = _rho(lc, sigma)
    gmf = margin_gmfs(q, delta)
    return float(
        rho ** 3 * a ** 3 * gmf[3]
        + 3.0 * rho ** 2 * a ** 2 * gmf[2]
        + 3.0 * rho * a * gmf[1]
        + gmf[0]
    )


def gaussian_margin_vf(q: int, delta: float) -> float:
    """Expected volume fraction of A_delta = margin_m0(q, delta)."""
    return margin_m0(q, delta)


def delta_for_vf(q: int, phi: float) -> float:
    """Invert margin_m0: find delta such that gaussian_margin_vf(q, delta) = phi."""
    if not 0.0 < phi < 1.0:
        raise ValueError("phi must be in (0, 1)")
    return float(
        optimize.brentq(lambda d: gaussian_margin_vf(q, d) - phi, 0.0, 10.0)
    )


# =============================================================================
# Coalescing lobes: Gaussian-space ball union
# Volume fraction is estimated by Monte Carlo (no closed-form M_j before contact).
# =============================================================================

def gaussian_lobe_vf_mc(
    centres: np.ndarray,
    radii: np.ndarray,
    n_samples: int = 2_000_000,
    seed: int = 0,
) -> np.ndarray:
    """P(min_a ||Z - c_a|| <= r) for each r in radii, estimated by MC.

    Z ~ N(0, I_d) where d = centres.shape[1].
    Before the contact radius r_c = ||c_a - c_b|| / 2, M_j(D_r) can in principle
    be derived analytically (non-central chi^2 CDF for a single ball). After
    r_c, ball intersections contribute through inclusion-exclusion, and the
    required Gaussian-weighted intersection volumes have not yet been derived.

    Args:
        centres:   Array of shape (K, d).
        radii:     Array of radii to evaluate.
        n_samples: Number of Monte Carlo samples.
        seed:      Random seed.

    Returns:
        Array of shape (len(radii),) with volume fraction estimates.
    """
    rng = np.random.default_rng(seed)
    radii = np.asarray(radii)
    counts = np.zeros(len(radii), dtype=np.int64)
    batch = 200_000
    for start in range(0, n_samples, batch):
        z = rng.standard_normal((min(batch, n_samples - start), centres.shape[1]))
        d2 = np.min(np.sum((z[:, None] - centres[None]) ** 2, axis=2), axis=1)
        counts += np.sum(d2[:, None] <= radii[None] ** 2, axis=0)
    return counts / n_samples
