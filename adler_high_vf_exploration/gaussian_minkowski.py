"""Low-order Gaussian Minkowski functionals for the margin hitting set."""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from math import comb


def normal_density(x: float) -> float:
    return float(np.exp(-0.5 * x * x) / np.sqrt(2.0 * np.pi))


def margin_m0(q: int, delta: float) -> float:
    """Gaussian measure of D_delta (the expected retained volume fraction)."""
    integrand = lambda t: normal_density(t) * ndtr(t - delta) ** (q - 1)
    value, _ = quad(integrand, -10.0, 10.0, epsabs=2e-11, epsrel=2e-11, limit=200)
    return float(q * value)


def margin_m1(q: int, delta: float) -> float:
    """Gaussian boundary measure of the q separated margin chambers."""
    integrand = (
        lambda t: normal_density(t)
        * normal_density(t - delta)
        * ndtr(t - delta) ** (q - 2)
    )
    value, _ = quad(integrand, -10.0, 10.0, epsabs=2e-11, epsrel=2e-11, limit=200)
    return float(np.sqrt(2.0) * q * (q - 1) * value)


def margin_stratum_measure(q: int, delta: float, codimension: int) -> float:
    """Total Gaussian measure of codimension-k boundary strata.

    A stratum has one winning coordinate, k runner-up coordinates exactly
    delta below it, and all remaining coordinates below those runners.
    """
    k = codimension
    if not 1 <= k <= min(3, q - 1):
        return 0.0
    integrand = (
        lambda t: normal_density(t + delta)
        * normal_density(t) ** k
        * ndtr(t) ** (q - 1 - k)
    )
    value, _ = quad(integrand, -10.0, 10.0, epsabs=2e-11, epsrel=2e-11, limit=300)
    multiplicity = q * comb(q - 1, k)
    surface_jacobian = np.sqrt(k + 1.0)
    return float(multiplicity * surface_jacobian * value)


def margin_gmfs(q: int, delta: float) -> np.ndarray:
    """General closed integral formulas for M0..M3 of the margin hitting set.

    The formulas follow from the Gaussian tube expansion of the polyhedral
    stratification. Facets, pairwise edges and triple strata contribute through
    codimensions 1, 2 and 3. Pairwise normals have angle pi/3; the codimension-3
    normal cone has solid angle 3*acos(1/3)-pi.
    """
    a = delta / np.sqrt(2.0)
    s1 = margin_stratum_measure(q, delta, 1)
    s2 = margin_stratum_measure(q, delta, 2)
    s3 = margin_stratum_measure(q, delta, 3)
    wedge_angle = np.pi / 3.0
    trihedral_solid_angle = 3.0 * np.arccos(1.0 / 3.0) - np.pi

    m0 = margin_m0(q, delta)
    m1 = s1
    m2 = a * s1 + wedge_angle * s2
    m3 = (
        (a * a - 1.0) * s1
        + (4.0 * a / np.sqrt(3.0)) * s2
        + 2.0 * trihedral_solid_angle * s3
    )
    return np.array([m0, m1, m2, m3], dtype=float)
