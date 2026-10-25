"""Random field generation on 3-D periodic and bounded grids.

All fields are zero-mean, unit-variance, and stationary. The spatial covariance
is C(r) ~ exp(-r^2 / lc^2) throughout. What changes between field families is
only the marginal distribution — the spatial correlation structure is identical.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GridSpec:
    """Uniform cubic grid description."""

    n: int
    length: float

    @property
    def dx(self) -> float:
        return self.length / self.n

    @property
    def volume(self) -> float:
        return self.length ** 3


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _spectral_amplitude(n: int, lc_voxels: float) -> np.ndarray:
    """Square-root of the spectral density for C(r) ~ exp(-r^2/lc^2).

    lc_voxels is the correlation length expressed in grid spacings (dx=1).
    The amplitude is normalised so that irfftn produces a unit-variance field
    after dividing by its empirical std.
    """
    kx = 2.0 * np.pi * np.fft.fftfreq(n)
    ky = 2.0 * np.pi * np.fft.fftfreq(n)
    kz = 2.0 * np.pi * np.fft.rfftfreq(n)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    k2 = KX ** 2 + KY ** 2 + KZ ** 2
    return np.exp(-lc_voxels ** 2 * k2 / 8.0)


def _generate_one(n: int, lc_voxels: float, rng: np.random.Generator) -> np.ndarray:
    amp = _spectral_amplitude(n, lc_voxels)
    noise = rng.standard_normal((n, n, n // 2 + 1)) + 1j * rng.standard_normal((n, n, n // 2 + 1))
    field = np.fft.irfftn(noise * amp, s=(n, n, n), axes=(0, 1, 2)).real
    field -= field.mean()
    field /= field.std()
    return field.astype(np.float32)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_gaussian_rf(grid: GridSpec, lc: float, rng: np.random.Generator) -> np.ndarray:
    """Single stationary correlated Gaussian RF, zero mean, unit variance.

    Args:
        grid: Grid specification.
        lc:   Correlation length in the same units as grid.length.
        rng:  NumPy random generator (use np.random.default_rng(seed)).

    Returns:
        float32 array of shape (n, n, n).
    """
    lc_voxels = lc / grid.dx
    return _generate_one(grid.n, lc_voxels, rng)


def generate_student_t_rf(
    grid: GridSpec, lc: float, nu: int, rng: np.random.Generator
) -> np.ndarray:
    """Worsley t_nu random field: T(x) = Z_0(x) / sqrt(chi2_field(x) / nu).

    Z_0 and each chi^2 component share the same correlation length lc. The
    spatially-varying denominator is what makes the Adler formula with t_nu
    M_j values exact for this construction. A scalar chi^2 denominator would
    merely rescale Z_0 per realisation and the formula would not apply.

    Valid range: nu >= 3 (nu=1 Cauchy has undefined lambda_2).

    Args:
        grid: Grid specification.
        lc:   Correlation length in grid.length units.
        nu:   Degrees of freedom (integer >= 3).
        rng:  NumPy random generator.

    Returns:
        float32 array of shape (n, n, n).
    """
    nu = int(nu)
    if nu < 3:
        raise ValueError("nu must be >= 3 for finite lambda_2")
    lc_voxels = lc / grid.dx
    Z0 = _generate_one(grid.n, lc_voxels, rng)
    chi2 = sum(_generate_one(grid.n, lc_voxels, rng) ** 2 for _ in range(nu))
    return (Z0 / np.sqrt(chi2 / nu)).astype(np.float32)


def generate_margin_fields(
    grid: GridSpec, lc: float, q: int, rng: np.random.Generator
) -> np.ndarray:
    """Stack of q independent Gaussian RFs for the margin tessellation.

    Args:
        grid: Grid specification.
        lc:   Correlation length.
        q:    Number of competing fields.
        rng:  NumPy random generator.

    Returns:
        float32 array of shape (q, n, n, n).
    """
    lc_voxels = lc / grid.dx
    return np.stack([_generate_one(grid.n, lc_voxels, rng) for _ in range(q)])


def generate_vector_field(
    grid: GridSpec,
    lc: float,
    n_components: int,
    rng: np.random.Generator,
    padding: int = 0,
) -> np.ndarray:
    """Stack of n_components independent Gaussian RFs, optionally cropped.

    When padding > 0 the field is generated on a (n + 2*padding)^3 periodic
    grid and then the central n^3 cube is returned. This avoids the periodic-
    boundary wrap-around artefacts that appear near the domain edges when the
    field is used as a non-periodic Gaussian-space map.

    Args:
        grid:        Grid specification for the *output* field.
        lc:          Correlation length in grid.length units.
        n_components: Number of independent field components.
        rng:         NumPy random generator.
        padding:     Number of voxels to pad on each side before cropping.

    Returns:
        float32 array of shape (n_components, n, n, n).
    """
    lc_voxels = lc / grid.dx
    if padding == 0:
        return np.stack([_generate_one(grid.n, lc_voxels, rng) for _ in range(n_components)])
    full_n = grid.n + 2 * padding
    crop = slice(padding, padding + grid.n)
    return np.stack([
        _generate_one(full_n, lc_voxels, rng)[crop, crop, crop]
        for _ in range(n_components)
    ])
