"""Random field generation on 3-D periodic and bounded grids.

Gaussian fields are generated via the gstools library (same engine SPAM uses)
through a thin wrapper `simulate_random_field` that mirrors SPAM's API. All
higher-level generators (`generate_gaussian_rf`, `generate_margin_fields`, etc.)
dispatch to it internally with the current project convention:

    C(r) = exp(-r^2 / lc^2)       (correlation function, our convention)
    gstools equivalent: Gaussian(len_scale=lc, rescale=1.0, var=1.0)

With `rescale=1.0` gstools uses exactly `C(r) = exp(-(r / len_scale)^2)` which
matches ours — one of the few models where our `lc` and gstools' `len_scale`
coincide directly. Other covariance models (Matern, Exponential, ...) can be
used via `simulate_random_field` with their own kwargs.

What changes between field families (chi², Student-t, margin) is only the
marginal distribution; the spatial correlation structure is identical.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import gstools
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

def _rng_to_seed(rng: np.random.Generator) -> int:
    """Draw a 32-bit seed from a numpy generator for gstools.

    gstools consumes 32-bit integer seeds; the parent numpy generator stays stateful
    so consecutive calls draw independent fields.
    """
    return int(rng.integers(0, 2**31 - 1))


def simulate_random_field(
    lengths: float | list[float] = 1.0,
    nNodes: int | list[int] = 100,
    covarianceModel: str = "Gaussian",
    covarianceParameters: dict[str, Any] | None = None,
    dim: int = 3,
    nRea: int = 1,
    seed: int | None = None,
) -> list[np.ndarray]:
    """Generate correlated-random-field realisations via gstools.

    API mirrors `spam.excursions.simulateRandomField` so results and parameter
    conventions match SPAM-based pipelines directly.

    Parameters
    ----------
    lengths: float or list of floats
        Physical length of the domain in each direction (len == dim). Scalar is
        broadcast across all directions.
    nNodes: int or list of ints
        Number of grid nodes in each direction (len == dim). Scalar is broadcast.
    covarianceModel: str
        Name of a gstools.covmodel class ("Gaussian", "Exponential", "Matern", ...).
    covarianceParameters: dict
        Kwargs forwarded to the covariance constructor. For `"Gaussian"`, pass
        `len_scale=lc` and keep `rescale=1.0` to match our convention
        `C(r) = exp(-r^2/lc^2)`. If `rescale` is not supplied it defaults to 1.0
        here (deviation from gstools' upstream default, matching SPAM's wrapper).
    dim: int
        Spatial dimension.
    nRea: int
        Number of independent realisations.
    seed: int or None
        Master seed. None uses a system random seed.

    Returns
    -------
    list of numpy arrays, each of shape `nNodes`
        The `nRea` realisations.
    """
    def _as_list(v, types):
        return [v] * dim if isinstance(v, types) else list(v)

    lengths = _as_list(lengths, (int, float))
    nNodes  = _as_list(nNodes,  int)
    if len(lengths) != dim:
        raise ValueError(f"lengths length {len(lengths)} != dim {dim}")
    if len(nNodes) != dim:
        raise ValueError(f"nNodes length {len(nNodes)} != dim {dim}")

    cov_kw = dict(covarianceParameters or {})
    cov_kw.setdefault("rescale", 1.0)
    model = getattr(gstools.covmodel, covarianceModel)(dim=dim, **cov_kw)

    # Nodes at 0, dx, 2dx, ..., (n-1)dx — matches the project's existing grid
    # convention (dx = L / n). SPAM uses linspace(0, L, n) giving dx = L/(n-1);
    # we keep L/n so experiment outputs stay compatible with prior runs.
    pos = [np.arange(n) * (l / n) for n, l in zip(nNodes, lengths)]

    srf = gstools.SRF(model)
    srf.set_pos(pos, "structured")
    master = gstools.random.MasterRNG(seed)

    for i in range(nRea):
        srf(seed=master(), store=f"field_{i}")
    return list(srf.all_fields)


def _generate_one(n: int, lc_voxels: float, rng: np.random.Generator) -> np.ndarray:
    """Single Gaussian RF realisation on an n^3 cubic grid.

    Internal shim used by the higher-level generators. `lc_voxels` is the
    correlation length expressed in grid spacings; the field is generated on a
    unit-spacing grid (lengths=n, nNodes=n) so that `len_scale = lc_voxels`
    maps directly into gstools' coordinate system.
    """
    seed = _rng_to_seed(rng)
    fields = simulate_random_field(
        lengths=float(n),
        nNodes=n,
        dim=3,
        nRea=1,
        covarianceModel="Gaussian",
        covarianceParameters={"len_scale": float(lc_voxels), "var": 1.0},
        seed=seed,
    )
    return fields[0].astype(np.float32)


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


def generate_chi2_rf(
    grid: GridSpec, lc: float, k: int, rng: np.random.Generator
) -> np.ndarray:
    """Chi-squared RF with k degrees of freedom: X(x) = sum_i f_i(x)^2.

    Each f_i is an independent unit-variance Gaussian RF with correlation length lc.
    The resulting field has mean k and variance 2k; it is NOT standardised so that
    the marginal stays chi2(k) and the Worsley M_j formulas apply directly.

    Args:
        grid: Grid specification.
        lc:   Correlation length in grid.length units.
        k:    Degrees of freedom (integer >= 1).
        rng:  NumPy random generator.

    Returns:
        float32 array of shape (n, n, n) with chi2(k) marginal (mean=k, var=2k).
    """
    k = int(k)
    if k < 1:
        raise ValueError("k must be >= 1")
    lc_voxels = lc / grid.dx
    return sum(_generate_one(grid.n, lc_voxels, rng) ** 2 for _ in range(k)).astype(np.float32)


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
