"""Numerical Gaussian tube coefficients for the margin hitting set.

For D_delta = {max(z)-secondmax(z) >= delta}, estimate

  gamma(Tube(D_delta, rho))
    = M0 + rho M1 + rho^2 M2/2 + rho^3 M3/6 + O(rho^4)

using common Sobol points in Gaussian value space.  Euclidean distance to the
nearest margin chamber is computed by projection onto its linear inequalities.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import ndtri
from scipy.stats import qmc

from gaussian_minkowski import margin_m0, margin_m1


@dataclass(frozen=True)
class TubeEstimate:
    q: int
    delta: float
    radii: np.ndarray
    probabilities: np.ndarray
    coefficients: np.ndarray
    residual_rms: float

    @property
    def gmfs(self) -> np.ndarray:
        # Polynomial coefficients use c_j rho^j, whereas tube coefficients are M_j/j!.
        return self.coefficients[:4] * np.array([1.0, 1.0, 2.0, 6.0])


def distance_to_margin_set(samples: np.ndarray, delta: float) -> np.ndarray:
    """Exact Euclidean distance to the union of winner-margin chambers.

    By permutation symmetry and the rearrangement inequality, a closest
    chamber can be chosen with the currently largest coordinate as winner.
    Projection onto that chamber reduces to a one-dimensional convex problem.
    """
    z = np.asarray(samples, dtype=float)
    if z.ndim != 2 or z.shape[1] < 2:
        raise ValueError("samples must have shape (n_samples, q), q >= 2")
    n, q = z.shape
    order = np.argsort(z, axis=1)
    winner = order[:, -1]
    zi = z[np.arange(n), winner]
    others = np.take_along_axis(z, order[:, :-1], axis=1)[:, ::-1]
    margin = zi - others[:, 0]
    distance = np.zeros(n, dtype=float)
    outside = margin < delta
    if not np.any(outside):
        return distance

    zi_o = zi[outside]
    b = others[outside] + delta
    cumulative = np.cumsum(b, axis=1)
    k_values = np.arange(1, q, dtype=float)
    candidates = (zi_o[:, None] + cumulative) / (k_values[None, :] + 1.0)

    # For active set size k, the next inactive b must be <= projected winner t.
    valid = np.ones_like(candidates, dtype=bool)
    valid[:, :-1] = b[:, 1:] <= candidates[:, :-1]
    selected_index = np.argmax(valid, axis=1)
    k_selected = selected_index + 1
    t = candidates[np.arange(len(zi_o)), selected_index]

    active = np.arange(q - 1)[None, :] < k_selected[:, None]
    squared = (t - zi_o) ** 2 + np.sum(np.where(active, (b - t[:, None]) ** 2, 0.0), axis=1)
    distance[outside] = np.sqrt(np.maximum(squared, 0.0))
    return distance


def estimate_tube_gmfs(
    q: int,
    delta: float,
    power: int = 20,
    max_radius: float = 0.10,
    n_radii: int = 21,
    fit_degree: int = 5,
    seed: int = 12345,
) -> TubeEstimate:
    """Estimate M0..M3, using higher polynomial terms to reduce truncation bias."""
    if fit_degree < 3:
        raise ValueError("fit_degree must be at least 3")
    sampler = qmc.Sobol(d=q, scramble=True, seed=seed)
    uniform = sampler.random_base2(power)
    # Avoid infinities if a generator ever returns an exact endpoint.
    normal = ndtri(np.clip(uniform, np.finfo(float).eps, 1.0 - np.finfo(float).eps))
    distances = distance_to_margin_set(normal, delta)
    radii = np.linspace(0.0, max_radius, n_radii)
    probabilities = np.array([np.mean(distances <= radius) for radius in radii])

    # Anchor known analytic M0 and M1, then fit M2 and M3 plus nuisance O(rho^4,5) terms.
    m0 = margin_m0(q, delta)
    m1 = margin_m1(q, delta)
    remainder = probabilities - m0 - m1 * radii
    design = np.column_stack([radii**degree for degree in range(2, fit_degree + 1)])
    higher, *_ = np.linalg.lstsq(design, remainder, rcond=None)
    coefficients = np.concatenate(([m0, m1], higher))
    fitted = sum(coefficients[j] * radii**j for j in range(len(coefficients)))
    residual = float(np.sqrt(np.mean((probabilities - fitted) ** 2)))
    return TubeEstimate(q, delta, radii, probabilities, coefficients, residual)


def q2_exact_gmfs(delta: float) -> np.ndarray:
    """Closed forms for the two separated half-spaces when q=2."""
    a = delta / np.sqrt(2.0)
    phi = np.exp(-0.5 * a * a) / np.sqrt(2.0 * np.pi)
    m0 = margin_m0(2, delta)
    return np.array([m0, 2.0 * phi, 2.0 * a * phi, 2.0 * (a * a - 1.0) * phi])

