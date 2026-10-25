"""Excursion set construction: thresholding, margin tessellation, Gaussian-space lobes."""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Scalar field excursions
# ---------------------------------------------------------------------------

def threshold(field: np.ndarray, kappa: float) -> np.ndarray:
    """Return boolean phase {field >= kappa}."""
    return field >= kappa


def threshold_at_vf(field: np.ndarray, phi: float) -> tuple[np.ndarray, float]:
    """Choose kappa so the excursion has volume fraction phi, return (phase, kappa).

    Args:
        field: Scalar field array.
        phi:   Target volume fraction in (0, 1).

    Returns:
        (phase array of bool, actual kappa used).
    """
    if not 0.0 < phi < 1.0:
        raise ValueError("phi must be in (0, 1)")
    kappa = float(np.quantile(field, 1.0 - phi))
    return field >= kappa, kappa


# ---------------------------------------------------------------------------
# Margin tessellation  A_delta = {max_i G_i - secondmax_i G_i >= delta}
# ---------------------------------------------------------------------------

def margin_excursion(
    fields: np.ndarray, delta: float
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the margin-tessellation excursion set and winner labels.

    The hitting set D_delta = {z : z_(q) - z_(q-1) >= delta} is the disjoint
    union of q polyhedral winner chambers. Each spatial voxel is retained when
    its winning field exceeds the runner-up by at least delta.

    Args:
        fields: float array of shape (q, n, n, n).
        delta:  Margin threshold (>= 0).

    Returns:
        phase:   bool array of shape (n, n, n), True = retained (in A_delta).
        winners: int16 array of shape (n, n, n), index of the winning field.
    """
    top2 = np.partition(fields, kth=fields.shape[0] - 2, axis=0)[-2:]
    margin = top2[1] - top2[0]
    winners = np.argmax(fields, axis=0).astype(np.int16)
    return margin >= delta, winners


# ---------------------------------------------------------------------------
# Coalescing Gaussian-space lobes
# ---------------------------------------------------------------------------

def tetrahedral_centres(rho: float) -> np.ndarray:
    """Four vertices of a regular tetrahedron inscribed in a sphere of radius rho.

    All pairwise distances equal rho * sqrt(8/3).
    The first ball-contact radius is r_c = rho * sqrt(2/3).

    Returns:
        Array of shape (4, 3).
    """
    signs = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)
    return rho * signs / np.sqrt(3.0)


def lobe_distances(
    vector_field: np.ndarray, centres: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Distance from each voxel's Gaussian-space point to its nearest centre.

    Args:
        vector_field: float array of shape (n_components, n, n, n).
        centres:      float array of shape (K, n_components).

    Returns:
        distances: float array of shape (n, n, n), min distance to any centre.
        nearest:   int8 array of shape (n, n, n), index of nearest centre.
    """
    spatial_shape = vector_field.shape[1:]
    best = np.full(spatial_shape, np.inf, dtype=np.float32)
    nearest = np.zeros(spatial_shape, dtype=np.int8)
    for label, centre in enumerate(centres):
        d2 = sum((vector_field[j] - centre[j]) ** 2 for j in range(len(centre)))
        update = d2 < best
        best[update] = d2[update]
        nearest[update] = label
    return np.sqrt(best, out=best), nearest


def lobe_excursion(
    vector_field: np.ndarray, centres: np.ndarray, radius: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return the ball-union excursion A_r = {min_a ||G(x) - c_a|| <= r}.

    Before the contact radius r_c = ||c_a - c_b|| / 2, different lobes are
    disjoint in Gaussian value space. At r_c they first touch; above r_c they
    overlap and inclusion-exclusion terms (not yet derived) enter the M_j
    formulas.

    Args:
        vector_field: float array of shape (n_components, n, n, n).
        centres:      float array of shape (K, n_components).
        radius:       Ball radius r in Gaussian-space units.

    Returns:
        phase:         bool array of shape (n, n, n).
        nearest_label: int8 array of shape (n, n, n).
    """
    distances, nearest = lobe_distances(vector_field, centres)
    return distances <= radius, nearest
