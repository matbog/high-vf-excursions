"""Empirical measurement of Lipschitz-Killing curvatures (LKCs) and morphology.

Quantities measured:
    L0  Euler characteristic — cubical complex formula (exact on the voxel grid).
    L1  Mean breadth — not measured here; requires a careful Crofton or
        signed-dihedral estimator and a resolution study before it can be used
        for Adler comparisons. Placeholder only.
    L2  Half surface area — marching-cubes estimate. Carries O(dx) resolution bias.
    L3  Volume — exact on the voxel grid.
"""

from __future__ import annotations

import numpy as np
from skimage import measure

from .fields import GridSpec


def measure_euler(phase: np.ndarray) -> float:
    """Euler characteristic via scikit-image cubical complex (6-connectivity)."""
    return float(measure.euler_number(phase, connectivity=1))


def measure_surface_area(phase: np.ndarray, dx: float) -> float:
    """Half-surface-area L2 via marching cubes.

    The phase is zero-padded before marching cubes so that the outer boundary
    of a phase that touches the domain face is properly closed.
    """
    padded = np.pad(phase.astype(np.float32), 1, constant_values=0.0)
    try:
        verts, faces, _, _ = measure.marching_cubes(padded, level=0.5, spacing=(dx, dx, dx))
    except ValueError:
        return 0.0
    tris = verts[faces]
    area = float(
        0.5 * np.linalg.norm(
            np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1
        ).sum()
    )
    return area / 2.0


def measure_percolation(phase: np.ndarray) -> dict[str, bool]:
    """Check whether any connected component spans from face 0 to face -1 along each axis.

    Returns:
        Dict with keys 'x', 'y', 'z' mapping to bool.
    """
    labels = measure.label(phase, connectivity=1)
    result: dict[str, bool] = {}
    for axis, key in enumerate(("x", "y", "z")):
        face0 = set(np.unique(np.take(labels, 0, axis=axis)))
        face1 = set(np.unique(np.take(labels, -1, axis=axis)))
        result[key] = bool(len((face0 & face1) - {0}) > 0)
    return result


def measure_all(phase: np.ndarray, grid: GridSpec) -> dict:
    """Measure all available morphological descriptors.

    Returns dict with keys:
        euler          float
        surface_area   float   (full surface area, not half)
        half_surface   float   (L2 = surface_area / 2)
        volume         float   (L3 in length^3 units)
        vf             float   (volume fraction)
        percolates_x   bool
        percolates_y   bool
        percolates_z   bool
        percolates_any bool
    """
    volume = float(phase.sum()) * grid.dx ** 3
    half_surf = measure_surface_area(phase, grid.dx)
    euler = measure_euler(phase)
    perc = measure_percolation(phase)
    return {
        "euler": euler,
        "surface_area": half_surf * 2.0,
        "half_surface": half_surf,
        "volume": volume,
        "vf": volume / grid.volume,
        "percolates_x": perc["x"],
        "percolates_y": perc["y"],
        "percolates_z": perc["z"],
        "percolates_any": perc["x"] or perc["y"] or perc["z"],
    }
