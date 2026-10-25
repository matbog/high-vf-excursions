"""Empirical measurement of Lipschitz-Killing curvatures (LKCs) and morphology.

Quantities measured:
    L0  Euler characteristic — cubical complex formula (exact on the voxel grid).
    L1  Integrated mean curvature / pi — marching-cubes + PyVista mean-curvature
        estimate. Carries O(dx) resolution bias; treat as indicative, not exact.
    L2  Half surface area — marching-cubes estimate. Carries O(dx) resolution bias.
    L3  Volume — exact on the voxel grid.
"""

from __future__ import annotations

import numpy as np
import pyvista as pv
from skimage import measure

from .fields import GridSpec


def measure_euler(phase: np.ndarray) -> float:
    """Euler characteristic via scikit-image cubical complex (6-connectivity)."""
    return float(measure.euler_number(phase, connectivity=1))


def measure_l1(phase: np.ndarray, dx: float) -> float:
    """Integrated mean curvature / pi via marching-cubes + PyVista mean curvature.

    L1 = (1/pi) * integral_surface H dA  where H = (kappa1+kappa2)/2.
    Carries O(dx) resolution bias from both the surface discretisation and the
    mean-curvature estimator. Treat as indicative; validate resolution dependence
    before using for quantitative Adler comparisons.
    """
    padded = np.pad(phase.astype(np.float32), 1, constant_values=0.0)
    try:
        verts, faces, _, _ = measure.marching_cubes(padded, level=0.5, spacing=(dx, dx, dx))
    except ValueError:
        return 0.0
    if len(faces) == 0:
        return 0.0
    tris = verts[faces]
    triangle_areas = 0.5 * np.linalg.norm(
        np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1
    )
    pv_faces = np.column_stack((np.full(len(faces), 3), faces)).astype(np.int64).ravel()
    surface = pv.PolyData(verts, pv_faces)
    mean_curv = np.asarray(surface.curvature(curv_type="mean"))
    vertex_areas = np.zeros(len(verts), dtype=float)
    for corner in range(3):
        np.add.at(vertex_areas, faces[:, corner], triangle_areas / 3.0)
    return float(np.sum(mean_curv * vertex_areas)) / np.pi


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


def measure_percolation_per_winner(
    phase: np.ndarray, winners: np.ndarray, q: int
) -> bool:
    """Check whether any single winner-type cluster spans the domain.

    Used for the margin tessellation model, where 'spanning' means a cluster
    of a single winner type (not the union of all types) forms a face-to-face
    connected path.  This is the correct metric for the 'disconnected inclusion'
    claim: the inclusion phase is disconnected if no individual winner type spans.

    Args:
        phase:   Binary inclusion mask (True = included voxel).
        winners: Integer array of the same shape; value in [0, q-1] is the
                 dominant field index at each voxel (ignored outside phase).
        q:       Number of winner types.

    Returns:
        True if any single winner-type connected component spans in any axis.
    """
    for label in range(q):
        single_type = phase & (winners == label)
        if not np.any(single_type):
            continue
        perc = measure_percolation(single_type)
        if perc["x"] or perc["y"] or perc["z"]:
            return True
    return False


def measure_all(phase: np.ndarray, grid: GridSpec, compute_l1: bool = False) -> dict:
    """Measure all available morphological descriptors.

    Args:
        phase:      Binary inclusion array.
        grid:       Grid specification (for physical units).
        compute_l1: If True, compute L1 via PyVista mean curvature (slow).
                    Default False to avoid the PyVista overhead in experiments
                    that don't need L1.

    Returns dict with keys:
        euler          float
        l1             float   (only if compute_l1=True; 0.0 otherwise)
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
    l1 = measure_l1(phase, grid.dx) if compute_l1 else 0.0
    euler = measure_euler(phase)
    perc = measure_percolation(phase)
    return {
        "euler": euler,
        "l1": l1,
        "surface_area": half_surf * 2.0,
        "half_surface": half_surf,
        "volume": volume,
        "vf": volume / grid.volume,
        "percolates_x": perc["x"],
        "percolates_y": perc["y"],
        "percolates_z": perc["z"],
        "percolates_any": perc["x"] or perc["y"] or perc["z"],
    }
