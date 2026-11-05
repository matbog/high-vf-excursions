"""Distinct voxel-union, ball-preimage and nearest-centre measurements.

The union scalar is r - min_a d_a at the sample points. Marching cubes of
that sampled scalar approximates its interface; it does not reconstruct the
original Gaussian components or detect every sub-grid ball crossing.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from src.fields import GridSpec
from src.measures import (measure_euler, measure_l1_from_field,
                          measure_surface_area_crofton,
                          measure_surface_area_from_field)


def measure_lobes(distances: np.ndarray, radius: float, grid: GridSpec,
                  *, include_mesh: bool = True) -> dict:
    """Measure ball distances without assuming disjointness.

    include_mesh=False skips costly mesh measurements and records them as NaN;
    voxel volume, topology, Crofton area and partition diagnostics are retained.
    The default preserves the full union_v2 measurement policy.
    """
    distances = np.asarray(distances)
    if (distances.ndim != 4 or distances.shape[0] < 1
            or distances.shape[1:] != (grid.n,) * 3
            or not np.all(np.isfinite(distances)) or np.any(distances < 0)
            or not np.isfinite(radius) or radius <= 0):
        raise ValueError("Require nonnegative finite ball distances and positive radius")
    nearest = distances.argmin(axis=0)
    union_signed = radius - distances.min(axis=0)
    phase = union_signed >= 0
    occupied = int(phase.sum())
    result = {
        "measured_phi": occupied / phase.size,
        "union_L3": occupied * grid.dx**3,
        "union_L1_field": measure_l1_from_field(union_signed, grid.dx) if include_mesh else float("nan"),
        "union_L2_crofton": measure_surface_area_crofton(phase, grid.dx),
        "union_L2_field": measure_surface_area_from_field(union_signed, grid.dx) if include_mesh else float("nan"),
    }
    for rank, connectivity in ((1, 6), (3, 26)):
        labels, count = ndimage.label(phase, ndimage.generate_binary_structure(3, rank))
        sizes = np.bincount(labels.ravel())
        spans = []
        for axis, name in enumerate("xyz"):
            common = np.intersect1d(np.take(labels, 0, axis=axis), np.take(labels, -1, axis=axis))
            spans.append(bool(np.any(common > 0)))
            result[f"union_spans_{connectivity}_{name}"] = spans[-1]
        result.update({
            f"union_euler_{connectivity}": int(measure_euler(phase, rank)),
            f"union_components_{connectivity}": int(count),
            f"union_largest_fraction_{connectivity}": float(sizes[1:].max(initial=0) / occupied) if occupied else 0.,
            f"union_spans_{connectivity}": any(spans),
        })
    sums = {key: 0. for key in (
        "sum_ball_euler_26", "sum_ball_L1_field", "sum_ball_L2_crofton",
        "sum_ball_L2_field", "sum_ball_L3", "partition_euler_26", "partition_L2_crofton")}
    if not include_mesh:
        sums["sum_ball_L1_field"] = sums["sum_ball_L2_field"] = float("nan")
    for i, distance in enumerate(distances):
        signed = radius - distance
        ball = signed >= 0
        partition = phase & (nearest == i)
        sums["sum_ball_euler_26"] += measure_euler(ball, 3)
        if include_mesh:
            sums["sum_ball_L1_field"] += measure_l1_from_field(signed, grid.dx)
        sums["sum_ball_L2_crofton"] += measure_surface_area_crofton(ball, grid.dx)
        if include_mesh:
            sums["sum_ball_L2_field"] += measure_surface_area_from_field(signed, grid.dx)
        sums["sum_ball_L3"] += ball.sum() * grid.dx**3
        sums["partition_euler_26"] += measure_euler(partition, 3)
        sums["partition_L2_crofton"] += measure_surface_area_crofton(partition, grid.dx)
    return result | sums
