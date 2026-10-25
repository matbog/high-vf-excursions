"""Euler characteristic of voxel phases on a three-dimensional torus."""

from __future__ import annotations

import numpy as np
from itertools import combinations, product


def periodic_euler_nd(phase: np.ndarray) -> int:
    """Euler characteristic of a periodic cubical complex in dimensions 1--3."""
    p = np.asarray(phase, dtype=bool)
    dimension = p.ndim
    if dimension not in (1, 2, 3):
        raise ValueError("only dimensions 1, 2, and 3 are supported")

    cell_counts: list[int] = []
    all_axes = set(range(dimension))
    for cell_dimension in range(dimension + 1):
        count = 0
        for parallel_axes in combinations(range(dimension), cell_dimension):
            perpendicular_axes = sorted(all_axes - set(parallel_axes))
            incident = np.zeros_like(p)
            for offsets in product((0, 1), repeat=len(perpendicular_axes)):
                shifted = p
                for axis, offset in zip(perpendicular_axes, offsets):
                    if offset:
                        shifted = np.roll(shifted, 1, axis=axis)
                incident |= shifted
            count += int(np.count_nonzero(incident))
        cell_counts.append(count)
    return sum((-1) ** k * count for k, count in enumerate(cell_counts))


def periodic_euler_binary(phase: np.ndarray) -> int:
    """Euler characteristic of a periodic cubical complex formed by occupied voxels."""
    if np.asarray(phase).ndim != 3:
        raise ValueError("phase must be three-dimensional")
    return periodic_euler_nd(phase)


def periodic_euler_label_aware(phase: np.ndarray, winners: np.ndarray, q: int) -> int:
    """Add ECs of separated winner chambers."""
    return sum(periodic_euler_nd(phase & (winners == winner)) for winner in range(q))
