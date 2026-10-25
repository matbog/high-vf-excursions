"""Lipschitz--Killing curvatures of bounded 3D cubical voxel sets."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
import pyvista as pv
from skimage import measure


@dataclass(frozen=True)
class CubicalLKCs:
    L0: float
    L1: float
    L2: float
    L3: float

    def as_array(self) -> np.ndarray:
        return np.array([self.L0, self.L1, self.L2, self.L3], dtype=float)


def cubical_cell_counts(phase: np.ndarray) -> tuple[int, int, int, int]:
    """Return unique vertex, edge, face and cube counts for a bounded voxel set."""
    p = np.asarray(phase, dtype=bool)
    if p.ndim != 3:
        raise ValueError("phase must be three-dimensional")

    n3 = int(np.count_nonzero(p))

    n2 = 0
    for axis in range(3):
        first = np.take(p, 0, axis=axis)
        last = np.take(p, -1, axis=axis)
        left = np.take(p, np.arange(p.shape[axis] - 1), axis=axis)
        right = np.take(p, np.arange(1, p.shape[axis]), axis=axis)
        n2 += int(np.count_nonzero(first)) + int(np.count_nonzero(last))
        n2 += int(np.count_nonzero(left | right))

    n1 = 0
    for parallel_axis in range(3):
        perpendicular = [axis for axis in range(3) if axis != parallel_axis]
        padding = [(0, 0)] * 3
        for axis in perpendicular:
            padding[axis] = (1, 1)
        padded = np.pad(p, padding, mode="constant", constant_values=False)
        incident_shape = list(p.shape)
        for axis in perpendicular:
            incident_shape[axis] += 1
        incident = np.zeros(incident_shape, dtype=bool)
        for offsets in product((0, 1), repeat=2):
            slices = [slice(None)] * 3
            for axis, offset in zip(perpendicular, offsets):
                slices[axis] = slice(offset, offset + incident_shape[axis])
            incident |= padded[tuple(slices)]
        n1 += int(np.count_nonzero(incident))

    padded = np.pad(p, 1, mode="constant", constant_values=False)
    incident_vertices = np.zeros(tuple(size + 1 for size in p.shape), dtype=bool)
    for offsets in product((0, 1), repeat=3):
        slices = tuple(slice(offset, offset + size + 1) for offset, size in zip(offsets, p.shape))
        incident_vertices |= padded[slices]
    n0 = int(np.count_nonzero(incident_vertices))
    return n0, n1, n2, n3


def cubical_lkcs(phase: np.ndarray, dx: float = 1.0) -> CubicalLKCs:
    """Intrinsic volumes in the convention used by Roubin's System 2.17."""
    n0, n1, n2, n3 = cubical_cell_counts(phase)
    return CubicalLKCs(
        L0=float(n0 - n1 + n2 - n3),
        L1=float((n1 - 2 * n2 + 3 * n3) * dx),
        L2=float((n2 - 3 * n3) * dx**2),
        L3=float(n3 * dx**3),
    )


def label_aware_lkcs(
    phase: np.ndarray, winners: np.ndarray, q: int, dx: float = 1.0
) -> CubicalLKCs:
    """Add LKCs of mutually separated winner chambers."""
    total = np.zeros(4, dtype=float)
    for label in range(q):
        total += cubical_lkcs(phase & (winners == label), dx).as_array()
    return CubicalLKCs(*total)


def mesh_l1_l2(phase: np.ndarray, dx: float = 1.0) -> tuple[float, float]:
    """Estimate smooth L1 and L2 from a closed marching-cubes surface.

    L2 is half surface area. L1 is the integrated mean curvature divided by pi.
    """
    p = np.asarray(phase, dtype=np.float32)
    if not np.any(p):
        return 0.0, 0.0
    padded = np.pad(p, 1, mode="constant", constant_values=0.0)
    vertices, faces, _, _ = measure.marching_cubes(
        padded, level=0.5, spacing=(dx, dx, dx), gradient_direction="descent"
    )
    triangles = vertices[faces]
    triangle_areas = 0.5 * np.linalg.norm(
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]), axis=1
    )
    l2 = float(triangle_areas.sum() / 2.0)

    pv_faces = np.column_stack((np.full(len(faces), 3), faces)).astype(np.int64).ravel()
    surface = pv.PolyData(vertices, pv_faces)
    mean_curvature = np.asarray(surface.curvature(curv_type="mean"))
    vertex_areas = np.zeros(len(vertices), dtype=float)
    for corner in range(3):
        np.add.at(vertex_areas, faces[:, corner], triangle_areas / 3.0)
    integral_mean_curvature = float(np.sum(mean_curvature * vertex_areas))
    l1 = abs(integral_mean_curvature) / np.pi
    return l1, l2


def label_aware_mesh_lkcs(
    phase: np.ndarray, winners: np.ndarray, q: int, dx: float = 1.0
) -> CubicalLKCs:
    """Cubical L0/L3 with isotropic mesh estimates for L1/L2."""
    cubical = label_aware_lkcs(phase, winners, q, dx)
    l1 = 0.0
    l2 = 0.0
    for label in range(q):
        component_l1, component_l2 = mesh_l1_l2(phase & (winners == label), dx)
        l1 += component_l1
        l2 += component_l2
    return CubicalLKCs(cubical.L0, l1, l2, cubical.L3)
