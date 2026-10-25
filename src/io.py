"""ParaView-readable VTI and VTP export."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyvista as pv
from skimage import measure

from .fields import GridSpec


def save_vti(
    path: str | Path,
    data: np.ndarray,
    grid: GridSpec,
    array_name: str,
) -> None:
    """Save a 3D array as a ParaView VTI ImageData file.

    Args:
        path:       Output file path (should end in .vti).
        data:       3D array of shape (n, n, n).
        grid:       Grid specification.
        array_name: Name of the data array in ParaView.
    """
    arr = np.asarray(data)
    if arr.ndim != 3:
        raise ValueError("data must be 3D")
    nx, ny, nz = arr.shape
    image = pv.ImageData(
        dimensions=(nx, ny, nz),
        spacing=(grid.dx, grid.dx, grid.dx),
        origin=(0.0, 0.0, 0.0),
    )
    image.point_data[array_name] = arr.flatten(order="F")
    image.save(str(path))


def save_vtp(path: str | Path, phase: np.ndarray, grid: GridSpec) -> None:
    """Save a marching-cubes surface of a binary phase as a ParaView VTP file.

    Args:
        path:  Output file path (should end in .vtp).
        phase: Boolean 3D array.
        grid:  Grid specification.
    """
    padded = np.pad(phase.astype(np.float32), 1, constant_values=0.0)
    verts, faces, _, _ = measure.marching_cubes(
        padded, level=0.5, spacing=(grid.dx, grid.dx, grid.dx)
    )
    pv_faces = np.hstack([np.full((len(faces), 1), 3), faces]).astype(np.int64).ravel()
    mesh = pv.PolyData(verts, pv_faces)
    mesh.save(str(path))
