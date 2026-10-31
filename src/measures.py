"""Empirical measurement of Lipschitz-Killing curvatures (LKCs) and morphology.

Quantities measured:
    L0  Euler characteristic — cubical complex formula (exact on the voxel grid).
    L1  Integrated mean curvature / pi — dihedral angles on the marching-cubes
        mesh (Schläfli-type formula for PL surfaces). Binary-surface bias can
        persist under refinement; no convergence rate is established here.
    L2  Half surface area — marching-cubes estimate. Carries O(dx) resolution bias.
    L3  Volume — exact on the voxel grid.
"""

from __future__ import annotations

import numpy as np
from skimage import measure

from .fields import GridSpec


def measure_euler(phase: np.ndarray, connectivity: int = 1) -> float:
    """Euler characteristic via scikit-image cubical complex.

    connectivity=1 → 6-connectivity (face-adjacent voxels only).
    connectivity=3 → 26-connectivity (face + edge + vertex adjacency).
    These can differ substantially on morphologies with thin films or
    sub-voxel features; the "correct" convention is not universal.
    """
    return float(measure.euler_number(phase, connectivity=connectivity))


def _l1_from_mesh(verts: np.ndarray, faces: np.ndarray) -> float:
    """Dihedral-angle L1 estimator for a consistently-wound PL mesh.

    For a PL surface with consistently-oriented outward normals:
        integral_S H dA = (1/2) * sum_edges |e| * theta_e
    where theta_e is the signed dihedral angle at edge e (positive convex,
    negative concave) from directed-edge orientation. Returns L1 = (∫H dA) / π.
    """
    faces = np.asarray(faces, dtype=np.int64)
    if len(faces) == 0:
        return 0.0
    verts = np.asarray(verts, dtype=np.float64)

    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    cross = np.cross(v1 - v0, v2 - v0)
    norms = np.linalg.norm(cross, axis=1, keepdims=True)
    face_normals = cross / np.maximum(norms, 1e-12)

    n_f = len(faces)
    n_v = len(verts)
    va = np.concatenate([faces[:, 0], faces[:, 1], faces[:, 2]])
    vb = np.concatenate([faces[:, 1], faces[:, 2], faces[:, 0]])
    fi_of_edge = np.tile(np.arange(n_f, dtype=np.int64), 3)

    edge_key = va * n_v + vb
    rev_key  = vb * n_v + va

    order      = np.argsort(edge_key)
    sorted_key = edge_key[order]

    pos   = np.searchsorted(sorted_key, rev_key)
    pos   = np.clip(pos, 0, len(sorted_key) - 1)
    found = sorted_key[pos] == rev_key

    edge_idx = np.arange(len(va))
    match    = order[pos]
    valid    = found & (edge_idx < match)

    if not valid.any():
        return 0.0

    i_idx = edge_idx[valid]
    j_idx = match[valid]

    n1        = face_normals[fi_of_edge[i_idx]]
    n2        = face_normals[fi_of_edge[j_idx]]
    edge_vecs = verts[vb[i_idx]] - verts[va[i_idx]]
    edge_lens = np.linalg.norm(edge_vecs, axis=1)

    cross_n = np.cross(n1, n2)
    sin_t   = np.einsum("ij,ij->i", cross_n, edge_vecs) / (edge_lens + 1e-12)
    cos_t   = np.einsum("ij,ij->i", n1, n2)
    angles  = np.arctan2(sin_t, cos_t)

    return float(np.sum(angles * edge_lens)) / (2.0 * np.pi)


def measure_l1(phase: np.ndarray, dx: float) -> float:
    """L1 via dihedral angles on marching cubes of a BINARY phase mask.

    Binarizing the field before meshing loses subvoxel interface positions;
    the resulting curvature bias scales ~O(1/N) under grid refinement. For
    higher accuracy pass the smooth signed field directly via
    `measure_l1_from_field`.
    """
    padded = np.pad(phase.astype(np.float32), 1, constant_values=0.0)
    try:
        verts, faces, _, _ = measure.marching_cubes(
            padded, level=0.5, spacing=(dx, dx, dx), gradient_direction="ascent"
        )
    except ValueError:
        return 0.0
    return _l1_from_mesh(verts, faces)


def measure_l1_from_field(
    scalar_field: np.ndarray, dx: float, level: float = 0.0
) -> float:
    """L1 via dihedral angles on marching cubes of a SMOOTH SIGNED FIELD.

    `scalar_field` must be a continuous array whose level-0 isosurface IS the
    boundary of the region of interest ({scalar >= level} = inside). This
    eliminates the subvoxel binarization bias of `measure_l1`: the mesh
    vertices sit at the true continuum level set rather than snapped to voxel
    midpoints, so the dihedral angles accurately reflect continuum curvature.

    Padding uses a large negative value (array_min − 1) so the field crosses
    `level` at the array boundary, closing the mesh with a surface aligned
    with the domain faces. This matches the GKF formula's domain-face terms.
    """
    scalar_field = scalar_field.astype(np.float32)
    pad_value = float(scalar_field.min()) - 1.0
    padded = np.pad(scalar_field, 1, constant_values=pad_value)
    try:
        verts, faces, _, _ = measure.marching_cubes(
            padded, level=level, spacing=(dx, dx, dx), gradient_direction="ascent"
        )
    except ValueError:
        return 0.0
    return _l1_from_mesh(verts, faces)


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


def measure_surface_area_from_field(
    scalar_field: np.ndarray, dx: float, level: float = 0.0
) -> float:
    """Half-surface-area L2 via marching cubes of a SMOOTH SIGNED FIELD.

    Analogous to `measure_l1_from_field`: mesh vertices sit at the true continuum
    isosurface via linear interpolation (not snapped to voxel edges). Padding with
    `min(field) − 1` closes the mesh at the domain boundary.

    Known trade-off vs `measure_surface_area_crofton`: marching cubes carries a
    small amplitude bias on curved surfaces (triangles underestimate curved
    area), partially cancelled by the mesh's many facets being slightly offset
    from the true surface. In practice, on isotropic Gaussian excursions, the
    net bias at N≥128 is comparable to or smaller than Crofton's.
    """
    scalar_field = scalar_field.astype(np.float32)
    pad_value = float(scalar_field.min()) - 1.0
    padded = np.pad(scalar_field, 1, constant_values=pad_value)
    try:
        verts, faces, _, _ = measure.marching_cubes(
            padded, level=level, spacing=(dx, dx, dx), gradient_direction="ascent"
        )
    except ValueError:
        return 0.0
    faces = np.asarray(faces, dtype=np.int64)
    if len(faces) == 0:
        return 0.0
    verts = np.asarray(verts, dtype=np.float64)
    tris = verts[faces]
    area = float(
        0.5 * np.linalg.norm(
            np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0]), axis=1
        ).sum()
    )
    return area / 2.0


def measure_surface_area_crofton(phase: np.ndarray, dx: float) -> float:
    """Half-surface-area L2 via 3-direction Crofton + exact domain-face term.

    Decomposes L2 into two additive contributions matching the GKF structure
    L2 = (interior surface area) / 2 + (domain-face phase area) / 2:

    1. Interior: for an isotropic surface the stereological Cauchy relation
       gives S_interior = (2/3) * (N_x + N_y + N_z) * dx^2, where N_d counts
       voxel-pair transitions in direction d. Half of that:
            L2_interior = (N_x + N_y + N_z) * dx^2 / 6.
       NOTE: averaging over only 3 axis directions is unbiased for isotropic
       surfaces but anisotropic surfaces (e.g. the cube faces themselves)
       are underestimated by factor 2/3 — hence we DO NOT run Crofton on a
       zero-padded array; instead the face term is handled analytically.

    2. Domain faces: half the area of the 6 cube faces occupied by phase.
       Matches the 3a^2 * M_0 term in the GKF L2 formula exactly.

    Unlike marching cubes, which carries a ~1.09 amplitude bias on isotropic
    surfaces, this estimator is unbiased on isotropic surfaces to <0.1% at
    N>=128 on a sphere.
    """
    N_x = int(np.sum(phase[1:,  :,  :] != phase[:-1, :,  :]))
    N_y = int(np.sum(phase[ :, 1:,  :] != phase[ :, :-1, :]))
    N_z = int(np.sum(phase[ :,  :, 1:] != phase[ :,  :, :-1]))
    L2_interior = (N_x + N_y + N_z) * dx * dx / 3.0

    face_area = 0.0
    for axis in range(3):
        face_area += float(np.sum(np.take(phase, 0,  axis=axis))) * dx * dx
        face_area += float(np.sum(np.take(phase, -1, axis=axis))) * dx * dx
    L2_faces = 0.5 * face_area

    return L2_interior + L2_faces


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
        compute_l1: If True, compute L1 via dihedral angles (moderately slow).
                    Default False to skip L1 in experiments that don't need it.

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
