"""Piecewise-affine margin regions on a shared tetrahedral mesh.

Each original Gaussian component is affine inside an original tetrahedron.
A winner region is clipped by every G_i-G_j >= delta inequality. Clipping the
sampled top-two gap is deliberately never used. Class unions are evaluated by
inclusion-exclusion of their disjoint winner regions, including all cross-class
intersections. No surface-Boolean approximation is needed for L0 or L3.
"""
from __future__ import annotations

from itertools import permutations

import numpy as np
import pyvista as pv
import vtk


def tetrahedral_grid(n: int, length: float) -> pv.UnstructuredGrid:
    """n cells and n+1 boundary-inclusive nodes per axis; six tets per cube."""
    if n < 1 or not np.isfinite(length) or length <= 0:
        raise ValueError("Require n >= 1 and finite length > 0.")
    axis = np.linspace(0.0, length, n + 1)
    points = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    origins = np.indices((n, n, n)).reshape(3, -1).T
    stride = np.array([(n + 1)**2, n + 1, 1])
    base = origins @ stride
    cells = []
    for order in permutations(range(3)):
        offsets = np.r_[0, np.cumsum(stride[list(order)])]
        tet = base[:, None] + offsets
        # Positive orientation, required for consistently outward boundary faces.
        xyz = points[tet[0]]
        if np.linalg.det((xyz[1:] - xyz[0]).T) < 0:
            tet[:, [1, 2]] = tet[:, [2, 1]]
        cells.append(tet)
    cells = np.concatenate(cells)
    packed = np.column_stack([np.full(len(cells), 4), cells])
    mesh = pv.UnstructuredGrid(packed, np.full(len(cells), pv.CellType.TETRA, np.uint8), points)
    mesh.cell_data["parent_tet"] = np.arange(mesh.n_cells, dtype=np.int64)
    return mesh


def clip_linear(mesh: pv.UnstructuredGrid, values: np.ndarray, level: float = 0.0):
    """Clip an affine nodal quantity with double-precision VTK geometry."""
    if mesh.n_cells == 0:
        return mesh
    values = np.asarray(values, dtype=np.float64) - level
    if values.min() >= 0:
        return mesh
    if values.max() <= 0:
        # Zero-volume tangencies are excluded; production thresholds must be regular.
        return pv.UnstructuredGrid()
    work = mesh.copy(deep=False)
    work.point_data["_clip_value"] = values
    work.set_active_scalars("_clip_value", preference="point")
    clip = vtk.vtkTableBasedClipDataSet()
    clip.SetInputData(work)
    clip.SetValue(0.0)
    clip.SetInsideOut(False)
    clip.SetMergeTolerance(1e-12)
    clip.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION)
    clip.Update()
    triangulate = vtk.vtkDataSetTriangleFilter()
    triangulate.SetInputConnection(clip.GetOutputPort())
    triangulate.SetTetrahedraOnly(True)
    triangulate.Update()
    result = pv.wrap(triangulate.GetOutput()).copy()
    if "_clip_value" in result.point_data:
        del result.point_data["_clip_value"]
    return result


def clip_winner(mesh, names: list[str], winner: int, delta: float):
    if not np.isfinite(delta) or delta <= 0:
        raise ValueError("Positive margins are required for disjoint winner regions.")
    result = mesh
    for j, name in enumerate(names):
        if j == winner:
            continue
        result = clip_linear(result, result.point_data[names[winner]] - result.point_data[name], delta)
        if result.n_cells == 0:
            break
    return result


def _unique_rows(array):
    return np.unique(np.sort(array, axis=1), axis=0)


def measure_region(mesh, length: float, tolerance: float = 0.0):
    """Measure a clipped volumetric complex; reject broken boundary topology.

    Zero tolerance merges exactly coincident vertices only. A positive fraction
    of domain length is optional, but can collapse valid thin regions. Volume
    and Euler use the same tetrahedral geometry.
    """
    if mesh.n_cells == 0:
        return {"L0": 0, "L3": 0.0, "components": 0}, mesh, np.empty(0)
    mesh = mesh.clean(tolerance=tolerance * length, remove_unused_points=True)
    if np.any(mesh.celltypes != pv.CellType.TETRA):
        raise RuntimeError("Expected tetrahedra after clipping.")
    tets = mesh.cells.reshape(-1, 5)[:, 1:]
    xyz = mesh.points[tets]
    volumes = np.abs(np.linalg.det(xyz[:, 1:] - xyz[:, :1])) / 6.0
    if np.any(volumes <= 0):
        # Clipping a vertex/edge can emit exactly flat tetrahedra. They do not
        # belong to the regular three-dimensional solid. Keep positive-volume
        # slivers; the boundary audit below must still pass after removing flats.
        mesh = mesh.extract_cells(volumes > 0).clean(tolerance=tolerance*length,
                                                     remove_unused_points=True)
        if mesh.n_cells == 0:
            return {"L0": 0, "L3": 0.0, "components": 0}, mesh, np.empty(0)
        tets = mesh.cells.reshape(-1, 5)[:, 1:]
        xyz = mesh.points[tets]
        volumes = np.abs(np.linalg.det(xyz[:, 1:] - xyz[:, :1])) / 6.0
    faces = np.concatenate([tets[:, f] for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))])
    face_keys, multiplicity = np.unique(np.sort(faces, axis=1), axis=0, return_counts=True)
    if np.any(multiplicity > 2):
        raise RuntimeError("Nonmanifold tetrahedral complex.")
    edges = _unique_rows(np.concatenate([tets[:, e] for e in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))]))
    euler = mesh.n_points - len(edges) + len(face_keys) - len(tets)
    boundary = face_keys[multiplicity == 1]
    boundary_edges, counts = np.unique(np.sort(np.concatenate(
        [boundary[:, e] for e in ((0, 1), (0, 2), (1, 2))]), axis=1), axis=0, return_counts=True)
    if np.any(counts != 2):
        raise RuntimeError("Clipped mesh has a nonmanifold/open boundary; refusing an Euler measurement.")
    boundary_euler = len(np.unique(boundary)) - len(boundary_edges) + len(boundary)
    if boundary_euler != 2 * euler:
        raise RuntimeError("Boundary/volume Euler consistency failed.")
    mesh.cell_data["piece_volume"] = volumes
    mesh = mesh.connectivity(label_regions=True)
    regions = np.asarray(mesh.cell_data["RegionId"])
    return {"L0": int(euler), "L3": float(volumes.sum()),
            "components": int(regions.max()) + 1}, mesh, np.asarray(mesh.cell_data["piece_volume"])


class DisjointSets:
    def __init__(self):
        self.parents = []

    def add(self):
        i = len(self.parents)
        self.parents.append(i)
        return i

    def find(self, i):
        while self.parents[i] != i:
            self.parents[i] = self.parents[self.parents[i]]
            i = self.parents[i]
        return i

    def join(self, a, b):
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parents[b] = a


def measure_union(base, class_names, deltas, length, *, tolerance=0.0,
                  export_directory=None, progress=None):
    """Exact inclusion-exclusion for the PL union, up to clipping roundoff.

    Raw connected regions are graph vertices. Positive-volume intersections
    connect them. Tangency-only connections at singular thresholds are outside
    this regular-threshold implementation. Exported raw volumes overlap across
    classes: they are a representation of the set union, not a disjoint mesh.
    """
    if not 1 <= len(class_names) <= 3 or len(deltas) != len(class_names):
        raise ValueError("Use one to three classes with matching margins.")
    sets = DisjointSets()
    lookup = {}  # (class, winner) -> original tetrahedron -> raw component
    boundary_flags = []
    signed_volumes = []
    terms = []
    raw_records = []
    totals = np.zeros(2)
    added_volumes = np.zeros(len(class_names))

    def register_raw(mesh, k, winner):
        local_to_global = {}
        cell_regions = np.asarray(mesh.cell_data["RegionId"])
        point_regions = np.asarray(mesh.point_data["RegionId"])
        for region in np.unique(cell_regions):
            global_id = sets.add()
            local_to_global[int(region)] = global_id
            points = mesh.points[point_regions == region]
            flags = []
            for axis in range(3):
                face_eps = max(tolerance, 32*np.finfo(float).eps)*length
                flags.extend([bool(np.any(points[:, axis] <= face_eps)),
                              bool(np.any(points[:, axis] >= length-face_eps))])
            boundary_flags.append(np.asarray(flags))
        parents = np.asarray(mesh.cell_data["parent_tet"])
        mapping = {}
        for parent, region in zip(parents, cell_regions):
            value = local_to_global[int(region)]
            if int(parent) in mapping and mapping[int(parent)] != value:
                raise RuntimeError("One convex winner piece in a parent tetrahedron split into components.")
            mapping[int(parent)] = value
        lookup[k, winner] = mapping

    def record(mesh, chosen):
        try:
            stats, mesh, volumes = measure_region(mesh, length, tolerance)
        except RuntimeError as error:
            if export_directory is not None:
                mesh.save(export_directory / "failed_mesh.vtu")
            raise RuntimeError(f"Intersection {chosen}: {error}") from error
        if not mesh.n_cells:
            return mesh
        sign = 1 if len(chosen) % 2 else -1
        totals[:] += sign * np.array([stats["L0"], stats["L3"]])
        added_volumes[chosen[-1][0]] += sign * stats["L3"]
        terms.append({"order": len(chosen), "classes": ",".join(str(k) for k, _ in chosen),
                      "winners": ",".join(str(w) for _, w in chosen), "sign": sign, **stats})
        if len(chosen) == 1:
            k, winner = chosen[0]
            register_raw(mesh, k, winner)
            raw_records.append({"class_index": k, "winner": winner, **stats})
            if export_directory is not None:
                exported = mesh.copy()
                exported.clear_data()
                exported.cell_data["class_index"] = np.full(mesh.n_cells, k)
                exported.cell_data["winner"] = np.full(mesh.n_cells, winner)
                exported.save(export_directory / f"class_{k}_winner_{winner}.vtu")
        primary = lookup[chosen[0]]
        parent_ids = np.asarray(mesh.cell_data["parent_tet"])
        for parent in np.unique(parent_ids):
            node = primary[int(parent)]
            for key in chosen[1:]:
                if int(parent) not in lookup[key]:
                    raise RuntimeError("Intersection has no corresponding raw region.")
                sets.join(node, lookup[key][int(parent)])
        node_ids = np.array([primary[int(p)] for p in parent_ids])
        for node in np.unique(node_ids):
            signed_volumes.append((int(node), sign * float(volumes[node_ids == node].sum())))
        return mesh

    # Register every raw region before intersection recursion uses its mapping.
    raw_meshes = {}
    for k, names in enumerate(class_names):
        for winner in range(len(names)):
            mesh = clip_winner(base, names, winner, deltas[k])
            if mesh.n_cells:
                raw_meshes[k, winner] = record(mesh, [(k, winner)])
        if progress:
            progress(f"clipped raw class {k+1}/{len(class_names)}")

    def intersections(mesh, chosen):
        for k in range(chosen[-1][0] + 1, len(class_names)):
            for winner in range(len(class_names[k])):
                if (k, winner) not in raw_meshes:
                    continue
                clipped = clip_winner(mesh, class_names[k], winner, deltas[k])
                if clipped.n_cells:
                    clipped = record(clipped, chosen + [(k, winner)])
                    if clipped.n_cells:
                        intersections(clipped, chosen + [(k, winner)])

    for key, mesh in raw_meshes.items():
        intersections(mesh, [key])
    components = {}
    for i, flags in enumerate(boundary_flags):
        root = sets.find(i)
        if root not in components:
            components[root] = {"volume": 0.0, "flags": np.zeros(6, bool), "raw_regions": 0}
        components[root]["flags"] |= flags
        components[root]["raw_regions"] += 1
    for node, volume in signed_volumes:
        components[sets.find(node)]["volume"] += volume
    component_rows = []
    for i, item in enumerate(components.values()):
        flags = item.pop("flags")
        if item["volume"] <= 0:
            raise RuntimeError("Nonpositive union component volume after inclusion-exclusion.")
        component_rows.append({"component": i, **item,
                               "touches_boundary": bool(flags.any()),
                               "spans": bool(np.any(flags[::2] & flags[1::2])),
                               "equivalent_diameter": float((6*item["volume"]/np.pi)**(1/3))})
    volume = float(totals[1])
    if not np.isclose(sum(r["volume"] for r in component_rows), volume, atol=1e-9*length**3, rtol=1e-9):
        raise RuntimeError("Component and inclusion-exclusion volumes disagree.")
    if volume < -1e-9*length**3 or volume > (1+1e-9)*length**3:
        raise RuntimeError("Union volume is outside the domain volume.")
    summary = {"L0": int(round(totals[0])), "L3": volume, "phi": volume / length**3,
               "components": len(component_rows), "spans": any(r["spans"] for r in component_rows),
               "largest_fraction": max((r["volume"] for r in component_rows), default=0.0)/volume if volume else 0.0,
               "nonempty_intersection_terms": len(terms)}
    return summary, component_rows, terms, raw_records, added_volumes
