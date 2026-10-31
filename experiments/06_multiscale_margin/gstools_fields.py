"""Experiment-local GSTools sampling on explicit boundary-inclusive nodes."""
from __future__ import annotations

import gstools as gs
import numpy as np


def sample_components(n, length, classes, seed, realization=0, mode_no=1000,
                      threads=1, progress=None):
    """Return nodal components; fixed streams permit paired grid refinement.

    Matches src.fields' Gaussian(len_scale=lc, rescale=1, var=1) and RandMeth
    engine. Explicit coordinates include both cube faces, unlike its existing
    voxel-oriented wrapper. No empirical centering or variance normalization.
    """
    if n < 1 or length <= 0 or mode_no < 1 or threads < 1:
        raise ValueError("Invalid grid, mode count, or thread count.")
    gs.config.NUM_THREADS = threads
    axis = np.linspace(0.0, length, n + 1)
    groups, diagnostics = [], []
    for k, item in enumerate(classes):
        model = gs.Gaussian(dim=3, len_scale=item["lc"], rescale=1.0, var=1.0, nugget=0.0)
        fields = []
        for component in range(item["q"]):
            rng = np.random.default_rng(np.random.SeedSequence([seed, realization, k, component]))
            component_seed = int(rng.integers(0, 2**31 - 1))
            srf = gs.SRF(model, generator="RandMeth", mode_no=mode_no, seed=component_seed)
            field = np.asarray(srf.structured((axis, axis, axis), store=False), dtype=np.float64)
            fields.append(field)
            record = {"realization": realization+1, "class_index": k, "component": component,
                      "seed": component_seed, "mode_no": mode_no,
                      "mean": float(field.mean()), "second_moment": float(np.mean(field**2)),
                      "lambda2_nominal": 2/item["lc"]**2}
            for axis_index, axis_name in enumerate("xyz"):
                record["fd_energy_"+axis_name] = float(np.mean((np.diff(field, axis=axis_index)/(length/n))**2))
            diagnostics.append(record)
        groups.append(fields)
        if progress:
            progress(f"generated GSTools class {k+1}/{len(classes)}")
    return groups, diagnostics


def attach_components(mesh, groups):
    names = []
    for k, group in enumerate(groups):
        class_names = []
        for j, field in enumerate(group):
            name = f"g_{k}_{j}"
            mesh.point_data[name] = np.asarray(field).ravel(order="C")
            class_names.append(name)
        names.append(class_names)
    return names


def center_samples(groups):
    """Evaluate the same Freudenthal PL fields at cube centres.

    The cube body diagonal belongs to all six tetrahedra. Its midpoint value is
    the average of the two diagonal endpoints, before winner competition.
    """
    return [[0.5*(f[:-1, :-1, :-1] + f[1:, 1:, 1:]) for f in group] for group in groups]
