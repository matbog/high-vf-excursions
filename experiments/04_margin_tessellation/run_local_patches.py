"""Small paired 3D reconstruction study of the saved q=32 Gaussian realization.

All reconstruction and direct-reference topology measurements remain digital.
Fixed patches are two original grid cells wide; no full-domain meshing is done.
"""
from __future__ import annotations

import argparse
import csv
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
import numpy as np
import pyvista as pv
from scipy import ndimage

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_local_films as films
saved = films.saved

SCHEMA = "local_margin_patches_v1"
METHODS = ("direct", "trilinear_components", "tetra_components", "trilinear_margin")
CORNERS = np.array(list(np.ndindex(2, 2, 2)))
COARSE_POINTS = np.array(list(np.ndindex(3, 3, 3)))


def read_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def grid_points(refinement):
    axis = np.linspace(0., 2., 2*refinement+1)
    return np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)


def interpolate(coarse, points, kind):
    """Trilinear or conforming Freudenthal tetrahedral interpolation.

    The tetrahedra in each cell follow the descending order of fractional
    coordinates, using the common 000--111 diagonal. All channels are retained.
    """
    coarse = np.asarray(coarse, dtype=float)
    points = np.asarray(points, dtype=float)
    shape = np.array(coarse.shape[1:])
    if (coarse.ndim != 4 or np.any(shape < 2) or points.ndim != 2
            or points.shape[1] != 3 or np.any(points < 0) or np.any(points > shape-1)):
        raise ValueError("Interpolation requires channels on a 3D grid and points inside it.")
    base = np.minimum(np.floor(points).astype(int), shape-2)
    fraction = points-base
    result = np.zeros((len(coarse), len(points)))
    if kind == "trilinear":
        for corner in CORNERS:
            weight = np.prod(np.where(corner, fraction, 1-fraction), axis=1)
            index = base+corner
            result += coarse[:, index[:, 0], index[:, 1], index[:, 2]]*weight
    elif kind == "tetra":
        order = np.argsort(-fraction, axis=1, kind="stable")
        sorted_fraction = np.take_along_axis(fraction, order, axis=1)
        weights = np.column_stack((1-sorted_fraction[:, 0],
                                   sorted_fraction[:, 0]-sorted_fraction[:, 1],
                                   sorted_fraction[:, 1]-sorted_fraction[:, 2],
                                   sorted_fraction[:, 2]))
        index = base.copy()
        for vertex in range(4):
            if vertex:
                index[np.arange(len(points)), order[:, vertex-1]] += 1
            result += coarse[:, index[:, 0], index[:, 1], index[:, 2]]*weights[:, vertex]
    else:
        raise ValueError(f"Unknown interpolation: {kind}")
    return result


def fraction_integral(values):
    """Normalized tensor trapezoid quadrature on a fixed closed patch.

    Endpoint weights are halved; a full patch always integrates to one.
    For threshold masks this is an approximate volume, not a clipped solid.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 3 or len(set(values.shape)) != 1 or values.shape[0] < 2:
        raise ValueError("Require a cubic array with at least two samples per axis.")
    weight = np.ones(values.shape[0])
    weight[[0, -1]] = .5
    return float(np.einsum("i,j,k,ijk->", weight, weight, weight, values, optimize=True)
                 / (values.shape[0]-1)**3)


def topology(mask, connectivity):
    """Bounded sample topology and component incidence on the six patch faces."""
    if not np.any(mask):
        return dict(L0=0, components=0, component_faces=[], face_pairs=0, interior_components=0)
    rank = {6: 1, 26: 3}[connectivity]
    labels, count = ndimage.label(mask, ndimage.generate_binary_structure(3, rank))
    faces = np.zeros(count+1, dtype=np.int64)
    for axis in range(3):
        for side, position in enumerate((0, -1)):
            faces[np.unique(np.take(labels, position, axis=axis))] |= 1 << (2*axis+side)
    signatures = sorted(int(v) for v in faces[1:])
    pairs = 0
    bit = 0
    for a in range(6):
        for b in range(a+1, 6):
            if any((value & (1 << a)) and (value & (1 << b)) for value in signatures):
                pairs |= 1 << bit
            bit += 1
    return dict(L0=int(saved.measures.measure_euler(mask, rank)), components=int(count),
                component_faces=signatures, face_pairs=pairs,
                interior_components=sum(value == 0 for value in signatures))


def signature(stats):
    return stats["L0"], stats["components"], tuple(stats["component_faces"])


def labelled_topology(phase, winners, connectivity):
    union = topology(phase, connectivity)
    labels = {int(w): topology(phase & (winners == w), connectivity)
              for w in np.unique(winners[phase])}
    return union, labels


def same_topology(a, b):
    """Compare per-label Euler, count and face incidence, plus union signature."""
    return (signature(a[0]) == signature(b[0]) and set(a[1]) == set(b[1])
            and all(signature(a[1][w]) == signature(b[1][w]) for w in a[1]))


def select_patches(margin, winners, film_rows, delta, seed=62713, candidates=8192):
    """Purposeful stress cases; not a representative or independent sample."""
    n = margin.shape[0]
    if n < 8:
        raise ValueError("Automatic selection needs N>=8; supply --patches for tiny scenes.")
    selected, used = [], []

    def add(lower, kind, rationale, **details):
        lower = np.asarray(lower, dtype=int)
        if np.any(lower < 1) or np.any(lower+2 > n-2):
            return False
        if any(np.all(np.abs(lower-previous) <= 2) for previous in used):
            return False
        selected.append(dict(patch_id=len(selected), kind=kind, lower=lower.tolist(),
                             cells_per_axis=2, selection=rationale, **details))
        used.append(lower)
        return True

    contacts = [r for r in film_rows if r["kind"] == "cross_label" and float(r["target_phi"]) == .85]
    if not contacts:
        raise ValueError("The film input must contain cross-label segments at target VF .85.")
    median = np.median([float(r["direct_gap_fraction"]) for r in contacts])
    orders = [(sorted(contacts, key=lambda r: -float(r["affine_phase_disagreement_fraction"])),
               "largest affine/direct disagreement in prior 85% edge sample"),
              (sorted(contacts, key=lambda r: abs(float(r["direct_gap_fraction"])-median)),
               "gap nearest the median in prior 85% edge sample")]
    for rows, why in orders:
        for row in rows:
            if add([int(row[k])-1 for k in ("i", "j", "k")], "interface", why,
                   film_segment_id=int(row["segment_id"])):
                break
        else:
            raise ValueError("Could not choose nonoverlapping interior interface patches.")

    rng = np.random.default_rng(seed)
    lowers = np.unique(rng.integers(1, n-3, size=(candidates, 3)), axis=0)
    points = lowers[:, None, :]+COARSE_POINTS[None, :, :]
    sampled_m = margin[points[:, :, 0], points[:, :, 1], points[:, :, 2]]
    sampled_w = winners[points[:, :, 0], points[:, :, 1], points[:, :, 2]]
    counts = np.array([len(np.unique(w[m >= delta])) for w, m in zip(sampled_w, sampled_m)])
    junctions = np.flatnonzero(counts >= 3)
    junctions = sorted(junctions, key=lambda i: (-int(counts[i]), float(sampled_m[i].min()), int(i)))
    for i in junctions:
        add(lowers[i], "junction_candidate", "most distinct occupied coarse labels at target VF .85",
            occupied_coarse_labels_at_85=int(counts[i]))
        if sum(r["kind"] == "junction_candidate" for r in selected) == 2:
            break
    same = np.all(sampled_w == sampled_w[:, :1], axis=1) & np.all(sampled_m >= delta, axis=1)
    controls = sorted(np.flatnonzero(same), key=lambda i: (float(sampled_m[i].min()-delta), int(i)))
    for i in controls:
        add(lowers[i], "same_label_near_threshold", "smallest positive coarse margin excess at target VF .85",
            minimum_coarse_excess_at_85=float(sampled_m[i].min()-delta))
        if sum(r["kind"] == "same_label_near_threshold" for r in selected) == 2:
            break
    if len(selected) != 6:
        raise ValueError("Could not select two patches per category; increase --candidates or supply --patches.")
    return selected, dict(candidate_draws=candidates, unique_candidate_boxes=len(lowers),
                          junction_candidates=len(junctions), same_label_occupied_candidates=len(controls),
                          selection_seed=seed, selection_delta=delta,
                          selection_target_phi=.85, sampling="purposeful, six disjoint patch boxes")


def validate_patches(patches, n):
    if not patches or len({r["patch_id"] for r in patches}) != len(patches):
        raise ValueError("Require nonempty patches with distinct IDs.")
    for row in patches:
        lower = np.asarray(row["lower"])
        if (lower.shape != (3,) or np.any(lower != np.floor(lower))
                or np.any(lower < 0) or np.any(lower+2 >= n)
                or row.get("cells_per_axis", 2) != 2
                or not isinstance(row["patch_id"], int) or row["patch_id"] < 0):
            raise ValueError("Each patch needs an integer ID and an in-domain two-cell box.")


def refine_direct(sampler, lower, old_margin, old_winner, refinement):
    size = 2*refinement+1
    ratio, remainder = divmod(size-1, len(old_margin)-1)
    if remainder or ratio < 1:
        raise ValueError("Direct evaluation requires nested grids.")
    margin = np.empty((size,)*3)
    winner = np.empty((size,)*3, dtype=np.int16)
    margin[::ratio, ::ratio, ::ratio] = old_margin
    winner[::ratio, ::ratio, ::ratio] = old_winner
    new = np.ones((size,)*3, dtype=bool)
    new[::ratio, ::ratio, ::ratio] = False
    points = grid_points(refinement)
    if np.any(new):
        values = sampler.evaluate(points[new.ravel()]+lower)
        m, w = films.margin_and_winner(values)
        margin[new], winner[new] = m, w
    return margin, winner, int(new.sum())


def reconstruct(coarse_fields, saved_margin, refinement):
    points = grid_points(refinement)
    shape = (2*refinement+1,)*3
    arrays = {}
    for interpolation in ("trilinear", "tetra"):
        m, w = films.margin_and_winner(interpolate(coarse_fields, points, interpolation))
        arrays[f"{interpolation}_components"] = (m.reshape(shape), w.reshape(shape))
    # A scalar margin has no label information. Assign the trilinear component
    # winner only for bookkeeping; union measurements do not depend on it.
    m = interpolate(saved_margin[None], points, "trilinear")[0].reshape(shape)
    arrays["trilinear_margin"] = (m, arrays["trilinear_components"][1])
    return arrays


def measure_level(patch, refinement, arrays, targets, deltas, dx, previous, volume_tolerance):
    rows, winner_rows, cache = [], [], {}
    patch_volume = (2*dx)**3
    for phi, delta in zip(targets, deltas):
        direct_m, direct_w = arrays["direct"]
        reference_phase = direct_m >= delta
        reference_labelled = np.where(reference_phase, direct_w, -1)
        references = {c: labelled_topology(reference_phase, direct_w, c) for c in (6, 26)}
        for method in METHODS:
            margin, winners = arrays[method]
            phase = margin >= delta
            volume_fraction = fraction_integral(phase)
            contacts = saved.cross_label_contacts(phase, winners)
            for connectivity in (6, 26):
                key = phi, method, connectivity
                stats = references[connectivity] if method == "direct" else labelled_topology(phase, winners, connectivity)
                union, labels = stats
                prior = previous.get(key)
                drift = abs(volume_fraction-prior["phi"]) if prior else None
                topology_stable = same_topology(stats, prior["topology"]) if prior else None
                stable = bool(topology_stable and drift <= volume_tolerance) if prior else False
                streak = (prior["streak"]+1 if stable else 0) if prior else 0
                cache[key] = dict(phi=volume_fraction, topology=stats, streak=streak)
                identity = dict(patch_id=patch["patch_id"], kind=patch["kind"], target_phi=phi,
                                delta=delta, method=method, refinement=refinement,
                                points_per_axis=2*refinement+1, spacing=dx/refinement,
                                connectivity=connectivity)
                rows.append(dict(**identity, patch_side=2*dx, patch_volume=patch_volume,
                    local_phi=volume_fraction, local_L3_quadrature=volume_fraction*patch_volume,
                    union_L0=union["L0"], union_components=union["components"],
                    union_component_faces=json.dumps(union["component_faces"]), union_face_pairs=union["face_pairs"],
                    per_winner_L0=sum(s["L0"] for s in labels.values()),
                    per_winner_components=sum(s["components"] for s in labels.values()),
                    present_winners=len(labels), cross_label_contacts=contacts[connectivity]["cross_label_contacts"],
                    occupied_neighbor_pairs=contacts[connectivity]["occupied_neighbor_pairs"],
                    phase_disagreement_phi=fraction_integral(phase != reference_phase),
                    labelled_disagreement_phi=fraction_integral(np.where(phase, winners, -1) != reference_labelled),
                    matches_direct_topology=same_topology(stats, references[connectivity]),
                    phi_change_from_previous=drift, topology_matches_previous=topology_stable,
                    stable_comparisons=streak,
                    two_step_stable=streak >= 2))
                for winner, item in sorted(labels.items()):
                    winner_rows.append(dict(**identity, winner=winner, **{k: json.dumps(v) if isinstance(v, list) else v
                                                                         for k, v in item.items()}))
    return rows, winner_rows, cache


def save_patch(output, patch, refinement, arrays, coarse_fields, coarse_margin, dx, origin):
    patch_id = patch["patch_id"]
    payload = dict(lower=np.array(patch["lower"]), refinement=np.array(refinement),
                   coarse_fields=coarse_fields, coarse_saved_margin=coarse_margin)
    mesh = pv.ImageData(dimensions=(2*refinement+1,)*3, spacing=(dx/refinement,)*3,
                        origin=tuple(np.asarray(origin)+np.array(patch["lower"])*dx))
    for method, (margin, winner) in arrays.items():
        payload[f"{method}_margin"] = margin
        payload[f"{method}_winner"] = winner.astype(np.int16)
        mesh.point_data[f"{method}_margin"] = margin.ravel(order="F")
        mesh.point_data[f"{method}_winner"] = winner.astype(np.int16).ravel(order="F")
    np.savez_compressed(output/f"patch_{patch_id:02d}.npz", **payload)
    mesh.save(output/f"patch_{patch_id:02d}.vti")


def plot_slices(output, patch, arrays, targets, deltas, dx):
    colors = np.vstack(([.94, .94, .94, 1.], plt.get_cmap("tab20")(np.linspace(0, 1, 20))))
    q = int(max(w.max() for _, w in arrays.values()))+1
    cmap = ListedColormap(np.vstack((colors[:1], np.tile(colors[1:], (int(np.ceil(q/20)), 1))[:q])))
    norm = BoundaryNorm(np.arange(-1.5, q+.5), cmap.N)
    for phi, delta in zip(targets, deltas):
        fig, axes = plt.subplots(3, 4, figsize=(12, 9), constrained_layout=True)
        for axis in range(3):
            for column, method in enumerate(METHODS):
                m, w = arrays[method]
                labelled = np.where(m >= delta, w, -1)
                section = np.take(labelled, len(m)//2, axis=axis).T
                axes[axis, column].imshow(section, origin="lower", extent=(0, 2*dx, 0, 2*dx),
                                         interpolation="nearest", cmap=cmap, norm=norm)
                axes[axis, column].set_title(f"{method}\n{'xyz'[axis]} = patch midplane", fontsize=9)
                axes[axis, column].set_xlabel("Local coordinate (domain units)")
        fig.suptitle(f"Patch {patch['patch_id']}: {patch['kind']}; target VF {phi:.0%}\n"
                     "Gray = matrix; colors = winner labels (palette repeats after 20)")
        fig.savefig(output/f"patch_{patch['patch_id']:02d}_phi{round(100*phi)}.png", dpi=130)
        plt.close(fig)


def write_summary(output, rows, complete):
    final_refinement = max(r["refinement"] for r in rows)
    finest = [r for r in rows if r["refinement"] == final_refinement and r["connectivity"] == 6]
    lines = ["# Local 3D patch reconstruction", "", f"Status: {'complete' if complete else 'partial'}.",
        "Six purposefully selected patches by default; one realization, paired thresholds.",
        "All topology is digital, on bounded patches. No global LKC or continuum-topology validation.",
        "Volume uses fixed-domain trapezoid quadrature; errors below are percentage points of patch volume.", "",
        "| Patch | Target VF | Direct local VF | Trilinear components error | Tetra components error | Scalar margin error |",
        "|---|---|---|---|---|---|"]
    for patch_id, phi in sorted({(r["patch_id"], r["target_phi"]) for r in finest}):
        group = {r["method"]: r for r in finest if r["patch_id"] == patch_id and r["target_phi"] == phi}
        direct = group["direct"]
        errors = " | ".join(f"{100*group[m]['phase_disagreement_phi']:.4f}" for m in METHODS[1:])
        lines.append(f"| {patch_id} | {phi:.0%} | {direct['local_phi']:.4%} | {errors} |")
    lines += ["", "Errors measure spatial phase disagreement, not just net volume bias.",
              "The scalar margin baseline uses trilinear-component winners for label bookkeeping only.",
              "Topology and refinement diagnostics are in measurements.csv and winners.csv.",
              "A stable digital signature is not certified continuum topology; small or tangential features can remain.",
              "Patch-face connections are local and do not imply spanning of the original domain.", ""]
    (output/"summary.md").write_text("\n".join(lines), encoding="utf-8")


def run(config, source, output, film_directory, patch_file=None):
    source, output, film_directory = Path(source).resolve(), Path(output).resolve(), Path(film_directory).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output is not empty; choose a new --output.")
    margin, winners, grid, scene = saved.load_scene(source, config["q"], config["source_phi"])
    first = pv.read(source/"field_g1.vti")
    if (tuple(first.dimensions) != margin.shape or not np.allclose(first.spacing, scene["spacing"])
            or not np.allclose(first.origin, scene["origin"])):
        raise ValueError("First field and scene geometry differ.")
    scene["input_sha256"]["field_g1.vti"] = saved.sha256(source/"field_g1.vti")
    extra_inputs = {}
    if patch_file is None:
        prior = json.loads((film_directory/"metadata.json").read_text())
        if (prior["scene"]["input_sha256"] != scene["input_sha256"]
                or any(prior["config"][k] != config[k] for k in ("q", "lc", "field_seed", "mode_no"))):
            raise ValueError("Film audit and patch study must use the same realization and parameters.")
        patches, selection = select_patches(margin, winners, read_csv(film_directory/"measurements.csv"),
            saved.adler.delta_for_vf(config["q"], .85), config["selection_seed"], config["candidates"])
        extra_inputs = {str(film_directory/name): saved.sha256(film_directory/name)
                        for name in ("metadata.json", "measurements.csv")}
    else:
        patch_file = Path(patch_file).resolve()
        patches = json.loads(patch_file.read_text())
        selection = dict(sampling="explicit patch file", path=str(patch_file))
        extra_inputs[str(patch_file)] = saved.sha256(patch_file)
    validate_patches(patches, grid.n)
    sampler = films.OriginalFields(config["q"], config["lc"], grid.dx, config["field_seed"],
                                   config["mode_no"], config["threads"])
    nodes = np.unique(np.concatenate([COARSE_POINTS+np.array(p["lower"]) for p in patches]), axis=0)
    validation = films.validate_nodes(sampler, nodes, margin, winners, source)
    print(f"Recovered saved values at {len(nodes)} patch nodes; {len(patches)} patches selected.", flush=True)
    output.mkdir(parents=True, exist_ok=True)
    saved.write_json(output/"patches.json", patches)
    sources = [Path(__file__), HERE/"run_local_films.py", HERE/"run_saved_sweep.py",
               saved.ROOT/"src/fields.py", saved.ROOT/"src/adler.py", saved.ROOT/"src/measures.py"]
    metadata = dict(schema=SCHEMA, status="running", config=config, source=str(source), scene=scene,
        selection=selection, seed_validation=validation, component_seeds=sampler.seeds,
        extra_input_sha256=extra_inputs, source_sha256={p.relative_to(saved.ROOT).as_posix(): saved.sha256(p) for p in sources},
        versions={v: version(v) for v in ("numpy", "scipy", "scikit-image", "gstools", "pyvista", "matplotlib")},
        independent_realizations=1, completed_patches=0, total_patches=len(patches),
        volume_convention="tensor trapezoid over fixed patch side 2*dx; sample topology under 6/26 adjacency",
        python=sys.version)
    saved.write_json(output/"metadata.json", metadata)
    rows, winner_rows = [], []
    start = time.perf_counter()
    deltas = [saved.adler.delta_for_vf(config["q"], phi) for phi in config["phis"]]
    try:
        for patch in patches:
            lower = np.array(patch["lower"])
            coarse = sampler.evaluate(COARSE_POINTS+lower).reshape((config["q"], 3, 3, 3))
            coarse_margin = margin[tuple(slice(int(v), int(v)+3) for v in lower)]
            dm, dw = films.margin_and_winner(coarse)
            previous = {}
            for refinement in config["refinements"]:
                print(f"Patch {patch['patch_id']} ({patch['kind']}), {2*refinement+1}^3 direct samples...", flush=True)
                dm, dw, added = refine_direct(sampler, lower, dm, dw, refinement)
                arrays = dict(direct=(dm, dw), **reconstruct(coarse, coarse_margin, refinement))
                measured, per_winner, previous = measure_level(patch, refinement, arrays, config["phis"],
                    deltas, grid.dx, previous, config["volume_tolerance"])
                rows.extend(measured)
                winner_rows.extend(per_winner)
                saved.write_csv(output/"measurements.csv", rows)
                saved.write_csv(output/"winners.csv", winner_rows)
                print(f"  {added} new locations evaluated; measurements saved.", flush=True)
            save_patch(output, patch, refinement, arrays, coarse, coarse_margin, grid.dx, scene["origin"])
            plot_slices(output, patch, arrays, config["phis"], deltas, grid.dx)
            metadata["completed_patches"] += 1
            saved.write_json(output/"metadata.json", metadata)
            write_summary(output, rows, complete=False)
        metadata.update(status="complete", elapsed_probe_and_report_seconds=time.perf_counter()-start)
        saved.write_json(output/"metadata.json", metadata)
        write_summary(output, rows, complete=True)
    except BaseException as exc:
        metadata.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                        error=f"{type(exc).__name__}: {exc}")
        saved.write_json(output/"metadata.json", metadata)
        raise
    print(f"Complete: {len(patches)} patches. {output}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE/"paraview/q32_phi75")
    parser.add_argument("--films", type=Path, default=HERE/"results_local_films_refined")
    parser.add_argument("--patches", type=Path, help="Optional fixed selection JSON, e.g. a previous patches.json")
    parser.add_argument("--output", type=Path, default=HERE/"results_local_patches")
    parser.add_argument("--q", type=int, default=32)
    parser.add_argument("--lc", type=float, default=5.)
    parser.add_argument("--source-phi", type=float, default=.75)
    parser.add_argument("--field-seed", type=int, default=321729)
    parser.add_argument("--mode-no", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--selection-seed", type=int, default=62713)
    parser.add_argument("--candidates", type=int, default=8192)
    parser.add_argument("--phis", nargs="+", type=float, default=[.75, .85])
    parser.add_argument("--refinements", nargs="+", type=int, default=[4, 8, 16])
    parser.add_argument("--volume-tolerance", type=float, default=.005,
                        help="Maximum local volume-fraction drift for a stable comparison")
    args = parser.parse_args()
    config = {k: v for k, v in vars(args).items() if k not in ("source", "films", "patches", "output")}
    config["phis"] = sorted(set(args.phis))
    config["refinements"] = sorted(set(args.refinements))
    if (args.q < 2 or args.q > 32767 or min(args.mode_no, args.threads, args.candidates) < 1
            or min(args.field_seed, args.selection_seed) < 0 or not np.isfinite(args.lc) or args.lc <= 0
            or not np.isfinite(args.volume_tolerance) or args.volume_tolerance <= 0
            or any(not np.isfinite(p) or not 0 < p < 1 for p in [args.source_phi, *args.phis])
            or any(r < 1 or r & (r-1) for r in config["refinements"])):
        parser.error("Require positive counts/scales, nonnegative seeds, 0<phi<1, and power-of-two refinements.")
    try:
        run(config, args.source, args.output, args.films, args.patches)
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    main()
