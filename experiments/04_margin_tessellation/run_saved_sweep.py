"""Rethreshold one saved margin/winner scene; no random-field generation.

Default: the local q=32, N=128 scene, target VFs .60/.70/.75/.80/.85.
Point samples retain the repository's voxel convention: domain side N*dx.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy import ndimage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from src import adler, measures
from src.fields import GridSpec

SCHEMA = "saved_margin_sweep_v1"


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_csv(path, rows, empty_fields=()):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else list(empty_fields))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def load_scene(source, q, source_phi):
    arrays, reference = {}, None
    hashes = {}
    for filename, name in (("margin.vti", "margin"), ("winners.vti", "winner_label"),
                           ("phase.vti", "phase")):
        path = source / filename
        mesh = pv.read(path)
        if not isinstance(mesh, pv.ImageData) or name not in mesh.point_data:
            raise ValueError(f"{path}: expected ImageData point array {name!r}.")
        geometry = (tuple(mesh.dimensions), tuple(mesh.spacing), tuple(mesh.origin),
                    tuple(np.asarray(mesh.direction_matrix).ravel()))
        if reference is None:
            reference = geometry
        elif geometry != reference:
            raise ValueError("Saved arrays must have identical grid geometry.")
        arrays[name] = np.asarray(mesh.point_data[name]).reshape(mesh.dimensions, order="F").copy()
        hashes[filename] = sha256(path)
    shape, spacing, origin, direction = reference
    if (len(set(shape)) != 1 or shape[0] < 2 or not np.allclose(spacing, spacing[0])
            or not np.isfinite(spacing).all() or min(spacing) <= 0
            or not np.allclose(np.array(direction).reshape(3, 3), np.eye(3))):
        raise ValueError("Require a cubic, axis-aligned grid with isotropic positive spacing.")
    margin, winners, baseline = arrays["margin"], arrays["winner_label"], arrays["phase"]
    if not np.isfinite(margin).all() or np.any(margin < 0):
        raise ValueError("Margins must be finite and nonnegative.")
    if (not np.isfinite(winners).all() or np.any(winners != np.floor(winners))
            or winners.min() < 0 or winners.max() >= q):
        raise ValueError("Winner labels must be integers in [0, q).")
    if not np.isin(baseline, [0, 1]).all():
        raise ValueError("The saved reference phase must be binary.")
    delta = adler.delta_for_vf(q, source_phi)
    mismatches = int(np.count_nonzero((baseline > 0) != (margin >= delta)))
    if mismatches:
        raise ValueError(f"Saved reference phase disagrees with q/source-phi at {mismatches} samples.")
    grid = GridSpec(shape[0], shape[0] * spacing[0])
    info = dict(dimensions=list(shape), spacing=list(spacing), origin=list(origin),
                domain_length=grid.length, point_coordinate_span=(shape[0]-1)*spacing[0],
                volume_convention="N^3 voxel samples of volume dx^3; nominal domain side N*dx",
                input_sha256=hashes, reference_phi=source_phi, reference_delta=delta,
                reference_measured_phi=float(baseline.mean()), reference_mask_mismatches=mismatches)
    return margin, winners.astype(np.int32), grid, info


def component_statistics(mask, dx, connectivity):
    """Bounded digital components, matching the experiment 05 conventions."""
    rank = {6: 1, 26: 3}[connectivity]
    labels, count = ndimage.label(mask, ndimage.generate_binary_structure(3, rank))
    sizes = np.bincount(labels.ravel(), minlength=count+1)
    sizes[0] = 0
    boundary = np.zeros(count+1, dtype=bool)
    spanning = np.zeros((3, count+1), dtype=bool)
    for axis in range(3):
        low = np.unique(np.take(labels, 0, axis=axis))
        high = np.unique(np.take(labels, -1, axis=axis))
        boundary[low] = boundary[high] = True
        spanning[axis, np.intersect1d(low, high)] = True
    boundary[0] = False
    spanning[:, 0] = False
    occupied = int(sizes.sum())
    stats = dict(L0=int(measures.measure_euler(mask, rank)), components=int(count),
                 voxels=occupied, volume=float(occupied*dx**3),
                 largest_volume=float(sizes.max()*dx**3),
                 largest_fraction=float(sizes.max()/occupied) if occupied else 0.,
                 spans=bool(spanning.any()),
                 spanning_volume=float(sizes[spanning.any(axis=0)].sum()*dx**3),
                 boundary_volume=float(sizes[boundary].sum()*dx**3))
    stats.update({f"spans_{axis}": bool(spanning[j].any()) for j, axis in enumerate("xyz")})
    rows = []
    for label in range(1, count+1):
        volume = float(sizes[label]*dx**3)
        rows.append(dict(component=label, voxels=int(sizes[label]), volume=volume,
                         equivalent_diameter=float((6*volume/np.pi)**(1/3)),
                         touches_boundary=bool(boundary[label]), spans=bool(spanning[:, label].any()),
                         **{f"spans_{axis}": bool(spanning[j, label]) for j, axis in enumerate("xyz")}))
    return stats, rows


def cross_label_contacts(phase, winners):
    """Count each undirected occupied pair once; 26 includes face pairs."""
    totals, cross = {6: 0, 26: 0}, {6: 0, 26: 0}
    for offset in np.ndindex(3, 3, 3):
        offset = tuple(v-1 for v in offset)
        if offset <= (0, 0, 0):
            continue
        left = tuple(slice(max(0, -d), n-max(0, d)) for d, n in zip(offset, phase.shape))
        right = tuple(slice(max(0, d), n-max(0, -d)) for d, n in zip(offset, phase.shape))
        occupied = phase[left] & phase[right]
        count = int(occupied.sum())
        contacts = int(np.count_nonzero(occupied & (winners[left] != winners[right])))
        totals[26] += count
        cross[26] += contacts
        if sum(abs(d) for d in offset) == 1:
            totals[6] += count
            cross[6] += contacts
    return {c: dict(occupied_neighbor_pairs=totals[c], cross_label_contacts=cross[c],
                    cross_label_fraction=cross[c]/totals[c] if totals[c] else 0.) for c in (6, 26)}


def size_quantile(components, probability):
    """Volume-weighted equivalent diameter; NaN is not needed for empty masks."""
    if not components:
        return 0.
    order = sorted(components, key=lambda row: row["equivalent_diameter"])
    cumulative = np.cumsum([r["volume"] for r in order])
    index = min(int(np.searchsorted(cumulative, probability*cumulative[-1])), len(order)-1)
    return order[index]["equivalent_diameter"]


def measure_threshold(margin, winners, grid, q, lc, target_phi):
    delta = adler.delta_for_vf(q, target_phi)
    phase = margin >= delta
    occupied = int(phase.sum())
    volume = float(occupied*grid.dx**3)
    theory = adler.expected_lkcs_margin(grid, lc, q, delta)
    contacts = cross_label_contacts(phase, winners)
    measurements, winner_rows, component_rows = [], [], []
    for connectivity in (6, 26):
        identity = dict(target_phi=target_phi, delta=delta, connectivity=connectivity)
        union, _ = component_statistics(phase, grid.dx, connectivity)
        stats, components = [], []
        for winner in range(q):
            current, detail = component_statistics(phase & (winners == winner), grid.dx, connectivity)
            stats.append(current)
            winner_rows.append(dict(**identity, winner=winner, **current))
            components.extend(dict(**identity, winner=winner, **r) for r in detail)
        if sum(r["voxels"] for r in stats) != occupied:
            raise RuntimeError("Per-winner volume accounting failed.")
        if sum(r["voxels"] for r in components) != occupied:
            raise RuntimeError("Component volume accounting failed.")
        component_rows.extend(components)
        fraction = lambda value: float(value/volume) if volume else 0.
        row = dict(**identity, N=grid.n, q=q, length=grid.length, lc=lc,
                   measured_phi=float(phase.mean()), measured_L3=volume,
                   theory_phi=theory["phi"], **{f"theory_L{i}": theory[f"L{i}"] for i in range(4)},
                   **{f"union_{k}": v for k, v in union.items()},
                   per_winner_L0=sum(r["L0"] for r in stats),
                   per_winner_components=sum(r["components"] for r in stats),
                   spanning_winners=sum(r["spans"] for r in stats),
                   largest_separate_component_fraction=fraction(max(r["largest_volume"] for r in stats)),
                   separate_spanning_volume_fraction=fraction(sum(r["spanning_volume"] for r in stats)),
                   separate_boundary_volume_fraction=fraction(sum(r["boundary_volume"] for r in stats)),
                   union_spanning_volume_fraction=fraction(union["spanning_volume"]),
                   euler_additivity_gap=union["L0"]-sum(r["L0"] for r in stats),
                   **contacts[connectivity],
                   **{f"diameter_volume_p{p}": size_quantile(components, p/100) for p in (10, 50, 90)})
        measurements.append(row)
    return measurements, winner_rows, component_rows


def write_summary(output, rows, complete):
    lines = ["# Saved q-scene threshold sweep", "",
             f"Status: {'complete' if complete else 'partial'}; {len(rows)//2} thresholds, one saved realization.",
             "Thresholds and winner labels are correlated observations, not independent trials.",
             "6 / 26 always denotes digital connectivity; no continuum claim follows from either.", "",
             "| Target VF | Measured VF | Delta | Spanning winners, 6 / 26 | Union spans, 6 / 26 | Largest separate / phase, 6 / 26 | Largest union / phase, 6 / 26 |",
             "|---|---|---|---|---|---|---|"]
    groups = {phi: [next(r for r in rows if r["target_phi"] == phi and r["connectivity"] == c)
                    for c in (6, 26)] for phi in sorted({r["target_phi"] for r in rows})}
    for phi, (a, b) in groups.items():
        lines.append(f"| {phi:.0%} | {a['measured_phi']:.3%} | {a['delta']:.6f} | "
                     f"{a['spanning_winners']} / {b['spanning_winners']} | {a['union_spans']} / {b['union_spans']} | "
                     f"{a['largest_separate_component_fraction']:.3%} / {b['largest_separate_component_fraction']:.3%} | "
                     f"{a['union_largest_fraction']:.3%} / {b['union_largest_fraction']:.3%} |")
    lines += ["", "| Target VF | Expected L0 | Per-winner L0, 6 / 26 | Union L0, 6 / 26 | Separate components, 6 / 26 | Cross-label pairs, 6 / 26 |",
              "|---|---|---|---|---|---|"]
    for phi, (a, b) in groups.items():
        lines.append(f"| {phi:.0%} | {a['theory_L0']:.2f} | {a['per_winner_L0']} / {b['per_winner_L0']} | "
                     f"{a['union_L0']} / {b['union_L0']} | {a['per_winner_components']} / {b['per_winner_components']} | "
                     f"{a['cross_label_contacts']} / {b['cross_label_contacts']} |")
    lines += ["", "## Interpretation", "",
              "The theory columns are ensemble GKF expectations under the declared q and lc;",
              "one saved scene cannot validate these expectations or supply uncertainty intervals.",
              "Winner masks are measured separately to expose cross-label voxel mergers.",
              "Their Euler sums may still have within-label discretization errors.",
              "Cross-label contacts at positive margin indicate unresolved continuum separation.",
              "This rethresholding does not reconstruct the underlying Gaussian components or their films.",
              "Point data use the existing N*dx voxel-volume convention, not the (N-1)*dx VTI point extent.",
              "The size statistics include boundary-cut and spanning objects. components.csv records",
              "their flags; equivalent diameter is volume-based and is not a calibrated sieve diameter.",
              "No spanning is a statement about these masks, not a probability or continuum guarantee.", ""]
    (output/"summary.md").write_text("\n".join(lines), encoding="utf-8")


def plot_summary(output, rows):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    first = [r for r in rows if r["connectivity"] == 6]
    x = [100*r["target_phi"] for r in first]
    axes[0, 0].plot(x, x, "k--", label="Target / expected")
    axes[0, 0].plot(x, [100*r["measured_phi"] for r in first], "o-", label="Saved scene")
    axes[0, 1].plot(x, [r["theory_L0"] for r in first], "k--", label="GKF expectation")
    for c, color, marker in ((6, "tab:blue", "o"), (26, "tab:orange", "s")):
        group = [r for r in rows if r["connectivity"] == c]
        label = f"{c} neighbours"
        axes[0, 1].plot(x, [r["per_winner_L0"] for r in group], marker+"-", color=color, label=label)
        axes[0, 2].plot(x, [r["union_L0"] for r in group], marker+"-", color=color, label=label)
        axes[1, 0].plot(x, [100*r["union_largest_fraction"] for r in group], marker+"-", color=color, label=f"Union, {c}")
        axes[1, 0].plot(x, [100*r["largest_separate_component_fraction"] for r in group], marker+"--", color=color, label=f"Separate, {c}")
        axes[1, 1].plot(x, [100*r["cross_label_fraction"] for r in group], marker+"-", color=color, label=label)
        for p, style in ((50, "-"), (90, "--")):
            axes[1, 2].plot(x, [r[f"diameter_volume_p{p}"] for r in group], marker+style, color=color, label=f"D{p}, {c}")
    settings = [
        ("Volume fraction", "Measured VF (%)"),
        ("Euler summed over winner labels", "L0"),
        ("Euler of the binary union", "L0"),
        ("Largest component / aggregate volume", "Percentage"),
        ("Cross-label / occupied neighbour pairs", "Percentage"),
        ("Volume-weighted separate-component sizes\nIncludes boundary-cut and spanning components", "Equivalent diameter (domain units)")]
    for ax, (title, ylabel) in zip(axes.flat, settings):
        ax.set(title=title, xlabel="Target VF (%)", ylabel=ylabel)
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.suptitle(f"Saved margin scene: q={first[0]['q']}, N={first[0]['N']}; one realization, both connectivities")
    fig.savefig(output/"threshold_sweep.png", dpi=160)
    plt.close(fig)


def run(source, output, q=32, lc=5., source_phi=.75, targets=(.60, .70, .75, .80, .85)):
    source, output = Path(source).resolve(), Path(output).resolve()
    targets = sorted(set(float(phi) for phi in targets))
    if q < 2 or not isinstance(q, int) or not np.isfinite(lc) or lc <= 0:
        raise ValueError("Require integer q>=2 and finite lc>0.")
    if not targets or any(not np.isfinite(phi) or not 0 < phi < 1 for phi in [source_phi, *targets]):
        raise ValueError("VF values must be finite and strictly between zero and one.")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output is not empty; choose a new --output.")
    margin, winners, grid, scene = load_scene(source, q, source_phi)
    output.mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__).resolve(), ROOT/"src/adler.py", ROOT/"src/measures.py", ROOT/"src/fields.py"]
    metadata = dict(schema=SCHEMA, status="running", expected_thresholds=len(targets), completed_thresholds=0,
                    config=dict(q=q, lc=lc, source_phi=source_phi, target_phis=targets),
                    source=str(source), scene=scene, independent_realizations=1,
                    theoretical_assumptions="i.i.d. unit Gaussian competitors, correlation exp(-r^2/lc^2); q/lc supplied by caller",
                    source_sha256={p.relative_to(ROOT).as_posix(): sha256(p) for p in sources},
                    versions={name: version(name) for name in ("numpy", "scipy", "scikit-image", "pyvista", "matplotlib")},
                    python=sys.version)
    metadata_path = output/"metadata.json"
    write_json(metadata_path, metadata)
    rows, winner_rows, component_rows = [], [], []
    start = time.perf_counter()
    try:
        for phi in targets:
            print(f"Target VF={phi:.0%}: rethresholding saved arrays and measuring...", flush=True)
            case_start = time.perf_counter()
            measured, per_winner, components = measure_threshold(margin, winners, grid, q, lc, phi)
            elapsed = time.perf_counter()-case_start
            rows.extend(dict(**r, threshold_measurement_seconds=elapsed) for r in measured)
            winner_rows.extend(per_winner)
            component_rows.extend(components)
            write_csv(output/"measurements.csv", rows)
            write_csv(output/"winners.csv", winner_rows)
            write_csv(output/"components.csv", component_rows,
                      ("target_phi", "delta", "connectivity", "winner", "component", "voxels",
                       "volume", "equivalent_diameter", "touches_boundary", "spans", "spans_x", "spans_y", "spans_z"))
            write_summary(output, rows, complete=False)
            metadata["completed_thresholds"] += 1
            write_json(metadata_path, metadata)
            a, b = measured
            print(f"  VF={a['measured_phi']:.4f}; spanning winners 6/26="
                  f"{a['spanning_winners']}/{b['spanning_winners']}; "
                  f"union spans={a['union_spans']}/{b['union_spans']}; {elapsed:.1f}s", flush=True)
        plot_summary(output, rows)
        write_summary(output, rows, complete=True)
        metadata.update(status="complete", elapsed_seconds=time.perf_counter()-start)
        write_json(metadata_path, metadata)
    except BaseException as exc:
        metadata.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                        error=f"{type(exc).__name__}: {exc}")
        write_json(metadata_path, metadata)
        raise
    print(f"Complete: {len(targets)} thresholds. Results: {output}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE/"paraview/q32_phi75")
    parser.add_argument("--output", type=Path, default=HERE/"results_saved_q32_sweep")
    parser.add_argument("--q", type=int, default=32)
    parser.add_argument("--lc", type=float, default=5.)
    parser.add_argument("--source-phi", type=float, default=.75)
    parser.add_argument("--target-phis", type=float, nargs="+", default=[.60, .70, .75, .80, .85])
    args = parser.parse_args()
    try:
        run(args.source, args.output, args.q, args.lc, args.source_phi, args.target_phis)
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    main()
