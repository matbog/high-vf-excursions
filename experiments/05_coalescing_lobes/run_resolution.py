"""Paired pre-contact resolution study, local to experiment 05.

Default: N=50,100,200; r/rc=.85,.98; one realization. Full measurements
include scalar meshes. --topology-only skips mesh curvature/surface.
No full study is run on import. Existing outputs require explicit --resume.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from importlib.metadata import version
from itertools import product
import json
import math
from pathlib import Path
import sys
import time

import gstools as gs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
from src import adler, excursions, measures
from src.fields import GridSpec
from measurement import measure_lobes


def component_seeds(seed: int, realization: int) -> list[int]:
    """Replicate both seed draws of the shared GSTools wrapper."""
    rng = np.random.default_rng(seed + realization - 1)
    return [int(gs.random.MasterRNG(int(rng.integers(0, 2**31 - 1)))())
            for _ in range(3)]


def sample_vector(n: int, config: dict, realization: int) -> np.ndarray:
    """Evaluate fixed RandMeth coefficients at fixed physical coordinates.

    Reference padding is a physical translation, independent of the current N.
    For reference_n=100, length=100 and reference_padding=20 the samples span
    [20,120-dx]. This matches the old N=100 cropped field, not its N-dependent
    padding policy on every grid. No per-realization normalization is applied.
    """
    axis = config["origin"] + np.arange(n, dtype=float)*(config["length"]/n)
    gs.config.NUM_THREADS = config["threads"]
    model = gs.Gaussian(dim=3, len_scale=config["lc"], var=1., rescale=1., nugget=0.)
    result = np.empty((3, n, n, n), dtype=np.float32)
    for i, seed in enumerate(component_seeds(config["seed"], realization)):
        srf = gs.SRF(model, generator="RandMeth", seed=seed, mode_no=config["mode_no"])
        for start in range(0, n, config["slab"]):
            end = min(start + config["slab"], n)
            result[i, start:end] = srf.structured((axis[start:end], axis, axis), store=False)
    return result


def ball_distances(vector: np.ndarray, centres: np.ndarray) -> np.ndarray:
    """Bound temporary memory: accumulate one component at a time."""
    # Match the original runner's float64 distance calculation after float32 sampling.
    distances = np.empty((len(centres), *vector.shape[1:]), dtype=np.float64)
    scratch = np.empty(vector.shape[1:], dtype=np.float64)
    for i, centre in enumerate(centres):
        distances[i].fill(0)
        for j in range(3):
            np.subtract(vector[j], centre[j], out=scratch, dtype=np.float64)
            np.square(scratch, out=scratch)
            distances[i] += scratch
        np.sqrt(distances[i], out=distances[i])
    return distances


def common_node_error(anchor: np.ndarray, other: np.ndarray) -> tuple[float, int]:
    """Compare shared points even when resolutions are not integer multiples."""
    a, b = anchor.shape[1], other.shape[1]
    common = math.gcd(a, b)
    error = np.max(np.abs(anchor[:, ::a//common, ::a//common, ::a//common]
                          - other[:, ::b//common, ::b//common, ::b//common]))
    return float(error), common**3


def cross_label_contacts(phase: np.ndarray, labels: np.ndarray) -> dict:
    """Count undirected occupied neighbour pairs joining distinct value-space balls."""
    crossing = {6: 0, 26: 0}
    occupied = {6: 0, 26: 0}
    for offset in product((-1, 0, 1), repeat=3):
        nonzero = [v for v in offset if v]
        if not nonzero or nonzero[0] < 0:
            continue
        left = tuple(slice(max(0, -v), n-max(0, v)) for n, v in zip(phase.shape, offset))
        right = tuple(slice(max(0, v), n-max(0, -v)) for n, v in zip(phase.shape, offset))
        pair = phase[left] & phase[right]
        count = int(pair.sum())
        cross = int(np.count_nonzero(pair & (labels[left] != labels[right])))
        occupied[26] += count
        crossing[26] += cross
        if sum(abs(v) for v in offset) == 1:
            occupied[6] += count
            crossing[6] += cross
    return {key: value for rank in (6, 26) for key, value in (
        (f"cross_label_contacts_{rank}", crossing[rank]),
        (f"occupied_neighbor_pairs_{rank}", occupied[rank]),
        (f"cross_label_fraction_{rank}", crossing[rank]/occupied[rank] if occupied[rank] else 0.))}


def measure_resolution(distances, radius, grid, include_mesh=True):
    result = measure_lobes(distances, radius, grid, include_mesh=include_mesh)
    phase = distances.min(axis=0) <= radius
    result.update(cross_label_contacts(phase, distances.argmin(axis=0)))
    result["sum_ball_euler_6"] = 0
    spans = {6: False, 26: False}
    for distance in distances:
        ball = distance <= radius
        result["sum_ball_euler_6"] += int(measures.measure_euler(ball, 1))
        for rank, connectivity in ((6, 1), (26, 3)):
            spans[rank] |= any(measures.measure_percolation(ball, connectivity).values())
    for rank in (6, 26):
        result[f"any_ball_spans_{rank}"] = spans[rank]
        result[f"euler_additivity_gap_{rank}"] = (
            result[f"union_euler_{rank}"] - result[f"sum_ball_euler_{rank}"])
    result["L1_additivity_gap"] = result["union_L1_field"] - result["sum_ball_L1_field"]
    result["L2_crofton_additivity_gap"] = result["union_L2_crofton"] - result["sum_ball_L2_crofton"]
    return result


def write_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    temporary.replace(path)


def write_csv(path, rows):
    if not rows:
        return
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def read_csv(path):
    if not path.exists():
        return []
    def number(v):
        return v == "True" if v in ("True", "False") else float(v)
    with path.open(newline="", encoding="utf-8") as stream:
        return [{k: number(v) for k, v in row.items()} for row in csv.DictReader(stream)]


def source_hashes():
    paths = [Path(__file__), HERE/"measurement.py", ROOT/"src/adler.py",
             ROOT/"src/measures.py", ROOT/"src/fields.py", ROOT/"src/excursions.py"]
    return {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def summary(output, rows, config, complete):
    """Plot individual paired paths. One realization is not an ensemble validation."""
    fig, axes = plt.subplots(len(config["radii"]), 4, squeeze=False,
                             figsize=(19, 4*len(config["radii"])), constrained_layout=True)
    text = ["# Paired pre-contact resolution study", "",
            f"Status: {'complete' if complete else 'partial'}. {len(rows)} completed cases.",
            "Curves join resolutions of the same realization. GKF curves are ensemble expectations.",
            "One realization cannot validate an expectation; refinement does not remove sampling uncertainty.",
            "", "| r/rc | N | completed realisations | VF | Euler 6 | Euler 26 | ball sum 26 | cross-label pairs 6 / 26 |",
            "|---|---|---|---|---|---|---|---|"]
    for index, ratio in enumerate(config["radii"]):
        group = [r for r in rows if r["radius_over_rc"] == ratio]
        for realization in range(1, config["realizations"]+1):
            cells = sorted([r for r in group if r["realization"] == realization], key=lambda r: r["N"])
            if not cells:
                continue
            xs = [r["N"] for r in cells]
            for key, style, label in (("union_euler_6", "s--", "union 6"),
                                      ("union_euler_26", "o-", "union 26"),
                                      ("sum_ball_euler_26", "x:", "ball sum 26")):
                axes[index, 0].plot(xs, [r[key] for r in cells], style,
                                    label=f"{label}, realization {realization}")
            for rank, style in ((6, "s--"), (26, "o-")):
                axes[index, 1].plot(xs, [r[f"cross_label_contacts_{rank}"] for r in cells], style,
                                    label=f"{rank} neighbours, realization {realization}")
            axes[index, 2].plot(xs, [r["measured_phi"] for r in cells], "o-", label=f"realization {realization}")
            if config["include_mesh"]:
                for key, style in (("union_L1_field", "o-"), ("sum_ball_L1_field", "x:")):
                    axes[index, 3].plot(xs, [r[key] for r in cells], style, label=f"{key}, realization {realization}")
            else:
                for rank, style in ((6, "s--"), (26, "o-")):
                    axes[index, 3].plot(xs, [r[f"euler_additivity_gap_{rank}"] for r in cells],
                                       style, label=f"{rank} neighbours, realization {realization}")
        if group:
            axes[index, 0].axhline(group[0]["theory_L0"], color="black", ls="--", label="GKF expectation")
            axes[index, 2].axhline(group[0]["theory_phi"], color="black", ls="--", label="GKF expectation")
            axes[index, 3].axhline(group[0]["theory_L1"] if config["include_mesh"] else 0,
                                   color="black", ls="--", label="GKF expectation" if config["include_mesh"] else "additivity")
        for n in config["grids"]:
            cells = [r for r in group if r["N"] == n]
            if cells:
                mean = lambda key: np.mean([r[key] for r in cells])
                text.append(f"| {ratio:.2f} | {n} | {len(cells)} | {mean('measured_phi'):.5f} | "
                            f"{mean('union_euler_6'):.2f} | {mean('union_euler_26'):.2f} | "
                            f"{mean('sum_ball_euler_26'):.2f} | "
                            f"{mean('cross_label_contacts_6'):.1f} / {mean('cross_label_contacts_26'):.1f} |")
        titles = ["Euler characteristic", "Cross-label neighbour pairs", "Volume fraction",
                  "Mesh L1" if config["include_mesh"] else "Union minus ball Euler"]
        for ax, title in zip(axes[index], titles):
            ax.set(title=f"r/rc={ratio:.2f}: {title}", xlabel="N")
            ax.set_xticks(config["grids"])
            ax.grid(alpha=.25)
            handles, _ = ax.get_legend_handles_labels()
            if handles:
                ax.legend(fontsize=6)
        axes[index, 1].set_yscale("symlog", linthresh=1)
    fig.savefig(output/"resolution.png", dpi=150)
    plt.close(fig)
    text += ["", "Cross-label contacts are digital contacts forbidden between distinct continuum lobes for r/rc<1.",
             "Raw pair counts need not decrease monotonically: the number of neighbour pairs grows with N.",
             "The CSV also records contact fractions, same-convention additivity gaps, per-ball spanning,",
             "and common-node sampling errors. Neither digital convention certifies continuum connectivity.",
             "Scalar-mesh interpolation and domain padding retain their usual limitations.", ""]
    (output/"summary.md").write_text("\n".join(text), encoding="utf-8")


def run(config, output, resume=False):
    hashes = source_hashes()
    versions = {k: version(k) for k in ("numpy", "scipy", "gstools", "scikit-image", "matplotlib")}
    occupied = output.exists() and any(output.iterdir())
    metadata_path = output/"metadata.json"
    if occupied:
        if not resume:
            raise ValueError(f"Output is not empty: {output}. Use --resume or another --output.")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("schema") != "lobes_paired_resolution_v1":
            raise ValueError("This output does not contain a paired-resolution checkpoint.")
        if metadata["config"] != config or metadata["source_sha256"] != hashes or metadata["versions"] != versions:
            raise ValueError("Resume requires identical parameters, sources and dependency versions; choose a new --output.")
        rows = read_csv(output/"measurements.csv")
        diagnostics = read_csv(output/"field_diagnostics.csv")
    else:
        if resume:
            raise ValueError("--resume requires an existing checkpoint.")
        output.mkdir(parents=True, exist_ok=True)
        rows, diagnostics = [], []
        metadata = {"schema": "lobes_paired_resolution_v1", "config": config, "versions": versions,
                    "source_sha256": hashes, "python": sys.version,
                    "sampler": "Fixed GSTools RandMeth streams; physical coordinates; float32 samples",
                    "component_seeds": {str(i): component_seeds(config["seed"], i)
                                        for i in range(1, config["realizations"]+1)}}
    done = {(int(r["realization"]), int(r["N"]), r["radius_over_rc"]) for r in rows}
    expected = {(i, n, ratio) for i in range(1, config["realizations"]+1)
                for n in config["grids"] for ratio in config["radii"]}
    if len(done) != len(rows) or not done <= expected:
        raise ValueError("Checkpoint contains duplicate or unexpected cases.")
    metadata.update(status="running", completed_cases=len(rows), expected_cases=len(expected))
    write_json(metadata_path, metadata)
    centres = excursions.tetrahedral_centres(config["rho"])
    rc = np.linalg.norm(centres[0]-centres[1])/2
    try:
        for realization in range(1, config["realizations"]+1):
            pending = [(n, ratio) for n in config["grids"] for ratio in config["radii"]
                       if (realization, n, ratio) not in done]
            if not pending:
                continue
            anchor_n = min(config["grids"])
            print(f"Realization {realization}: sampling reference grid N={anchor_n}", flush=True)
            start = time.perf_counter()
            anchor = sample_vector(anchor_n, config, realization)
            anchor_seconds = time.perf_counter()-start
            for n in config["grids"]:
                ratios = [ratio for nn, ratio in pending if nn == n]
                if not ratios:
                    continue
                print(f"  N={n}: sampling components", flush=True)
                start = time.perf_counter()
                vector = anchor if n == anchor_n else sample_vector(n, config, realization)
                sample_seconds = anchor_seconds if n == anchor_n else time.perf_counter()-start
                error, common = common_node_error(anchor, vector)
                if error > 2e-5:
                    raise RuntimeError(f"Paired sampling failed at N={n}: common-node error {error:g}")
                diagnostics = [r for r in diagnostics if (int(r["realization"]), int(r["N"])) != (realization, n)]
                for c in range(3):
                    diagnostics.append(dict(realization=realization, N=n, component=c,
                        seed=component_seeds(config["seed"], realization)[c],
                        mean=float(vector[c].mean()), second_moment=float(np.mean(vector[c]**2)),
                        common_node_max_abs=error, common_nodes=common, sampling_seconds=sample_seconds))
                write_csv(output/"field_diagnostics.csv", diagnostics)
                distances = ball_distances(vector, centres)
                del vector
                grid = GridSpec(n=n, length=config["length"])
                for ratio in ratios:
                    print(f"  N={n}, r/rc={ratio}: measuring (mesh={config['include_mesh']})", flush=True)
                    start = time.perf_counter()
                    radius = float(ratio*rc)
                    theory = adler.expected_lkcs_lobes_disjoint(grid, config["lc"], centres, radius)
                    measured = measure_resolution(distances, radius, grid, config["include_mesh"])
                    row = dict(realization=realization, N=n, dx=grid.dx, radius_over_rc=ratio, radius=radius,
                               theory_phi=theory["phi"], **{f"theory_L{i}": theory[f"L{i}"] for i in range(4)},
                               **measured, measurement_seconds=time.perf_counter()-start)
                    if not np.isclose(row["union_L3"], row["sum_ball_L3"], rtol=1e-12):
                        raise RuntimeError("Pre-contact voxel ball volumes should add exactly.")
                    rows.append(row)
                    done.add((realization, n, ratio))
                    write_csv(output/"measurements.csv", rows)
                    metadata["completed_cases"] = len(rows)
                    write_json(metadata_path, metadata)
                    summary(output, rows, config, complete=False)
                    print(f"    saved: VF={row['measured_phi']:.5f}, "
                          f"L0(6/26)={row['union_euler_6']}/{row['union_euler_26']}, "
                          f"contacts(6/26)={row['cross_label_contacts_6']}/{row['cross_label_contacts_26']}, "
                          f"{row['measurement_seconds']:.1f}s", flush=True)
                del distances
            del anchor
        summary(output, rows, config, complete=True)
        metadata.update(status="complete", completed_cases=len(rows))
        metadata.pop("error", None)
        write_json(metadata_path, metadata)
        print(f"Complete: {len(rows)} cases in {output}", flush=True)
    except BaseException as exc:
        metadata.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                        error=f"{type(exc).__name__}: {exc}", completed_cases=len(rows))
        write_json(metadata_path, metadata)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grids", type=int, nargs="+", default=[50, 100, 200])
    parser.add_argument("--radii", type=float, nargs="+", default=[.85, .98], help="Strictly pre-contact ratios")
    parser.add_argument("--realizations", type=int, default=1)
    parser.add_argument("--seed", type=int, default=9417)
    parser.add_argument("--length", type=float, default=100.)
    parser.add_argument("--lc", type=float, default=10.)
    parser.add_argument("--rho", type=float, default=1.5)
    parser.add_argument("--reference-n", type=int, default=100)
    parser.add_argument("--reference-padding", type=int, default=20)
    parser.add_argument("--mode-no", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--slab", type=int, default=16, help="Sampling slab depth; does not alter the field")
    parser.add_argument("--topology-only", action="store_true", help="Skip scalar meshes; retain Euler/VF/Crofton/contacts")
    parser.add_argument("--output", type=Path, default=HERE/"results_resolution")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if (any(n < 4 for n in args.grids) or args.realizations < 1 or args.seed < 0
            or min(args.reference_n, args.mode_no, args.threads, args.slab) < 1
            or args.reference_padding < 0
            or any(not np.isfinite(v) or v <= 0 for v in (args.length, args.lc, args.rho))
            or any(not np.isfinite(v) or not 0 < v < 1 for v in args.radii)):
        parser.error("Require grids>=4, positive scales/counts, nonnegative seed/padding and 0<r/rc<1.")
    config = {k: getattr(args, k) for k in ("realizations", "seed", "length", "lc", "rho",
              "reference_n", "reference_padding", "mode_no", "threads", "slab")}
    config.update(grids=sorted(set(args.grids)), radii=sorted(set(args.radii)),
                  origin=args.reference_padding*args.length/args.reference_n,
                  include_mesh=not args.topology_only)
    try:
        run(config, args.output, args.resume)
    except ValueError as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    main()
