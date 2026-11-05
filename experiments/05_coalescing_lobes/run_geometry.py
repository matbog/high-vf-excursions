"""Experiment 05: pre-contact tetrahedron geometry sweep, L0/L3 and grain diagnostics.

Default N=100, rho=1.5/2/2.5, r/rc=.85/.95/.98, four realizations.
One Gaussian vector realization is reused for every geometry and radius.
Only voxel topology/volume are measured. No scalar meshes or full-grid exports.
"""
from __future__ import annotations

import argparse
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
from scipy import ndimage

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
import run_resolution as resolution
from src import adler, excursions, measures
from src.fields import GridSpec

SCHEMA = "lobes_geometry_v1"


def component_statistics(mask, dx, connectivity):
    """Bounded digital components; diameter is volume-equivalent, not sieve size."""
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
    summary = dict(euler=int(measures.measure_euler(mask, rank)), components=int(count),
                   occupied_voxels=occupied, volume=float(occupied*dx**3),
                   largest_volume=float(sizes.max()*dx**3),
                   largest_fraction=float(sizes.max()/occupied) if occupied else 0.,
                   spans=bool(spanning.any()),
                   boundary_volume=float(sizes[boundary].sum()*dx**3),
                   spanning_volume=float(sizes[spanning.any(axis=0)].sum()*dx**3))
    summary.update({f"spans_{axis}": bool(spanning[j].any()) for j, axis in enumerate("xyz")})
    rows = []
    for label in range(1, count+1):
        volume = float(sizes[label]*dx**3)
        rows.append(dict(component=label, voxels=int(sizes[label]), volume=volume,
                         equivalent_diameter=float((6*volume/np.pi)**(1/3)),
                         touches_boundary=bool(boundary[label]), spans=bool(spanning[:, label].any()),
                         **{f"spans_{axis}": bool(spanning[j, label]) for j, axis in enumerate("xyz")}))
    return summary, rows


def measure_geometry(distances, radius, grid):
    """Measure the union and four disjoint ball masks, with both conventions."""
    masks = distances <= radius
    occupancy = masks.sum(axis=0)
    if np.any(occupancy > 1):
        raise ValueError("Pre-contact measurement requires disjoint ball memberships at each sample.")
    phase = occupancy > 0
    volume = float(phase.sum()*grid.dx**3)
    result = dict(measured_phi=volume/grid.volume, union_L3=volume)
    result.update(resolution.cross_label_contacts(phase, distances.argmin(axis=0)))
    ball_rows, component_rows = [], []
    for connectivity in (6, 26):
        union, _ = component_statistics(phase, grid.dx, connectivity)
        for key in ("euler", "components", "largest_fraction", "spans", "spans_x", "spans_y", "spans_z"):
            result[f"union_{key}_{connectivity}"] = union[key]
        balls = []
        for ball, mask in enumerate(masks):
            stats, components = component_statistics(mask, grid.dx, connectivity)
            balls.append(stats)
            ball_rows.append(dict(ball=ball, connectivity=connectivity,
                                  vf=stats["volume"]/grid.volume, **stats))
            component_rows.extend(dict(ball=ball, connectivity=connectivity, **r) for r in components)
        total = sum(r["volume"] for r in balls)
        if not np.isclose(total, volume, rtol=1e-12, atol=1e-10):
            raise RuntimeError("Pre-contact volume accounting failed.")
        result.update({
            f"sum_ball_L0_{connectivity}": sum(r["euler"] for r in balls),
            f"sum_ball_components_{connectivity}": sum(r["components"] for r in balls),
            f"sum_ball_L3_{connectivity}": total,
            f"any_ball_spans_{connectivity}": any(r["spans"] for r in balls),
            f"spanning_balls_{connectivity}": sum(r["spans"] for r in balls),
            f"largest_separate_component_fraction_{connectivity}": (
                max(r["largest_volume"] for r in balls)/volume if volume else 0.),
            f"max_within_ball_largest_fraction_{connectivity}": max(r["largest_fraction"] for r in balls),
            f"separate_spanning_volume_fraction_{connectivity}": (
                sum(r["spanning_volume"] for r in balls)/volume if volume else 0.),
            f"separate_boundary_volume_fraction_{connectivity}": (
                sum(r["boundary_volume"] for r in balls)/volume if volume else 0.),
            f"euler_additivity_gap_{connectivity}": union["euler"]-sum(r["euler"] for r in balls),
        })
    return result, ball_rows, component_rows


def theory_table(config):
    grid = GridSpec(config["n"], config["length"])
    rows = []
    for rho in config["rhos"]:
        centres = excursions.tetrahedral_centres(rho)
        rc = float(np.linalg.norm(centres[0]-centres[1])/2)
        for ratio in config["radii"]:
            theory = adler.expected_lkcs_lobes_disjoint(grid, config["lc"], centres, ratio*rc)
            rows.append(dict(rho=rho, radius_over_rc=ratio, radius=ratio*rc, contact_radius=rc,
                             theory_phi=theory["phi"],
                             **{f"theory_L{i}": theory[f"L{i}"] for i in range(4)}))
    return rows


def expectations_key(row):
    return float(row["rho"]), float(row["radius_over_rc"])


def case_key(row):
    return int(row["realization"]), *expectations_key(row)


def case_paths(config, output):
    return {(r, rho, ratio): output/"cases"/f"r{r:04d}_rho{i:03d}_ratio{j:03d}.json"
            for r in range(1, config["realizations"]+1)
            for i, rho in enumerate(config["rhos"]) for j, ratio in enumerate(config["radii"])}


def load_cases(paths):
    records = {}
    for key, path in paths.items():
        if path.exists():
            record = json.loads(path.read_text(encoding="utf-8"))
            if (record.get("status") != "complete" or case_key(record["measurement"]) != key
                    or len(record["balls"]) != 8):
                raise ValueError(f"Invalid case checkpoint: {path}")
            records[key] = record
    if paths:
        folder = next(iter(paths.values())).parent
        if folder.exists() and set(folder.glob("*.json")) - set(paths.values()):
            raise ValueError("Unexpected case checkpoint files; choose a clean output directory.")
    return records


def aggregate(records):
    rows = [r["measurement"] for _, r in sorted(records.items())]
    out = []
    for rho, ratio in sorted({expectations_key(r) for r in rows}):
        group = [r for r in rows if expectations_key(r) == (rho, ratio)]
        count = len(group)
        mean = lambda key: float(np.mean([r[key] for r in group]))
        sem = lambda key: float(np.std([r[key] for r in group], ddof=1)/np.sqrt(count)) if count > 1 else float("nan")
        item = dict(rho=rho, radius_over_rc=ratio, realizations=count,
                    theory_phi=group[0]["theory_phi"], theory_L0=group[0]["theory_L0"],
                    measured_phi_mean=mean("measured_phi"), measured_phi_sem=sem("measured_phi"))
        for rank in (6, 26):
            for key in (f"sum_ball_L0_{rank}", f"any_ball_spans_{rank}",
                        f"largest_separate_component_fraction_{rank}", f"euler_additivity_gap_{rank}",
                        f"cross_label_fraction_{rank}", f"separate_spanning_volume_fraction_{rank}"):
                item[f"{key}_mean"] = mean(key)
            item[f"sum_ball_L0_{rank}_sem"] = sem(f"sum_ball_L0_{rank}")
        out.append(item)
    return rows, out


def write_outputs(output, records, config, theory, complete):
    rows, stats = aggregate(records)
    balls = [b for _, record in sorted(records.items()) for b in record["balls"]]
    components = [c for _, record in sorted(records.items()) for c in record["components"]]
    for filename, data in (("measurements.csv", rows), ("balls.csv", balls),
                           ("components.csv", components), ("aggregate.csv", stats)):
        resolution.write_csv(output/filename, data)
    lines = ["# Pre-contact tetrahedron geometry sweep", "",
             f"Status: {'complete' if complete else 'partial'}; {len(rows)} completed 3D cases.",
             "rho is the value-space centre radius; r/rc is the morphology control after rho is fixed.",
             "The four balls share one Gaussian vector: they are not four independent realizations.", "",
             "| rho | r/rc | repetitions | expected VF | measured VF | expected L0 | ball sum L0, 6 / 26 | any ball spans, 6 / 26 |",
             "|---|---|---|---|---|---|---|---|"]
    for r in stats:
        lines.append(f"| {r['rho']:g} | {r['radius_over_rc']:g} | {r['realizations']} | "
                     f"{r['theory_phi']:.5f} | {r['measured_phi_mean']:.5f} | {r['theory_L0']:.2f} | "
                     f"{r['sum_ball_L0_6_mean']:.2f} / {r['sum_ball_L0_26_mean']:.2f} | "
                     f"{r['any_ball_spans_6_mean']:.2f} / {r['any_ball_spans_26_mean']:.2f} |")
    if config["theory_only"]:
        lines += ["", "Theory-only run: no 3D fields were generated. See theory.csv."]
    lines += ["", "Spanning frequencies use independent vector realizations as trials.",
              "Zero spans in four trials is weak evidence, not proof of a disconnected continuum morphology.",
              "Euler is not a component count. Ball sums are continuum-additive before contact,",
              "but their digital estimators can still miss within-ball topology.",
              "Components include boundary-cut and spanning objects; flags permit explicit filtering.",
              "Equivalent diameter is volume-based, not a sieve diameter.",
              "A single N=100 sweep is screening, not a resolution or continuum validation.",
              "aggregate.csv uses standard errors across realizations; SEM is NaN with one realization.", ""]
    (output/"summary.md").write_text("\n".join(lines), encoding="utf-8")


def plot_summary(output, records, config):
    _, stats = aggregate(records)
    if not stats:
        return
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    colors = plt.colormaps["tab10"]
    for i, rho in enumerate(config["rhos"]):
        group = sorted([r for r in stats if r["rho"] == rho], key=lambda r: r["radius_over_rc"])
        if not group:
            continue
        color = colors(i % 10)
        ratios = [r["radius_over_rc"] for r in group]
        phis = [r["theory_phi"] for r in group]
        axes[0, 0].plot(ratios, phis, "--", color=color)
        axes[0, 0].plot(ratios, [r["measured_phi_mean"] for r in group], "o-", color=color, label=f"rho={rho:g}")
        axes[0, 1].plot(phis, [r["theory_L0"] for r in group], "-", color=color, alpha=.4)
        for rank, marker in ((6, "s--"), (26, "o:")):
            label = f"rho={rho:g}, {rank} neighbours"
            for ax, key in ((axes[0, 1], f"sum_ball_L0_{rank}_mean"),
                            (axes[0, 2], f"any_ball_spans_{rank}_mean"),
                            (axes[1, 0], f"largest_separate_component_fraction_{rank}_mean"),
                            (axes[1, 1], f"euler_additivity_gap_{rank}_mean")):
                ax.plot(phis, [r[key] for r in group], marker, color=color, label=label)
        selected = [c for record in records.values() for c in record["components"]
                    if c["rho"] == rho and c["radius_over_rc"] == max(config["radii"]) and c["connectivity"] == 6]
        if selected:
            selected.sort(key=lambda c: c["equivalent_diameter"])
            weights = np.array([c["volume"] for c in selected])
            axes[1, 2].step([c["equivalent_diameter"] for c in selected],
                            np.cumsum(weights)/weights.sum(), where="post", color=color, label=f"rho={rho:g}")
    axes[0, 0].set(title="VF: measured means (solid), theory (dashed)", xlabel="r/rc", ylabel="VF")
    axes[0, 1].set(title="Ball-sum Euler; pale lines: GKF expectation", xlabel="Expected VF", ylabel="L0")
    axes[0, 2].set(title="At least one ball spans", xlabel="Expected VF", ylabel="Frequency across realizations", ylim=(-.05, 1.05))
    axes[1, 0].set(title="Largest separate component / total phase", xlabel="Expected VF", ylabel="Volume fraction", ylim=(0, 1))
    axes[1, 1].set(title="Union Euler minus ball-sum Euler", xlabel="Expected VF", ylabel="Additivity gap")
    axes[1, 1].axhline(0, color="black", lw=.6)
    axes[1, 2].set(title=f"Volume-weighted sizes at r/rc={max(config['radii']):g}\n6-neighbour, includes boundary/spanning components",
                   xlabel="Equivalent diameter", ylabel="Cumulative volume fraction", ylim=(0, 1))
    for ax in axes.flat:
        ax.grid(alpha=.25)
        ax.legend(fontsize=6)
    fig.suptitle(f"Geometry screening: N={config['n']}, up to {config['realizations']} independent vector realizations")
    fig.savefig(output/"geometry_summary.png", dpi=160)
    plt.close(fig)


def run(config, output, resume=False):
    hashes = resolution.source_hashes()
    hashes[Path(__file__).relative_to(ROOT).as_posix()] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    versions = {name: version(name) for name in ("numpy", "scipy", "gstools", "scikit-image", "matplotlib")}
    occupied = output.exists() and any(output.iterdir())
    metadata_path = output/"metadata.json"
    if occupied:
        if not resume:
            raise ValueError(f"Output is not empty: {output}. Use --resume or a new --output.")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if (metadata.get("schema") != SCHEMA or metadata.get("config") != config
                or metadata.get("source_sha256") != hashes or metadata.get("versions") != versions):
            raise ValueError("Resume requires matching schema, parameters, source hashes and dependency versions.")
    else:
        if resume:
            raise ValueError("--resume requires existing checkpoints.")
        output.mkdir(parents=True, exist_ok=True)
        metadata = dict(schema=SCHEMA, config=config, source_sha256=hashes, versions=versions,
                        python=sys.version, sampler="GSTools RandMeth; shared fields across all rho and radius values",
                        component_seeds={str(r): resolution.component_seeds(config["seed"], r)
                                         for r in range(1, config["realizations"]+1)})
    paths = {} if config["theory_only"] else case_paths(config, output)
    records = load_cases(paths)
    diagnostics_path = output/"field_diagnostics.csv"
    diagnostics = resolution.read_csv(diagnostics_path)
    theory = theory_table(config)
    predictions = {expectations_key(row): row for row in theory}
    resolution.write_csv(output/"theory.csv", theory)
    metadata.update(status="running", expected_cases=len(paths), completed_cases=len(records))
    resolution.write_json(metadata_path, metadata)
    write_outputs(output, records, config, theory, complete=False)
    try:
        for realization in range(1, config["realizations"]+1):
            pending = [key for key in paths if key[0] == realization and key not in records]
            if not pending:
                continue
            print(f"Realization {realization}/{config['realizations']}: generating one vector field for {len(pending)} cases", flush=True)
            start = time.perf_counter()
            vector = resolution.sample_vector(config["n"], config, realization)
            elapsed = time.perf_counter()-start
            diagnostics = [r for r in diagnostics if int(r["realization"]) != realization]
            for c, seed in enumerate(resolution.component_seeds(config["seed"], realization)):
                diagnostics.append(dict(realization=realization, component=c, seed=seed,
                    mean=float(vector[c].mean(dtype=np.float64)),
                    second_moment=float(np.mean(vector[c]**2, dtype=np.float64)), sampling_seconds=elapsed))
            resolution.write_csv(diagnostics_path, diagnostics)
            print(f"  field ready in {elapsed:.1f}s; reused for every centre/radius setting", flush=True)
            grid = GridSpec(config["n"], config["length"])
            for rho in config["rhos"]:
                selected = [key for key in pending if key[1] == rho]
                if not selected:
                    continue
                centres = excursions.tetrahedral_centres(rho)
                distances = resolution.ball_distances(vector, centres)
                for key in selected:
                    ratio = key[2]
                    prediction = predictions[(rho, ratio)]
                    start = time.perf_counter()
                    measured, balls, components = measure_geometry(distances, prediction["radius"], grid)
                    identity = dict(realization=realization, N=config["n"], rho=rho, radius_over_rc=ratio)
                    row = dict(**identity, **{k: v for k, v in prediction.items() if k not in identity},
                               **measured, measurement_seconds=time.perf_counter()-start)
                    payload = dict(status="complete", measurement=row,
                        balls=[dict(**identity, theory_phi_ball=prediction["theory_phi"]/4,
                                    theory_L0_ball=prediction["theory_L0"]/4, **r) for r in balls],
                        components=[dict(**identity, **r) for r in components])
                    paths[key].parent.mkdir(exist_ok=True)
                    # The complete case is the atomic source of truth; CSVs can be rebuilt on resume.
                    resolution.write_json(paths[key], payload)
                    records[key] = payload
                    metadata["completed_cases"] = len(records)
                    resolution.write_json(metadata_path, metadata)
                    write_outputs(output, records, config, theory, complete=False)
                    print(f"  rho={rho:g}, r/rc={ratio:g}: VF={row['measured_phi']:.4f}, "
                          f"ball L0(6/26)={row['sum_ball_L0_6']}/{row['sum_ball_L0_26']}, "
                          f"ball spans(6/26)={row['any_ball_spans_6']}/{row['any_ball_spans_26']}; saved", flush=True)
                del distances
            del vector
            plot_summary(output, records, config)
        write_outputs(output, records, config, theory, complete=True)
        plot_summary(output, records, config)
        metadata.update(status="complete", completed_cases=len(records))
        metadata.pop("error", None)
        resolution.write_json(metadata_path, metadata)
        print(f"Complete: {len(records)} 3D cases. Outputs: {output}", flush=True)
    except BaseException as exc:
        metadata.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                        completed_cases=len(records), error=f"{type(exc).__name__}: {exc}")
        resolution.write_json(metadata_path, metadata)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--rhos", type=float, nargs="+", default=[1.5, 2., 2.5])
    parser.add_argument("--radii", type=float, nargs="+", default=[.85, .95, .98])
    parser.add_argument("--realizations", type=int, default=4)
    parser.add_argument("--seed", type=int, default=9417)
    parser.add_argument("--length", type=float, default=100.)
    parser.add_argument("--lc", type=float, default=10.)
    parser.add_argument("--reference-n", type=int, default=100)
    parser.add_argument("--reference-padding", type=int, default=20)
    parser.add_argument("--mode-no", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--slab", type=int, default=16)
    parser.add_argument("--theory-only", action="store_true")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if (args.n < 4 or min(args.realizations, args.reference_n, args.mode_no, args.threads, args.slab) < 1
            or args.reference_padding < 0 or args.seed < 0
            or any(not np.isfinite(v) or v <= 0 for v in [args.length, args.lc, *args.rhos])
            or any(not np.isfinite(v) or not 0 < v < 1 for v in args.radii)):
        parser.error("Require n>=4, positive counts/scales, nonnegative seed/padding, and 0<r/rc<1.")
    config = {key: getattr(args, key) for key in ("n", "realizations", "seed", "length", "lc",
                "reference_n", "reference_padding", "mode_no", "threads", "slab", "theory_only")}
    config.update(rhos=sorted(set(args.rhos)), radii=sorted(set(args.radii)),
                  origin=args.reference_padding*args.length/args.reference_n)
    output = args.output or HERE/("results_geometry_theory" if args.theory_only else "results_geometry")
    try:
        run(config, output, args.resume)
    except ValueError as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    main()
