"""Experiment 07: multiscale margin unions for aggregate grading.

Run from any directory: python run.py --help
Only this experiment's output directory is written. Existing outputs require
--overwrite. Gaussian draws are reused over the total-VF sweep.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import scipy
from scipy import ndimage
import skimage
from skimage.measure import euler_number

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from src import adler
from model import (GaussianSampler, class_probabilities, component_table,
                   cube_lkcs, intersection_operator, sample_margin, topology,
                   union_expectations)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config.json")
    parser.add_argument("--output", type=Path, default=HERE / "results")
    parser.add_argument("--n", type=int, help="Grid side; physical length stays fixed.")
    parser.add_argument("--realizations", type=int)
    parser.add_argument("--lcs", type=float, nargs="+", help="One Lc per class, coarse first.")
    parser.add_argument("--q", type=int, help="Override competitor count for every class.")
    parser.add_argument("--phis", type=float, nargs="+", help="Total aggregate VF values.")
    parser.add_argument("--theory-only", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def load_config(args):
    config = json.loads(args.config.read_text(encoding="utf-8"))
    for key in ("n", "realizations"):
        if getattr(args, key) is not None:
            config[key] = getattr(args, key)
    if args.phis is not None:
        config["target_phis"] = sorted(set(args.phis))
    classes = config["classes"]
    if args.lcs is not None:
        if len(args.lcs) != len(classes):
            raise ValueError("--lcs needs one value per configured class.")
        for item, lc in zip(classes, args.lcs):
            item["lc"] = lc
    if args.q is not None:
        for item in classes:
            item["q"] = args.q
    if (not isinstance(config["n"], int) or config["n"] < 8
            or not isinstance(config["realizations"], int) or config["realizations"] < 1
            or not np.isfinite(config["length"]) or config["length"] <= 0
            or not np.isfinite(config["padding_length"]) or config["padding_length"] < 0
            or not isinstance(config["seed"], int) or config["seed"] < 0):
        raise ValueError("Invalid grid, padding, seed, or realization count.")
    if not 0 <= config["reference_class"] < len(classes):
        raise ValueError("reference_class must select one configured class.")
    if not classes or len(classes) > 255 or len({c["name"] for c in classes}) != len(classes):
        raise ValueError("Use 1-255 uniquely named classes.")
    for item in classes:
        if (not np.isfinite(item["lc"]) or item["lc"] <= 0
                or not isinstance(item["q"], int) or item["q"] < 2):
            raise ValueError("Each class requires lc > 0 and integer q >= 2.")
    lcs = [c["lc"] for c in classes]
    if any(a < b for a, b in zip(lcs, lcs[1:])):
        raise ValueError("List classes in descending Lc order for coarse-first ownership.")
    if not config["target_phis"]:
        raise ValueError("Specify at least one target VF.")
    config["target_phis"] = sorted(set(config["target_phis"]))
    for phi in config["target_phis"]:
        class_probabilities(phi, [c["share"] for c in classes])
    return config


def write_csv(path, rows):
    if not rows:
        return
    # Some measurement rows also carry winner or ownership diagnostics.
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def add_theory(row, lkcs, volume):
    row.update({f"theory_L{j}": float(value) for j, value in enumerate(lkcs)})
    row["theory_phi"] = float(lkcs[3] / volume)
    return row


def make_theory(config, spectral):
    classes = config["classes"]
    length = config["length"]
    table, plans = [], {}
    for phi in config["target_phis"]:
        probabilities = class_probabilities(phi, [c["share"] for c in classes])
        deltas = [adler.delta_for_vf(c["q"], p) for c, p in zip(classes, probabilities)]
        gmfs = [adler.margin_gmfs(c["q"], d) for c, d in zip(classes, deltas)]
        operators = [intersection_operator(m, s["rho"]) for m, s in zip(gmfs, spectral)]
        stages = union_expectations(operators, length)
        approximate_euler = 0.0
        for k, (item, delta, probability, stage) in enumerate(zip(classes, deltas, probabilities, stages)):
            previous_phi = phi * sum(c["share"] for c in classes[:k])
            approximate_euler += stage["raw"][0] * (1 - previous_phi)
            base = {"target_phi": phi, "class_index": k, "class": item["name"],
                    "lc": item["lc"], "q": item["q"], "delta": float(delta),
                    "raw_target_phi": float(probability),
                    "added_target_phi": phi * item["share"],
                    "lambda2_spectral": spectral[k]["lambda2"],
                    "lambda2_nominal": 2 / item["lc"]**2}
            for kind in ("raw", "intersection", "union"):
                row = add_theory({**base, "kind": kind}, stage[kind], length**3)
                row["paper_approximation_L0"] = float(approximate_euler) if kind == "union" else ""
                table.append(row)
        ref = config["reference_class"]
        ref_delta = adler.delta_for_vf(classes[ref]["q"], phi)
        ref_lkcs = intersection_operator(adler.margin_gmfs(classes[ref]["q"], ref_delta),
                                         spectral[ref]["rho"]) @ cube_lkcs(length)
        plans[phi] = {"deltas": deltas, "stages": stages, "reference_delta": ref_delta,
                      "reference_lkcs": ref_lkcs, "approximate_euler": approximate_euler}
        table.append(add_theory({"target_phi": phi, "kind": "reference",
                                "class_index": ref, "class": classes[ref]["name"],
                                "lc": classes[ref]["lc"], "q": classes[ref]["q"],
                                "delta": ref_delta}, ref_lkcs, length**3))
    return table, plans


def measure(mask, common, expected, length, *, owner=None, dx=None, n_classes=0):
    row = add_theory(dict(common), expected, length**3)
    row["measured_phi"] = float(mask.mean())
    row["measured_L3"] = float(mask.sum() * (length / mask.shape[0])**3)
    result6, labels, sizes = topology(mask, 6, return_labels=True)
    row.update(result6)
    row.update(topology(mask, 26))
    components = []
    if owner is not None:
        components = [{**common, **r} for r in component_table(labels, sizes, owner, dx, n_classes)]
        occupied = int(mask.sum())
        mixed_voxels = sum(r["voxels"] for r in components if r["mixed_classes"])
        row["mixed_class_component_volume_fraction"] = mixed_voxels / occupied if occupied else 0.0
    return row, components


def winner_diagnostics(mask, winners, q, dx, common):
    """Only for a raw single-scale phase; never discard cross-scale union links."""
    e6 = e26 = 0
    grains = []
    for winner in range(q):
        one = mask & (winners == winner)
        e6 += int(euler_number(one, connectivity=1))
        e26 += int(euler_number(one, connectivity=3))
        labels, _ = ndimage.label(one)
        sizes = np.bincount(labels.ravel())
        boundary = set()
        for axis in range(3):
            boundary.update(np.unique(np.take(labels, 0, axis=axis)))
            boundary.update(np.unique(np.take(labels, -1, axis=axis)))
        for label in range(1, len(sizes)):
            volume = sizes[label] * dx**3
            grains.append({**common, "winner": winner, "component": label,
                           "volume": float(volume),
                           "equivalent_diameter": float((6 * volume / np.pi)**(1/3)),
                           "touches_boundary": label in boundary})
    return {"winner_euler_6": e6, "winner_euler_26": e26,
            "winner_components_6": len(grains)}, grains


def plot_results(config, rows, class_rows, components, samples, out):
    phis = config["target_phis"]
    last = len(config["classes"]) - 1
    full = [r for r in rows if r["kind"] == "union" and r["class_index"] == last]
    reference = [r for r in rows if r["kind"] == "reference"]

    def curve(selected, key):
        values = [[r[key] for r in selected if r["target_phi"] == p] for p in phis]
        mean = np.array([np.mean(v) for v in values])
        sem = np.array([np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0 for v in values])
        return mean, sem

    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    ax = axes[0, 0]
    mean, sem = curve(full, "measured_phi")
    ax.errorbar(phis, mean, sem, fmt="o-", label="Union mean ± SEM")
    ax.plot(phis, phis, "k--", label="Exact expected VF")
    ax.set(title="Total aggregate volume fraction", ylabel="Measured VF")
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    for conn, marker in ((6, "o"), (26, "s")):
        mean, sem = curve(full, f"euler_{conn}")
        ax.errorbar(phis, mean, sem, fmt=marker+"-", label=f"Union, {conn}-connectivity")
    ax.plot(phis, curve(full, "theory_L0")[0], "k--", label="Conditional GKF")
    ax.plot(phis, curve(full, "paper_approximation_L0")[0], ":", color="gray", label="Paper approximation")
    ax.axhline(0, color="gray", lw=.6)
    ax.set(title="Euler characteristic", ylabel="L0")
    ax.legend(fontsize=8)
    ax = axes[0, 2]
    for selected, label, style in ((full, "Union", "o-"), (reference, "Single-scale reference", "s--")):
        ax.plot(phis, curve(selected, "spans_6")[0], style, label=label)
    ax.plot(phis, curve(full, "spans_26")[0], "^:", label="Union, 26-connectivity")
    ax.set(title="Observed spanning frequency", ylabel="Fraction of realizations", ylim=(-.04, 1.04))
    ax.legend(fontsize=8)
    ax = axes[1, 0]
    for selected, label, style in ((full, "Union", "o-"), (reference, "Single-scale reference", "s--")):
        mean, sem = curve(selected, "largest_fraction_6")
        ax.errorbar(phis, mean, sem, fmt=style, label=label)
    ax.set(title="Largest connected component (6)", ylabel="Fraction of aggregate volume", ylim=(0, 1.04))
    ax.legend(fontsize=8)
    ax = axes[1, 1]
    for k, item in enumerate(config["classes"]):
        selected = [r for r in class_rows if r["class_index"] == k]
        mean, sem = curve(selected, "measured_added_phi")
        line = ax.errorbar(phis, mean, sem, fmt="o-", label=item["name"])
        ax.plot(phis, np.asarray(phis)*item["share"], "--", color=line[0].get_color())
    ax.set(title="Added class volumes (dashed = target)", ylabel="Fraction of domain")
    ax.legend(fontsize=8)
    ax = axes[1, 2]
    for phi in phis:
        selected = [r for r in components if r["target_phi"] == phi]
        selected.sort(key=lambda r: r["equivalent_diameter"])
        volumes = np.array([r["volume"] for r in selected])
        if volumes.size and volumes.sum() > 0:
            ax.step([r["equivalent_diameter"] for r in selected], np.cumsum(volumes)/volumes.sum(),
                    where="post", label=f"VF {phi:.0%}")
    ax.set(title="Actual union component sizes (6)", xlabel=f"Equivalent diameter [{config['length_unit']}]",
           ylabel="Cumulative aggregate volume", xscale="log")
    ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.25)
        if ax is not axes[1, 2]:
            ax.set_xlabel("Target total aggregate VF")
    fig.suptitle(f"Multiscale margin union — N={config['n']}, {config['realizations']} realizations\n"
                 "Component diameter is not a sieve diameter; boundary-cut and spanning components included")
    fig.savefig(out / "multiscale_summary.png", dpi=150)
    plt.close(fig)

    cmap = ListedColormap(["#f2f0e8", "#5c6f82", "#c18449", "#ddbd6a"])
    if len(config["classes"]) != 3:
        cmap = plt.get_cmap("tab20", len(config["classes"]) + 1)
    fig, axes = plt.subplots(len(phis), 3, figsize=(10, 3.0*len(phis)), squeeze=False,
                             constrained_layout=True)
    for i, phi in enumerate(phis):
        for axis in range(3):
            ax = axes[i, axis]
            ax.imshow(samples[phi][axis], origin="lower", cmap=cmap, interpolation="nearest",
                      vmin=0, vmax=len(config["classes"]), extent=(0, config["length"], 0, config["length"]))
            ax.set_title(f"VF {phi:.0%} — normal to {'xyz'[axis]}", fontsize=9)
            ax.set_xlabel(config["length_unit"])
    fig.suptitle("Realization 1: matrix and coarse-first class ownership\n"
                 "Colours identify added material classes, not disconnected grains")
    fig.savefig(out / "class_slices.png", dpi=150)
    plt.close(fig)


def write_summary(config, rows, out, elapsed):
    last = len(config["classes"]) - 1
    full = [r for r in rows if r["kind"] == "union" and r["class_index"] == last]
    lines = ["# Multiscale margin pilot results", "", f"Grid: {config['n']}³; domain side: {config['length']} {config['length_unit']}.",
             f"Realizations: {config['realizations']}; elapsed: {elapsed:.1f} s.", "",
             "| Target VF | Mean VF | Mean L0 (6 / 26) | GKF L0 | Spanning (6 / 26) | Mean largest fraction (6) |",
             "| --- | --- | --- | --- | --- | --- |"]
    for phi in config["target_phis"]:
        selected = [r for r in full if r["target_phi"] == phi]
        mean = lambda key: float(np.mean([r[key] for r in selected]))
        span = lambda key: sum(bool(r[key]) for r in selected)
        lines.append(f"| {phi:.0%} | {mean('measured_phi'):.3%} | {mean('euler_6'):.1f} / {mean('euler_26'):.1f} | "
                     f"{mean('theory_L0'):.1f} | {span('spans_6')}/{len(selected)} / {span('spans_26')}/{len(selected)} | "
                     f"{mean('largest_fraction_6'):.3%} |")
    lines += ["", "This is a finite-resolution experiment, not validation of a concrete mixture.",
              "The sand class is a scale proxy. A grading curve requires calibration of measured grain sizes.",
              "Coarse-first ownership recovers class-volume accounting but does not remove physical contacts.",
              "GKF predictions inherit the margin GMFs' validation status; 6/26 disagreement diagnoses discretization sensitivity.",
              "Spanning is measured on the actual union with bounded face-to-face connectivity.",
              "", "See metadata.json for the effective configuration and source hashes."]
    (out / "summary.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


def main():
    started = time.perf_counter()
    args = parse_args()
    config = load_config(args)
    out = args.output.resolve()
    if out.exists() and any(out.iterdir()) and not args.overwrite:
        raise SystemExit(f"Output directory is not empty: {out}. Choose a new --output or use --overwrite.")
    out.mkdir(parents=True, exist_ok=True)
    classes = config["classes"]
    n, length = config["n"], config["length"]
    dx = length / n
    # Store tiny spectral metadata, not all class filters simultaneously.
    spectral = []
    for item in classes:
        sampler = GaussianSampler(n, length, item["lc"], config["padding_length"])
        spectral.append({"full_n": sampler.full_n, "lambda2": sampler.lambda2,
                         "rho": sampler.rho, "forward_difference_variance": sampler.fd_variance})
        del sampler
    table, plans = make_theory(config, spectral)
    write_csv(out / "theory.csv", table)
    source_paths = [HERE / "run.py", HERE / "model.py", ROOT / "src" / "adler.py"]
    metadata = {"config": config, "spectral": spectral, "mode": "theory_only" if args.theory_only else "simulation",
                "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                "skimage": skimage.__version__, "generator": "real-white FFT; ensemble-normalized; cropped; Nyquist removed",
                "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}}
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    print(f"Experiment 07: N={n}, L={length}, Lcs={[c['lc'] for c in classes]}, output={out}", flush=True)
    if min(c["lc"] for c in classes) / dx < 4:
        print("Pilot resolution: the smallest Lc is below 4 voxels. Do not treat this as continuum validation.", flush=True)
    if args.theory_only:
        print("Saved analytical expectations; no fields generated.")
        return
    rows, class_rows, field_rows, components, grains = [], [], [], [], []
    samples = {}
    for realization in range(config["realizations"]):
        margins, winners = [], []
        for k, item in enumerate(classes):
            sampler = GaussianSampler(n, length, item["lc"], config["padding_length"])
            margin, winner, diagnostics = sample_margin(sampler, item["q"], config["seed"], realization, k)
            margins.append(margin)
            winners.append(winner)
            field_rows.extend({**r, "class": item["name"]} for r in diagnostics)
            del sampler
            print(f"  realization {realization+1}: generated {item['name']}", flush=True)
        for phi in config["target_phis"]:
            plan = plans[phi]
            union = np.zeros((n,)*3, dtype=bool)
            owner = np.zeros((n,)*3, dtype=np.uint8)
            for k, item in enumerate(classes):
                raw = margins[k] >= plan["deltas"][k]
                common = {"realization": realization+1, "target_phi": phi, "class_index": k,
                          "class": item["name"], "lc": item["lc"], "delta": plan["deltas"][k]}
                raw_row, _ = measure(raw, {**common, "kind": "raw"}, plan["stages"][k]["raw"], length)
                diag, raw_grains = winner_diagnostics(raw, winners[k], item["q"], dx, common)
                raw_row.update(diag)
                rows.append(raw_row)
                grains.extend(raw_grains)
                if k:
                    inter_row, _ = measure(union & raw, {**common, "kind": "intersection"},
                                           plan["stages"][k]["intersection"], length)
                    rows.append(inter_row)
                added = raw & ~union
                owner[added] = k+1
                class_rows.append({**common, "target_added_phi": phi*item["share"],
                                   "raw_target_phi": float(adler.margin_m0(item["q"], plan["deltas"][k])),
                                   "measured_raw_phi": float(raw.mean()), "measured_added_phi": float(added.mean()),
                                   "measured_overlap_phi": float(np.mean(raw & union))})
                union |= raw
                union_row, component_rows = measure(
                    union, {**common, "kind": "union"}, plan["stages"][k]["union"], length,
                    owner=owner if k == len(classes)-1 else None, dx=dx, n_classes=len(classes))
                if k == len(classes)-1:
                    union_row["paper_approximation_L0"] = float(plan["approximate_euler"])
                    components.extend(component_rows)
                rows.append(union_row)
            ref = config["reference_class"]
            ref_mask = margins[ref] >= plan["reference_delta"]
            reference, _ = measure(ref_mask, {"realization": realization+1, "target_phi": phi,
                                              "class_index": ref, "kind": "reference", "class": classes[ref]["name"]},
                                   plan["reference_lkcs"], length)
            rows.append(reference)
            print(f"    VF {phi:.0%}: measured {union_row['measured_phi']:.3%}; "
                  f"L0(6/26)={union_row['euler_6']}/{union_row['euler_26']}; "
                  f"theory={union_row['theory_L0']:.1f}; spans6={union_row['spans_6']}; "
                  f"largest6={union_row['largest_fraction_6']:.1%}", flush=True)
            if realization == 0:
                samples[phi] = [np.take(owner, n//2, axis=axis).copy() for axis in range(3)]
                if config.get("save_volumes", False):
                    np.savez_compressed(out / f"sample_vf_{phi:.4f}.npz", phase=union, class_owner=owner,
                                        length=length, dx=dx, target_phi=phi,
                                        deltas=plan["deltas"], class_names=[c["name"] for c in classes])
        # Checkpoint completed realizations so a long run retains usable data.
        write_csv(out / "measurements.csv", rows)
        write_csv(out / "class_volumes.csv", class_rows)
        write_csv(out / "field_diagnostics.csv", field_rows)
        write_csv(out / "union_components.csv", components)
        write_csv(out / "raw_winner_components.csv", grains)
    plot_results(config, rows, class_rows, components, samples, out)
    elapsed = time.perf_counter() - started
    metadata["elapsed_seconds"] = elapsed
    metadata["completed_realizations"] = config["realizations"]
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    write_summary(config, rows, out, elapsed)
    print(f"Done in {elapsed:.1f}s. See {out / 'summary.md'}", flush=True)


if __name__ == "__main__":
    main()
