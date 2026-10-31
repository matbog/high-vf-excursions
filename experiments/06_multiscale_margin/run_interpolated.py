"""Experiment 07 only: GSTools components, interpolated margins, volumetric LKCs.

The configured N=100 ensemble is expensive. For a small first run use:
  python run_interpolated.py --n 6 --q 3 --realizations 1 --phis .6 .8 --output results_pl_smoke
N denotes cells per axis; N+1 nodal samples include both domain boundaries.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time

# This process runs alongside other research jobs: keep numerical libraries at
# one thread. These settings do not change another process or the environment.
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import gstools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import scipy
import vtk

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from src import adler
from gstools_fields import attach_components, center_samples, sample_components
from interpolated_geometry import measure_union, tetrahedral_grid
from model import class_probabilities, intersection_operator, margin_from_fields, topology, union_expectations


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def read_config():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=HERE / "config_interpolated.json")
    parser.add_argument("--output", type=Path, default=HERE / "results_interpolated")
    parser.add_argument("--n", type=int, help="Cells per axis, with N+1 boundary-inclusive nodes.")
    parser.add_argument("--q", type=int)
    parser.add_argument("--realizations", type=int)
    parser.add_argument("--mode-no", type=int)
    parser.add_argument("--phis", nargs="+", type=float)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-meshes", action="store_true")
    parser.add_argument("--no-fields", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    for key in ("n", "realizations", "mode_no"):
        if getattr(args, key) is not None:
            config[key] = getattr(args, key)
    if args.phis is not None:
        config["target_phis"] = args.phis
    if args.q is not None:
        for item in config["classes"]:
            item["q"] = args.q
    if args.no_meshes:
        config["save_meshes"] = False
    if args.no_fields:
        config["save_fields"] = False
    if any(not isinstance(config[key], int) or config[key] < 1 for key in ("n", "realizations", "mode_no")):
        raise ValueError("N, realizations and mode_no must be positive integers.")
    if not isinstance(config["seed"], int) or config["seed"] < 0:
        raise ValueError("Require a nonnegative integer seed.")
    if not np.isfinite(config["length"]) or config["length"] <= 0:
        raise ValueError("Require a finite positive domain length.")
    if not 0 <= config["merge_tolerance"] < 1e-5:
        raise ValueError("merge_tolerance must be zero or a small positive fraction of domain length.")
    classes = config["classes"]
    if not 1 <= len(classes) <= 3:
        raise ValueError("The interpolated runner supports one to three classes.")
    if len({c["name"] for c in classes}) != len(classes):
        raise ValueError("Class names must be unique.")
    if any(not np.isfinite(c["lc"]) or c["lc"] <= 0 or not isinstance(c["q"], int) or c["q"] < 2 for c in classes):
        raise ValueError("Require positive Lc and integer q >= 2.")
    if any(a["lc"] < b["lc"] for a, b in zip(classes, classes[1:])):
        raise ValueError("List classes coarse first.")
    config["target_phis"] = sorted(set(config["target_phis"]))
    if not config["target_phis"]:
        raise ValueError("Specify at least one target VF.")
    for phi in config["target_phis"]:
        class_probabilities(phi, [c["share"] for c in classes])
    return args, config


def theoretical_plan(config):
    plans, rows = {}, []
    for phi in config["target_phis"]:
        probabilities = class_probabilities(phi, [c["share"] for c in config["classes"]])
        deltas = [adler.delta_for_vf(c["q"], p) for c, p in zip(config["classes"], probabilities)]
        operators = [intersection_operator(adler.margin_gmfs(c["q"], d), 1/(c["lc"]*np.sqrt(np.pi)))
                     for c, d in zip(config["classes"], deltas)]
        stages = union_expectations(operators, config["length"])
        plans[phi] = {"deltas": deltas, "lkcs": stages[-1]["union"]}
        for k, (item, p, d, stage) in enumerate(zip(config["classes"], probabilities, deltas, stages)):
            for kind, lkcs in stage.items():
                rows.append({"target_phi": phi, "class_index": k, "class": item["name"],
                             "kind": kind, "delta": d, "raw_target_phi": p,
                             **{f"L{j}": float(v) for j, v in enumerate(lkcs)}})
    return plans, rows


def voxel_comparison(center_groups, deltas):
    mask = np.zeros(center_groups[0][0].shape, bool)
    for fields, delta in zip(center_groups, deltas):
        margin, _ = margin_from_fields(fields)
        mask |= margin >= delta
    stats = {"voxel_phi": float(mask.mean())}
    for conn in (6, 26):
        stats.update({"voxel_"+k: v for k, v in topology(mask, conn).items()})
    return stats


def plot_and_summarize(config, rows, out):
    phis = config["target_phis"]
    def mean(key):
        return np.array([np.mean([r[key] for r in rows if r["target_phi"] == p]) for p in phis])
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    axes[0].plot(phis, phis, "k--", label="Continuum target")
    axes[0].plot(phis, mean("phi"), "o-", label="Interpolated volume")
    axes[0].plot(phis, mean("voxel_phi"), "s:", label="Same PL fields, voxel centres")
    axes[0].set(title="Volume fraction", ylabel="Measured VF")
    axes[1].plot(phis, mean("theory_L0"), "k--", label="Continuum GKF")
    axes[1].plot(phis, mean("L0"), "o-", label="Interpolated complex")
    for conn in (6, 26):
        axes[1].plot(phis, mean(f"voxel_euler_{conn}"), ":", label=f"Voxels, connectivity {conn}")
    axes[1].set(title="Euler characteristic", ylabel="L0")
    axes[2].plot(phis, mean("largest_fraction"), "o-", label="Interpolated union")
    axes[2].plot(phis, mean("voxel_largest_fraction_6"), "s:", label="Voxels, connectivity 6")
    axes[2].set(title="Largest connected component", ylabel="Fraction of aggregate volume", ylim=(0, 1.03))
    for ax in axes:
        ax.set_xlabel("Target aggregate VF")
        ax.grid(alpha=.25)
        ax.legend(fontsize=7)
    fig.suptitle(f"GSTools + component interpolation: {config['n']} cells/axis, {config['realizations']} realizations\n"
                 "Finite-grid PL geometry; continuum expectations are reference values")
    fig.savefig(out / "comparison.png", dpi=150)
    plt.close(fig)
    lines = ["# Interpolated prototype (deferred)", "",
             f"Grid: {config['n']} cells and {config['n']+1} nodes per axis; full domain side {config['length']} {config['length_unit']}.",
             f"GSTools RandMeth; {config['mode_no']} modes/component; {config['realizations']} realizations; one thread.", "",
             "| Target VF | Interpolated VF | Centre-voxel VF | L0 interpolated | L0 voxels 6 / 26 | GKF L0 | Largest interpolated | Spanning |",
             "| --- | --- | --- | --- | --- | --- | --- | --- |"]
    for phi in phis:
        selected = [r for r in rows if r["target_phi"] == phi]
        m = lambda key: float(np.mean([r[key] for r in selected]))
        lines.append(f"| {phi:.0%} | {m('phi'):.3%} | {m('voxel_phi'):.3%} | {m('L0'):.1f} | "
                     f"{m('voxel_euler_6'):.1f} / {m('voxel_euler_26'):.1f} | {m('theory_L0'):.1f} | "
                     f"{m('largest_fraction'):.3%} | {sum(r['spans'] for r in selected)}/{len(selected)} |")
    lines += ["", "All reported mesh terms passed closed-boundary and Euler-consistency audits.",
              "The geometric calculations describe the piecewise-linear component interpolant, not the exact smooth GSTools field.",
              "Coarse interpolation reduces interior variance and can miss oscillations; the continuum VF and GKF values need not match at this resolution.",
              "The voxel comparison uses the same interpolated components at cell centres, so backend and seed are held fixed.",
              "", "Open the exported .vtm files in ParaView to display the union of raw class/winner volumes.",
              "These blocks overlap across classes. Their volumes must not be summed directly; the CSV measurements include intersection corrections.",
              "No single conforming Boolean-union mesh is exported. L1/L2 are not measured in this implementation.",
              "", "See INTERPOLATED.md for methods, limitations and refinement commands."]
    (out / "summary.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


def main():
    args, config = read_config()
    out = args.output.resolve()
    # Keep every new result inside this experiment directory.
    if not out.is_relative_to(HERE) or out == HERE:
        raise ValueError("Choose an output subdirectory inside experiments/06_multiscale_margin.")
    if out.exists() and any(out.iterdir()):
        if not args.overwrite:
            raise SystemExit(f"Output is not empty: {out}; use a new --output or --overwrite.")
        existing = out/"metadata.json"
        if (not existing.is_file() or json.loads(existing.read_text()).get("geometry")
                != "component-wise affine on six tetrahedra/cube"):
            raise ValueError("--overwrite only replaces this interpolated runner's own result directories.")
        # out is resolved and verified above to be a strict child of HERE.
        # Remove old meshes too, so an absent winner cannot leave a stale block.
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    vtk.vtkMultiThreader.SetGlobalMaximumNumberOfThreads(1)
    vtk.vtkSMPTools.Initialize(1)
    started = time.perf_counter()
    sources = [Path(__file__), HERE/"interpolated_geometry.py", HERE/"gstools_fields.py", HERE/"model.py", ROOT/"src/adler.py"]
    metadata = {"config": config, "status": "running", "completed_realizations": 0,
                "backend": "GSTools RandMeth", "geometry": "component-wise affine on six tetrahedra/cube",
                "n_cells_axis": config["n"], "n_nodes_axis": config["n"]+1,
                "python": platform.python_version(), "gstools": gstools.__version__,
                "numpy": np.__version__, "scipy": scipy.__version__, "pyvista": pv.__version__,
                "vtk": vtk.vtkVersion.GetVTKVersion(), "threads": 1,
                "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
    def save_metadata():
        (out/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    save_metadata()
    plans, theory = theoretical_plan(config)
    write_csv(out/"theory.csv", theory)
    rows, components, terms, raw_rows, field_rows, class_rows = [], [], [], [], [], []
    def checkpoint():
        for name, table in (("measurements", rows), ("union_components", components),
                            ("intersection_terms", terms), ("raw_regions", raw_rows),
                            ("field_diagnostics", field_rows), ("class_volumes", class_rows)):
            write_csv(out/(name+".csv"), table)
        save_metadata()
    log = lambda message: print(message, flush=True)
    log(f"Experiment 07 interpolated: N={config['n']} cells, {config['n']+1} nodes, one thread, output={out}")
    try:
        for realization in range(config["realizations"]):
            groups, diagnostics = sample_components(config["n"], config["length"], config["classes"],
                                                    config["seed"], realization, config["mode_no"], progress=log)
            field_rows.extend(diagnostics)
            base = tetrahedral_grid(config["n"], config["length"])
            names = attach_components(base, groups)
            center_groups = center_samples(groups)
            if config["save_fields"]:
                payload = {name: field for group_names, fields in zip(names, groups) for name, field in zip(group_names, fields)}
                np.savez_compressed(out/f"fields_r{realization+1}.npz", length=config["length"], **payload)
            for phi in config["target_phis"]:
                common = {"realization": realization+1, "target_phi": phi}
                plan = plans[phi]
                mesh_dir = None
                if config["save_meshes"]:
                    mesh_dir = out/f"meshes_r{realization+1}_vf_{phi:.4f}"
                    mesh_dir.mkdir(exist_ok=True)
                log(f"  realization {realization+1}, VF {phi:.0%}: clipping components")
                stats, comps, intersections, raw, added = measure_union(
                    base, names, plan["deltas"], config["length"], tolerance=config["merge_tolerance"],
                    export_directory=mesh_dir, progress=lambda m: log("    "+m))
                row = {**common, **stats, **voxel_comparison(center_groups, plan["deltas"]),
                       "theory_L0": float(plan["lkcs"][0]), "theory_L3": float(plan["lkcs"][3])}
                rows.append(row)
                components.extend({**common, **r} for r in comps)
                terms.extend({**common, **r} for r in intersections)
                raw_rows.extend({**common, **r} for r in raw)
                class_rows.extend({**common, "class_index": k, "class": item["name"],
                                   "target_added_phi": phi*item["share"],
                                   "measured_added_phi": float(added[k]/config["length"]**3)}
                                  for k, item in enumerate(config["classes"]))
                if mesh_dir is not None:
                    # A lightweight multiblock manifest references the measured raw
                    # pieces. Cross-class overlap is intentional and documented.
                    files = sorted(mesh_dir.glob("class_*_winner_*.vtu"))
                    xml = ['<?xml version="1.0"?>', '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian">', '<vtkMultiBlockDataSet>']
                    xml += [f'<DataSet index="{i}" name="{p.stem}" file="{mesh_dir.name}/{p.name}"/>' for i, p in enumerate(files)]
                    xml += ['</vtkMultiBlockDataSet>', '</VTKFile>']
                    (out/(mesh_dir.name+".vtm")).write_text("\n".join(xml)+"\n", encoding="utf-8")
                log(f"    PL VF={stats['phi']:.3%}, L0={stats['L0']}, components={stats['components']}, "
                    f"largest={stats['largest_fraction']:.2%}; voxel L0={row['voxel_euler_6']}/{row['voxel_euler_26']}")
                checkpoint()
            metadata["completed_realizations"] = realization+1
            checkpoint()
        plot_and_summarize(config, rows, out)
        metadata["status"] = "complete"
        metadata["elapsed_seconds"] = time.perf_counter()-started
        save_metadata()
        log(f"Done in {metadata['elapsed_seconds']:.1f}s. See {out/'summary.md'}")
    except Exception as error:
        metadata["status"] = "failed"
        metadata["error"] = f"{type(error).__name__}: {error}"
        checkpoint()
        raise


if __name__ == "__main__":
    main()
