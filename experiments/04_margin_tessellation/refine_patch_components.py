"""Refine interpolation inputs while reusing the completed 3D patch reference.

The direct comparison grid stays fixed. This isolates component interpolation
error and does not improve or certify the reference's digital topology.
"""
from __future__ import annotations

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_local_patches as patches


def component_reconstruction(coarse, evaluation_refinement, input_refinement):
    if coarse.shape[1:] != (2*input_refinement+1,)*3:
        raise ValueError("Component input grid does not match its refinement.")
    points = patches.grid_points(evaluation_refinement)*input_refinement
    shape = (2*evaluation_refinement+1,)*3
    arrays = {}
    for kind in ("trilinear", "tetra"):
        m, w = patches.films.margin_and_winner(patches.interpolate(coarse, points, kind))
        arrays[f"{kind}_components"] = (m.reshape(shape), w.reshape(shape))
    return arrays


def report(output, rows):
    lines = ["# Refining the component input grids", "",
             "The direct 33^3 reference grid is unchanged in the default study.",
             "Numbers below are phase-disagreement percentage points of local patch volume.",
             "Remaining reference sampling and topology uncertainties are unchanged.", "",
             "| Patch | Target VF | Input refinement | Trilinear error | Tetrahedral error |",
             "|---|---|---|---|---|"]
    six = [r for r in rows if int(r["connectivity"]) == 6]
    for key in sorted({(int(r["patch_id"]), float(r["target_phi"]), int(r["input_refinement"])) for r in six}):
        pid, phi, factor = key
        selected = {r["method"]: r for r in six if (int(r["patch_id"]), float(r["target_phi"]),
                                                   int(r["input_refinement"])) == key}
        lines.append(f"| {pid} | {phi:.0%} | {factor} | "
                     f"{100*float(selected['trilinear_components']['phase_disagreement_phi']):.5f} | "
                     f"{100*float(selected['tetra_components']['phase_disagreement_phi']):.5f} |")
    lines += ["", "Input refinement 1 uses the original 3^3 component nodes per patch;",
              "2 and 4 use 5^3 and 9^3 nodes on the same fixed box.",
              "Topology comparisons use Euler, component count and component/face incidence",
              "under both digital adjacencies; agreement is not a continuum certificate.", ""]
    (output/"summary.md").write_text("\n".join(lines), encoding="utf-8")
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), constrained_layout=True)
    for ax, pid in zip(axes.flat, sorted({int(r["patch_id"]) for r in six})):
        for phi, color in ((.75, "tab:blue"), (.85, "tab:orange")):
            for method, style in (("trilinear_components", "o-"), ("tetra_components", "s--")):
                group = sorted([r for r in six if int(r["patch_id"]) == pid
                                and float(r["target_phi"]) == phi and r["method"] == method],
                               key=lambda r: int(r["input_refinement"]))
                if group:
                    ax.plot([int(r["input_refinement"]) for r in group],
                            [100*float(r["phase_disagreement_phi"]) for r in group], style,
                            color=color, label=f"{phi:.0%}, {method.split('_')[0]}")
        ax.set(title=f"Patch {pid}", xlabel="Component input refinement factor", ylabel="Phase disagreement (pp)", xticks=[1, 2, 4])
        ax.grid(alpha=.25)
        ax.legend(fontsize=7)
    fig.suptitle("Refining component values; fixed direct comparison grid")
    fig.savefig(output/"input_refinement.png", dpi=150)
    plt.close(fig)


def run(source, output, factors=(2, 4)):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output is not empty; choose a new --output.")
    metadata = json.loads((source/"metadata.json").read_text())
    if metadata["status"] != "complete":
        raise ValueError("The source patch study must be complete.")
    config = metadata["config"]
    dx = metadata["scene"]["spacing"][0]
    evaluation = max(config["refinements"])
    if any(f <= 1 or f > evaluation or f & (f-1) for f in factors):
        raise ValueError("Input factors must be powers of two, above one and no finer than the reference.")
    selection = json.loads((source/"patches.json").read_text())
    original = patches.read_csv(source/"measurements.csv")
    methods = ("trilinear_components", "tetra_components")
    # Source hashes make a silently changed evaluator or geometry routine an error.
    for name, digest in metadata["source_sha256"].items():
        if patches.saved.sha256(patches.saved.ROOT/name) != digest:
            raise ValueError(f"Source study code has changed: {name}")
    sampler = patches.films.OriginalFields(config["q"], config["lc"], dx, config["field_seed"],
                                         config["mode_no"], config["threads"])
    references = {(int(r["patch_id"]), float(r["target_phi"]), int(r["connectivity"])): r
                  for r in original if r["method"] == "direct" and int(r["refinement"]) == evaluation}

    def extract(row, factor):
        row = dict(row)
        for key in ("phi_change_from_previous", "topology_matches_previous", "stable_comparisons", "two_step_stable"):
            row.pop(key)
        ref = references[int(row["patch_id"]), float(row["target_phi"]), int(row["connectivity"])]
        row.update(input_refinement=factor, input_spacing=dx/factor,
                   reference_two_step_stable=ref["two_step_stable"],
                   reference_topology_matches_previous=ref["topology_matches_previous"])
        return row

    rows = [extract(r, 1) for r in original if int(r["refinement"]) == evaluation and r["method"] in methods]
    output.mkdir(parents=True, exist_ok=True)
    inputs = [source/name for name in ("metadata.json", "measurements.csv", "patches.json")]
    inputs += [source/f"patch_{p['patch_id']:02d}.npz" for p in selection]
    sources = [Path(__file__), HERE/"run_local_patches.py", HERE/"run_local_films.py"]
    info = dict(schema="local_patch_component_refinement_v1", status="running", source=str(source),
                input_factors=[1, *factors], evaluation_refinement=evaluation, config=config,
                independent_realizations=1, completed_patches=0, total_patches=len(selection),
                input_sha256={str(p): patches.saved.sha256(p) for p in inputs},
                source_sha256={p.relative_to(patches.saved.ROOT).as_posix(): patches.saved.sha256(p) for p in sources},
                versions={v: version(v) for v in ("numpy", "scipy", "gstools", "scikit-image")},
                reference="unchanged direct samples from source; no additional direct-reference refinement")
    patches.saved.write_json(output/"metadata.json", info)
    start = time.perf_counter()
    deltas = [patches.saved.adler.delta_for_vf(config["q"], phi) for phi in config["phis"]]
    try:
        for patch in selection:
            pid = patch["patch_id"]
            with np.load(source/f"patch_{pid:02d}.npz") as stored:
                np.testing.assert_array_equal(sampler.evaluate(patches.COARSE_POINTS+patch["lower"]).reshape(config["q"], 3, 3, 3),
                                              stored["coarse_fields"])
                direct = stored["direct_margin"].copy(), stored["direct_winner"].copy()
                baseline = stored["trilinear_margin_margin"].copy(), stored["trilinear_margin_winner"].copy()
            payload = {}
            for factor in factors:
                print(f"Patch {pid}: component values on {(2*factor+1)}^3 nodes; direct reference fixed.", flush=True)
                coarse = sampler.evaluate(patches.grid_points(factor)+patch["lower"]).reshape((config["q"],)+(2*factor+1,)*3)
                arrays = dict(direct=direct, trilinear_margin=baseline,
                              **component_reconstruction(coarse, evaluation, factor))
                measured, _, _ = patches.measure_level(patch, evaluation, arrays, config["phis"], deltas,
                                                       dx, {}, config["volume_tolerance"])
                rows.extend(extract(r, factor) for r in measured if r["method"] in methods)
                payload[f"input_{factor}_fields"] = coarse
                for method in methods:
                    payload[f"input_{factor}_{method}_margin"] = arrays[method][0]
                    payload[f"input_{factor}_{method}_winner"] = arrays[method][1].astype(np.int16)
            np.savez_compressed(output/f"patch_{pid:02d}.npz", **payload)
            patches.saved.write_csv(output/"measurements.csv", rows)
            info["completed_patches"] += 1
            patches.saved.write_json(output/"metadata.json", info)
        report(output, rows)
        info.update(status="complete", elapsed_seconds=time.perf_counter()-start)
        patches.saved.write_json(output/"metadata.json", info)
    except BaseException as exc:
        info.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                    error=f"{type(exc).__name__}: {exc}")
        patches.saved.write_json(output/"metadata.json", info)
        raise
    print(f"Complete: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE/"results_local_patches")
    parser.add_argument("--output", type=Path, default=HERE/"results_local_patch_components_refined")
    parser.add_argument("--input-factors", type=int, nargs="+", default=[2, 4])
    args = parser.parse_args()
    try:
        run(args.source, args.output, sorted(set(args.input_factors)))
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    main()
