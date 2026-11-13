"""Refine the film pilot using saved modes and adaptive hidden-winner splits.

The original pilot is retained. This follow-up uses the same fields and
lines, doubles the line resolution, and repeats normal-profile sampling
from the improved crossing pools. Direct spectral evaluation bypasses SRF
position bookkeeping; tests verify agreement with the original evaluator.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import brentq
from gstools.field.generator import _summate

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_film_widths as base


class StoredFields(base.Fields):
    def __init__(self, path):
        with np.load(path) as data:
            self.k = data["k"]
            self.cos = data["cosine"]
            self.sin = data["sine"]
            self.seeds = data["component_seeds"].tolist()
        self.fields = [None]*len(self.seeds)

    def evaluate(self, points, labels=None):
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        labels = range(len(self.fields)) if labels is None else labels
        return np.array([_summate(self.k[i], self.cos[i], self.sin[i], points.T, num_threads=1)
                         for i in labels])


def adaptive_crossings(fields, origins, axis, t, values):
    winner = np.argmax(values, axis=0)
    roots, splits = [], 0

    def visit(line, low, i, high, j, depth):
        nonlocal splits
        point = origins[line].copy()

        def difference(s):
            point[axis] = s
            return fields.value_gradient(point, i)[0]-fields.value_gradient(point, j)[0]

        root = brentq(difference, low, high, xtol=2e-11)
        point[axis] = root
        all_values = fields.evaluate(point)[:, 0]
        other = np.delete(all_values, (i, j))
        clearance = float(all_values[i]-max(other)) if len(other) else float("inf")
        if clearance < -1e-8:
            if depth >= 30 or min(root-low, high-root) < 1e-10:
                raise RuntimeError("Adaptive hidden-winner bracket failed to separate")
            splits += 1
            k = int(np.argmax(all_values))
            visit(line, low, i, root, k, depth+1)
            visit(line, root, k, high, j, depth+1)
            return
        gradient = fields.value_gradient(point, i)[1]-fields.value_gradient(point, j)[1]
        norm = np.linalg.norm(gradient)
        if norm < 1e-10:
            raise RuntimeError("Degenerate tie gradient")
        normal = gradient/norm
        roots.append(dict(axis=axis, line=int(line), coordinate=float(root), winner_i=int(i), winner_j=int(j),
                          x=float(point[0]), y=float(point[1]), z=float(point[2]),
                          nx=float(normal[0]), ny=float(normal[1]), nz=float(normal[2]),
                          gradient_norm=float(norm), area_weight=float(base.axis_weight(normal)),
                          third_clearance=clearance, tie_residual=float(abs(all_values[i]-all_values[j]))))

    for line, cell in np.argwhere(winner[:, :-1] != winner[:, 1:]):
        visit(line, t[cell], int(winner[line, cell]), t[cell+1], int(winner[line, cell+1]), 0)
    return roots, splits


def run(source, output, subdivisions):
    source, output = Path(source), Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Choose an empty output directory")
    parent = json.loads((source/"metadata.json").read_text())
    if parent["status"] != "complete":
        raise ValueError("Source study must be complete")
    input_names = ["metadata.json"]+[f"spectral_modes_{r}.npz" for r in range(parent["config"]["realizations"])]
    for name in input_names[1:]:
        if hashlib.sha256((source/name).read_bytes()).hexdigest() != parent["output_sha256"][name]:
            raise ValueError(f"Changed input: {name}")
    config = dict(parent["config"], subdivisions=subdivisions)
    if subdivisions <= parent["config"]["subdivisions"]:
        raise ValueError("Follow-up must refine the source line spacing")
    refinement_factor, remainder = divmod(subdivisions, parent["config"]["subdivisions"])
    if remainder or refinement_factor & (refinement_factor-1):
        raise ValueError("Use a power-of-two multiple of the source subdivisions")
    output.mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__), HERE/"test_film_width_refinement.py", HERE/"run_film_widths.py",
               base.ROOT/"theory/margin_film_widths.py", base.ROOT/"theory/check_margin_tubes.py",
               base.ROOT/"src/adler.py", base.ROOT/"src/fields.py"]
    metadata = dict(schema="normal_film_widths_refined_v1", status="running", config=config,
                    executed_utc=datetime.now(timezone.utc).isoformat(), source=str(source),
                    versions=dict(python=base.platform.python_version(), numpy=np.__version__,
                                  scipy=base.scipy.__version__, gstools=base.gs.__version__),
                    input_sha256={name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in input_names},
                    source_sha256={p.relative_to(base.ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    save_metadata = lambda: (output/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    save_metadata()
    start = time.perf_counter()
    crossings, convergence, profiles = [], [], []
    try:
        for realization in range(config["realizations"]):
            fields = StoredFields(source/f"spectral_modes_{realization}.npz")
            maximum_fields = fields.fields
            rng = np.random.default_rng(config["seed"]+10000+realization)
            origins = [rng.uniform(0, config["length"], size=(config["lines_per_axis"], 3)) for _ in range(3)]
            cells = int(np.ceil(config["length"]/(config["lc"]/parent["config"]["subdivisions"])))
            cells += cells % 2
            cells *= refinement_factor
            t = np.linspace(0, config["length"], cells+1)
            pools = {q: [] for q in config["qs"]}
            for axis in range(3):
                points = np.repeat(origins[axis][:, None, :], len(t), axis=1)
                points[:, :, axis] = t
                fields.fields = maximum_fields
                values = fields.evaluate(points.reshape(-1, 3)).reshape(len(maximum_fields), len(origins[axis]), len(t))
                for q in config["qs"]:
                    fields.fields = maximum_fields[:q]
                    coarse, splits_coarse = adaptive_crossings(fields, origins[axis], axis, t[::2], values[:q, :, ::2])
                    fine, splits_fine = adaptive_crossings(fields, origins[axis], axis, t, values[:q])
                    convergence.append(dict(q=q, realization=realization, axis=axis,
                                            coarse_spacing=float(t[2]-t[0]), fine_spacing=float(t[1]-t[0]),
                                            adaptive_splits_coarse=splits_coarse, adaptive_splits_fine=splits_fine,
                                            **base.match_crossings(coarse, fine)))
                    pools[q].extend(fine)
                print(f"Refined realization={realization}, axis={axis}: crossings complete", flush=True)
            for q in config["qs"]:
                fields.fields = maximum_fields[:q]
                pool = pools[q]
                for row in pool:
                    row.update(q=q, realization=realization, crossing_id=len(crossings))
                    crossings.append(row)
                chosen = rng.choice(len(pool), min(config["profiles_per_group"], len(pool)), replace=False)
                specs = [(vf, factor, base.delta_for_vf(q, vf)*factor) for vf in config["vfs"] for factor in config["factors"]]
                payload = {}
                for index in chosen:
                    crossing = pool[index]
                    rows, raw = base.normal_profile(fields, crossing, [d for _, _, d in specs], config["lc"], config["profile_points"])
                    for row, (vf, factor, _) in zip(rows, specs):
                        profiles.append(dict(q=q, realization=realization, crossing_id=crossing["crossing_id"],
                                             target_vf=vf, margin_factor=factor,
                                             crossing_weight=crossing["area_weight"],
                                             subsample_inverse_probability=len(pool)/len(chosen),
                                             area_weight=crossing["area_weight"]*len(pool)/len(chosen), **row))
                    payload[f"t_{crossing['crossing_id']}"] = raw["t"]
                    payload[f"values_{crossing['crossing_id']}"] = raw["values"]
                np.savez_compressed(output/f"profiles_q{q}_r{realization}.npz", **payload)
                print(f"Refined q={q}, realization={realization}: {len(pool)} ties and {len(chosen)} profiles complete", flush=True)
                base.write_csv(output/"crossings.csv", crossings)
                base.write_csv(output/"crossing_refinement.csv", convergence)
                base.write_csv(output/"profiles.csv", profiles)
        groups, widths = base.summarize(crossings, profiles, config)
        base.write_csv(output/"gradient_summary.csv", groups)
        base.write_csv(output/"width_summary.csv", widths)
        base.plot_results(output, crossings, profiles, config)
        metadata.update(status="complete", elapsed_seconds=time.perf_counter()-start,
                        output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                       for p in sorted(output.iterdir()) if p.suffix in (".csv", ".npz", ".png")})
        save_metadata()
    except BaseException as exc:
        metadata.update(status="failed", error=repr(exc))
        save_metadata()
        raise
    print(f"Refinement complete in {metadata['elapsed_seconds']:.1f} seconds", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE/"results_film_widths")
    parser.add_argument("--output", type=Path, default=HERE/"results_film_widths_refined")
    parser.add_argument("--subdivisions", type=int, default=32)
    args = parser.parse_args()
    run(args.source, args.output, args.subdivisions)
