"""Small spatial check of the area-weighted normal film-width law.

Uniform axis-parallel lines locate winning ties. Pooling equal numbers of
lines in three directions and weighting crossings by 1/sum(abs(normal))
estimates interface-area averages. Direct normal profiles then test the
small-margin approximation. No 3D grid or mesh is allocated.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import gstools as gs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.stats import chi

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from src.adler import delta_for_vf
from theory import margin_film_widths as law


class Fields:
    """GSTools RandMeth in physical coordinates, with analytic mode gradients.

    The gradient uses private spectral arrays, checked against direct SRF
    evaluation and finite differences. Saved versions and mode arrays make
    this dependency explicit. No empirical realization normalization is used.
    """
    def __init__(self, q, lc, seed, modes=1000):
        gs.config.NUM_THREADS = 1
        rng = np.random.default_rng(seed)
        self.seeds = rng.integers(0, 2**31-1, size=q).tolist()
        model = gs.Gaussian(dim=3, len_scale=lc, rescale=1., var=1.)
        self.fields = [gs.SRF(model, generator="RandMeth", seed=s, mode_no=modes)
                       for s in self.seeds]
        self.k = np.stack([f.generator._cov_sample for f in self.fields])
        self.cos = np.stack([f.generator._z_1 for f in self.fields])/np.sqrt(modes)
        self.sin = np.stack([f.generator._z_2 for f in self.fields])/np.sqrt(modes)

    def evaluate(self, points, labels=None):
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        labels = range(len(self.fields)) if labels is None else labels
        return np.array([self.fields[i](points.T, mesh_type="unstructured", store=False)
                         for i in labels])

    def value_gradient(self, point, label):
        phase = np.asarray(point)@self.k[label]
        cosine, sine = np.cos(phase), np.sin(phase)
        value = self.cos[label]@cosine+self.sin[label]@sine
        gradient = self.k[label]@(-self.cos[label]*sine+self.sin[label]*cosine)
        return float(value), gradient

    def save(self, path):
        np.savez_compressed(path, k=self.k, cosine=self.cos, sine=self.sin,
                            component_seeds=np.array(self.seeds))


def gap_and_winner(values):
    top = np.partition(values, -2, axis=0)[-2:]
    return top[1]-top[0], np.argmax(values, axis=0)


def axis_weight(normal):
    return 1/np.sum(np.abs(normal))


def locate_crossings(fields, origins, axis, t, values):
    """Bracket winner changes; reject roots hidden below a third winner."""
    winner = np.argmax(values, axis=0)
    roots, rejected = [], 0
    for line, cell in np.argwhere(winner[:, :-1] != winner[:, 1:]):
        i, j = int(winner[line, cell]), int(winner[line, cell+1])
        point = origins[line].copy()

        def difference(s):
            point[axis] = s
            return fields.value_gradient(point, i)[0]-fields.value_gradient(point, j)[0]

        root = brentq(difference, t[cell], t[cell+1], xtol=2e-11)
        point[axis] = root
        all_values = fields.evaluate(point)[:, 0]
        others = np.delete(all_values, (i, j))
        clearance = float(all_values[i]-np.max(others)) if len(others) else float("inf")
        if clearance < -1e-8:
            rejected += 1
            continue
        gradient = fields.value_gradient(point, i)[1]-fields.value_gradient(point, j)[1]
        norm = np.linalg.norm(gradient)
        if norm < 1e-10:
            raise RuntimeError("Near-degenerate winning tie; cannot define its normal")
        normal = gradient/norm
        roots.append(dict(axis=axis, line=int(line), coordinate=float(root),
                          winner_i=i, winner_j=j, x=float(point[0]), y=float(point[1]), z=float(point[2]),
                          nx=float(normal[0]), ny=float(normal[1]), nz=float(normal[2]),
                          gradient_norm=float(norm), area_weight=float(axis_weight(normal)),
                          third_clearance=clearance, tie_residual=float(abs(all_values[i]-all_values[j]))))
    return roots, rejected


def match_crossings(coarse, fine, tolerance=1e-7):
    def key(row):
        return row["axis"], row["line"], tuple(sorted((row["winner_i"], row["winner_j"])))
    matched = 0
    for row in coarse:
        if any(key(other) == key(row) and abs(row["coordinate"]-other["coordinate"]) < tolerance
               for other in fine):
            matched += 1
    return dict(coarse_count=len(coarse), fine_count=len(fine), matched_count=matched,
                coarse_unmatched=len(coarse)-matched, fine_unmatched=len(fine)-matched)


def normal_profile(fields, crossing, deltas, lc, points=129):
    """Measure the connected matrix interval containing the winning tie.

    Two nested probe grids locate brackets. Brent solves the original
    component margin inside each bracket. Stability is diagnostic, not
    certified isolation of all roots. Censored profiles remain recorded.
    """
    centre = np.array([crossing[k] for k in ("x", "y", "z")])
    normal = np.array([crossing[k] for k in ("nx", "ny", "nz")])
    extent = min(lc, max(.05*lc, 4*max(deltas)/crossing["gradient_norm"]))
    t = np.linspace(-extent, extent, points)
    values = fields.evaluate(centre+t[:, None]*normal)
    margin, _ = gap_and_winner(values)
    rows = []
    for delta in deltas:
        boundaries = []
        for stride in (2, 1):
            probe_t, probe_margin = t[::stride], margin[::stride]
            mid = len(probe_t)//2
            endpoints, labels, runners = [], [], []
            for side in (-1, 1):
                indexes = range(mid+side, len(probe_t), side) if side > 0 else range(mid-1, -1, -1)
                endpoint, winner, runner = None, None, None
                for index in indexes:
                    if probe_margin[index] >= delta:
                        bracket = sorted((probe_t[index-side], probe_t[index]))

                        def excess(s):
                            m, _ = gap_and_winner(fields.evaluate(centre+s*normal))
                            return float(m[0]-delta)

                        endpoint = float(brentq(excess, *bracket, xtol=2e-11))
                        endpoint_values = fields.evaluate(centre+endpoint*normal)[:, 0]
                        top_two = np.argsort(endpoint_values)[-2:]
                        runner, winner = int(top_two[0]), int(top_two[1])
                        break
                endpoints.append(endpoint)
                labels.append(winner)
                runners.append(runner)
            boundaries.append((endpoints, labels, runners))
        coarse, fine = boundaries
        censored = any(x is None for x in fine[0])
        stable = not censored and all(x is not None for x in coarse[0]) and max(
            abs(a-b) for a, b in zip(coarse[0], fine[0])) < 1e-7*lc
        width = None if censored else fine[0][1]-fine[0][0]
        pair_preserved = not censored and set(fine[1]) == {crossing["winner_i"], crossing["winner_j"]}
        pair_controls = pair_preserved and all(
            {winner, runner} == {crossing["winner_i"], crossing["winner_j"]}
            for winner, runner in zip(fine[1], fine[2]))
        predicted = 2*delta/crossing["gradient_norm"]
        rows.append(dict(delta=delta, extent=extent, profile_points=points,
                         censored=censored, stable=stable, pair_preserved=pair_preserved,
                         pair_controls_boundaries=pair_controls,
                         left=fine[0][0], right=fine[0][1], left_winner=fine[1][0], right_winner=fine[1][1],
                         left_runner=fine[2][0], right_runner=fine[2][1],
                         measured_width=width, linearized_width=predicted,
                         relative_error=None if censored else (width-predicted)/predicted))
    return rows, dict(t=t, values=values)


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"No rows to save: {path}")
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(crossings, profiles, config):
    groups, width_groups = [], []
    for q, realization in sorted({(r["q"], r["realization"]) for r in crossings}):
        rows = [r for r in crossings if (r["q"], r["realization"]) == (q, realization)]
        g = np.array([r["gradient_norm"]*config["lc"]/2 for r in rows])
        weights = np.array([r["area_weight"] for r in rows])
        moments = law.weighted_quantile(1/g, weights, [.1, .5, .9])
        density = float(weights.sum()/(config["lines_per_axis"]*config["length"]))
        groups.append(dict(q=q, realization=realization, crossings=len(rows),
                           estimated_interface_area_per_volume=density,
                           theory_interface_area_per_volume=law.interface_area_density(q, config["lc"]),
                           weighted_gradient_cdf_distance=law.weighted_cdf_distance(g, weights, lambda x: chi.cdf(x, df=4)),
                           mean_normalized_linear_width=float(np.average(1/g, weights=weights)),
                           theory_mean_normalized_width=law.width_moment(1, 1., 1.),
                           normalized_width_p10=float(moments[0]), normalized_width_median=float(moments[1]),
                           normalized_width_p90=float(moments[2])))
    for q, vf, factor in sorted({(r["q"], r["target_vf"], r["margin_factor"]) for r in profiles}):
        rows = [r for r in profiles if (r["q"], r["target_vf"], r["margin_factor"]) == (q, vf, factor)]
        valid = [r for r in rows if not r["censored"] and r["stable"]]
        paired = [r for r in valid if r["pair_controls_boundaries"]]
        weights = np.array([r["area_weight"] for r in valid])
        errors = np.array([abs(r["relative_error"]) for r in valid])
        med = law.weighted_quantile(errors, weights, [.5, .9]) if len(valid) else [None, None]
        paired_median = (float(law.weighted_quantile([abs(r["relative_error"]) for r in paired],
                                                    [r["area_weight"] for r in paired], [.5])[0]) if paired else None)
        cdf = lambda x: law.width_cdf(x, rows[0]["delta"], config["lc"])
        width_groups.append(dict(q=q, target_vf=vf, margin_factor=factor, delta=rows[0]["delta"],
                                 profiles=len(rows), stable_uncensored=len(valid),
                                 censored=sum(r["censored"] for r in rows),
                                 pair_preserved=sum(r["pair_preserved"] for r in valid),
                                 pair_controls_boundaries=len(paired),
                                 median_absolute_relative_error=None if med[0] is None else float(med[0]),
                                 p90_absolute_relative_error=None if med[1] is None else float(med[1]),
                                 pair_controlled_median_error=paired_median,
                                 measured_cdf_distance=law.weighted_cdf_distance(
                                     [r["measured_width"] for r in valid], weights, cdf) if valid else None,
                                 paired_linear_cdf_distance=law.weighted_cdf_distance(
                                     [r["linearized_width"] for r in valid], weights, cdf) if valid else None,
                                 measured_mean_width=float(np.average([r["measured_width"] for r in valid], weights=weights)) if valid else None,
                                 theory_mean_width=law.width_moment(1, rows[0]["delta"], config["lc"])))
    return groups, width_groups


def plot_results(output, crossings, profiles, config):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), constrained_layout=True)
    x = np.linspace(.05, 2, 300)
    axes[0].plot(x, law.width_cdf(x, 1., 1.), "k--", label="Area law: 1/chi(4)")
    for q in config["qs"]:
        rows = [r for r in crossings if r["q"] == q]
        values = np.array([2/(r["gradient_norm"]*config["lc"]) for r in rows])
        weights = np.array([r["area_weight"] for r in rows])
        order = np.argsort(values)
        axes[0].step(values[order], np.cumsum(weights[order])/weights.sum(), where="post", label=f"q={q}, spatial ties")
        selected = [r for r in profiles if r["q"] == q and r["margin_factor"] == 1 and r["stable"] and not r["censored"]]
        axes[1].scatter([r["linearized_width"] for r in selected], [r["measured_width"] for r in selected],
                        s=10, alpha=.4, label=f"q={q}")
    limit = max(axes[1].get_xlim()[1], axes[1].get_ylim()[1])
    axes[1].plot([0, limit], [0, limit], "k--")
    for q in config["qs"]:
        for vf in config["vfs"]:
            delta = delta_for_vf(q, vf)
            widths = np.linspace(.01, config["lc"]*.3, 300)
            axes[2].plot(widths, law.width_cdf(widths, delta, config["lc"]), label=f"q={q}, VF={vf:.0%}")
    axes[0].set(xlim=(0, 2), ylim=(0, 1), xlabel="Linearized width / (delta Lc)", ylabel="Area-weighted CDF")
    axes[1].set(xlabel="Linearized normal width", ylabel="Measured normal width", title="Finite margins, VF 75% / 85%")
    axes[2].set(xlabel="Width threshold (physical units)", ylabel="Predicted interface-area fraction below threshold",
                title="Small-margin approximation")
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend(fontsize=7)
    fig.savefig(output/"film_widths.png", dpi=160)
    plt.close(fig)


def run(config, output):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Choose a new, empty output directory")
    output.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    sources = [Path(__file__), HERE/"test_film_widths.py", ROOT/"theory/margin_film_widths.py",
               ROOT/"theory/check_margin_tubes.py", ROOT/"src/adler.py", ROOT/"src/fields.py"]
    metadata = dict(schema="normal_film_widths_v1", status="running", config=config,
                    executed_utc=datetime.now(timezone.utc).isoformat(),
                    versions=dict(python=platform.python_version(), numpy=np.__version__,
                                  scipy=scipy.__version__, gstools=gs.__version__),
                    source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    write_metadata = lambda: (output/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    write_metadata()
    crossings, convergence, profiles, diagnostics = [], [], [], []
    try:
        for realization in range(config["realizations"]):
            fields = Fields(max(config["qs"]), config["lc"], config["seed"]+realization, config["modes"])
            probe = np.array([1.23, 4.56, 7.89])
            observed = fields.evaluate(probe)[:, 0]
            reconstructed = np.array([fields.value_gradient(probe, i)[0] for i in range(len(fields.fields))])
            if np.max(abs(observed-reconstructed)) > 1e-10:
                raise RuntimeError("Private spectral arrays do not reproduce direct GSTools values")
            fields.save(output/f"spectral_modes_{realization}.npz")
            rng = np.random.default_rng(config["seed"]+10000+realization)
            lines = [rng.uniform(0, config["length"], size=(config["lines_per_axis"], 3)) for _ in range(3)]
            cells = int(np.ceil(config["length"]/(config["lc"]/config["subdivisions"])))
            cells += cells % 2
            t = np.linspace(0, config["length"], cells+1)
            pools = {q: [] for q in config["qs"]}
            original_fields = fields.fields
            for axis in range(3):
                points = np.repeat(lines[axis][:, None, :], len(t), axis=1)
                points[:, :, axis] = t
                fields.fields = original_fields
                values = fields.evaluate(points.reshape(-1, 3)).reshape(len(original_fields), len(lines[axis]), len(t))
                for q in config["qs"]:
                    fields.fields = original_fields[:q]
                    coarse, rejected_coarse = locate_crossings(fields, lines[axis], axis, t[::2], values[:q, :, ::2])
                    fine, rejected_fine = locate_crossings(fields, lines[axis], axis, t, values[:q])
                    convergence.append(dict(q=q, realization=realization, axis=axis,
                                            coarse_spacing=float(t[2]-t[0]), fine_spacing=float(t[1]-t[0]),
                                            rejected_coarse=rejected_coarse, rejected_fine=rejected_fine,
                                            **match_crossings(coarse, fine)))
                    pools[q].extend(fine)
            for q in config["qs"]:
                fields.fields = original_fields[:q]
                pool = pools[q]
                if not pool:
                    raise RuntimeError("No winning ties found; increase the line sample")
                for row in pool:
                    row.update(q=q, realization=realization, crossing_id=len(crossings))
                    crossings.append(row)
                chosen = rng.choice(len(pool), min(config["profiles_per_group"], len(pool)), replace=False)
                profiles_payload = {}
                specifications = [(vf, factor, delta_for_vf(q, vf)*factor)
                                  for vf in config["vfs"] for factor in config["factors"]]
                print(f"q={q}, realization={realization}: {len(pool)} ties, {len(chosen)} normal profiles", flush=True)
                for count, index in enumerate(chosen):
                    crossing = pool[index]
                    rows, payload = normal_profile(fields, crossing, [d for _, _, d in specifications],
                                                   config["lc"], config["profile_points"])
                    for row, (vf, factor, _) in zip(rows, specifications):
                        profiles.append(dict(q=q, realization=realization, crossing_id=crossing["crossing_id"],
                                             target_vf=vf, margin_factor=factor,
                                             crossing_weight=crossing["area_weight"],
                                             subsample_inverse_probability=len(pool)/len(chosen),
                                             area_weight=crossing["area_weight"]*len(pool)/len(chosen), **row))
                    profiles_payload[f"t_{crossing['crossing_id']}"] = payload["t"]
                    profiles_payload[f"values_{crossing['crossing_id']}"] = payload["values"]
                    if (count+1) % 16 == 0:
                        print(f"  {count+1}/{len(chosen)} profiles", flush=True)
                np.savez_compressed(output/f"profiles_q{q}_r{realization}.npz", **profiles_payload)
                # Deterministic source/spectrum diagnostics, not per-sample normalization.
                point = np.array([1.23, 4.56, 7.89])
                direct = fields.evaluate(point)[:, 0]
                analytic = np.array([fields.value_gradient(point, i)[0] for i in range(q)])
                moment = np.mean(fields.k[:q]**2, axis=2)
                diagnostics.append(dict(q=q, realization=realization,
                                        spectral_value_max_error=float(np.max(abs(direct-analytic))),
                                        ideal_component_gradient_variance=2/config["lc"]**2,
                                        spectrum_gradient_variance_min=float(moment.min()),
                                        spectrum_gradient_variance_max=float(moment.max())))
                write_csv(output/"crossings.csv", crossings)
                write_csv(output/"crossing_refinement.csv", convergence)
                write_csv(output/"profiles.csv", profiles)
                write_csv(output/"field_diagnostics.csv", diagnostics)
            fields.fields = original_fields
        groups, width_groups = summarize(crossings, profiles, config)
        write_csv(output/"gradient_summary.csv", groups)
        write_csv(output/"width_summary.csv", width_groups)
        predictions = []
        for q in config["qs"]:
            for vf in config["vfs"]:
                delta = delta_for_vf(q, vf)
                predictions.append(dict(q=q, target_vf=vf, lc=config["lc"], delta=delta,
                                        mean_width=law.width_moment(1, delta, config["lc"]),
                                        p10=float(law.width_quantile(.1, delta, config["lc"])),
                                        median=float(law.width_quantile(.5, delta, config["lc"])),
                                        p90=float(law.width_quantile(.9, delta, config["lc"])),
                                        fraction_below_dx_0p5=float(law.width_cdf(.5, delta, config["lc"])),
                                        fraction_below_dx_0p25=float(law.width_cdf(.25, delta, config["lc"]))))
        write_csv(output/"predictions.csv", predictions)
        plot_results(output, crossings, profiles, config)
        metadata.update(status="complete", elapsed_seconds=time.perf_counter()-start,
                        output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                       for p in sorted(output.iterdir()) if p.suffix in (".csv", ".npz", ".png")})
        write_metadata()
    except BaseException as exc:
        metadata.update(status="failed", error=repr(exc))
        write_metadata()
        raise
    print(f"Complete in {metadata['elapsed_seconds']:.1f} seconds: {output}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE/"results_film_widths")
    parser.add_argument("--qs", type=int, nargs="+", default=[16, 32])
    parser.add_argument("--realizations", type=int, default=2)
    parser.add_argument("--lines-per-axis", type=int, default=24)
    parser.add_argument("--profiles-per-group", type=int, default=64)
    parser.add_argument("--subdivisions", type=int, default=16, help="Fine line samples per Lc")
    parser.add_argument("--profile-points", type=int, default=129)
    parser.add_argument("--modes", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=481703)
    parser.add_argument("--lc", type=float, default=5.)
    parser.add_argument("--length", type=float, default=64.)
    args = parser.parse_args()
    config = {k: v for k, v in vars(args).items() if k != "output"}
    config.update(qs=sorted(set(args.qs)), vfs=[.75, .85], factors=[1., .5, .25])
    if (min(config["qs"]) < 2 or min(args.realizations, args.lines_per_axis, args.profiles_per_group,
                                     args.subdivisions, args.modes) < 1 or args.seed < 0
            or args.profile_points < 9 or (args.profile_points-1) % 4
            or not np.isfinite(args.lc) or args.lc <= 0
            or not np.isfinite(args.length) or args.length <= 0):
        parser.error("Require q>=2, positive sizes/counts, nonnegative seed, profile-points=4k+1>=9")
    run(config, args.output)


if __name__ == "__main__":
    main()
