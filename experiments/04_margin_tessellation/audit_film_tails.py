"""Bounded normal-profile audit: linear, curved pair, and full competition.

Reuses the refined film pilot's interface points and weights. All unresolved
or >10% error profiles, plus random controls, receive longer/nested scans.
No field generation, 3D grid, mesh, or particle-size analysis is performed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refine_film_widths as previous
base = previous.base


class BudgetReached(RuntimeError):
    pass


class Profile:
    def __init__(self, fields, crossing, t, values, max_points):
        self.fields, self.crossing = fields, crossing
        self.centre = np.array([crossing[k] for k in ("x", "y", "z")])
        self.normal = np.array([crossing[k] for k in ("nx", "ny", "nz")])
        self.cache = {float(s): v.copy() for s, v in zip(t, values.T)}
        self.max_points = max_points
        self.initial_points = len(self.cache)

    def evaluate(self, t):
        t = np.atleast_1d(t).astype(float)
        missing = list(dict.fromkeys(float(s) for s in t if float(s) not in self.cache))
        if len(self.cache)+len(missing) > self.max_points:
            raise BudgetReached("Profile evaluation cap reached")
        if missing:
            values = self.fields.evaluate(self.centre+np.array(missing)[:, None]*self.normal)
            self.cache.update({s: v.copy() for s, v in zip(missing, values.T)})
        return np.stack([self.cache[float(s)] for s in t], axis=1)

    def margin(self, t, pair=False):
        values = self.evaluate(t)
        if pair:
            return abs(values[self.crossing["winner_i"]]-values[self.crossing["winner_j"]])
        return base.gap_and_winner(values)[0]


def boundaries(profile, grid, delta, pair):
    margin = profile.margin(grid, pair)
    mid = len(grid)//2
    if margin[mid] >= delta:
        raise ValueError("The central point must be inside the matrix")
    ends = []
    for side in (-1, 1):
        candidates = np.arange(mid+side, len(grid) if side > 0 else -1, side)
        hits = candidates[margin[candidates] >= delta]
        endpoint = None
        if len(hits):
            k = hits[0]
            lo, hi = sorted((grid[k-side], grid[k]))
            endpoint = float(brentq(lambda s: float(profile.margin(s, pair)[0]-delta), lo, hi, xtol=2e-11))
        ends.append(endpoint)
    return ends


def same_endpoints(a, b, tolerance):
    return all((x is None and y is None) or
               (x is not None and y is not None and abs(x-y) <= tolerance)
               for x, y in zip(a, b))


def measure_profile(profile, source_rows, lc, max_extent_lc=2., max_refinements=5):
    """Expand the window without coarsening, then halve spacing until stable.

    Two consecutive refinement comparisons are required. The endpoint search
    remains numerical, not certified root isolation. Unclosed sides have
    conditional lower bounds at the searched endpoints, and infinite upper
    bounds. Unstable widths receive [0, infinity] in population summaries.
    """
    original_extent = max(r["extent"] for r in source_rows)
    cap = max_extent_lc*lc
    if cap < original_extent:
        raise ValueError("Extent cap must cover the original profile")
    step, half = original_extent/128, 128
    maximum_half = int(np.floor(cap/step+1e-10))
    previous_result = None
    streak = {(j, kind): 0 for j in range(len(source_rows)) for kind in ("full", "pair")}
    log, last, refinements = [], None, 0
    stop = "refinement_cap"
    while True:
        grid = np.arange(-half, half+1)*step
        try:
            profile.evaluate(grid)
            result = {}
            for j, row in enumerate(source_rows):
                for kind in ("full", "pair"):
                    result[j, kind] = boundaries(profile, grid, row["delta"], pair=kind == "pair")
        except BudgetReached:
            stop = "point_cap"
            break
        for key, endpoints in result.items():
            # An extent change alone is not a resolution check.
            comparable = previous_result is not None and last["half_extent"] == half*step and last["step"] > step
            match = comparable and same_endpoints(endpoints, previous_result[key], 1e-7*lc)
            streak[key] = streak[key]+1 if match else 0
        last = dict(result=result, half_extent=half*step, step=step, streak=dict(streak))
        censored = sum(any(x is None for x in ends) for ends in result.values())
        log.append(dict(iteration=len(log), grid_points=len(grid), unique_evaluations=len(profile.cache),
                        half_extent=half*step, step=step, unclosed_objects=censored,
                        stable_objects=sum(s >= 2 for s in streak.values())))
        if censored and half < maximum_half:
            previous_result = result
            half = min(2*half, maximum_half)
            continue
        if all(s >= 2 for s in streak.values()):
            stop = "stable" if not censored else "extent_cap"
            break
        if refinements >= max_refinements:
            break
        previous_result = result
        half *= 2
        maximum_half *= 2
        step /= 2
        refinements += 1
    if last is None:
        raise RuntimeError("Point cap too small for the first complete measurement")
    rows = []
    for j, old in enumerate(source_rows):
        full, pair = last["result"][j, "full"], last["result"][j, "pair"]
        full_closed, pair_closed = all(x is not None for x in full), all(x is not None for x in pair)
        full_stable = last["streak"][j, "full"] >= 2
        pair_stable = last["streak"][j, "pair"] >= 2
        full_width = full[1]-full[0] if full_closed else None
        pair_width = pair[1]-pair[0] if pair_closed else None
        lower = sum(abs(x) if x is not None else last["half_extent"] for x in full) if full_stable else 0.
        upper = full_width if full_stable and full_closed else float("inf")
        if full_stable and full_closed:
            lower = full_width
        curvature = competition = None
        if full_stable and full_closed and pair_stable and pair_closed:
            curvature = (pair_width-old["linearized_width"])/old["linearized_width"]
            competition = (full_width-pair_width)/old["linearized_width"]
        winners, runners = [], []
        for x in full:
            if x is None:
                winners.append(None)
                runners.append(None)
            else:
                top = np.argsort(profile.evaluate(x)[:, 0])[-2:]
                runners.append(int(top[0]))
                winners.append(int(top[1]))
        original_pair = {profile.crossing["winner_i"], profile.crossing["winner_j"]}
        row = dict(old)
        row.update(measured_width=full_width, stable=full_stable and full_closed, censored=not full_closed,
                   relative_error=None if full_width is None else (full_width-old["linearized_width"])/old["linearized_width"],
                   left=full[0], right=full[1], left_winner=winners[0], right_winner=winners[1],
                   left_runner=runners[0], right_runner=runners[1],
                   pair_preserved=full_closed and set(winners) == original_pair,
                   pair_controls_boundaries=full_closed and all({w, r} == original_pair for w, r in zip(winners, runners)),
                   extent=last["half_extent"], profile_points=len(profile.cache),
                   audit_selected=True, original_width=old["measured_width"], original_stable=old["stable"],
                   original_censored=old["censored"], width_lower=lower, width_upper=upper,
                   pair_width=pair_width, pair_stable=pair_stable and pair_closed,
                   pair_left=pair[0], pair_right=pair[1], curvature_contribution=curvature,
                   competition_contribution=competition, audit_stop=stop,
                   third_clearance_over_delta=profile.crossing["third_clearance"]/old["delta"])
        rows.append(row)
    return rows, log


def load_csv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    integers = {"q", "realization", "crossing_id", "axis", "line", "winner_i", "winner_j",
                "profile_points", "left_winner", "right_winner", "left_runner", "right_runner"}
    for row in rows:
        for key, value in row.items():
            if value == "":
                row[key] = None
            elif value in ("True", "False"):
                row[key] = value == "True"
            else:
                row[key] = int(value) if key in integers else float(value)
    return rows


def select_profiles(rows, threshold=.1, controls=8, seed=47211):
    groups = {}
    for row in rows:
        groups.setdefault(row["crossing_id"], []).append(row)
    reasons = {}
    for key, group in groups.items():
        flags = []
        if any(r["censored"] for r in group): flags.append("censored")
        if any(not r["censored"] and not r["stable"] for r in group): flags.append("unstable")
        if any(r["relative_error"] is not None and abs(r["relative_error"]) > threshold for r in group): flags.append("large_error")
        if flags: reasons[key] = ";".join(flags)
    rng = np.random.default_rng(seed)
    for q, realization in sorted({(r["q"], r["realization"]) for r in rows}):
        candidates = sorted(k for k, group in groups.items()
                            if k not in reasons and (group[0]["q"], group[0]["realization"]) == (q, realization))
        for key in rng.choice(candidates, min(controls, len(candidates)), replace=False):
            reasons[int(key)] = "control"
    return groups, reasons


def inherited(row):
    known = row["stable"] and not row["censored"]
    out = dict(row)
    out.update(audit_selected=False, audit_reason="not_selected", original_width=row["measured_width"],
               original_stable=row["stable"], original_censored=row["censored"],
               width_lower=row["measured_width"] if known else 0.,
               width_upper=row["measured_width"] if known else float("inf"),
               pair_width=None, pair_stable=False, pair_left=None, pair_right=None,
               curvature_contribution=None, competition_contribution=None,
               audit_stop="inherited", third_clearance_over_delta=None)
    return out


def summarize(rows):
    summaries, cdf_rows = [], []
    for q, vf, factor in sorted({(r["q"], r["target_vf"], r["margin_factor"]) for r in rows}):
        group = [r for r in rows if (r["q"], r["target_vf"], r["margin_factor"]) == (q, vf, factor)]
        known = [r for r in group if r["stable"] and not r["censored"]]
        w = np.array([r["area_weight"] for r in group])
        lower = np.array([r["width_lower"] for r in group])
        upper = np.array([r["width_upper"] for r in group])
        errors = base.law.weighted_quantile([abs(r["relative_error"]) for r in known],
                                            [r["area_weight"] for r in known], [.5, .9]) if known else [None, None]
        decomposed = [r for r in group if r["competition_contribution"] is not None and abs(r["relative_error"]) > .1]
        summaries.append(dict(q=q, target_vf=vf, margin_factor=factor, profiles=len(group),
                              rechecked=sum(r["audit_selected"] for r in group), stable_uncensored=len(known),
                              censored=sum(r["censored"] for r in group),
                              unstable_uncensored=sum(not r["stable"] and not r["censored"] for r in group),
                              mean_width_lower=float(np.average(lower, weights=w)),
                              mean_width_upper=float(np.average(upper, weights=w)),
                              median_absolute_relative_error=errors[0], p90_absolute_relative_error=errors[1],
                              large_errors_decomposed=len(decomposed),
                              competition_dominant=sum(abs(r["competition_contribution"]) > abs(r["curvature_contribution"]) for r in decomposed)))
        for cut in (.1, .25, .5, 1., 2.):
            cdf_rows.append(dict(q=q, target_vf=vf, margin_factor=factor, width_threshold=cut,
                                 area_cdf_lower=float(np.sum(w[upper <= cut])/sum(w)),
                                 area_cdf_upper=float(np.sum(w[lower <= cut])/sum(w))))
    return summaries, cdf_rows


def plot_results(output, rows):
    fig, axes = base.plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    for q in sorted({r["q"] for r in rows}):
        data = [r for r in rows if r["q"] == q and r["margin_factor"] == 1 and r["competition_contribution"] is not None]
        axes[0].scatter([abs(r["curvature_contribution"]) for r in data],
                        [abs(r["competition_contribution"]) for r in data], s=14, alpha=.6, label=f"q={q}")
        axes[1].scatter([r["linearized_width"] for r in data], [r["measured_width"] for r in data],
                        s=14, alpha=.6, label=f"q={q}, audited")
    axes[0].set(xscale="symlog", yscale="symlog", xlabel="Absolute curvature contribution / linear width",
                ylabel="Absolute competition contribution / linear width", title="Selected profiles: baseline margins")
    limit = max(axes[1].get_xlim()[1], axes[1].get_ylim()[1])
    axes[1].plot([0, limit], [0, limit], "k--")
    axes[1].set(xlabel="Linearized width", ylabel="Full-field normal width", title="Longer windows and finer probes")
    for ax in axes:
        ax.grid(alpha=.2)
        ax.legend()
    fig.savefig(output/"tail_audit.png", dpi=160)
    base.plt.close(fig)


def run(source, modes_source, output, max_points=8193, max_extent_lc=2., controls=8):
    source, modes_source, output = map(Path, (source, modes_source, output))
    if output.exists() and any(output.iterdir()):
        raise ValueError("Choose a new empty output directory")
    parent = json.loads((source/"metadata.json").read_text())
    mode_metadata = json.loads((modes_source/"metadata.json").read_text())
    if parent["status"] != "complete" or mode_metadata["status"] != "complete":
        raise ValueError("Both source studies must be complete")
    input_files = [source/"metadata.json", modes_source/"metadata.json"]
    for name, digest in parent["source_sha256"].items():
        if hashlib.sha256((base.ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Changed parent-study code: {name}")
    for folder, report, names in ((source, parent, ["profiles.csv", "crossings.csv"]+
                                 [f"profiles_q{q}_r{r}.npz" for q in parent["config"]["qs"] for r in range(parent["config"]["realizations"])]),
                                (modes_source, mode_metadata, [f"spectral_modes_{r}.npz" for r in range(parent["config"]["realizations"])])):
        for name in names:
            path = folder/name
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != report["output_sha256"][name]:
                raise ValueError(f"Changed saved evidence: {path}")
            if folder == modes_source and digest != parent["input_sha256"][name]:
                raise ValueError("Mode source differs from the source used by the refined pilot")
            input_files.append(path)
    old_rows = load_csv(source/"profiles.csv")
    crossing = {r["crossing_id"]: r for r in load_csv(source/"crossings.csv")}
    groups, reasons = select_profiles(old_rows, controls=controls)
    output.mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__), HERE/"test_film_tail_audit.py", HERE/"refine_film_widths.py", HERE/"run_film_widths.py",
               base.ROOT/"theory/margin_film_widths.py", base.ROOT/"theory/check_margin_tubes.py",
               base.ROOT/"src/adler.py", base.ROOT/"src/fields.py"]
    metadata = dict(schema="normal_film_tail_audit_v1", status="running", executed_utc=datetime.now(timezone.utc).isoformat(),
                    source=str(source.resolve()), modes_source=str(modes_source.resolve()),
                    config=dict(max_points=max_points, max_extent_lc=max_extent_lc, max_refinements=5,
                                controls_per_q_realization=controls, large_error_threshold=.1, selection_seed=47211),
                    versions=dict(python=base.platform.python_version(), numpy=np.__version__, scipy=base.scipy.__version__, gstools=base.gs.__version__),
                    selected_profiles=len(reasons), original_profiles=len(groups),
                    source_sha256={p.relative_to(base.ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
                    input_sha256={p.relative_to(base.ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in input_files})
    save_metadata = lambda: (output/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    save_metadata()
    selection = [dict(crossing_id=k, q=groups[k][0]["q"], realization=groups[k][0]["realization"], reason=reasons[k]) for k in sorted(reasons)]
    base.write_csv(output/"selection.csv", selection)
    updated, logs = {}, []
    start = time.perf_counter()
    try:
        for realization in range(parent["config"]["realizations"]):
            fields = previous.StoredFields(modes_source/f"spectral_modes_{realization}.npz")
            for q in parent["config"]["qs"]:
                fields.fields = [None]*q
                payload = {}
                selected = [k for k in sorted(reasons) if (groups[k][0]["q"], groups[k][0]["realization"]) == (q, realization)]
                with np.load(source/f"profiles_q{q}_r{realization}.npz") as raw:
                    for count, key in enumerate(selected):
                        t, values = raw[f"t_{key}"], raw[f"values_{key}"]
                        # Check this stream against stored profile nodes before extending.
                        c = crossing[key]
                        centre = np.array([c[k] for k in ("x", "y", "z")])
                        normal = np.array([c[k] for k in ("nx", "ny", "nz")])
                        if np.max(abs(fields.evaluate(centre+t[[0, len(t)//2, -1], None]*normal)-values[:, [0, len(t)//2, -1]])) > 1e-10:
                            raise ValueError("Stored spectrum does not reproduce the source profile")
                        profile = Profile(fields, c, t, values, max_points)
                        rows, iterations = measure_profile(profile, groups[key], parent["config"]["lc"], max_extent_lc)
                        for row in rows:
                            row["audit_reason"] = reasons[key]
                            updated[key, row["target_vf"], row["margin_factor"]] = row
                        logs.extend(dict(crossing_id=key, q=q, realization=realization, **entry) for entry in iterations)
                        positions = np.array(sorted(profile.cache))
                        payload[f"t_{key}"] = positions
                        payload[f"values_{key}"] = np.stack([profile.cache[x] for x in positions], axis=1)
                        if (count+1) % 16 == 0:
                            print(f"q={q}, realization={realization}: {count+1}/{len(selected)} audited", flush=True)
                if payload:
                    np.savez_compressed(output/f"profiles_q{q}_r{realization}.npz", **payload)
                print(f"q={q}, realization={realization}: {len(selected)} profiles complete", flush=True)
        combined = []
        for old in old_rows:
            row = inherited(old)
            row.update(updated.get((old["crossing_id"], old["target_vf"], old["margin_factor"]), {}))
            combined.append(row)
        base.write_csv(output/"profiles.csv", combined)
        base.write_csv(output/"refinement.csv", logs)
        summary, cdf = summarize(combined)
        base.write_csv(output/"summary.csv", summary)
        base.write_csv(output/"cdf_bounds.csv", cdf)
        plot_results(output, combined)
        metadata.update(status="complete", elapsed_seconds=time.perf_counter()-start,
                        output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir()) if p.suffix in (".csv", ".npz", ".png")})
        save_metadata()
    except BaseException as exc:
        metadata.update(status="failed", error=repr(exc))
        save_metadata()
        raise
    print(f"Audit complete: {len(reasons)} profiles, {metadata['elapsed_seconds']:.1f} seconds", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE/"results_film_widths_refined")
    parser.add_argument("--modes-source", type=Path, default=HERE/"results_film_widths")
    parser.add_argument("--output", type=Path, default=HERE/"results_film_tails")
    parser.add_argument("--max-points", type=int, default=8193)
    parser.add_argument("--max-extent-lc", type=float, default=2.)
    parser.add_argument("--controls", type=int, default=8)
    args = parser.parse_args()
    if args.max_points < 1025 or not np.isfinite(args.max_extent_lc) or args.max_extent_lc < 1 or args.controls < 0:
        parser.error("Require max-points>=1025, max-extent-lc>=1, controls>=0")
    run(args.source, args.modes_source, args.output, args.max_points, args.max_extent_lc, args.controls)
