"""Local directional film audit of the original saved GSTools realization.

Samples face-neighbour segments at target VF .75 and .85. Compares linear
margin interpolation, exact restrictions of affine Gaussian components,
and nested direct RandMeth samples. This is not a global geometry/LKC study.
"""
from __future__ import annotations

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import sys
import time

import gstools as gs
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_saved_sweep as saved

SCHEMA = "local_margin_films_v1"


def margin_and_winner(values):
    top = np.partition(values, -2, axis=0)[-2:]
    return top[1]-top[0], np.argmax(values, axis=0)


class OriginalFields:
    """Same double seed draw and voxel coordinates as src.fields._generate_one."""
    def __init__(self, q, lc, dx, seed, mode_no=1000, threads=1, batch=1024):
        gs.config.NUM_THREADS = threads
        rng = np.random.default_rng(seed)
        self.seeds = [int(gs.random.MasterRNG(int(rng.integers(0, 2**31-1)))()) for _ in range(q)]
        model = gs.Gaussian(dim=3, len_scale=lc/dx, var=1., rescale=1.)
        self.fields = [gs.SRF(model, generator="RandMeth", seed=s, mode_no=mode_no) for s in self.seeds]
        self.batch = batch

    def evaluate(self, points):
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        values = np.empty((len(self.fields), len(points)))
        for begin in range(0, len(points), self.batch):
            end = min(begin+self.batch, len(points))
            for c, field in enumerate(self.fields):
                values[c, begin:end] = field(points[begin:end].T, mesh_type="unstructured", store=False)
        return values


def select_segments(margin, winners, delta, count, kind, rng):
    """Uniform without replacement over all eligible positive-axis face pairs."""
    phase = margin >= delta
    shapes, counts = [], []

    def eligible(axis):
        low = [slice(None)]*3
        high = [slice(None)]*3
        low[axis], high[axis] = slice(None, -1), slice(1, None)
        low, high = tuple(low), tuple(high)
        different = winners[low] != winners[high]
        return phase[low] & phase[high] & (different if kind == "cross_label" else ~different)

    for axis in range(3):
        mask = eligible(axis)
        shapes.append(mask.shape)
        counts.append(int(mask.sum()))
    total = sum(counts)
    chosen = np.sort(rng.choice(total, size=min(count, total), replace=False))
    rows, offset = [], 0
    for axis, population in enumerate(counts):
        positions = chosen[(chosen >= offset) & (chosen < offset+population)]-offset
        if len(positions):
            flat = np.flatnonzero(eligible(axis))[positions]
            for start in np.array(np.unravel_index(flat, shapes[axis])).T:
                end = start.copy()
                end[axis] += 1
                rows.append(dict(kind=kind, axis=axis, i=int(start[0]), j=int(start[1]), k=int(start[2]),
                                 winner_start=int(winners[tuple(start)]), winner_end=int(winners[tuple(end)]),
                                 margin_start=float(margin[tuple(start)]), margin_end=float(margin[tuple(end)])))
        offset += population
    return rows, dict(kind=kind, pool_total=total, pool_per_axis=counts, selected=len(rows))


def validate_nodes(sampler, nodes, margin, winners, source, tolerance=2e-6):
    """Round regenerated fields to float32 before comparing with saved arrays."""
    mesh = pv.read(source/"field_g1.vti")
    if (tuple(mesh.dimensions) != margin.shape or "G_1" not in mesh.point_data
            or not np.allclose(mesh.direction_matrix, np.eye(3))):
        raise ValueError("Invalid saved first-field geometry or scalar array.")
    first = np.asarray(mesh.point_data["G_1"]).reshape(mesh.dimensions, order="F")
    values = sampler.evaluate(nodes).astype(np.float32)
    reconstructed_margin, reconstructed_winner = margin_and_winner(values)
    index = tuple(nodes.T)
    diagnostics = dict(nodes=len(nodes),
                       first_field_max_abs=float(np.max(np.abs(values[0]-first[index]))),
                       margin_max_abs=float(np.max(np.abs(reconstructed_margin-margin[index]))),
                       winner_mismatches=int(np.count_nonzero(reconstructed_winner != winners[index])),
                       first_field_exact_mismatches=int(np.count_nonzero(values[0] != first[index])),
                       margin_exact_mismatches=int(np.count_nonzero(reconstructed_margin != margin[index])))
    if (diagnostics["first_field_max_abs"] > tolerance or diagnostics["margin_max_abs"] > tolerance
            or diagnostics["winner_mismatches"]):
        raise ValueError(f"Original realization was not reproduced: {diagnostics}")
    return diagnostics


def complement(occupied):
    cursor, gaps = 0., []
    for left, right in sorted(occupied):
        left, right = max(0., float(left)), min(1., float(right))
        if right <= left:
            continue
        if left > cursor:
            gaps.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1.:
        gaps.append((cursor, 1.))
    return gaps


def affine_intervals(endpoints, delta):
    """Exact 1D winner regions for affine components; includes all competitors."""
    regions = []
    for winner in range(len(endpoints)):
        low, high = 0., 1.
        for other in range(len(endpoints)):
            if winner == other:
                continue
            difference = endpoints[winner, 0]-endpoints[other, 0]
            intercept = difference-delta
            slope = (endpoints[winner, 1]-endpoints[other, 1])-difference
            if slope == 0:
                if intercept < 0:
                    high = -1.
                    break
            elif slope > 0:
                low = max(low, -intercept/slope)
            else:
                high = min(high, -intercept/slope)
            if high <= low:
                break
        if high > low:
            regions.append((float(low), float(high), winner))
    regions.sort()
    return regions, complement([(lo, hi) for lo, hi, _ in regions])


def sampled_gaps(t, signed_margin):
    """Negative intervals of the piecewise-linear sampled signed margin."""
    intervals = []
    for a, b, fa, fb in zip(t[:-1], t[1:], signed_margin[:-1], signed_margin[1:]):
        if fa >= 0 and fb >= 0:
            continue
        if fa < 0 and fb < 0:
            left, right = a, b
        else:
            root = a-fa*(b-a)/(fb-fa)
            left, right = (a, root) if fa < 0 else (root, b)
        if right <= left:
            continue
        if intervals and abs(intervals[-1][1]-left) < 1e-13:
            intervals[-1] = (intervals[-1][0], float(right))
        else:
            intervals.append((float(left), float(right)))
    return intervals


def total_width(gaps):
    return float(sum(right-left for left, right in gaps))


def boundary_shift(a, b):
    if len(a) != len(b):
        return None
    return float(np.max(np.abs(np.array(a)-np.array(b)))) if a else 0.


def symmetric_difference(a, b):
    overlap = sum(max(0., min(hi, hj)-max(lo, lj)) for lo, hi in a for lj, hj in b)
    return max(0., total_width(a)+total_width(b)-2*overlap)


def unresolved_switches(margin, winners, delta):
    return int(np.count_nonzero((winners[:-1] != winners[1:]) &
                               (margin[:-1] >= delta) & (margin[1:] >= delta)))


def refine_segments(sampler, segments, endpoints, levels, tolerance, output=None):
    """Nested probes with two stable successive comparisons; finite diagnostic."""
    states, convergence = [], []
    for row, fields in zip(segments, endpoints):
        m, w = margin_and_winner(fields)
        states.append(dict(t=np.array([0., 1.]), margin=m, winners=w, gaps=None,
                           streak=0, converged=False, last_shift=None))
    for size in levels:
        active = [i for i, state in enumerate(states) if not state["converged"]]
        if not active:
            break
        t = np.linspace(0., 1., size)
        positions, lengths = [], []
        for index in active:
            row, state = segments[index], states[index]
            new_t = t[1:-1] if len(state["t"]) == 2 else t[1::2]
            start = np.array([row["i"], row["j"], row["k"]], dtype=float)
            points = np.broadcast_to(start, (len(new_t), 3)).copy()
            points[:, row["axis"]] += new_t
            positions.append(points)
            lengths.append(len(new_t))
        print(f"Direct evaluation: {len(active)} segments, nested grid {size} points...", flush=True)
        values = sampler.evaluate(np.concatenate(positions))
        new_margins, new_winners = margin_and_winner(values)
        offset = 0
        for index, length in zip(active, lengths):
            state, row = states[index], segments[index]
            m, w = np.empty(size), np.empty(size, dtype=int)
            if len(state["t"]) == 2:
                m[[0, -1]], w[[0, -1]] = state["margin"], state["winners"]
                m[1:-1], w[1:-1] = new_margins[offset:offset+length], new_winners[offset:offset+length]
            else:
                m[::2], w[::2] = state["margin"], state["winners"]
                m[1::2], w[1::2] = new_margins[offset:offset+length], new_winners[offset:offset+length]
            offset += length
            gaps = sampled_gaps(t, m-row["delta"])
            shift = None if state["gaps"] is None else boundary_shift(state["gaps"], gaps)
            switches = unresolved_switches(m, w, row["delta"])
            stable = shift is not None and shift <= tolerance and switches == 0
            if row["kind"] == "cross_label" and not gaps:
                stable = False
            streak = state["streak"]+1 if stable else 0
            converged = size >= 65 and streak >= 2
            state.update(t=t.copy(), margin=m, winners=w, gaps=gaps, streak=streak,
                         converged=converged, last_shift=shift, unresolved_switches=switches)
            convergence.append(dict(segment_id=row["segment_id"], points=size, gap_count=len(gaps),
                                    gap_fraction=total_width(gaps), boundary_shift_over_dx=shift,
                                    unresolved_switches=switches, stable_comparisons=streak,
                                    converged=converged))
        if output is not None:
            saved.write_csv(output/"convergence.csv", convergence)
        print(f"  {sum(state['converged'] for state in states)}/{len(states)} meet the stability criterion.", flush=True)
    return states, convergence


def analyze_segments(segments, endpoints, states, dx):
    rows, intervals = [], []
    for segment, fields, state in zip(segments, endpoints, states):
        regions, affine = affine_intervals(fields, segment["delta"])
        naive = sampled_gaps(np.array([0., 1.]),
                             np.array([segment["margin_start"], segment["margin_end"]])-segment["delta"])
        direct = state["gaps"]
        endpoint_labels = {segment["winner_start"], segment["winner_end"]}
        row = dict(**segment, dx=dx, reference_points=len(state["t"]), reference_stable=state["converged"],
                   reference_last_boundary_shift_over_dx=state["last_shift"],
                   reference_unresolved_switches=state["unresolved_switches"],
                   direct_gap_count=len(direct), direct_gap_fraction=total_width(direct),
                   direct_gap_length=dx*total_width(direct),
                   direct_max_gap_length=dx*max((hi-lo for lo, hi in direct), default=0.),
                   scalar_gap_count=len(naive), scalar_gap_fraction=total_width(naive),
                   affine_gap_count=len(affine), affine_gap_fraction=total_width(affine),
                   affine_gap_length=dx*total_width(affine),
                   affine_abs_length_error_over_dx=abs(total_width(affine)-total_width(direct)),
                   affine_phase_disagreement_fraction=symmetric_difference(affine, direct),
                   affine_max_boundary_error_over_dx=boundary_shift(affine, direct),
                   direct_additional_winner=any(int(w) not in endpoint_labels and m >= segment["delta"]
                                                for w, m in zip(state["winners"], state["margin"])),
                   affine_additional_winner=any(w not in endpoint_labels for _, _, w in regions))
        rows.append(row)
        for method, gaps in (("direct_sampled", direct), ("affine_components", affine), ("linear_margin", naive)):
            intervals.extend(dict(segment_id=segment["segment_id"], method=method, gap=j,
                                  t_start=lo, t_end=hi, length_over_dx=hi-lo, length=dx*(hi-lo))
                             for j, (lo, hi) in enumerate(gaps))
    return rows, intervals


def summarize(rows):
    groups = []
    for phi, kind in sorted({(r["target_phi"], r["kind"]) for r in rows}):
        group = [r for r in rows if r["target_phi"] == phi and r["kind"] == kind]
        values = lambda name: [r[name] for r in group]
        groups.append(dict(target_phi=phi, kind=kind, segments=len(group),
            stable_references=sum(values("reference_stable")),
            direct_gap_segments=sum(r["direct_gap_count"] > 0 for r in group),
            affine_gap_segments=sum(r["affine_gap_count"] > 0 for r in group),
            scalar_gap_segments=sum(r["scalar_gap_count"] > 0 for r in group),
            direct_gap_fraction_p10=float(np.quantile(values("direct_gap_fraction"), .1)),
            direct_gap_fraction_median=float(np.median(values("direct_gap_fraction"))),
            direct_gap_fraction_p90=float(np.quantile(values("direct_gap_fraction"), .9)),
            affine_abs_length_error_over_dx_median=float(np.median(values("affine_abs_length_error_over_dx"))),
            affine_abs_length_error_over_dx_max=float(np.max(values("affine_abs_length_error_over_dx"))),
            affine_phase_disagreement_median=float(np.median(values("affine_phase_disagreement_fraction"))),
            affine_phase_disagreement_max=float(np.max(values("affine_phase_disagreement_fraction"))),
            additional_direct_winner_segments=sum(values("direct_additional_winner"))))
    return groups


def write_report(output, rows, groups, diagnostics):
    lines = ["# Local directional film check", "",
        "One saved realization; 75%/85% refer to target VF. Segments are sampled",
        "conditionally on occupied face-neighbour endpoints, separately by label relation.",
        "These are directional gap lengths, not normal film thicknesses or global statistics.", "",
        f"Seed validation: {diagnostics['nodes']} nodes; first-field max error "
        f"{diagnostics['first_field_max_abs']:g}, margin max error {diagnostics['margin_max_abs']:g}, "
        f"winner mismatches {diagnostics['winner_mismatches']}.", "",
        "| Target VF | Pair type | Segments | Stable references | Gaps: direct / affine / scalar | Direct gap / dx, median | Affine width error / dx, median / max |",
        "|---|---|---|---|---|---|---|"]
    for r in groups:
        lines.append(f"| {r['target_phi']:.0%} | {r['kind']} | {r['segments']} | {r['stable_references']} | "
                     f"{r['direct_gap_segments']} / {r['affine_gap_segments']} / {r['scalar_gap_segments']} | "
                     f"{r['direct_gap_fraction_median']:.6f} | "
                     f"{r['affine_abs_length_error_over_dx_median']:.6f} / {r['affine_abs_length_error_over_dx_max']:.6f} |")
    lines += ["", "The affine-component calculation clips all winner inequalities exactly along",
              "each segment; it is not a linear interpolation of the already-derived margin.",
              "The direct reference evaluates the finite-mode GSTools realization between nodes.",
              "Nested probes and interpolated threshold crossings provide a numerical stability",
              "check, not certified root isolation. Features narrower than the probes may remain.",
              "Reference stability requires two successive equal-gap-count comparisons with",
              "boundary motion below the configured tolerance, at least 65 samples, and no",
              "unresolved winner switches between occupied consecutive probes.",
              "Unstable references remain flagged in measurements.csv and are included in the",
              "summary; they must not be treated as validated reference widths.",
              "The scalar interpolation has no gap when both endpoint margins exceed delta.",
              "All intervals and interpolation disagreements are measured along the segment.",
              "Same-label controls do not establish 3D connectivity. No global VF or L0 is",
              "re-estimated here; the Adler hitting set and underlying field model are unchanged.", ""]
    (output/"summary.md").write_text("\n".join(lines), encoding="utf-8")


def plot_results(output, rows, states, endpoints):
    groups = summarize(rows)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for phi, color in zip(sorted({r["target_phi"] for r in rows}), ("tab:blue", "tab:orange")):
        group = [r for r in rows if r["target_phi"] == phi and r["kind"] == "cross_label"]
        if not group:
            continue
        widths = np.sort([r["direct_gap_fraction"] for r in group])
        axes[0].step(widths, np.arange(1, len(widths)+1)/len(widths), where="post", color=color, label=f"VF {phi:.0%}")
        axes[1].scatter([r["direct_gap_fraction"] for r in group], [r["affine_gap_fraction"] for r in group],
                        s=14, alpha=.6, color=color, label=f"VF {phi:.0%}")
        axes[2].scatter([r["direct_gap_fraction"] for r in group],
                        [r["affine_phase_disagreement_fraction"] for r in group],
                        s=14, alpha=.6, color=color, label=f"VF {phi:.0%}")
    axes[0].set(title="Cross-label directional gaps", xlabel="Total matrix length / dx", ylabel="Empirical cumulative fraction")
    axes[1].plot([0, 1], [0, 1], "k--", lw=.8)
    axes[1].set(title="Affine components versus direct probes", xlabel="Direct gap length / dx", ylabel="Affine gap length / dx")
    axes[2].set(title="Location and width disagreement", xlabel="Direct gap length / dx", ylabel="Disagreeing segment fraction")
    for ax in axes:
        ax.grid(alpha=.25)
        ax.legend(fontsize=8)
    fig.savefig(output/"film_comparison.png", dpi=160)
    plt.close(fig)
    chosen = []
    for group in groups:
        candidates = [i for i, r in enumerate(rows) if r["target_phi"] == group["target_phi"] and r["kind"] == group["kind"]]
        chosen.append(max(candidates, key=lambda i: rows[i]["affine_phase_disagreement_fraction"]))
    fig, axes = plt.subplots(len(chosen), 1, figsize=(9, 3*len(chosen)), squeeze=False, constrained_layout=True)
    example_rows = []
    for ax, index in zip(axes.flat, chosen):
        row, state, fields = rows[index], states[index], endpoints[index]
        t = np.linspace(0, 1, 1025)
        affine_margin, _ = margin_and_winner(fields[:, :1]*(1-t)+fields[:, 1:] * t)
        scalar_margin = row["margin_start"]*(1-t)+row["margin_end"]*t
        ax.plot(state["t"], state["margin"], "-", label="Direct Gaussian evaluation", lw=2)
        ax.plot(t, affine_margin, "--", label="Affine Gaussian components")
        ax.plot(t, scalar_margin, ":", label="Linear saved margin")
        ax.axhline(row["delta"], color="black", lw=.8, label="Threshold delta")
        for lo, hi in state["gaps"]:
            ax.axvspan(lo, hi, color="gray", alpha=.2)
        ax.set(title=f"VF {row['target_phi']:.0%}, {row['kind']}, segment {row['segment_id']} (largest disagreement in group)",
               xlabel="Distance along segment / dx", ylabel="Margin (field-value units)")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)
        example_rows.append(dict(segment_id=row["segment_id"], target_phi=row["target_phi"], kind=row["kind"],
                                 selection="largest affine/direct disagreement in group"))
    fig.savefig(output/"profile_examples.png", dpi=160)
    plt.close(fig)
    saved.write_csv(output/"profile_examples.csv", example_rows)


def run(config, source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output is not empty; choose a new --output.")
    margin, winners, grid, scene = saved.load_scene(source, config["q"], config["source_phi"])
    first_mesh = pv.read(source/"field_g1.vti")
    if (tuple(first_mesh.dimensions) != margin.shape or not np.allclose(first_mesh.spacing, scene["spacing"])
            or not np.allclose(first_mesh.origin, scene["origin"])):
        raise ValueError("Saved first field does not share the scene geometry.")
    scene["input_sha256"]["field_g1.vti"] = saved.sha256(source/"field_g1.vti")
    sampler = OriginalFields(config["q"], config["lc"], grid.dx, config["field_seed"],
                             config["mode_no"], config["threads"])
    rng = np.random.default_rng(config["selection_seed"])
    segments, pools = [], []
    for phi in config["phis"]:
        delta = saved.adler.delta_for_vf(config["q"], phi)
        for kind, count in (("cross_label", config["contacts"]), ("same_label", config["controls"])):
            selected, pool = select_segments(margin, winners, delta, count, kind, rng)
            pools.append(dict(target_phi=phi, **pool))
            for row in selected:
                segments.append(dict(segment_id=len(segments), target_phi=phi, delta=delta, **row))
    if not segments:
        raise ValueError("No eligible segments were found.")
    starts = np.array([[s["i"], s["j"], s["k"]] for s in segments])
    ends = starts.copy()
    ends[np.arange(len(ends)), [s["axis"] for s in segments]] += 1
    random_nodes = rng.integers(0, grid.n, size=(config["validation_nodes"], 3))
    nodes = np.unique(np.concatenate([starts, ends, random_nodes]), axis=0)
    print(f"Verifying original field streams at {len(nodes)} nodes...", flush=True)
    diagnostics = validate_nodes(sampler, nodes, margin, winners, source)
    raw = sampler.evaluate(np.stack([starts, ends], axis=1).reshape(-1, 3))
    endpoints = raw.reshape(config["q"], len(segments), 2).transpose(1, 0, 2)
    for row, fields in zip(segments, endpoints):
        endpoint_margin, endpoint_winner = margin_and_winner(fields)
        if np.any(endpoint_margin < row["delta"]) or tuple(endpoint_winner) != (row["winner_start"], row["winner_end"]):
            raise ValueError("Double-precision endpoint membership differs from selected saved endpoints.")
    output.mkdir(parents=True, exist_ok=True)
    sources = [Path(__file__).resolve(), HERE/"run_saved_sweep.py", saved.ROOT/"src/fields.py",
               saved.ROOT/"src/adler.py"]
    metadata = dict(schema=SCHEMA, status="running", config=config, source=str(source), scene=scene,
                    component_seeds=sampler.seeds, seed_validation=diagnostics, independent_realizations=1,
                    segment_count=len(segments), pools=pools, coordinate_convention="original voxel-index coordinates; physical length = index length * dx",
                    source_sha256={p.relative_to(saved.ROOT).as_posix(): saved.sha256(p) for p in sources},
                    versions={v: version(v) for v in ("numpy", "scipy", "gstools", "pyvista", "matplotlib")},
                    python=sys.version)
    saved.write_json(output/"metadata.json", metadata)
    saved.write_csv(output/"segments.csv", segments)
    saved.write_csv(output/"sampling_pools.csv", pools)
    start = time.perf_counter()
    try:
        levels = [2**power+1 for power in range(4, int(np.log2(config["max_points"]-1))+1)]
        states, _ = refine_segments(sampler, segments, endpoints, levels, config["tolerance"], output)
        rows, intervals = analyze_segments(segments, endpoints, states, grid.dx)
        saved.write_csv(output/"measurements.csv", rows)
        saved.write_csv(output/"gap_intervals.csv", intervals,
                        ("segment_id", "method", "gap", "t_start", "t_end", "length_over_dx", "length"))
        groups = summarize(rows)
        saved.write_csv(output/"aggregate.csv", groups)
        payload = {"endpoint_fields": endpoints}
        for row, state in zip(segments, states):
            i = row["segment_id"]
            payload.update({f"t_{i}": state["t"], f"margin_{i}": state["margin"], f"winner_{i}": state["winners"]})
        np.savez_compressed(output/"profiles.npz", **payload)
        write_report(output, rows, groups, diagnostics)
        plot_results(output, rows, states, endpoints)
        metadata.update(status="complete", stable_references=sum(s["converged"] for s in states),
                        elapsed_probe_and_report_seconds=time.perf_counter()-start)
        saved.write_json(output/"metadata.json", metadata)
    except BaseException as exc:
        metadata.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed",
                        error=f"{type(exc).__name__}: {exc}")
        saved.write_json(output/"metadata.json", metadata)
        raise
    print(f"Complete: {len(segments)} segments; {metadata['stable_references']} stable references. {output}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=HERE/"paraview/q32_phi75")
    parser.add_argument("--output", type=Path, default=HERE/"results_local_films")
    parser.add_argument("--q", type=int, default=32)
    parser.add_argument("--lc", type=float, default=5.)
    parser.add_argument("--source-phi", type=float, default=.75)
    parser.add_argument("--field-seed", type=int, default=321729)
    parser.add_argument("--selection-seed", type=int, default=119731)
    parser.add_argument("--mode-no", type=int, default=1000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--contacts", type=int, default=100)
    parser.add_argument("--controls", type=int, default=100)
    parser.add_argument("--validation-nodes", type=int, default=256)
    parser.add_argument("--phis", type=float, nargs="+", default=[.75, .85])
    parser.add_argument("--max-points", type=int, default=257)
    parser.add_argument("--tolerance", type=float, default=1e-4, help="Boundary shift in units of dx")
    args = parser.parse_args()
    values = vars(args)
    config = {k: v for k, v in values.items() if k not in ("source", "output")}
    config["phis"] = sorted(set(args.phis))
    if (args.q < 2 or min(args.contacts, args.controls, args.validation_nodes, args.mode_no, args.threads) < 1
            or min(args.field_seed, args.selection_seed) < 0
            or not np.isfinite(args.lc) or args.lc <= 0
            or not np.isfinite(args.tolerance) or args.tolerance <= 0
            or args.max_points < 65 or (args.max_points-1) & (args.max_points-2)
            or any(not np.isfinite(p) or not 0 < p < 1 for p in [args.source_phi, *args.phis])):
        parser.error("Require positive counts/scales, 0<phi<1, nonnegative seeds, and max-points=2^k+1>=65.")
    try:
        run(config, args.source, args.output)
    except (ValueError, FileNotFoundError) as exc:
        parser.exit(2, f"{exc}\n")


if __name__ == "__main__":
    main()
