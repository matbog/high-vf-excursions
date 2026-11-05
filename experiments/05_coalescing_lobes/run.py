"""Experiment 05: Coalescing Gaussian-space lobes, union measurement schema v2.

A 3D Gaussian vector field G(x) = (G1(x), G2(x), G3(x)) maps each voxel to a
point in auxiliary Gaussian space R^3. Four fixed target centres c_a are placed
at the vertices of a regular tetrahedron of radius rho. A voxel is included when
its Gaussian-space point lies within distance r of any centre:

    A_r = {x : min_a ||G(x) - c_a|| <= r}

The single parameter r controls morphology continuously. Lobe a and lobe b first
touch in value space at the contact radius:

    r_c = ||c_a - c_b|| / 2 = rho * sqrt(2/3)

Before r_c the lobes are disjoint value-space balls (GKF M_j are additive).
After r_c they overlap and inclusion-exclusion terms enter the M_j formulas —
this contribution has not yet been analytically derived.

Measure the binary union under both connectivities. Keep ball-preimage sums
and nearest-centre partition diagnostics in explicitly separate columns.
The historical results/ files use an older schema and are never overwritten.

Run:
    mamba run -n rf_explo python run.py

Outputs (by default in results_union_v2/):
    coalescing_lobes_sweep.csv      — all realisations, all radii
    coalescing_lobes_slices.png     — winner-coloured central slices
    coalescing_lobes_diagnostics.png — 4-panel diagnostic

Seeds: SEED=9417.
"""

from __future__ import annotations

import csv
import argparse
import hashlib
import json
import sys
from importlib.metadata import version
from pathlib import Path
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, excursions, fields
from src.adler import expected_lkcs_lobes_disjoint
from src.fields import GridSpec
from measurement import measure_lobes

RESULTS = Path(__file__).parent / "results_union_v2"

SEED = 9417
N = 100
LENGTH = 100.0
LC = 10.0
PADDING = 20
RHO = 1.5
N_REALIZATIONS = 4
RADIUS_OVER_RC = [0.45, 0.65, 0.85, 0.98, 1.02, 1.15, 1.35, 1.60]
MC_SAMPLES = 2_000_000


Row = SimpleNamespace


def _sample_std(values) -> float:
    """No dispersion estimate is available for the one-realization smoke run."""
    return float(np.std(values, ddof=1)) if len(values) > 1 else 0.0


def run() -> None:
    if RESULTS.exists() and any(RESULTS.iterdir()):
        raise SystemExit(f"Output directory is not empty: {RESULTS}. Choose a new --output.")
    RESULTS.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[2]
    sources = [Path(__file__).resolve(), Path(__file__).with_name("measurement.py").resolve(),
               root / "src/adler.py", root / "src/measures.py", root / "src/fields.py"]
    metadata = {
        "measurement_schema": "union_v2", "status": "running", "seed": SEED,
        "n": N, "length": LENGTH, "lc": LC, "rho": RHO, "padding": PADDING,
        "realizations": N_REALIZATIONS, "radius_over_rc": RADIUS_OVER_RC,
        "mc_samples": MC_SAMPLES, "python": sys.version,
        "versions": {name: version(name) for name in ("numpy", "scipy", "gstools", "scikit-image", "matplotlib")},
        "generator": "GSTools Gaussian rescale=1; shared src.fields wrapper",
        "source_sha256": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    }
    (RESULTS / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    grid = GridSpec(n=N, length=LENGTH)
    centres = excursions.tetrahedral_centres(RHO)
    rc = float(np.linalg.norm(centres[0] - centres[1]) / 2.0)
    ratios = np.array(RADIUS_OVER_RC)
    radii = rc * ratios

    print(f"Coalescing lobes: n={N}, lc={LC}, rho={RHO}, r_c={rc:.6f}")
    print(f"r/r_c sweep: {RADIUS_OVER_RC}")
    print(f"NOTE: GKF L0..L3 theory applies ONLY for r/r_c < 1 (balls disjoint).\n")

    print("Estimating Gaussian-space volume fractions by MC...")
    theory_phi = adler.gaussian_lobe_vf_mc(centres, radii, n_samples=MC_SAMPLES,
                                            seed=SEED + 100_000)

    rows: list[Row] = []
    distances_all, nearest_all = [], []

    for realization in range(N_REALIZATIONS):
        vf = fields.generate_vector_field(grid, LC, 3, np.random.default_rng(SEED + realization),
                                          padding=PADDING)
        distances, nearest = excursions.lobe_distances(vf, centres)
        distances_all.append(distances)
        nearest_all.append(nearest)

        # Per-ball distance field: d_a(x) = ||G(x) - c_a||. The per-ball signed
        # field at radius r is  s_a(x) = r - d_a(x)  (positive inside ball a).
        # For disjoint balls (r < r_c) the per-ball preimages partition A_r,
        # so per-ball LKC sum = LKC(A_r).
        per_ball_dists = np.stack([
            np.linalg.norm(vf - c[:, None, None, None], axis=0) for c in centres
        ])   # shape (K, N, N, N)

        for radius, ratio, phi in zip(radii, ratios, theory_phi):
            measured = measure_lobes(per_ball_dists, float(radius), grid)

            # Full LKC theory (only valid for r/r_c < 1)
            if ratio < 1.0:
                th = expected_lkcs_lobes_disjoint(grid, LC, centres, float(radius))
                th_L0, th_L1, th_L2, th_L3 = th["L0"], th["L1"], th["L2"], th["L3"]
            else:
                th_L0 = th_L1 = th_L2 = th_L3 = float("nan")

            rows.append(Row(
                realization=realization + 1,
                radius=float(radius),
                radius_over_rc=float(ratio),
                theory_phi=float(phi),
                theory_L0=th_L0, theory_L1=th_L1, theory_L2=th_L2, theory_L3=th_L3,
                **measured,
            ))
            def _th(v: float) -> str:
                return f"{v:7.1f}" if (ratio < 1.0 and not np.isnan(v)) else "    n/a"
            print(f"  r={realization+1} r/r_c={ratio:.2f} phi={measured['measured_phi']:.3f}  "
                  f"union L0(6/26)={measured['union_euler_6']}/{measured['union_euler_26']} "
                  f"theory={_th(th_L0)}")

    with (RESULTS / "coalescing_lobes_sweep.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(vars(rows[0])))
        w.writeheader()
        w.writerows(vars(r) for r in rows)
    print("\n  saved coalescing_lobes_sweep.csv")

    _plot_slices(distances_all, nearest_all, radii, ratios, theory_phi)
    _plot_diagnostics(rows, ratios, theory_phi)
    _plot_lkcs(rows, ratios)
    metadata["status"] = "complete"
    (RESULTS / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")


def _plot_slices(distances_all, nearest_all, radii, ratios, theory_phi) -> None:
    cmap = plt.colormaps["tab10"].with_extremes(bad="white")
    mid = N // 2
    fig, axes = plt.subplots(N_REALIZATIONS, len(radii), squeeze=False,
                             figsize=(2.7 * len(radii), 2.6 * N_REALIZATIONS),
                             constrained_layout=True)
    for row in range(N_REALIZATIONS):
        for col, (radius, ratio, phi) in enumerate(zip(radii, ratios, theory_phi)):
            phase = distances_all[row][:, :, mid] <= radius
            image = np.where(phase.T, nearest_all[row][:, :, mid].T.astype(float), np.nan)
            ax = axes[row, col]
            ax.imshow(image, origin="lower", cmap=cmap, vmin=0, vmax=9,
                      interpolation="nearest")
            ax.set(xticks=[], yticks=[])
            if row == 0:
                ax.set_title(rf"$r/r_c={ratio:.2f}$" + "\n" + rf"$\phi={phi:.2f}$")
            ax.text(0.03, 0.04, rf"$\phi_{{slice}}={phase.mean():.2f}$",
                    transform=ax.transAxes, fontsize=7,
                    bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
        axes[row, 0].set_ylabel(f"realisation {row + 1}")
    fig.suptitle(rf"Tetrahedral Gaussian lobes (value-space contact at $r/r_c=1$)",
                 fontsize=13)
    fig.savefig(RESULTS / "coalescing_lobes_slices.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  saved coalescing_lobes_slices.png")


def _plot_diagnostics(rows: list[Row], ratios, theory_phi) -> None:
    def mean_by_ratio(attr: str) -> np.ndarray:
        return np.array([
            np.mean([getattr(r, attr) for r in rows if np.isclose(r.radius_over_rc, ratio)])
            for ratio in ratios
        ])
    def std_by_ratio(attr: str) -> np.ndarray:
        return np.array([
            _sample_std([getattr(r, attr) for r in rows if np.isclose(r.radius_over_rc, ratio)])
            for ratio in ratios
        ])
    def theory_vals(attr: str) -> np.ndarray:
        return np.array([
            getattr(next(r for r in rows if np.isclose(r.radius_over_rc, ratio)), attr)
            for ratio in ratios
        ])

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax in axes.flat:
        ax.axvline(1.0, color="crimson", ls="--", lw=1.2, label=r"value-space contact $r_c$")
        ax.set_xlabel(r"$r/r_c$")
        ax.grid(alpha=0.3)

    axes[0, 0].plot(ratios, theory_phi, "k-", label="Gaussian-space MC")
    axes[0, 0].errorbar(ratios, mean_by_ratio("measured_phi"),
                         yerr=std_by_ratio("measured_phi"), fmt="o", ms=4, label="3D RF")
    axes[0, 0].set_ylabel("volume fraction")
    axes[0, 0].set_title("Volume calibration")
    axes[0, 0].legend(fontsize=8)

    # L0 (Euler) with theory curve for r/r_c < 1
    for rank, style in ((6, "s--"), (26, "o-")):
        key = f"union_euler_{rank}"
        axes[0, 1].errorbar(ratios, mean_by_ratio(key), yerr=std_by_ratio(key),
                             fmt=style, ms=4, label=f"binary union ({rank}-conn)")
    th_L0 = theory_vals("theory_L0")
    mask  = ~np.isnan(th_L0)
    axes[0, 1].plot(ratios[mask], th_L0[mask], "-", color="#ff6e7c", lw=2, label="GKF theory (disjoint)")
    axes[0, 1].axhline(0, color="k", lw=0.8)
    axes[0, 1].set_ylabel("Euler characteristic")
    axes[0, 1].set_title("L0 — Euler\n(theory valid only before contact)")
    axes[0, 1].legend(fontsize=8, loc="best")

    for rank in (6, 26):
        axes[1, 0].plot(ratios, mean_by_ratio(f"union_largest_fraction_{rank}"), "o-", ms=4, label=f"{rank}-conn")
    axes[1, 0].legend(fontsize=8)
    axes[1, 0].set_ylabel("largest component / phase")
    axes[1, 0].set_title("Component coalescence")

    for rank in (6, 26):
        axes[1, 1].plot(ratios, mean_by_ratio(f"union_spans_{rank}"), "o-", ms=4, label=f"{rank}-conn")
    axes[1, 1].legend(fontsize=8)
    axes[1, 1].set_ylim(-0.05, 1.05)
    axes[1, 1].set_ylabel("spanning frequency")
    axes[1, 1].set_title("Binary-union spanning (finite grid)")

    fig.suptitle(rf"Coalescing lobes  ($N={N}^3$, $L_c={LC}$, {N_REALIZATIONS} realisations)")
    fig.savefig(RESULTS / "coalescing_lobes_diagnostics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  saved coalescing_lobes_diagnostics.png")


def _plot_lkcs(rows: list[Row], ratios) -> None:
    """Union measurements throughout; pre-contact per-ball controls also shown."""
    def mean(attr):
        return np.array([
            np.mean([getattr(r, attr) for r in rows if np.isclose(r.radius_over_rc, ratio)])
            for ratio in ratios
        ])
    def std(attr):
        return np.array([
            _sample_std([getattr(r, attr) for r in rows if np.isclose(r.radius_over_rc, ratio)])
            for ratio in ratios
        ])
    def th(attr):
        return np.array([
            getattr(next(r for r in rows if np.isclose(r.radius_over_rc, ratio)), attr)
            for ratio in ratios
        ])

    specs = [
        ("union_euler_26", "theory_L0", "L0 — union Euler (26-conn)", "sum_ball_euler_26"),
        ("union_L1_field", "theory_L1", "L1 — union signed-field mesh", "sum_ball_L1_field"),
        ("union_L2_crofton", "theory_L2", "L2 — union half-area (Crofton)", "sum_ball_L2_crofton"),
        ("union_L2_field", "theory_L2", "L2 — union half-area (field mesh)", "sum_ball_L2_field"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.6), constrained_layout=True)
    for ax, (meas, th_key, title, control) in zip(axes, specs):
        m = mean(meas); s = std(meas); t = th(th_key)
        ok = ~np.isnan(t)
        ax.errorbar(ratios, m, yerr=s, fmt="o-", color="#3a8cb4", ecolor="#3a8cb4",
                    capsize=3, ms=4, lw=1.4, label="binary/scalar union")
        ax.plot(ratios[ok], mean(control)[ok], "x:", color="0.4",
                label="per-ball sum (pre-contact only)")
        ax.plot(ratios[ok], t[ok], "-", color="#ff6e7c", lw=2.0,
                label="GKF theory (r/r_c < 1)")
        ax.axvline(1.0, color="crimson", ls="--", lw=1.0,
                   label=r"value-space contact $r_c$")
        ax.set_xlabel(r"$r/r_c$", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3, lw=0.5)
        ax.legend(fontsize=8, loc="best")

    fig.suptitle(
        f"Coalescing lobes — LKC measurement vs GKF theory  "
        f"(N={N}, LC={LC}, {N_REALIZATIONS} realisations; theory valid for r/r_c < 1)",
        fontsize=12,
    )
    fig.savefig(RESULTS / "coalescing_lobes_lkcs.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("  saved coalescing_lobes_lkcs.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=RESULTS)
    parser.add_argument("--n", type=int, default=N)
    parser.add_argument("--realizations", type=int, default=N_REALIZATIONS)
    parser.add_argument("--radii", type=float, nargs="+", default=RADIUS_OVER_RC, help="Ratios r/r_c")
    parser.add_argument("--mc-samples", type=int, default=MC_SAMPLES)
    args = parser.parse_args()
    if args.n < 4 or args.realizations < 1 or args.mc_samples < 1 or any(not np.isfinite(r) or r <= 0 for r in args.radii):
        parser.error("Require n>=4, positive realizations/sample count and positive finite radii")
    RESULTS, N, N_REALIZATIONS = args.output, args.n, args.realizations
    RADIUS_OVER_RC, MC_SAMPLES = sorted(set(args.radii)), args.mc_samples
    run()
    print(f"Done. Results in {RESULTS}")
