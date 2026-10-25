"""Experiment 03: Gaussian margin tessellation — high VF, partially analytic.

The margin tessellation A_delta = {max_i G_i(x) - second_max_i G_i(x) >= delta}
retains the winner-cell cores and removes a matrix film around tie interfaces.
By increasing q (number of competitors) the winner cells become smaller and
more numerous, enabling disconnected volume fractions up to 90%+.

The hitting set D_delta in R^q is a disjoint union of polyhedral winner chambers,
so the morphology is a Gaussian-related excursion and the GKF applies.
Closed M_j formulas exist (theory/gaussian_minkowski_functionals.md) and have
been numerically validated, but no published formal stratified-tube proof exists yet.

Run:
    mamba run -n rf_explo python run.py

Outputs (all in results/):
    margin_sweep.csv                  — all realisations, all (q, phi) combinations
    margin_sweep.png                  — spanning frequency and Euler vs phi per q
    margin_euler_vs_theory.png        — measured Euler vs Adler+GMF prediction

Seeds: SEED=1729.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, excursions, fields, measures
from src.fields import GridSpec

RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

SEED = 1729
N = 64
LENGTH = 100.0
LC = 5.0
Q_VALUES = [2, 4, 8, 16]
TARGET_PHIS = [0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90]
N_REALIZATIONS = 5


@dataclass
class Row:
    realization: int
    q: int
    target_phi: float
    delta: float
    theory_phi: float
    measured_phi: float
    euler: int
    n_components: int
    largest_fraction_phase: float
    spans_any: bool
    theory_euler: float


def run_sweep() -> list[Row]:
    grid = GridSpec(n=N, length=LENGTH)
    rows: list[Row] = []
    max_q = max(Q_VALUES)

    print(f"Margin tessellation sweep: n={N}, lc={LC}, q={Q_VALUES}, seed={SEED}")
    print(f"NOTE: GMF formulas are numerically validated but not yet formally proved.\n")

    for realization in range(N_REALIZATIONS):
        rng = np.random.default_rng(SEED + realization)
        all_fields = fields.generate_margin_fields(grid, LC, max_q, rng)
        for q in Q_VALUES:
            fs = all_fields[:q]
            for phi in TARGET_PHIS:
                delta = adler.delta_for_vf(q, phi)
                theory_phi = adler.gaussian_margin_vf(q, delta)
                phase, winners = excursions.margin_excursion(fs, delta)
                m = measures.measure_all(phase, grid)
                theory_euler = adler.margin_expected_euler(grid, LC, q, delta)
                rows.append(Row(
                    realization=realization + 1,
                    q=q,
                    target_phi=phi,
                    delta=delta,
                    theory_phi=theory_phi,
                    measured_phi=m["vf"],
                    euler=int(m["euler"]),
                    n_components=0,   # heavy: skip full component count here
                    largest_fraction_phase=0.0,
                    spans_any=m["percolates_any"],
                    theory_euler=theory_euler,
                ))
                print(f"  r={realization+1} q={q:2d} phi={phi:.2f} "
                      f"measured={m['vf']:.3f} spans={m['percolates_any']} "
                      f"euler_meas={m['euler']:.0f} euler_theory={theory_euler:.1f}")

    return rows


def save_and_plot(rows: list[Row]) -> None:
    with (RESULTS / "margin_sweep.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(Row.__annotations__.keys()))
        w.writeheader()
        w.writerows(asdict(r) for r in rows)
    print("\n  saved margin_sweep.csv")

    # --- spanning frequency and Euler vs phi per q ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0, 0.5, len(Q_VALUES)))

    for q, color in zip(Q_VALUES, colors):
        q_rows = [r for r in rows if r.q == q]
        phis = sorted(set(r.target_phi for r in q_rows))
        span_freq, euler_mean, euler_std, euler_theory = [], [], [], []
        for phi in phis:
            phi_rows = [r for r in q_rows if r.target_phi == phi]
            span_freq.append(sum(r.spans_any for r in phi_rows) / len(phi_rows))
            euler_vals = [r.euler for r in phi_rows]
            euler_mean.append(np.mean(euler_vals))
            euler_std.append(np.std(euler_vals, ddof=1))
            euler_theory.append(phi_rows[0].theory_euler)

        axes[0].plot(phis, span_freq, "o-", color=color, label=f"q={q}")
        axes[1].errorbar(phis, euler_mean, yerr=euler_std, fmt="o", color=color,
                         capsize=3, label=f"q={q} measured")
        axes[1].plot(phis, euler_theory, "--", color=color, alpha=0.7,
                     label=f"q={q} Adler+GMF")

    axes[0].set_xlabel("volume fraction φ")
    axes[0].set_ylabel("spanning frequency")
    axes[0].set_title("Percolation vs volume fraction")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    axes[1].set_xlabel("volume fraction φ")
    axes[1].set_ylabel("Euler characteristic")
    axes[1].set_title("Euler characteristic vs φ\n(dashed = Adler+GMF, not yet formally proved)")
    axes[1].axhline(0, color="k", lw=0.8, ls=":")
    axes[1].legend(fontsize=7, ncol=2)
    axes[1].grid(alpha=0.3)
    fig.suptitle(f"Margin tessellation  (n={N}, lc={LC}, {N_REALIZATIONS} realisations)")
    fig.savefig(RESULTS / "margin_sweep.png", dpi=180)
    plt.close(fig)
    print("  saved margin_sweep.png")


if __name__ == "__main__":
    rows = run_sweep()
    save_and_plot(rows)
    print("Done. Results in experiments/03_margin_tessellation/results/")
