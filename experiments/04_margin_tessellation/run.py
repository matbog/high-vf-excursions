"""Experiment 04: Gaussian margin tessellation — high VF, partially analytic.

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
    margin_sweep.csv          — all realisations, all (q, phi) combinations
    margin_sweep.png          — 2-panel overview: spanning freq + L0 fill_between vs phi
    margin_main.png           — 4 rows × 5 cols: L0, L2, spanning, ratio, slices per q
    margin_euler_ratio.png    — measured/theory Euler ratio vs phi per q

Seeds: SEED=1729.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.gridspec as mgridspec
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, excursions, fields, measures
from src.adler import expected_lkcs_margin
from src.fields import GridSpec
from src.measures import (
    measure_l1_from_field,
    measure_percolation_per_winner,
    measure_surface_area_crofton,
    measure_surface_area_from_field,
)

RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

SEED = 1729
# Physical domain fixed, voxel resolution varied for a built-in convergence
# study. Previously split between exp04 (single resolution) and exp06
# (convergence); now merged. N=128 is the main validation resolution, with
# N=64 and N=256 bracketing it for the convergence story.
LENGTH = 64.0
LC = 5.0
N_VALUES = [64, 128, 256]
N_MAIN   = 128     # resolution used for the "main validation" figures
Q_VALUES = [2, 4, 8, 16]
TARGET_PHIS = [0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.85, 0.90]
N_REALIZATIONS = 5

# Style constants — match experiments 01–03
C_THEORY       = "#ff6e7c"
C_NUM          = "#3a8cb4"
Q_COLORS       = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759"]
PHI_VIS        = [0.50, 0.70, 0.90]
PHI_VIS_COLORS = ["#d62728", "#2ca02c", "#9467bd"]

COL_TITLES = [
    "L0 — Euler characteristic",
    "L1 — int. mean curvature / π",
    "L2 — half surface area",
    "L3 — volume fraction φ",
]


def _style_ax(ax: plt.Axes) -> None:
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_xlabel("volume fraction φ", fontsize=9)


@dataclass
class Row:
    realization: int
    N: int                 # grid resolution (added for the merged convergence study)
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
    l1_winner: float
    l2_winner: float       # Crofton (voxel-based) — primary L2 estimator
    l2_field: float        # MC on signed field — comparison L2 estimator
    theory_l1: float
    theory_l2: float


def run_sweep() -> list[Row]:
    rows: list[Row] = []
    max_q = max(Q_VALUES)

    print(f"Margin tessellation sweep: N={N_VALUES}, lc={LC}, q={Q_VALUES}, seed={SEED}")
    print(f"NOTE: GMF formulas are numerically validated but not yet formally proved.\n")

    for N in N_VALUES:
        grid = GridSpec(n=N, length=LENGTH)
        print(f"\n=== N = {N}  (dx={grid.dx:.3f}, LC_vox={LC/grid.dx:.1f}) ===")

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
                    theory_lkcs = expected_lkcs_margin(grid, LC, q, delta)
                    # Per-winner sums: disjoint winner types → additive for all LKCs.
                    # L0 uses 26-conn, L2 uses Crofton (primary) and signed-field MC
                    # (comparison), L1 uses marching cubes on the smooth signed field
                    # s_i(x) = G_i(x) − max_{j≠i} G_j(x) − δ.
                    sorted_top2 = np.partition(fs, kth=q - 2, axis=0)[-2:]
                    max_field   = sorted_top2[1]
                    sec_field   = sorted_top2[0]
                    euler_per_winner = 0
                    l1_per_winner    = 0.0
                    l2_crofton       = 0.0
                    l2_from_field    = 0.0
                    for label in range(q):
                        mask = phase & (winners == label)
                        euler_per_winner += int(measures.measure_euler(mask, connectivity=3))
                        l2_crofton       += float(measure_surface_area_crofton(mask, grid.dx))
                        is_winner        = (winners == label)
                        max_other        = np.where(is_winner, sec_field, max_field)
                        s_i              = fs[label] - max_other - delta
                        l1_per_winner    += float(measure_l1_from_field(s_i, grid.dx, level=0.0))
                        l2_from_field    += float(measure_surface_area_from_field(s_i, grid.dx, level=0.0))
                    spans_winner = measure_percolation_per_winner(phase, winners, q)
                    rows.append(Row(
                        realization=realization + 1,
                        N=N,
                        q=q,
                        target_phi=phi,
                        delta=delta,
                        theory_phi=theory_phi,
                        measured_phi=m["vf"],
                        euler=euler_per_winner,
                        n_components=0,
                        largest_fraction_phase=0.0,
                        spans_any=spans_winner,
                        theory_euler=theory_lkcs["L0"],
                        l1_winner=l1_per_winner,
                        l2_winner=l2_crofton,
                        l2_field=l2_from_field,
                        theory_l1=theory_lkcs["L1"],
                        theory_l2=theory_lkcs["L2"],
                    ))
                    print(f"  N={N:3d} r={realization+1} q={q:2d} phi={phi:.2f} "
                          f"L0={euler_per_winner:.0f}/{theory_lkcs['L0']:.1f} "
                          f"L1={l1_per_winner:.0f}/{theory_lkcs['L1']:.1f} "
                          f"L2_Cr={l2_crofton:.0f} L2_fld={l2_from_field:.0f} "
                          f"/{theory_lkcs['L2']:.1f}")

    return rows


def save_csv(rows: list[Row]) -> None:
    with (RESULTS / "margin_sweep.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(Row.__annotations__.keys()))
        w.writeheader()
        w.writerows(asdict(r) for r in rows)
    print("\n  saved margin_sweep.csv")


def make_overview_figure(rows: list[Row]) -> None:
    """2-panel overview: spanning freq (top) + L0 fill_between (bottom). Matches percolation.png."""
    fig, (ax_span, ax_L0) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, constrained_layout=True
    )

    for q, color in zip(Q_VALUES, Q_COLORS):
        q_rows  = [r for r in rows if r.q == q]
        phis    = sorted(set(r.target_phi for r in q_rows))
        l0_mu, l0_sd, l0_th, span_freq = [], [], [], []
        for phi in phis:
            pr = [r for r in q_rows if r.target_phi == phi]
            vals = [r.euler for r in pr]
            l0_mu.append(np.mean(vals))
            l0_sd.append(np.std(vals, ddof=1) if len(vals) > 1 else 0.0)
            l0_th.append(pr[0].theory_euler)
            span_freq.append(sum(r.spans_any for r in pr) / len(pr))

        phis_arr  = np.array(phis)
        l0_mu_arr = np.array(l0_mu)
        l0_sd_arr = np.array(l0_sd)

        ax_span.plot(phis, span_freq, "o-", color=color, lw=1.5, ms=5, label=f"q = {q}")

        ax_L0.fill_between(phis_arr, l0_mu_arr - l0_sd_arr, l0_mu_arr + l0_sd_arr,
                           color=color, alpha=0.15, lw=0)
        ax_L0.plot(phis_arr, l0_mu_arr, color=color, lw=1.5, alpha=0.8)
        ax_L0.plot(phis_arr, l0_th,     color=color, lw=2.0, ls="--", label=f"q = {q}")

    ax_span.axhline(0.5, color="0.3", lw=1.0, ls="--", label="p = 0.5")
    ax_span.set_ylabel("per-winner spanning frequency", fontsize=10)
    ax_span.set_title("Per-winner spanning frequency vs φ", fontsize=10)
    ax_span.legend(fontsize=8)
    ax_span.grid(alpha=0.25, lw=0.5)
    ax_span.set_ylim(-0.05, 1.05)

    ax_L0.axhline(0, color="0.3", lw=1.0, ls="--", label="E[L₀] = 0")
    ax_L0.set_xlabel("volume fraction φ", fontsize=10)
    ax_L0.set_ylabel("Euler characteristic L0", fontsize=10)
    ax_L0.set_title("L0 — Euler characteristic  (solid = MC mean ± std,  dashed = GKF theory)",
                    fontsize=10)
    ax_L0.legend(fontsize=8)
    ax_L0.grid(alpha=0.25, lw=0.5)

    fig.suptitle(
        f"Margin tessellation — overview  (n={N_MAIN}, lc={LC}, {N_REALIZATIONS} realisations)",
        fontsize=13,
    )
    fig.savefig(RESULTS / "margin_sweep.png", dpi=180)
    plt.close(fig)
    print("  saved margin_sweep.png")


def make_main_figure(rows: list[Row]) -> None:
    """4 rows × 5 cols: L0, L2, spanning, Euler ratio, slices. Matches mc_gaussian.png style.

    All rows should be at a single resolution (callers pass the N=N_MAIN subset).
    """
    grid     = GridSpec(n=N_MAIN, length=LENGTH)
    rng      = np.random.default_rng(SEED)
    all_flds = fields.generate_margin_fields(grid, LC, max(Q_VALUES), rng)
    cmap_w   = plt.cm.tab20

    n_rows = len(Q_VALUES)
    fig    = plt.figure(figsize=(21, 4 * n_rows))
    outer  = mgridspec.GridSpec(
        n_rows, 5, figure=fig,
        width_ratios=[1, 1, 1, 1, 1.1],
        left=0.06, right=0.985, top=0.91, bottom=0.08,
        wspace=0.38, hspace=0.42,
    )

    for row_idx, q in enumerate(Q_VALUES):
        q_rows = [r for r in rows if r.q == q]
        phis   = sorted(set(r.target_phi for r in q_rows))
        first_row = (row_idx == 0)

        l0_mean, l0_std, l0_th = [], [], []
        l1_mean, l1_std, l1_th = [], [], []
        l2_mean, l2_std, l2_th = [], [], []
        phi_mean, phi_std, phi_th = [], [], []
        for phi in phis:
            pr = [r for r in q_rows if r.target_phi == phi]
            l0_vals  = [r.euler      for r in pr]
            l1_vals  = [r.l1_winner  for r in pr]
            l2_vals  = [r.l2_winner  for r in pr]
            phi_vals = [r.measured_phi for r in pr]
            l0_mean.append(np.mean(l0_vals));  l0_std.append(np.std(l0_vals, ddof=1) if len(l0_vals) > 1 else 0.0)
            l0_th.append(pr[0].theory_euler)
            l1_mean.append(np.mean(l1_vals));  l1_std.append(np.std(l1_vals, ddof=1) if len(l1_vals) > 1 else 0.0)
            l1_th.append(pr[0].theory_l1)
            l2_mean.append(np.mean(l2_vals));  l2_std.append(np.std(l2_vals, ddof=1) if len(l2_vals) > 1 else 0.0)
            l2_th.append(pr[0].theory_l2)
            phi_mean.append(np.mean(phi_vals)); phi_std.append(np.std(phi_vals, ddof=1) if len(phi_vals) > 1 else 0.0)
            phi_th.append(pr[0].theory_phi)

        def _vlines(ax: plt.Axes) -> None:
            for pv, pc in zip(PHI_VIS, PHI_VIS_COLORS):
                ax.axvline(pv, color=pc, lw=1.2, ls=":", alpha=0.85, zorder=0)

        lkc_cols = [
            (l0_mean, l0_std, l0_th, True),   # L0: has zero-crossing
            (l1_mean, l1_std, l1_th, False),
            (l2_mean, l2_std, l2_th, False),
            (phi_mean, phi_std, phi_th, False),
        ]
        for col_idx, (mu, sd, th, has_zero) in enumerate(lkc_cols):
            ax = fig.add_subplot(outer[row_idx, col_idx])
            _vlines(ax)
            ax.errorbar(phis, mu, yerr=sd,
                        fmt="o", color=C_NUM, ecolor=C_NUM, capsize=3, ms=3, lw=1.0,
                        label=f"Numerical (n={N_REALIZATIONS})" if (first_row and col_idx == 0) else None)
            ax.plot(phis, th, color=C_THEORY, lw=2.0,
                    label="GKF theory" if (first_row and col_idx == 0) else None)
            if has_zero:
                ax.axhline(0, color="0.4", lw=0.8, ls=":")
            _style_ax(ax)
            if first_row:
                ax.set_title(COL_TITLES[col_idx], fontsize=10)
            if col_idx == 0:
                ax.set_ylabel(f"q = {q}", fontsize=11, fontweight="bold")
                if first_row:
                    ax.legend(fontsize=8, loc="best")

        # Col 4: 3 stacked mid-plane slices
        fs       = all_flds[:q]
        w_colors = cmap_w(np.linspace(0, 0.95, q))
        inner    = mgridspec.GridSpecFromSubplotSpec(
            3, 1, subplot_spec=outer[row_idx, 4], hspace=0.30
        )
        top_thumb = None
        for k, (phi_v, phi_color) in enumerate(zip(PHI_VIS, PHI_VIS_COLORS)):
            ax_img = fig.add_subplot(inner[k])
            if k == 0:
                top_thumb = ax_img
            delta  = adler.delta_for_vf(q, phi_v)
            phase, winners = excursions.margin_excursion(fs, delta)
            ps = phase[:, :, N_MAIN // 2]
            ws = winners[:, :, N_MAIN // 2]
            img = np.full((*ps.shape, 4), [0.92, 0.92, 0.92, 1.0])
            for label in range(q):
                mask = ps & (ws == label)
                img[mask] = [*w_colors[label][:3], 1.0]
            ax_img.imshow(img, origin="lower", interpolation="nearest")
            ax_img.set_xticks([])
            ax_img.set_yticks([])
            for sp in ax_img.spines.values():
                sp.set_edgecolor(phi_color)
                sp.set_linewidth(2.0)
            ax_img.set_title(f"φ = {phi_v:.2f}", fontsize=8, color=phi_color,
                             fontweight="bold", pad=3)

        if first_row and top_thumb is not None:
            pos = top_thumb.get_position()
            fig.text(
                (pos.x0 + pos.x1) / 2, pos.y1 + 0.012,
                "mid-plane slice  (realisation 1)",
                ha="center", va="bottom", fontsize=8, color="0.45",
            )

    fig.suptitle(
        f"Margin tessellation — Monte Carlo validation  "
        f"(n={N_MAIN}, lc={LC}, {N_REALIZATIONS} realisations)",
        fontsize=13,
    )
    fig.savefig(RESULTS / "margin_main.png", dpi=180)
    plt.close(fig)
    print("  saved margin_main.png")


def make_ratio_figure(rows: list[Row]) -> None:
    """Percent-deviation plot: 100 * (mean - theory) / theory vs phi, per q and LKC.

    Dots are the mean % deviation across realisations. The lighter error whiskers
    show the realisation-to-realisation sampling noise translated into % of theory:
    ± std(measured) / |theory| * 100%. Points whose |theory| is below a small
    threshold (where ratios would blow up) are skipped.
    """
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.5), sharey=True,
                             constrained_layout=True)
    lkc_specs = [
        ("euler",      "theory_euler", "L0 — Euler  (26-conn)"),
        ("l1_winner",  "theory_l1",    "L1 — mean curvature  (signed-field MC)"),
        ("l2_winner",  "theory_l2",    "L2 — surface area  (Crofton)"),
        ("l2_field",   "theory_l2",    "L2 — surface area  (signed-field MC)"),
    ]

    for ax, (meas_key, th_key, title) in zip(axes, lkc_specs):
        for q, color in zip(Q_VALUES, Q_COLORS):
            q_rows = [r for r in rows if r.q == q]
            phis   = sorted(set(r.target_phi for r in q_rows))
            pct_mean, pct_noise, phi_plot = [], [], []
            for phi in phis:
                pr   = [r for r in q_rows if r.target_phi == phi]
                vals = np.array([getattr(r, meas_key) for r in pr], dtype=float)
                th   = float(pr[0].__dict__[th_key])
                if abs(th) < 10.0 or len(vals) < 2:
                    continue
                pct_mean.append(100.0 * (vals.mean() - th) / th)
                pct_noise.append(100.0 * vals.std(ddof=1) / abs(th))
                phi_plot.append(phi)
            ax.errorbar(phi_plot, pct_mean, yerr=pct_noise, fmt="o-",
                        color=color, ecolor=color, capsize=3, ms=5, lw=1.6,
                        label=f"q = {q}")

        ax.axhline(0.0, color=C_THEORY, lw=2.0, ls="--", label="theory (0%)")
        ax.axhspan(-5.0, 5.0, color="0.90", alpha=0.5, zorder=0)
        ax.set_xlabel("volume fraction φ", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.25, lw=0.5)
    axes[0].set_ylabel("% deviation  (mean − theory) / theory × 100", fontsize=10)
    axes[-1].legend(fontsize=9, loc="best")

    fig.suptitle(
        "LKC measurement bias vs GKF theory  "
        "(whiskers = realisation std / |theory| × 100%; grey band = ±5%)",
        fontsize=12,
    )
    fig.savefig(RESULTS / "margin_euler_ratio.png", dpi=180)
    plt.close(fig)
    print("  saved margin_euler_ratio.png  (now a % deviation plot for L0/L1/L2)")


def run_analytical_checks() -> None:
    """Verify two identities before the main sweep.

    1. M1 derivative identity:  M1(q, delta) = -sqrt(2) * dM0/ddelta
    2. q=2 Gaussian identity:   margin_gmfs(2, delta) ≈ 2 * gaussian_Mj(delta/sqrt(2))
    """
    print("\nAnalytical checks:")
    eps = 1e-5

    # --- Check 1: M1 = -sqrt(2) * dM0/ddelta ---
    print("  [1] M1 = -sqrt(2) * dM0/ddelta")
    test_cases = [(2, 0.5), (4, 0.3), (8, 0.2), (16, 0.1)]
    max_err = 0.0
    for q, delta in test_cases:
        gmf = adler.margin_gmfs(q, delta)
        m0_plus = adler.margin_m0(q, delta + eps)
        m0_minus = adler.margin_m0(q, delta - eps)
        deriv_approx = (m0_plus - m0_minus) / (2.0 * eps)
        m1_from_deriv = -np.sqrt(2.0) * deriv_approx
        rel = abs(gmf[1] - m1_from_deriv) / (abs(m1_from_deriv) + 1e-15)
        max_err = max(max_err, rel)
        print(f"    q={q:2d} delta={delta:.2f}  M1={gmf[1]:.6f}  from_deriv={m1_from_deriv:.6f}  rel_err={rel:.2e}")
    status = "PASS" if max_err < 1e-4 else "FAIL"
    print(f"    => {status}  (max rel error {max_err:.2e})\n")

    # --- Check 2: q=2 identity ---
    # margin_gmfs(2, delta) should equal 2 * [M0, M1, M2, M3] of a Gaussian at delta/sqrt(2)
    print("  [2] q=2 identity: margin_gmfs(2, delta) == 2 * gaussian_Mj(delta/sqrt(2))")
    max_err2 = 0.0
    for delta in [0.2, 0.5, 1.0, 1.5]:
        gmf2 = adler.margin_gmfs(2, delta)
        u = delta / np.sqrt(2.0)
        g_mj = np.array([
            float(adler.gaussian_m0(np.array([u]))[0]),
            float(adler.gaussian_m1(np.array([u]))[0]),
            float(adler.gaussian_m2(np.array([u]))[0]),
            float(adler.gaussian_m3(np.array([u]))[0]),
        ])
        ref = 2.0 * g_mj
        rel = np.abs(gmf2 - ref) / (np.abs(ref) + 1e-15)
        max_err2 = max(max_err2, float(np.max(rel)))
        print(f"    delta={delta:.2f}  rel_err per M_j: {rel}")
    status2 = "PASS" if max_err2 < 1e-4 else "FAIL"
    print(f"    => {status2}  (max rel error {max_err2:.2e})\n")


def make_convergence_figure(rows: list[Row]) -> None:
    """Convergence of measured/theory LKC ratios vs N, per q, averaged over φ.

    Shows how the N=64 → 128 → 256 refinement reduces estimator biases for
    each LKC. Replaces the old exp06 convergence figure now that exp04 is
    the single home for margin-tessellation validation.
    """
    nq = len(Q_VALUES)
    fig, axes = plt.subplots(1, 4, figsize=(19, 4.5), sharex=True,
                             constrained_layout=True)
    lkc_specs = [
        ("euler",      "theory_euler", "L0 — Euler (26-conn)"),
        ("l1_winner",  "theory_l1",    "L1 — mean curvature (signed-field MC)"),
        ("l2_winner",  "theory_l2",    "L2 — surface area (Crofton)"),
        ("l2_field",   "theory_l2",    "L2 — surface area (signed-field MC)"),
    ]

    for ax, (meas_key, th_key, title) in zip(axes, lkc_specs):
        for q, color in zip(Q_VALUES, Q_COLORS):
            xs, pct, noise = [], [], []
            for N in N_VALUES:
                vals, ths = [], []
                for r in rows:
                    if r.N != N or r.q != q:
                        continue
                    th = float(r.__dict__[th_key])
                    if abs(th) < 50.0:
                        continue
                    vals.append(float(getattr(r, meas_key)) / th)
                    ths.append(th)
                if len(vals) < 2:
                    continue
                vals = np.array(vals)
                xs.append(N)
                pct.append(100.0 * (vals.mean() - 1.0))
                noise.append(100.0 * vals.std(ddof=1))
            if xs:
                ax.errorbar(xs, pct, yerr=noise, fmt="o-", color=color,
                            ecolor=color, capsize=3, ms=5, lw=1.6, label=f"q = {q}")
        ax.axhline(0.0, color=C_THEORY, lw=2.0, ls="--", label="theory (0%)")
        ax.axhspan(-5.0, 5.0, color="0.90", alpha=0.5, zorder=0)
        ax.set_xscale("log", base=2)
        ax.set_xticks(N_VALUES); ax.set_xticklabels([str(n) for n in N_VALUES])
        ax.set_xlabel("grid resolution N", fontsize=11)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.25, lw=0.5)
    axes[0].set_ylabel("% deviation  (mean − theory) / theory × 100", fontsize=10)
    axes[-1].legend(fontsize=9, loc="best")

    fig.suptitle(
        f"Resolution convergence  (LENGTH={LENGTH}, LC={LC}, φ-averaged, "
        f"{N_REALIZATIONS} realisations per N; grey band = ±5%)",
        fontsize=12,
    )
    fig.savefig(RESULTS / "margin_convergence.png", dpi=180)
    plt.close(fig)
    print("  saved margin_convergence.png")


if __name__ == "__main__":
    run_analytical_checks()
    rows = run_sweep()
    save_csv(rows)
    # Single-resolution validation figures use the main resolution only.
    rows_main = [r for r in rows if r.N == N_MAIN]
    make_overview_figure(rows_main)
    make_main_figure(rows_main)
    make_ratio_figure(rows_main)
    # Convergence figure uses all N.
    make_convergence_figure(rows)
    print("Done. Results in experiments/04_margin_tessellation/results/")
