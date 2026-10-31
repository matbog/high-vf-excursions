"""Experiment 01: Gaussian random field baseline with full Adler validation.

Runs for lc in LCS.  VTK files go to results/lc_XX/; combined figures
go to results/.

Run:
    mamba run -n rf_explo python run.py

Outputs:
    results/mc_gaussian.png         — 3-row x 5-col, LKC error bars vs theory + mid-plane slices
    results/percolation.png         — percolation probability P(kappa) + phi(kappa) per lc
    results/lc_XX/threshold_sweep.csv
    results/lc_XX/random_field.vti
    results/lc_XX/excursion.vti     — binary excursion at KAPPA_VTI
    results/lc_XX/excursion_surface.vtp
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.gridspec as mgridspec
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, fields, io, measures
from src.fields import GridSpec

BASE_RESULTS  = Path(__file__).parent / "results"
BASE_RESULTS.mkdir(exist_ok=True)

SEED_BASE      = 42
N_REALIZATIONS = 10
N              = 100
LENGTH         = 100.0
LCS            = [4, 7, 12]
KAPPA_VTI      = 0.5           # κ ≈ 31% VF — disconnected inclusions
THRESHOLDS     = np.linspace(-2.5, 2.5, 25)

KAPPA_VIS  = [-1.5, 0.0, 1.5]
VIS_COLORS = ["#d62728", "#2ca02c", "#9467bd"]
LC_COLORS  = ["#4e79a7", "#f28e2b", "#59a14f"]   # one per lc in LCS order

C_THEORY = "#ff6e7c"
C_NUM    = "#3a8cb4"

COL_TITLES = [
    "L0 — Euler characteristic",
    "L1 — integrated mean curvature / π\n(binary-surface bias)",
    "L2 — half surface area",
    "L3 — volume fraction φ",
]
COL_KEYS = ["L0", "L1", "L2", "phi"]


def _style_ax(ax: plt.Axes) -> None:
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_xlabel("threshold κ", fontsize=9)


def run_one(lc: float) -> dict:
    """Run N_REALIZATIONS for one lc; return dict of arrays, theory and field slice."""
    results = BASE_RESULTS / f"lc_{lc:02.0f}"
    results.mkdir(exist_ok=True)

    grid   = GridSpec(n=N, length=LENGTH)
    theory = adler.expected_lkcs_gaussian(grid, lc, THRESHOLDS)

    L0_all   = np.zeros((N_REALIZATIONS, len(THRESHOLDS)))
    L1_all   = np.zeros_like(L0_all)
    L2_all   = np.zeros_like(L0_all)
    phi_all  = np.zeros_like(L0_all)
    perc_all = np.zeros((N_REALIZATIONS, len(THRESHOLDS)), dtype=bool)
    field_slice = None

    print(f"  lc={lc}  running {N_REALIZATIONS} realisations...")
    for r in range(N_REALIZATIONS):
        rng   = np.random.default_rng(SEED_BASE + r)
        field = fields.generate_gaussian_rf(grid, lc, rng)
        print(f"    realisation {r + 1}/{N_REALIZATIONS}", end="\r")
        for j, kappa in enumerate(THRESHOLDS):
            phase = field >= kappa
            m = measures.measure_all(phase, grid, compute_l1=True)
            L0_all[r, j]   = m["euler"]
            L1_all[r, j]   = m["l1"]
            L2_all[r, j]   = m["half_surface"]
            phi_all[r, j]  = m["vf"]
            perc_all[r, j] = m["percolates_any"]
        if r == 0:
            field_slice = field[:, :, N // 2]
            phase_vti   = field >= KAPPA_VTI
            io.save_vti(results / "random_field.vti",      field,                      grid, "random_field")
            io.save_vti(results / "excursion.vti",         phase_vti.astype(np.uint8), grid, "excursion")
            io.save_vtp(results / "excursion_surface.vtp", phase_vti, grid)
    print()

    with (results / "threshold_sweep.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "kappa",
            "L0_mean", "L0_std", "L1_mean", "L1_std",
            "L2_mean", "L2_std", "phi_mean", "phi_std",
            "perc_prob",
            "theory_L0", "theory_L1", "theory_L2", "theory_phi",
        ])
        for j, kappa in enumerate(THRESHOLDS):
            w.writerow([
                kappa,
                L0_all[:, j].mean(),  L0_all[:, j].std(ddof=1),
                L1_all[:, j].mean(),  L1_all[:, j].std(ddof=1),
                L2_all[:, j].mean(),  L2_all[:, j].std(ddof=1),
                phi_all[:, j].mean(), phi_all[:, j].std(ddof=1),
                perc_all[:, j].mean(),
                theory["L0"][j], theory["L1"][j], theory["L2"][j], theory["phi"][j],
            ])

    print(f"  lc={lc}  VTK + CSV saved to {results}/")
    return {
        "lc": lc, "theory": theory, "field_slice": field_slice,
        "L0": L0_all, "L1": L1_all, "L2": L2_all, "phi": phi_all,
        "perc": perc_all,
    }


def _interp_half_crossing(x: np.ndarray, p: np.ndarray) -> float | None:
    """Return x where percolation probability p drops through 0.5 (linear interp)."""
    for i in range(len(p) - 1):
        if p[i] >= 0.5 >= p[i + 1]:
            t = (p[i] - 0.5) / (p[i] - p[i + 1])
            return float(x[i] + t * (x[i + 1] - x[i]))
    return None


def make_percolation_figure(all_data: list[dict], fname: str) -> None:
    """φ (top) and L0 (bottom) vs κ with shared x-axis.

    Vertical dashed lines at κ_c (upper E[L0]=0 crossing, Roubin & Colliat 2016)
    span both panels; φ_c is annotated on the top panel.
    """
    grid = GridSpec(n=N, length=LENGTH)
    fig, (ax_phi, ax_L0) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, constrained_layout=True
    )

    for data, color in zip(all_data, LC_COLORS):
        lc     = data["lc"]
        phi_mu = data["phi"].mean(0)
        phi_sd = data["phi"].std(0, ddof=1)
        L0_mu  = data["L0"].mean(0)
        L0_sd  = data["L0"].std(0, ddof=1)
        t_phi  = data["theory"]["phi"]
        t_L0   = data["theory"]["L0"]

        kp_hi, phi_hi = adler.find_percolation_threshold(grid, lc, u_search=(0.5,  5.0))
        kp_lo, phi_lo = adler.find_percolation_threshold(grid, lc, u_search=(-2.5, -0.1))

        # Volume fraction (top)
        ax_phi.fill_between(THRESHOLDS, phi_mu - phi_sd, phi_mu + phi_sd,
                            color=color, alpha=0.15, lw=0)
        ax_phi.plot(THRESHOLDS, phi_mu, color=color, lw=1.5, alpha=0.8)
        ax_phi.plot(THRESHOLDS, t_phi,  color=color, lw=2.0, ls="--",
                    label=f"lc = {lc}")

        # Euler characteristic (bottom)
        ax_L0.fill_between(THRESHOLDS, L0_mu - L0_sd, L0_mu + L0_sd,
                           color=color, alpha=0.15, lw=0)
        ax_L0.plot(THRESHOLDS, L0_mu, color=color, lw=1.5, alpha=0.8)
        ax_L0.plot(THRESHOLDS, t_L0,  color=color, lw=2.0, ls="--",
                   label=f"lc = {lc}")

        # Upper crossing: κ_c+, φ_c (disconnected inclusions onset)
        if np.isfinite(kp_hi):
            ax_phi.axvline(kp_hi, color=color, lw=1.2, ls="--", alpha=0.75)
            ax_L0.axvline(kp_hi, color=color, lw=1.2, ls="--", alpha=0.75)
            if np.isfinite(phi_hi):
                ax_phi.annotate(
                    f"φ_c = {phi_hi:.3f}",
                    xy=(kp_hi, phi_hi), xytext=(kp_hi + 0.12, phi_hi + 0.03),
                    fontsize=7.5, color=color,
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.8),
                )

        # Lower crossing: κ_c−, φ = 1−φ_c (matrix becomes disconnected)
        if np.isfinite(kp_lo):
            ax_phi.axvline(kp_lo, color=color, lw=1.2, ls=":", alpha=0.75)
            ax_L0.axvline(kp_lo, color=color, lw=1.2, ls=":", alpha=0.75)
            if np.isfinite(phi_lo):
                ax_phi.annotate(
                    f"1−φ_c = {phi_lo:.3f}",
                    xy=(kp_lo, phi_lo), xytext=(kp_lo - 0.12, phi_lo - 0.07),
                    fontsize=7.5, color=color, ha="right",
                    arrowprops=dict(arrowstyle="-", color=color, lw=0.8),
                )

    ax_L0.axhline(0, color="0.3", lw=1.0, ls="--", label="E[L₀] = 0")

    ax_phi.set_ylabel("volume fraction φ", fontsize=10)
    ax_phi.set_title("L3 — volume fraction  (solid = MC mean ± std,  dashed = Adler theory)",
                     fontsize=10)
    ax_phi.legend(fontsize=8)
    ax_phi.grid(alpha=0.25, lw=0.5)

    ax_L0.set_xlabel("threshold κ", fontsize=10)
    ax_L0.set_ylabel("Euler characteristic L0", fontsize=10)
    ax_L0.set_title("L0 — Euler characteristic\n"
                    "(-- = upper κ_c,  ··· = lower κ_c,  both where E[L₀] = 0)",
                    fontsize=10)
    ax_L0.legend(fontsize=8)
    ax_L0.grid(alpha=0.25, lw=0.5)

    fig.suptitle(
        f"Gaussian RF — volume fraction & Euler characteristic  "
        f"(n={N}, {N_REALIZATIONS} realisations)",
        fontsize=13,
    )
    fig.savefig(BASE_RESULTS / fname, dpi=180)
    plt.close(fig)
    print(f"  saved results/{fname}")


def make_figure(all_data: list[dict], fname: str) -> None:
    n_rows = len(all_data)
    fig    = plt.figure(figsize=(21, 4 * n_rows))
    outer  = mgridspec.GridSpec(
        n_rows, 5, figure=fig,
        width_ratios=[1, 1, 1, 1, 1.1],
        left=0.06, right=0.985, top=0.91, bottom=0.08,
        wspace=0.38, hspace=0.42,
    )

    for row, data in enumerate(all_data):
        lc          = data["lc"]
        theory      = data["theory"]
        arrs        = [data["L0"], data["L1"], data["L2"], data["phi"]]
        field_slice = data["field_slice"]
        first_row   = (row == 0)

        for col, (arr, key, title) in enumerate(zip(arrs, COL_KEYS, COL_TITLES)):
            ax = fig.add_subplot(outer[row, col])
            mu, sd = arr.mean(0), arr.std(0, ddof=1)

            for kv, vc in zip(KAPPA_VIS, VIS_COLORS):
                ax.axvline(kv, color=vc, lw=1.2, ls=":", alpha=0.85, zorder=0)

            ax.errorbar(THRESHOLDS, mu, yerr=sd,
                        fmt="o", color=C_NUM, ecolor=C_NUM, capsize=3, ms=3, lw=1.0,
                        label=f"Numerical experiment (n={N_REALIZATIONS})" if first_row and col == 0 else None)
            ax.plot(THRESHOLDS, theory[key], color=C_THEORY, lw=2.0,
                    label="Minkowski Functionals" if first_row and col == 0 else None)

            _style_ax(ax)
            if first_row:
                ax.set_title(title, fontsize=10)
            if col == 0:
                ax.set_ylabel(f"lc = {lc}", fontsize=11, fontweight="bold")
            if first_row and col == 0:
                ax.legend(fontsize=8, loc="best")

        # 5th column: 3 stacked mid-plane slice thumbnails
        inner = mgridspec.GridSpecFromSubplotSpec(
            3, 1, subplot_spec=outer[row, 4], hspace=0.30
        )
        top_thumb = None
        for k, (kv, vc) in enumerate(zip(KAPPA_VIS, VIS_COLORS)):
            ax_img = fig.add_subplot(inner[k])
            if k == 0:
                top_thumb = ax_img
            ax_img.imshow(
                (field_slice >= kv), cmap="gray_r", origin="lower",
                interpolation="nearest", vmin=0, vmax=1,
            )
            ax_img.set_xticks([])
            ax_img.set_yticks([])
            for spine in ax_img.spines.values():
                spine.set_edgecolor(vc)
                spine.set_linewidth(2.0)
            ax_img.set_title(f"κ = {kv:+.1f}", fontsize=8, color=vc,
                             fontweight="bold", pad=3)

        if first_row and top_thumb is not None:
            pos = top_thumb.get_position()
            fig.text(
                (pos.x0 + pos.x1) / 2, pos.y1 + 0.012,
                "mid-plane slice  (realisation 1)",
                ha="center", va="bottom", fontsize=8, color="0.45",
            )

    fig.suptitle(
        f"Gaussian RF — Monte Carlo validation  (n={N}, {N_REALIZATIONS} realisations)",
        fontsize=13,
    )
    fig.savefig(BASE_RESULTS / fname, dpi=180)
    plt.close(fig)
    print(f"  saved results/{fname}")


if __name__ == "__main__":
    all_data = []
    for lc in LCS:
        all_data.append(run_one(lc))
    make_figure(all_data, "mc_gaussian.png")
    make_percolation_figure(all_data, "percolation.png")
    print("Done.")
