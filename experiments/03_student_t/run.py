"""Experiment 03: Student-t random fields — Adler validation and percolation.

T(x) = Z_0(x) / sqrt(chi2_field(x) / nu): polynomial tails raise the EC-zero
disconnected-VF threshold above the Gaussian ~16% ceiling. Adler formula remains
exact (Worsley 1994). Sweeps nu in {inf (Gaussian), 5, 3}.

Run:
    mamba run -n rf_explo python run.py

Outputs (results/):
    mc_student_t.png           — 3-row x 5-col, LKC error bars vs theory + slices
    percolation.png            — 2-panel: phi(kappa) / L0(kappa) with EC=0 crossings
    percolation_diagram.png    — theory-only phi_c+ vs nu
    <label>/threshold_sweep.csv
    <label>/random_field.vti
    <label>/excursion.vti
    <label>/excursion_surface.vtp
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.gridspec as mgridspec
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats, optimize

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, fields, io, measures
from src.fields import GridSpec

BASE_RESULTS = Path(__file__).parent / "results"
BASE_RESULTS.mkdir(exist_ok=True)

SEED_BASE      = 42
N_REALIZATIONS = 10
N              = 100
LENGTH         = 100.0
LC             = 8.0

# Three cases for the MC figure: Gaussian + two t_nu
MC_CASES = [
    (None, "gaussian"),   # nu = inf
    (5,    "t_5"),
    (3,    "t_3"),
]
# Wider nu range for the theory-only percolation diagram
NU_DIAGRAM: list[int | None] = (
    [None]
    + sorted(set(int(round(n)) for n in np.logspace(np.log10(3.0), np.log10(200.0), 30)),
             reverse=True)
)

THRESHOLDS = np.linspace(-3.0, 3.0, 35)
KAPPA_VIS  = [-1.5, 0.0, 1.5]
KAPPA_VTI  = 0.5

CASE_COLORS = ["#4e79a7", "#f28e2b", "#59a14f"]
VIS_COLORS  = ["#d62728", "#2ca02c", "#9467bd"]
C_THEORY    = "#ff6e7c"
C_NUM       = "#3a8cb4"

COL_TITLES = [
    "L0 — Euler characteristic",
    "L1 — integrated mean curvature / π",
    "L2 — half surface area",
    "L3 — volume fraction φ",
]
COL_KEYS = ["L0", "L1", "L2", "phi"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _theory(grid: GridSpec, nu: int | None, kappas: np.ndarray) -> dict:
    if nu is None:
        return adler.expected_lkcs_gaussian(grid, LC, kappas)
    return adler.expected_lkcs_t(grid, LC, nu, kappas)


def _phi_at(nu: int | None, kappa: float) -> float:
    if nu is None:
        return float(stats.norm.sf(kappa))
    return float(stats.t.sf(kappa, df=nu))


def _ec_zero_bracket(
    grid: GridSpec, nu: int | None, lo: float = 0.01, hi: float = 5.0, n: int = 400
) -> tuple[float, float]:
    """Scan for the upper E[L0]=0 bracket (negative -> positive as u increases)."""
    scan = np.linspace(lo, hi, n)
    L0 = _theory(grid, nu, scan)["L0"]
    for i in range(len(L0) - 1):
        if L0[i] <= 0 < L0[i + 1]:
            return float(scan[i]), float(scan[i + 1])
    return float("nan"), float("nan")


def _style_ax(ax: plt.Axes) -> None:
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_xlabel("threshold κ", fontsize=9)


# ---------------------------------------------------------------------------
# Monte Carlo sweep (one nu / label)
# ---------------------------------------------------------------------------

def run_one(nu: int | None, label: str) -> dict:
    results = BASE_RESULTS / label
    results.mkdir(exist_ok=True)

    grid   = GridSpec(n=N, length=LENGTH)
    theory = _theory(grid, nu, THRESHOLDS)

    L0_all   = np.zeros((N_REALIZATIONS, len(THRESHOLDS)))
    L1_all   = np.zeros_like(L0_all)
    L2_all   = np.zeros_like(L0_all)
    phi_all  = np.zeros_like(L0_all)
    perc_all = np.zeros((N_REALIZATIONS, len(THRESHOLDS)), dtype=bool)
    field_slice = None

    print(f"  {label}  running {N_REALIZATIONS} realisations...")
    for r in range(N_REALIZATIONS):
        rng   = np.random.default_rng(SEED_BASE + r)
        field = (fields.generate_gaussian_rf(grid, LC, rng) if nu is None
                 else fields.generate_student_t_rf(grid, LC, nu, rng))
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

    print(f"  {label}  CSV + VTK saved to {results}/")
    return {
        "nu": nu, "label": label, "theory": theory, "field_slice": field_slice,
        "L0": L0_all, "L1": L1_all, "L2": L2_all, "phi": phi_all, "perc": perc_all,
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figure(all_data: list[dict], fname: str) -> None:
    """3-row × 5-col: L0 / L1 / L2 / phi error bars vs theory + mid-plane slices."""
    n_rows = len(all_data)
    fig    = plt.figure(figsize=(21, 4 * n_rows))
    outer  = mgridspec.GridSpec(
        n_rows, 5, figure=fig,
        width_ratios=[1, 1, 1, 1, 1.1],
        left=0.06, right=0.985, top=0.91, bottom=0.08,
        wspace=0.38, hspace=0.42,
    )

    for row, data in enumerate(all_data):
        label       = data["label"]
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
                        label=f"Numerical (n={N_REALIZATIONS})" if first_row and col == 0 else None)
            ax.plot(THRESHOLDS, theory[key], color=C_THEORY, lw=2.0,
                    label="Adler theory (Worsley M_j)" if first_row and col == 0 else None)

            _style_ax(ax)
            if first_row:
                ax.set_title(title, fontsize=10)
            if col == 0:
                ax.set_ylabel(label, fontsize=11, fontweight="bold")
            if first_row and col == 0:
                ax.legend(fontsize=8, loc="best")

        inner     = mgridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[row, 4], hspace=0.30)
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
            ax_img.set_title(f"κ = {kv:+.1f}", fontsize=8, color=vc, fontweight="bold", pad=3)

        if first_row and top_thumb is not None:
            pos = top_thumb.get_position()
            fig.text(
                (pos.x0 + pos.x1) / 2, pos.y1 + 0.012,
                "mid-plane slice  (realisation 1)",
                ha="center", va="bottom", fontsize=8, color="0.45",
            )

    cases_str = " / ".join(d["label"] for d in all_data)
    fig.suptitle(
        f"Student-t RF — Monte Carlo validation  ({cases_str},  n={N},  {N_REALIZATIONS} realisations)",
        fontsize=12,
    )
    fig.savefig(BASE_RESULTS / fname, dpi=180)
    plt.close(fig)
    print(f"  saved results/{fname}")


def make_percolation_figure(all_data: list[dict], fname: str) -> None:
    """2-panel: phi(κ) and L0(κ) with upper/lower EC=0 crossings marked."""
    grid = GridSpec(n=N, length=LENGTH)
    fig, (ax_phi, ax_L0) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, constrained_layout=True
    )

    for data, color in zip(all_data, CASE_COLORS):
        nu     = data["nu"]
        label  = data["label"]
        phi_mu = data["phi"].mean(0)
        phi_sd = data["phi"].std(0, ddof=1)
        L0_mu  = data["L0"].mean(0)
        L0_sd  = data["L0"].std(0, ddof=1)
        t_phi  = data["theory"]["phi"]
        t_L0   = data["theory"]["L0"]

        # Find upper EC=0 crossing (positive kappa side)
        bracket = _ec_zero_bracket(grid, nu)
        if np.isfinite(bracket[0]):
            kp_hi, phi_hi = adler.find_percolation_threshold(
                grid, LC, nu=nu, u_search=bracket)
        else:
            kp_hi, phi_hi = float("nan"), float("nan")

        # Lower crossing: symmetric for t_nu and Gaussian
        if np.isfinite(kp_hi):
            kp_lo  = -kp_hi
            phi_lo = 1.0 - phi_hi
        else:
            kp_lo, phi_lo = float("nan"), float("nan")

        ax_phi.fill_between(THRESHOLDS, phi_mu - phi_sd, phi_mu + phi_sd,
                            color=color, alpha=0.15, lw=0)
        ax_phi.plot(THRESHOLDS, phi_mu, color=color, lw=1.5, alpha=0.8)
        ax_phi.plot(THRESHOLDS, t_phi,  color=color, lw=2.0, ls="--", label=label)

        ax_L0.fill_between(THRESHOLDS, L0_mu - L0_sd, L0_mu + L0_sd,
                           color=color, alpha=0.15, lw=0)
        ax_L0.plot(THRESHOLDS, L0_mu, color=color, lw=1.5, alpha=0.8)
        ax_L0.plot(THRESHOLDS, t_L0,  color=color, lw=2.0, ls="--", label=label)

        for ax in (ax_phi, ax_L0):
            if np.isfinite(kp_hi):
                ax.axvline(kp_hi, color=color, lw=1.2, ls="--", alpha=0.70)
            if np.isfinite(kp_lo):
                ax.axvline(kp_lo, color=color, lw=1.2, ls=":",  alpha=0.70)

        if np.isfinite(kp_hi) and np.isfinite(phi_hi):
            ax_phi.annotate(
                f"φ_c = {phi_hi:.3f}",
                xy=(kp_hi, phi_hi), xytext=(kp_hi + 0.15, phi_hi + 0.03),
                fontsize=7.5, color=color,
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
    ax_L0.set_title("L0 — Euler characteristic  (-- = upper κ_c,  ··· = lower κ_c)", fontsize=10)
    ax_L0.legend(fontsize=8)
    ax_L0.grid(alpha=0.25, lw=0.5)

    cases_str = " / ".join(d["label"] for d in all_data)
    fig.suptitle(
        f"Student-t RF — percolation analysis  ({cases_str},  n={N},  {N_REALIZATIONS} realisations)",
        fontsize=12,
    )
    fig.savefig(BASE_RESULTS / fname, dpi=180)
    plt.close(fig)
    print(f"  saved results/{fname}")


def run_percolation_diagram() -> None:
    """Theory-only phi_c+ vs nu (EC=0 criterion), saved as percolation_diagram.png/.csv."""
    print("Percolation diagram phi_c+ vs nu (theory only)...")
    grid = GridSpec(n=N, length=LENGTH)

    kappa_list, phi_list, nu_plot = [], [], []
    rows = []
    for nu in NU_DIAGRAM:
        bracket = _ec_zero_bracket(grid, nu)
        if not np.isfinite(bracket[0]):
            print(f"  nu={'inf' if nu is None else nu}: no upper crossing found — skip")
            continue
        kc, phi_c = adler.find_percolation_threshold(grid, LC, nu=nu, u_search=bracket)
        lbl = "inf" if nu is None else str(nu)
        print(f"  nu={lbl:>6}: kappa_c+ = {kc:.4f}, phi_c+ = {phi_c:.4f}")
        kappa_list.append(kc)
        phi_list.append(phi_c)
        nu_plot.append(1000.0 if nu is None else float(nu))
        rows.append([lbl, kc, phi_c])

    nu_arr  = np.array(nu_plot)
    phi_arr = np.array(phi_list)
    kappa_arr = np.array(kappa_list)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    ax1.semilogx(nu_arr, phi_arr * 100, "o-", lw=2)
    ax1.axhline(phi_arr[0] * 100, color="gray", ls="--",
                label=f"Gaussian: {phi_arr[0]*100:.1f}%")
    ax1.set_xlabel("degrees of freedom ν")
    ax1.set_ylabel("φ_c⁺  [%]  (EC-zero criterion)")
    ax1.set_title("EC-zero VF threshold vs ν  (not the spanning threshold)")
    ax1.set_xticks([3, 5, 10, 30, 100, 1000])
    ax1.set_xticklabels(["3", "5", "10", "30", "100", "∞"])
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.semilogx(nu_arr, kappa_arr, "s-", lw=2, color="C1")
    ax2.set_xlabel("degrees of freedom ν")
    ax2.set_ylabel("κ_c⁺")
    ax2.set_title("EC-zero threshold κ_c⁺ vs ν")
    ax2.set_xticks([3, 5, 10, 30, 100, 1000])
    ax2.set_xticklabels(["3", "5", "10", "30", "100", "∞"])
    ax2.grid(alpha=0.3)

    fig.suptitle(f"Student-t RF — Adler EC criterion  (n={N}, lc={LC})", fontsize=13)
    fig.savefig(BASE_RESULTS / "percolation_diagram.png", dpi=180)
    plt.close(fig)

    with (BASE_RESULTS / "percolation_diagram.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nu", "kappa_c_plus", "phi_c_plus"])
        w.writerows(rows)
    print("  saved results/percolation_diagram.png/.csv")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    all_data = []
    for nu, label in MC_CASES:
        all_data.append(run_one(nu, label))
    make_figure(all_data, "mc_student_t.png")
    make_percolation_figure(all_data, "percolation.png")
    run_percolation_diagram()
    print("Done. Results in experiments/03_student_t/results/")
