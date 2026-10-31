"""Experiment 02: Chi-squared random field (k=1) — Adler validation and percolation.

X(x) = f(x)^2  where f is a unit-variance Gaussian RF with covariance C(r)=exp(-r²/lc²).
The marginal is chi2(1); the hitting-set M_j follow from Worsley (1994).

Run:
    mamba run -n rf_explo python run.py

Outputs:
    results/mc_chi2.png           — 3-row x 5-col, LKC error bars vs theory + mid-plane slices
    results/percolation.png       — 3-panel: phi(kappa) / L0(kappa) / P(kappa)
    results/lc_XX/threshold_sweep.csv
    results/lc_XX/random_field.vti
    results/lc_XX/excursion.vti
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

BASE_RESULTS = Path(__file__).parent / "results"
BASE_RESULTS.mkdir(exist_ok=True)

SEED_BASE      = 42
N_REALIZATIONS = 10
N              = 100
LENGTH         = 100.0
K              = 1                          # chi2 dof — k=1: X = f², both Gaussian tails
LCS            = [4, 7, 12]                # correlation lengths to sweep (like exp01)
THRESHOLDS_RAW = np.concatenate([
    np.linspace(0.10, 3.0, 35),   # dense: EC crossings ~0.44–0.80, spanning ~1.6–1.8
    np.geomspace(4.0, 20.0, 5),   # sparse: trivial tail
])

LC_COLORS  = ["#4e79a7", "#f28e2b", "#59a14f"]
VIS_COLORS = ["#d62728", "#2ca02c", "#9467bd"]

C_THEORY = "#ff6e7c"
C_NUM    = "#3a8cb4"

COL_TITLES = [
    "L0 — Euler characteristic",
    "L1 — integrated mean curvature / π\n(binary-surface bias)",
    "L2 — half surface area",
    "L3 — volume fraction φ",
]
COL_KEYS = ["L0", "L1", "L2", "phi"]


def _chi2_crossings(lc: float, grid: GridSpec) -> tuple[tuple, tuple]:
    """Return Brent brackets for both E[L0]=0 crossings of chi2(k=1) at given lc.

    Scans log-spaced grid to locate sign changes, returns tight brackets
    (lo, hi) for the lower crossing (positive→negative) and upper crossing
    (negative→positive).  Either may be (nan, nan) if not found.
    """
    scan  = np.geomspace(1e-4, 50.0, 400)
    L0    = adler.expected_lkcs_chi2(grid, lc, K, scan)["L0"]
    signs = np.sign(L0)
    changes = np.where(np.diff(signs) != 0)[0]

    lower_bracket = (float("nan"), float("nan"))
    upper_bracket = (float("nan"), float("nan"))

    for idx in changes:
        lo, hi = float(scan[idx]), float(scan[idx + 1])
        if signs[idx] > 0:      # positive → negative: lower crossing
            lower_bracket = (lo, hi)
        else:                   # negative → positive: upper crossing
            upper_bracket = (lo, hi)

    return lower_bracket, upper_bracket


def _kappa_vis() -> list[float]:
    """Three thresholds for chi2(1) corresponding to φ ≈ 0.55, 0.35, 0.12."""
    from scipy import stats
    return [float(stats.chi2.ppf(1 - phi, df=K)) for phi in [0.55, 0.35, 0.12]]


def _style_ax(ax: plt.Axes) -> None:
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_xlabel("threshold κ", fontsize=9)


def run_one(lc: float) -> dict:
    results = BASE_RESULTS / f"lc_{int(lc):02d}"
    results.mkdir(exist_ok=True)

    grid      = GridSpec(n=N, length=LENGTH)
    theory    = adler.expected_lkcs_chi2(grid, lc, K, THRESHOLDS_RAW)
    kappa_vis = _kappa_vis()
    kappa_vti = kappa_vis[2]   # ~12% VF — near upper percolation threshold

    L0_all   = np.zeros((N_REALIZATIONS, len(THRESHOLDS_RAW)))
    L1_all   = np.zeros_like(L0_all)
    L2_all   = np.zeros_like(L0_all)
    phi_all  = np.zeros_like(L0_all)
    perc_all = np.zeros((N_REALIZATIONS, len(THRESHOLDS_RAW)), dtype=bool)
    field_slice = None

    print(f"  lc={lc}  running {N_REALIZATIONS} realisations...")
    for r in range(N_REALIZATIONS):
        rng   = np.random.default_rng(SEED_BASE + r)
        field = fields.generate_chi2_rf(grid, lc, K, rng)
        print(f"    realisation {r + 1}/{N_REALIZATIONS}", end="\r")
        for j, kappa in enumerate(THRESHOLDS_RAW):
            phase = field >= kappa
            m = measures.measure_all(phase, grid, compute_l1=True)
            L0_all[r, j]   = m["euler"]
            L1_all[r, j]   = m["l1"]
            L2_all[r, j]   = m["half_surface"]
            phi_all[r, j]  = m["vf"]
            perc_all[r, j] = m["percolates_any"]
        if r == 0:
            field_slice = field[:, :, N // 2]
            phase_vti   = field >= kappa_vti
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
        for j, kappa in enumerate(THRESHOLDS_RAW):
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
        "kappa_vis": kappa_vis,
        "L0": L0_all, "L1": L1_all, "L2": L2_all, "phi": phi_all,
        "perc": perc_all,
    }


def make_percolation_figure(all_data: list[dict], fname: str) -> None:
    grid = GridSpec(n=N, length=LENGTH)
    fig, (ax_phi, ax_L0, ax_p) = plt.subplots(
        3, 1, figsize=(9, 9), sharex=True, constrained_layout=True
    )

    for data, color in zip(all_data, LC_COLORS):
        lc     = data["lc"]
        phi_mu = data["phi"].mean(0)
        phi_sd = data["phi"].std(0, ddof=1)
        L0_mu  = data["L0"].mean(0)
        L0_sd  = data["L0"].std(0, ddof=1)
        perc_p = data["perc"].mean(0)
        t_phi  = data["theory"]["phi"]
        t_L0   = data["theory"]["L0"]

        lb, ub = _chi2_crossings(lc, grid)
        kp_lo  = adler.find_percolation_threshold(grid, lc, k_chi2=K, u_search=lb)[0] if np.isfinite(lb[0]) else float("nan")
        phi_lo = float(adler.chi2_m0(np.asarray([kp_lo]), K)[0]) if np.isfinite(kp_lo) else float("nan")
        kp_hi  = adler.find_percolation_threshold(grid, lc, k_chi2=K, u_search=ub)[0] if np.isfinite(ub[0]) else float("nan")
        phi_hi = float(adler.chi2_m0(np.asarray([kp_hi]), K)[0]) if np.isfinite(kp_hi) else float("nan")

        ax_phi.fill_between(THRESHOLDS_RAW, phi_mu - phi_sd, phi_mu + phi_sd,
                            color=color, alpha=0.15, lw=0)
        ax_phi.plot(THRESHOLDS_RAW, phi_mu, color=color, lw=1.5, alpha=0.8)
        ax_phi.plot(THRESHOLDS_RAW, t_phi,  color=color, lw=2.0, ls="--", label=f"lc = {lc}")

        ax_L0.fill_between(THRESHOLDS_RAW, L0_mu - L0_sd, L0_mu + L0_sd,
                           color=color, alpha=0.15, lw=0)
        ax_L0.plot(THRESHOLDS_RAW, L0_mu, color=color, lw=1.5, alpha=0.8)
        ax_L0.plot(THRESHOLDS_RAW, t_L0,  color=color, lw=2.0, ls="--", label=f"lc = {lc}")

        ax_p.plot(THRESHOLDS_RAW, perc_p, color=color, lw=2.0, marker="o",
                  ms=3.5, label=f"lc = {lc}  (n={N_REALIZATIONS})")

        for ax in (ax_phi, ax_L0, ax_p):
            if np.isfinite(kp_hi):
                ax.axvline(kp_hi, color=color, lw=1.2, ls="--", alpha=0.70)
            if np.isfinite(kp_lo):
                ax.axvline(kp_lo, color=color, lw=1.2, ls=":",  alpha=0.70)

        if np.isfinite(kp_hi) and np.isfinite(phi_hi):
            ax_phi.annotate(
                f"φ_c = {phi_hi:.3f}",
                xy=(kp_hi, phi_hi), xytext=(kp_hi + 0.3, phi_hi + 0.03),
                fontsize=7.5, color=color,
                arrowprops=dict(arrowstyle="-", color=color, lw=0.8),
            )
        if np.isfinite(kp_lo) and np.isfinite(phi_lo):
            ax_phi.annotate(
                f"1−φ_c = {phi_lo:.3f}",
                xy=(kp_lo, phi_lo), xytext=(kp_lo - 0.3, phi_lo - 0.07),
                fontsize=7.5, color=color, ha="right",
                arrowprops=dict(arrowstyle="-", color=color, lw=0.8),
            )

    ax_L0.axhline(0, color="0.3", lw=1.0, ls="--", label="E[L₀] = 0")
    ax_p.axhline(0.5, color="0.5", lw=0.8, ls=":", label="P = 0.5")

    ax_phi.set_ylabel("volume fraction φ", fontsize=10)
    ax_phi.set_title("L3 — volume fraction  (solid = MC mean ± std,  dashed = Adler theory)",
                     fontsize=10)
    ax_phi.legend(fontsize=8)
    ax_phi.grid(alpha=0.25, lw=0.5)

    ax_L0.set_yscale("symlog", linthresh=10)
    ax_L0.set_ylabel("Euler characteristic L0  (symlog)", fontsize=10)
    ax_L0.set_title("L0 — Euler characteristic\n"
                    "(-- = upper κ_c,  ··· = lower κ_c,  both where E[L₀] = 0)",
                    fontsize=10)
    ax_L0.legend(fontsize=8)
    ax_L0.grid(alpha=0.25, lw=0.5)

    ax_p.set_xlabel("threshold κ", fontsize=10)
    ax_p.set_ylabel("percolation probability", fontsize=10)
    ax_p.set_ylim(-0.05, 1.10)
    ax_p.set_title(f"Percolation probability  (empirical, {N_REALIZATIONS} realisations)",
                   fontsize=10)
    ax_p.legend(fontsize=8)
    ax_p.grid(alpha=0.25, lw=0.5)

    lc_str = "/".join(str(lc) for lc in LCS)
    fig.suptitle(
        f"Chi² RF (k=1) — percolation analysis  (lc={lc_str}, n={N}, {N_REALIZATIONS} realisations)",
        fontsize=13,
    )
    fig.savefig(BASE_RESULTS / fname, dpi=180)
    plt.close(fig)
    print(f"  saved results/{fname}")


def make_figure(all_data: list[dict], fname: str) -> None:
    n_rows   = len(all_data)
    kappa_vis = _kappa_vis()    # same for all rows (chi2(1) marginal)
    fig       = plt.figure(figsize=(21, 4 * n_rows))
    outer     = mgridspec.GridSpec(
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

            for kv, vc in zip(kappa_vis, VIS_COLORS):
                ax.axvline(kv, color=vc, lw=1.2, ls=":", alpha=0.85, zorder=0)

            ax.errorbar(THRESHOLDS_RAW, mu, yerr=sd,
                        fmt="o", color=C_NUM, ecolor=C_NUM, capsize=3, ms=3, lw=1.0,
                        label=f"Numerical (n={N_REALIZATIONS})" if first_row and col == 0 else None)
            ax.plot(THRESHOLDS_RAW, theory[key], color=C_THEORY, lw=2.0,
                    label="Adler theory (Worsley M_j)" if first_row and col == 0 else None)

            _style_ax(ax)
            if first_row:
                ax.set_title(title, fontsize=10)
            if col == 0:
                ax.set_ylabel(f"lc = {lc}", fontsize=11, fontweight="bold")
            if first_row and col == 0:
                ax.legend(fontsize=8, loc="best")

        inner     = mgridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=outer[row, 4], hspace=0.30)
        top_thumb = None
        for idx, (kv, vc) in enumerate(zip(kappa_vis, VIS_COLORS)):
            ax_img = fig.add_subplot(inner[idx])
            if idx == 0:
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
            ax_img.set_title(f"κ = {kv:.2f}", fontsize=8, color=vc, fontweight="bold", pad=3)

        if first_row and top_thumb is not None:
            pos = top_thumb.get_position()
            fig.text(
                (pos.x0 + pos.x1) / 2, pos.y1 + 0.012,
                "mid-plane slice  (realisation 1)",
                ha="center", va="bottom", fontsize=8, color="0.45",
            )

    lc_str = "/".join(str(lc) for lc in LCS)
    fig.suptitle(
        f"Chi² RF (k=1) — Monte Carlo validation  (lc={lc_str}, n={N}, {N_REALIZATIONS} realisations)",
        fontsize=13,
    )
    fig.savefig(BASE_RESULTS / fname, dpi=180)
    plt.close(fig)
    print(f"  saved results/{fname}")


def make_phi_c_vs_k_figure(fname: str) -> None:
    """Analytic φ_c⁺ as a function of chi² degrees of freedom k.

    Pure Adler calculation — no Monte Carlo. Shows the central qualitative
    claim: as k increases the marginal becomes less skewed, both EC
    zero-crossings tighten toward the Gaussian value, and the achievable
    "disconnected inclusion" VF ceiling drops.
    """
    from scipy import stats
    grid = GridSpec(n=N, length=LENGTH)
    k_scan = np.arange(1, 21)
    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)

    for lc, color in zip(LCS, LC_COLORS):
        phi_hi = np.full_like(k_scan, np.nan, dtype=float)
        for i, k in enumerate(k_scan):
            scan = np.geomspace(1e-4, 50.0, 400)
            L0   = adler.expected_lkcs_chi2(grid, lc, int(k), scan)["L0"]
            changes = np.where(np.diff(np.sign(L0)) != 0)[0]
            for idx in changes:
                if L0[idx] < 0 and L0[idx + 1] > 0:   # upper (negative → positive)
                    kp, pc = adler.find_percolation_threshold(
                        grid, lc, k_chi2=int(k),
                        u_search=(float(scan[idx]), float(scan[idx + 1])),
                    )
                    phi_hi[i] = pc
                    break
        ax.plot(k_scan, phi_hi, color=color, lw=2.0, marker="o", ms=4,
                label=f"lc = {lc}  (L/lc = {LENGTH/lc:.1f})")

    # Bulk (infinite-domain) reference: M3 = 0 at u_c = (2k-1 + sqrt(8k-7))/2.
    # phi_c+ = chi2(k).sf(u_c) — independent of lc, since the bulk term is the
    # only surviving one in the limit L/lc -> infinity.
    u_c_bulk   = (2.0 * k_scan - 1.0 + np.sqrt(8.0 * k_scan - 7.0)) / 2.0
    phi_c_bulk = np.array([float(stats.chi2.sf(u, df=int(k)))
                           for u, k in zip(u_c_bulk, k_scan)])
    ax.plot(k_scan, phi_c_bulk, color="black", lw=2.0, ls="--",
            label="bulk limit  (L/lc → ∞)")

    ax.axhline(0.159, color="0.4", lw=0.8, ls=":", label="Gaussian bulk (k→∞): 15.9%")
    ax.set_xlabel("chi² degrees of freedom  k", fontsize=10)
    ax.set_ylabel("φ_c⁺  (upper EC-zero volume fraction)", fontsize=10)
    ax.set_title(f"Chi² RF — disconnected-phase VF ceiling vs k  (analytic, N={N}, L={LENGTH:.0f})",
                 fontsize=11)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=9)
    fig.savefig(BASE_RESULTS / fname, dpi=180)
    plt.close(fig)
    print(f"  saved results/{fname}")


if __name__ == "__main__":
    all_data = []
    for lc in LCS:
        all_data.append(run_one(lc))
    make_figure(all_data, "mc_chi2.png")
    make_percolation_figure(all_data, "percolation.png")
    make_phi_c_vs_k_figure("phi_c_vs_k.png")
    print("Done.")
