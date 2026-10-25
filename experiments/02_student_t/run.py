"""Experiment 02: Student-t random fields — heavier tails, higher VF ceiling.

The Worsley t_nu field T(x) = Z_0(x) / sqrt(chi2_field(x) / nu) has polynomial
tails, which pushes the disconnected volume fraction ceiling above the Gaussian ~16%.
The Adler formula remains exact (Worsley 1994), with t_nu M_j replacing Gaussian M_j.
Only M_j changes; the spatial covariance structure (and thus lambda_2) is identical.

Run:
    mamba run -n rf_explo python run.py

Outputs (all in results/):
    student_t_sweep.png/.csv    — threshold sweep for Gaussian and t_nu, measured vs theory
    mc_student_t.png/.csv       — Monte Carlo comparison for Gaussian, t_5, t_3
    percolation_diagram.png     — phi_c+ vs nu (max disconnected VF from Adler EC criterion)

Seeds: SEED=42 (single realisation), MC seeds per case.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t as t_dist

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, fields, measures
from src.fields import GridSpec

RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

SEED = 42
N = 80
LENGTH = 100.0
LC = 8.0
NU_VALUES = [None, 30, 10, 5, 3]   # None = Gaussian
MC_CASES = [(None, "Gaussian"), (5, "t_5"), (3, "t_3")]
N_RUNS_MC = 20


def run_single_sweep() -> None:
    print("Single-realisation threshold sweep for Gaussian and t_nu...")
    grid = GridSpec(n=N, length=LENGTH)
    thresholds = np.linspace(-3.0, 3.0, 25)

    fig, axes = plt.subplots(2, len(NU_VALUES), figsize=(4 * len(NU_VALUES), 8),
                             constrained_layout=True)

    rows = []
    for col, nu in enumerate(NU_VALUES):
        label = "Gaussian" if nu is None else f"t_{nu}"
        rng = np.random.default_rng(SEED)
        field = (fields.generate_gaussian_rf(grid, LC, rng) if nu is None
                 else fields.generate_student_t_rf(grid, LC, nu, rng))
        print(f"  {label}...", end=" ")
        L0_measured, phi_measured = [], []
        for kappa in thresholds:
            m = measures.measure_all(field >= kappa, grid)
            L0_measured.append(m["euler"])
            phi_measured.append(m["vf"])

        theory = (adler.expected_lkcs_gaussian(grid, LC, thresholds) if nu is None
                  else adler.expected_lkcs_t(grid, LC, nu, thresholds))

        rows.extend(
            [label, float(kappa), float(L0_measured[i]), float(phi_measured[i]),
             float(theory["L0"][i]), float(theory["phi"][i])]
            for i, kappa in enumerate(thresholds)
        )

        ax_ec = axes[0, col]
        ax_phi = axes[1, col]
        ax_ec.plot(thresholds, L0_measured, "o", ms=3, label="measured")
        ax_ec.plot(thresholds, theory["L0"], "-", lw=2, label="Adler")
        ax_ec.axhline(0, color="k", lw=0.8, ls="--")
        ax_ec.set_title(label)
        ax_ec.set_xlabel("κ")
        ax_ec.set_ylabel("L0" if col == 0 else "")
        ax_ec.legend(fontsize=7)
        ax_ec.grid(alpha=0.3)
        ax_phi.plot(thresholds, phi_measured, "o", ms=3, label="measured")
        ax_phi.plot(thresholds, theory["phi"], "-", lw=2, label="theory")
        ax_phi.set_title(label)
        ax_phi.set_xlabel("κ")
        ax_phi.set_ylabel("Φ" if col == 0 else "")
        ax_phi.legend(fontsize=7)
        ax_phi.grid(alpha=0.3)
        print("done")

    axes[0, 0].set_ylabel("L0 Euler characteristic")
    axes[1, 0].set_ylabel("Volume fraction Φ")
    fig.suptitle(f"Gaussian vs t_ν threshold sweep  (n={N}, lc={LC}, seed={SEED})")
    fig.savefig(RESULTS / "student_t_sweep.png", dpi=180)
    plt.close(fig)

    with (RESULTS / "student_t_sweep.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["field", "kappa", "L0_measured", "phi_measured", "L0_theory", "phi_theory"])
        w.writerows(rows)
    print("  saved student_t_sweep.png/.csv")


def run_monte_carlo() -> None:
    print(f"Monte Carlo ({N_RUNS_MC} runs per case)...")
    grid = GridSpec(n=N, length=LENGTH)
    thresholds = np.linspace(-2.5, 2.5, 17)

    fig, axes = plt.subplots(2, len(MC_CASES), figsize=(5 * len(MC_CASES), 8),
                             constrained_layout=True)
    rows = []

    for col, (nu, label) in enumerate(MC_CASES):
        print(f"  {label}...")
        L0_all = np.zeros((N_RUNS_MC, len(thresholds)))
        phi_all = np.zeros_like(L0_all)
        for r in range(N_RUNS_MC):
            rng = np.random.default_rng(col * 1000 + r)
            field = (fields.generate_gaussian_rf(grid, LC, rng) if nu is None
                     else fields.generate_student_t_rf(grid, LC, nu, rng))
            for j, kappa in enumerate(thresholds):
                m = measures.measure_all(field >= kappa, grid)
                L0_all[r, j] = m["euler"]
                phi_all[r, j] = m["vf"]
            print(f"    run {r + 1}/{N_RUNS_MC}", end="\r")
        print()

        theory = (adler.expected_lkcs_gaussian(grid, LC, thresholds) if nu is None
                  else adler.expected_lkcs_t(grid, LC, nu, thresholds))

        for j, kappa in enumerate(thresholds):
            rows.append([label, float(kappa),
                         float(L0_all[:, j].mean()), float(L0_all[:, j].std(ddof=1)),
                         float(phi_all[:, j].mean()), float(phi_all[:, j].std(ddof=1)),
                         float(theory["L0"][j]), float(theory["phi"][j])])

        ax_ec = axes[0, col]
        ax_phi = axes[1, col]
        ax_ec.errorbar(thresholds, L0_all.mean(0), yerr=L0_all.std(0, ddof=1),
                       fmt="o", capsize=3, ms=3, label=f"MC mean±std (n={N_RUNS_MC})")
        ax_ec.plot(thresholds, theory["L0"], "-", lw=2, label="Adler")
        ax_ec.axhline(0, color="k", lw=0.8, ls="--")
        ax_ec.set_title(label)
        ax_ec.set_xlabel("κ")
        ax_ec.set_ylabel("L0" if col == 0 else "")
        ax_ec.legend(fontsize=7)
        ax_ec.grid(alpha=0.3)
        ax_phi.errorbar(thresholds, phi_all.mean(0), yerr=phi_all.std(0, ddof=1),
                        fmt="o", capsize=3, ms=3)
        ax_phi.plot(thresholds, theory["phi"], "-", lw=2, label="theory")
        ax_phi.set_title(label)
        ax_phi.set_xlabel("κ")
        ax_phi.set_ylabel("Φ" if col == 0 else "")
        ax_phi.legend(fontsize=7)
        ax_phi.grid(alpha=0.3)

    fig.suptitle(f"Student-t MC validation  (n={N}, lc={LC}, {N_RUNS_MC} runs)")
    fig.savefig(RESULTS / "mc_student_t.png", dpi=180)
    plt.close(fig)

    with (RESULTS / "mc_student_t.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["field", "kappa", "L0_mean", "L0_std", "phi_mean", "phi_std",
                    "L0_theory", "phi_theory"])
        w.writerows(rows)
    print("  saved mc_student_t.png/.csv")


def run_percolation_diagram() -> None:
    print("Percolation diagram phi_c+ vs nu...")
    grid = GridSpec(n=N, length=LENGTH)
    nu_float = np.logspace(np.log10(3.0), np.log10(200.0), 30)[::-1]
    nu_vals: list[int | None] = [None] + sorted(set(int(round(n)) for n in nu_float), reverse=True)

    kappa_list, phi_list, nu_plot = [], [], []
    rows = []
    for nu in nu_vals:
        kc, phi_c = adler.find_percolation_threshold(grid, LC, nu)
        lbl = "inf" if nu is None else str(nu)
        print(f"  nu={lbl:>6}: kappa_c+ = {kc:.4f}, phi_c+ = {phi_c:.4f}")
        kappa_list.append(kc)
        phi_list.append(phi_c)
        nu_plot.append(1000.0 if nu is None else float(nu))
        rows.append([lbl, kc, phi_c])

    nu_arr = np.array(nu_plot)
    phi_arr = np.array(phi_list)
    kappa_arr = np.array(kappa_list)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    ax1.semilogx(nu_arr, phi_arr * 100, "o-", lw=2)
    ax1.axhline(phi_arr[0] * 100, color="gray", ls="--",
                label=f"Gaussian: {phi_arr[0]*100:.1f}%")
    ax1.set_xlabel("degrees of freedom ν")
    ax1.set_ylabel("Φ_c⁺  [%]")
    ax1.set_title("Max disconnected VF ceiling vs ν")
    ax1.set_xticks([3, 5, 10, 30, 100, 1000])
    ax1.set_xticklabels(["3", "5", "10", "30", "100", "∞"])
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax2.semilogx(nu_arr, kappa_arr, "s-", lw=2, color="C1")
    ax2.set_xlabel("degrees of freedom ν")
    ax2.set_ylabel("κ_c⁺")
    ax2.set_title("Percolation threshold κ_c⁺ vs ν")
    ax2.set_xticks([3, 5, 10, 30, 100, 1000])
    ax2.set_xticklabels(["3", "5", "10", "30", "100", "∞"])
    ax2.grid(alpha=0.3)
    fig.suptitle(f"Adler EC criterion  (n={N}, lc={LC})")
    fig.savefig(RESULTS / "percolation_diagram.png", dpi=180)
    plt.close(fig)

    with (RESULTS / "percolation_diagram.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nu", "kappa_c_plus", "phi_c_plus"])
        w.writerows(rows)
    print("  saved percolation_diagram.png/.csv")


if __name__ == "__main__":
    run_single_sweep()
    run_monte_carlo()
    run_percolation_diagram()
    print("Done. Results in experiments/02_student_t/results/")
