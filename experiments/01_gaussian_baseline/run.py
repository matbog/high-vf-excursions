"""Experiment 01: Gaussian random field baseline with full Adler validation.

Demonstrates the complete Adler-Taylor machinery for a scalar Gaussian RF
with Gaussian covariance. The Adler formula gives exact closed-form predictions
for E[L0], E[L2], E[L3] as functions of threshold kappa. This experiment
validates those predictions against single-realisation and Monte Carlo measurements.

Run:
    mamba run -n rf_explo python run.py

Outputs (all in results/):
    threshold_sweep.png/.csv  — single realisation vs Adler theory
    mc_gaussian.png/.csv      — Monte Carlo mean +/- std vs Adler theory
    random_field.vti          — ParaView scalar field
    excursion.vti             — ParaView binary excursion (kappa=1.0)
    excursion_surface.vtp     — ParaView surface mesh

Seeds: SEED=42 (single realisation), SEED_MC=10 (Monte Carlo).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, fields, io, measures
from src.fields import GridSpec

RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

SEED = 42
SEED_MC = 10
N_RUNS_MC = 50
N = 100
LENGTH = 100.0
LC = 8.0
KAPPA_VTI = 1.0


def run_single_sweep() -> None:
    print("Single-realisation threshold sweep...")
    grid = GridSpec(n=N, length=LENGTH)
    rng = np.random.default_rng(SEED)
    field = fields.generate_gaussian_rf(grid, LC, rng)

    thresholds = np.linspace(-2.5, 2.5, 25)
    rows = []
    m_L0, m_L2, m_L3, m_phi = [], [], [], []

    for kappa in thresholds:
        phase = field >= kappa
        m = measures.measure_all(phase, grid)
        m_L0.append(m["euler"])
        m_L2.append(m["half_surface"])
        m_L3.append(m["volume"])
        m_phi.append(m["vf"])

    theory = adler.expected_lkcs_gaussian(grid, LC, thresholds)

    with (RESULTS / "threshold_sweep.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kappa", "measured_L0", "measured_L2", "measured_L3", "measured_phi",
                    "theory_L0", "theory_L2", "theory_L3", "theory_phi"])
        for i, kappa in enumerate(thresholds):
            w.writerow([kappa, m_L0[i], m_L2[i], m_L3[i], m_phi[i],
                        theory["L0"][i], theory["L2"][i], theory["L3"][i], theory["phi"][i]])

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
    panels = [
        (axes[0, 0], "Volume fraction", m_phi, theory["phi"]),
        (axes[0, 1], "L3 volume", m_L3, theory["L3"]),
        (axes[1, 0], "L2 = half surface area", m_L2, theory["L2"]),
        (axes[1, 1], "L0 = Euler characteristic", m_L0, theory["L0"]),
    ]
    for ax, title, measured, expected in panels:
        ax.plot(thresholds, measured, "o", ms=4, label="measured")
        ax.plot(thresholds, expected, "-", lw=2, label="Adler theory")
        ax.set_xlabel("threshold κ")
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle(f"Gaussian RF — single realisation  (n={N}, lc={LC}, seed={SEED})")
    fig.savefig(RESULTS / "threshold_sweep.png", dpi=180)
    plt.close(fig)
    print(f"  saved threshold_sweep.png/.csv")

    phase = field >= KAPPA_VTI
    io.save_vti(RESULTS / "random_field.vti", field, grid, "random_field")
    io.save_vti(RESULTS / "excursion.vti", phase.astype(np.uint8), grid, "excursion")
    io.save_vtp(RESULTS / "excursion_surface.vtp", phase, grid)
    print(f"  saved VTI/VTP files (kappa={KAPPA_VTI})")


def run_monte_carlo() -> None:
    print(f"Monte Carlo ({N_RUNS_MC} runs)...")
    grid = GridSpec(n=N, length=LENGTH)
    thresholds = np.linspace(-2.0, 2.0, 17)
    L0_all = np.zeros((N_RUNS_MC, len(thresholds)))
    L2_all = np.zeros_like(L0_all)
    phi_all = np.zeros_like(L0_all)

    for r in range(N_RUNS_MC):
        rng = np.random.default_rng(SEED_MC + r)
        field = fields.generate_gaussian_rf(grid, LC, rng)
        for j, kappa in enumerate(thresholds):
            m = measures.measure_all(field >= kappa, grid)
            L0_all[r, j] = m["euler"]
            L2_all[r, j] = m["half_surface"]
            phi_all[r, j] = m["vf"]
        print(f"  run {r + 1}/{N_RUNS_MC}", end="\r")

    print()
    theory = adler.expected_lkcs_gaussian(grid, LC, thresholds)

    with (RESULTS / "mc_gaussian.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["kappa", "L0_mean", "L0_std", "L2_mean", "L2_std", "phi_mean", "phi_std",
                    "theory_L0", "theory_L2", "theory_phi"])
        for j, kappa in enumerate(thresholds):
            w.writerow([kappa,
                        L0_all[:, j].mean(), L0_all[:, j].std(ddof=1),
                        L2_all[:, j].mean(), L2_all[:, j].std(ddof=1),
                        phi_all[:, j].mean(), phi_all[:, j].std(ddof=1),
                        theory["L0"][j], theory["L2"][j], theory["phi"][j]])

    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True)
    for ax, key, label, t_key in [
        (axes[0], L0_all, "L0 Euler characteristic", "L0"),
        (axes[1], L2_all, "L2 half surface area", "L2"),
        (axes[2], phi_all, "Volume fraction", "phi"),
    ]:
        ax.errorbar(thresholds, key.mean(0), yerr=key.std(0, ddof=1),
                    fmt="o", capsize=3, ms=4, label=f"MC mean ± std (n={N_RUNS_MC})")
        ax.plot(thresholds, theory[t_key], "-", lw=2, label="Adler theory")
        ax.set_xlabel("threshold κ")
        ax.set_title(label)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle(f"Gaussian RF — Monte Carlo validation  (n={N}, lc={LC})")
    fig.savefig(RESULTS / "mc_gaussian.png", dpi=180)
    plt.close(fig)
    print("  saved mc_gaussian.png/.csv")


if __name__ == "__main__":
    run_single_sweep()
    run_monte_carlo()
    print("Done. Results in experiments/01_gaussian_baseline/results/")
