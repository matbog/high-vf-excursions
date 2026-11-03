"""Exp 04 ceiling push — spanning phase diagram + resolution convergence.

Complements the main validation sweep (run.py). This script does *not*
compare to GKF theory; it pushes q and phi higher to answer two
materials-goal questions:

  Phase A — spanning phase diagram at N = N_A (default 128)
    For each (q, phi) cell, measure the fraction of realizations in
    which any single winner-type cluster spans the domain face-to-face
    on any axis. Draw a 2D map of "disconnected" vs "spanning" regimes
    with the matrix-film thickness delta/dx (voxel units) as a side
    diagnostic.

  Phase B — resolution convergence at high q
    Fix a handful of (q, phi) points and sweep N in {64, 128, 256} to
    see how the spanning fraction and film thickness scale with
    resolution. Tests whether high-q, high-phi results survive
    refinement.

Primary output:
  results/ceiling_phase_a.csv       per-realization measurements
  results/ceiling_phase_a.png       heatmap of spanning + film contours
  results/ceiling_phase_a_vf.png    measured VF vs target per (q)
  results/ceiling_phase_b.csv       convergence table
  results/ceiling_phase_b.png       spanning + film vs N per (q, phi)

Run (overnight-safe, ~3-4 h expected):
  mamba run -n rf_explo python run_ceiling.py            # both phases
  mamba run -n rf_explo python run_ceiling.py --phase a  # only A
  mamba run -n rf_explo python run_ceiling.py --phase b  # only B
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parents[2]))
from src import adler, excursions, fields, measures  # noqa: E402
from src.fields import GridSpec  # noqa: E402

RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

LENGTH = 64.0
LC     = 5.0
SEED   = 1729

# Phase A: spanning phase diagram (N fixed)
Q_A        = [2, 4, 8, 16, 32, 64]
PHI_A      = [0.70, 0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98]
N_A        = 128
N_REAL_A   = 10

# Phase B: resolution convergence
Q_B        = [16, 32]
PHI_B      = [0.85, 0.90, 0.95]
N_B        = [64, 128, 256]
N_REAL_B   = 5


# ---------------------------------------------------------------------------
# Core per-realization measurement
# ---------------------------------------------------------------------------

def measure_cell(fs: np.ndarray, delta: float, grid: GridSpec) -> dict:
    """One realization measurement: phase, winners, spanning, film thickness."""
    phase, winners = excursions.margin_excursion(fs, delta)
    vf = float(phase.mean())
    spans = measures.measure_percolation_per_winner(phase, winners, fs.shape[0])
    return {
        "vf":              vf,
        "spans_any":       bool(spans),
        "delta_physical":  float(delta),
        "delta_voxels":    float(delta / grid.dx),
    }


def sweep_one_realization(grid: GridSpec, q: int, phi_values: list[float],
                           rng: np.random.Generator) -> list[dict]:
    """Generate q fields once, apply each phi target."""
    fs = fields.generate_margin_fields(grid, LC, q, rng)
    out = []
    for phi in phi_values:
        delta = adler.delta_for_vf(q, phi)
        m = measure_cell(fs, delta, grid)
        m["q"]          = q
        m["N"]          = grid.n
        m["phi_target"] = phi
        out.append(m)
    return out


# ---------------------------------------------------------------------------
# Phase A — spanning phase diagram at fixed N
# ---------------------------------------------------------------------------

def run_phase_a() -> list[dict]:
    print(f"=== Phase A: spanning phase diagram, N={N_A} ===")
    grid = GridSpec(n=N_A, length=LENGTH)
    all_rows = []
    t0_total = time.time()
    for q in Q_A:
        t_q = time.time()
        print(f"  q={q}, {N_REAL_A} realizations × {len(PHI_A)} phi targets")
        for r in range(N_REAL_A):
            rng = np.random.default_rng(SEED + 10000 * q + r)
            rows = sweep_one_realization(grid, q, PHI_A, rng)
            for row in rows:
                row["realization"] = r
                all_rows.append(row)
            # Checkpoint CSV after each realization
            write_phase_a_csv(all_rows)
            if r == 0:
                dt = time.time() - t_q
                print(f"    first realization: {dt:.1f}s  (ETA for q={q}: "
                      f"{dt * N_REAL_A:.0f}s)")
        print(f"  q={q} done in {time.time() - t_q:.1f}s")
    print(f"Phase A total: {time.time() - t0_total:.1f}s")
    return all_rows


def write_phase_a_csv(rows: list[dict]) -> None:
    with (RESULTS / "ceiling_phase_a.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["q", "N", "phi_target", "realization", "vf",
                    "spans_any", "delta_physical", "delta_voxels"])
        for row in rows:
            w.writerow([row["q"], row["N"], row["phi_target"],
                        row["realization"], row["vf"],
                        int(row["spans_any"]),
                        row["delta_physical"], row["delta_voxels"]])


def summarize_phase_a(rows: list[dict]) -> dict:
    """Aggregate per (q, phi): mean VF, spanning fraction, mean film thickness."""
    out = {}
    for q in Q_A:
        out[q] = {}
        for phi in PHI_A:
            cell = [r for r in rows if r["q"] == q and r["phi_target"] == phi]
            if not cell:
                continue
            out[q][phi] = {
                "vf_mean":     float(np.mean([r["vf"] for r in cell])),
                "vf_std":      float(np.std([r["vf"] for r in cell], ddof=1)),
                "spans_frac":  float(np.mean([r["spans_any"] for r in cell])),
                "delta_vox":   float(np.mean([r["delta_voxels"] for r in cell])),
                "n_real":      len(cell),
            }
    return out


def plot_phase_a(summary: dict) -> None:
    q_values = sorted(summary.keys())
    phi_values = sorted(set(phi for q in q_values for phi in summary[q].keys()))

    spans = np.full((len(q_values), len(phi_values)), np.nan)
    films = np.full_like(spans, np.nan)
    for i, q in enumerate(q_values):
        for j, phi in enumerate(phi_values):
            cell = summary[q].get(phi)
            if cell:
                spans[i, j] = cell["spans_frac"]
                films[i, j] = cell["delta_vox"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

    # Panel 1: spanning fraction heatmap
    im1 = ax1.pcolormesh(phi_values, q_values, spans, cmap="RdYlGn_r",
                          vmin=0, vmax=1, shading="auto")
    cbar1 = fig.colorbar(im1, ax=ax1, shrink=0.9)
    cbar1.set_label("fraction of realizations with any per-winner spanning",
                     fontsize=10)
    ax1.set_xlabel("target φ", fontsize=11)
    ax1.set_ylabel("q (competing fields)", fontsize=11)
    ax1.set_yscale("log")
    ax1.set_yticks(q_values)
    ax1.set_yticklabels([str(q) for q in q_values])
    ax1.set_title(f"Spanning phase diagram  (N={N_A}, {N_REAL_A} realizations/cell)",
                  fontsize=12)
    # overlay numeric values
    for i, q in enumerate(q_values):
        for j, phi in enumerate(phi_values):
            if not np.isnan(spans[i, j]):
                ax1.text(phi, q, f"{spans[i, j]:.1f}",
                         ha="center", va="center", fontsize=8, color="black")

    # Panel 2: film thickness in voxels
    im2 = ax2.pcolormesh(phi_values, q_values, films, cmap="viridis",
                          norm=plt.matplotlib.colors.LogNorm(vmin=0.1, vmax=10),
                          shading="auto")
    cbar2 = fig.colorbar(im2, ax=ax2, shrink=0.9)
    cbar2.set_label("matrix film thickness δ / dx  (voxel units)", fontsize=10)
    # 1-voxel contour
    ax2.contour(phi_values, q_values, films, levels=[0.5, 1.0, 2.0],
                colors="white", linewidths=1.2)
    ax2.set_xlabel("target φ", fontsize=11)
    ax2.set_ylabel("q", fontsize=11)
    ax2.set_yscale("log")
    ax2.set_yticks(q_values)
    ax2.set_yticklabels([str(q) for q in q_values])
    ax2.set_title("Matrix film thickness (contours at 0.5, 1, 2 voxels)",
                  fontsize=12)
    for i, q in enumerate(q_values):
        for j, phi in enumerate(phi_values):
            if not np.isnan(films[i, j]):
                ax2.text(phi, q, f"{films[i, j]:.2f}",
                         ha="center", va="center", fontsize=8, color="white"
                         if films[i, j] < 1 else "black")

    fig.suptitle(f"Exp 04 ceiling push — Phase A  "
                 f"(L={LENGTH}, lc={LC}, L/lc={LENGTH/LC:.1f})",
                 fontsize=13)
    fig.savefig(RESULTS / "ceiling_phase_a.png", dpi=150)
    plt.close(fig)
    print("  saved ceiling_phase_a.png")

    # Secondary plot: VF realisation accuracy
    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for i, q in enumerate(q_values):
        xs = [phi for phi in phi_values if phi in summary[q]]
        ys = [summary[q][phi]["vf_mean"] for phi in xs]
        es = [summary[q][phi]["vf_std"]  for phi in xs]
        ax.errorbar(xs, ys, yerr=es, marker="o", lw=1.5, ms=5,
                    capsize=3, color=cmap(i / max(1, len(q_values) - 1)),
                    label=f"q={q}")
    phi_line = np.linspace(min(phi_values), max(phi_values), 50)
    ax.plot(phi_line, phi_line, "k--", lw=1, alpha=0.6, label="y = x")
    ax.set_xlabel("target φ", fontsize=11)
    ax.set_ylabel("measured φ (mean ± std)", fontsize=11)
    ax.set_title(f"VF calibration across q (N={N_A})", fontsize=12)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.savefig(RESULTS / "ceiling_phase_a_vf.png", dpi=150)
    plt.close(fig)
    print("  saved ceiling_phase_a_vf.png")


# ---------------------------------------------------------------------------
# Phase B — resolution convergence at high q
# ---------------------------------------------------------------------------

def run_phase_b() -> list[dict]:
    print(f"=== Phase B: resolution convergence at high q ===")
    all_rows = []
    t0_total = time.time()
    for N in N_B:
        grid = GridSpec(n=N, length=LENGTH)
        for q in Q_B:
            t_q = time.time()
            print(f"  N={N}, q={q}, {N_REAL_B} realizations × "
                  f"{len(PHI_B)} phi targets")
            for r in range(N_REAL_B):
                rng = np.random.default_rng(SEED + 10000 * q + r + 100 * N)
                rows = sweep_one_realization(grid, q, PHI_B, rng)
                for row in rows:
                    row["realization"] = r
                    all_rows.append(row)
                write_phase_b_csv(all_rows)
            print(f"    q={q}, N={N}: {time.time() - t_q:.1f}s")
    print(f"Phase B total: {time.time() - t0_total:.1f}s")
    return all_rows


def write_phase_b_csv(rows: list[dict]) -> None:
    with (RESULTS / "ceiling_phase_b.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["q", "N", "phi_target", "realization", "vf",
                    "spans_any", "delta_physical", "delta_voxels"])
        for row in rows:
            w.writerow([row["q"], row["N"], row["phi_target"],
                        row["realization"], row["vf"],
                        int(row["spans_any"]),
                        row["delta_physical"], row["delta_voxels"]])


def plot_phase_b(rows: list[dict]) -> None:
    fig, axes = plt.subplots(2, len(Q_B), figsize=(6 * len(Q_B), 10),
                              sharex=True, constrained_layout=True)
    if len(Q_B) == 1:
        axes = axes[:, None]

    cmap = plt.get_cmap("plasma")

    for col, q in enumerate(Q_B):
        ax_span = axes[0, col]
        ax_film = axes[1, col]
        for j, phi in enumerate(PHI_B):
            spans_vs_N = []
            films_vs_N = []
            for N in N_B:
                cell = [r for r in rows if r["q"] == q and r["N"] == N
                        and r["phi_target"] == phi]
                spans_vs_N.append(np.mean([c["spans_any"] for c in cell])
                                   if cell else np.nan)
                films_vs_N.append(np.mean([c["delta_voxels"] for c in cell])
                                   if cell else np.nan)
            c = cmap(j / max(1, len(PHI_B) - 1))
            ax_span.plot(N_B, spans_vs_N, "o-", color=c, lw=1.8, ms=6,
                          label=f"φ = {phi}")
            ax_film.plot(N_B, films_vs_N, "o-", color=c, lw=1.8, ms=6,
                          label=f"φ = {phi}")

        ax_span.set_xscale("log")
        ax_span.set_xticks(N_B)
        ax_span.set_xticklabels([str(N) for N in N_B])
        ax_span.set_ylim(-0.05, 1.05)
        ax_span.set_ylabel("spanning fraction")
        ax_span.set_title(f"q = {q}", fontsize=11)
        ax_span.grid(alpha=0.3, which="both")
        ax_span.legend(fontsize=9)

        ax_film.set_xscale("log")
        ax_film.set_xticks(N_B)
        ax_film.set_xticklabels([str(N) for N in N_B])
        ax_film.axhline(1, color="gray", lw=0.8, ls="--", alpha=0.7)
        ax_film.text(N_B[-1] * 0.98, 1.05, "1 voxel",
                      ha="right", fontsize=8, color="gray")
        ax_film.set_ylabel("film thickness δ/dx (voxels)")
        ax_film.set_xlabel("N (grid)")
        ax_film.grid(alpha=0.3, which="both")

    fig.suptitle(f"Exp 04 Phase B — resolution convergence at high q  "
                 f"(L={LENGTH}, lc={LC}, {N_REAL_B} realizations)",
                 fontsize=13)
    fig.savefig(RESULTS / "ceiling_phase_b.png", dpi=150)
    plt.close(fig)
    print("  saved ceiling_phase_b.png")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["a", "b", "both"], default="both")
    args = parser.parse_args()

    if args.phase in ("a", "both"):
        rows_a = run_phase_a()
        summary_a = summarize_phase_a(rows_a)
        plot_phase_a(summary_a)

    if args.phase in ("b", "both"):
        rows_b = run_phase_b()
        plot_phase_b(rows_b)

    print("Done.")
