"""Experiment 04: Coalescing Gaussian-space lobes.

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

Key finding: physical percolation occurred before r/r_c = 1 (i.e., before the
value-space lobes touched), with a large negative Euler characteristic dip
near r_c suggesting a tunnel-rich transition.

Run:
    mamba run -n rf_explo python run.py

Outputs (all in results/):
    coalescing_lobes_sweep.csv      — all realisations, all radii
    coalescing_lobes_slices.png     — winner-coloured central slices
    coalescing_lobes_diagnostics.png — 4-panel diagnostic

Seeds: SEED=9417.
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

SEED = 9417
N = 100
LENGTH = 100.0
LC = 10.0
PADDING = 20
RHO = 1.5
N_REALIZATIONS = 4
RADIUS_OVER_RC = [0.45, 0.65, 0.85, 0.98, 1.02, 1.15, 1.35, 1.60]
MC_SAMPLES = 2_000_000


@dataclass
class Row:
    realization: int
    radius: float
    radius_over_rc: float
    theory_phi: float
    measured_phi: float
    euler: int
    largest_fraction_phase: float
    spans_any: bool


def run() -> None:
    grid = GridSpec(n=N, length=LENGTH)
    centres = excursions.tetrahedral_centres(RHO)
    rc = float(np.linalg.norm(centres[0] - centres[1]) / 2.0)
    ratios = np.array(RADIUS_OVER_RC)
    radii = rc * ratios

    print(f"Coalescing lobes: n={N}, lc={LC}, rho={RHO}, r_c={rc:.6f}")
    print(f"r/r_c sweep: {RADIUS_OVER_RC}")
    print(f"NOTE: After r/r_c > 1 the analytic M_j require ball-intersection")
    print(f"      inclusion-exclusion (not yet derived).\n")

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

        for radius, ratio, phi in zip(radii, ratios, theory_phi):
            phase = distances <= radius
            m = measures.measure_all(phase, grid)
            occupied = int(phase.sum())
            from skimage import measure as sk_measure
            from scipy import ndimage
            labels, _ = ndimage.label(phase, ndimage.generate_binary_structure(3, 1))
            sizes = np.bincount(labels.ravel())
            largest = int(sizes[1:].max(initial=0)) if occupied > 0 else 0
            rows.append(Row(
                realization=realization + 1,
                radius=float(radius),
                radius_over_rc=float(ratio),
                theory_phi=float(phi),
                measured_phi=float(phase.mean()),
                euler=int(m["euler"]),
                largest_fraction_phase=largest / occupied if occupied else 0.0,
                spans_any=m["percolates_any"],
            ))
            print(f"  r={realization+1} r/r_c={ratio:.2f} phi={phase.mean():.3f} "
                  f"euler={m['euler']:6.0f} spans={m['percolates_any']}")

    with (RESULTS / "coalescing_lobes_sweep.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(Row.__annotations__.keys()))
        w.writeheader()
        w.writerows(asdict(r) for r in rows)
    print("\n  saved coalescing_lobes_sweep.csv")

    _plot_slices(distances_all, nearest_all, radii, ratios, theory_phi)
    _plot_diagnostics(rows, ratios, theory_phi)


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
    fig.suptitle(rf"Tetrahedral Gaussian lobes (contact at $r/r_c=1$, dashed)",
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
            np.std([getattr(r, attr) for r in rows if np.isclose(r.radius_over_rc, ratio)], ddof=1)
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

    axes[0, 1].errorbar(ratios, mean_by_ratio("euler"), yerr=std_by_ratio("euler"),
                         fmt="o-", ms=4)
    axes[0, 1].axhline(0, color="k", lw=0.8)
    axes[0, 1].set_ylabel("Euler characteristic")
    axes[0, 1].set_title("Euler characteristic\n(no analytic prediction after contact)")

    axes[1, 0].plot(ratios, mean_by_ratio("largest_fraction_phase"), "o-", ms=4)
    axes[1, 0].set_ylabel("largest component / phase")
    axes[1, 0].set_title("Component coalescence")

    axes[1, 1].plot(ratios, mean_by_ratio("spans_any"), "o-", ms=4)
    axes[1, 1].set_ylim(-0.05, 1.05)
    axes[1, 1].set_ylabel("spanning frequency")
    axes[1, 1].set_title("Physical percolation\n(occurs before value-space contact)")

    fig.suptitle(rf"Coalescing lobes  ($N={N}^3$, $L_c={LC}$, {N_REALIZATIONS} realisations)")
    fig.savefig(RESULTS / "coalescing_lobes_diagnostics.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  saved coalescing_lobes_diagnostics.png")


if __name__ == "__main__":
    run()
    print("Done. Results in experiments/04_coalescing_lobes/results/")
