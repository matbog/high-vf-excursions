"""Numerical pilot for coalescing Gaussian-space ball excursions."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from skimage.measure import euler_number

from margin_tessellation import generate_gaussian_field


@dataclass
class Result:
    realization: int
    radius: float
    radius_over_rc: float
    theory_phi_mc: float
    measured_phi: float
    euler: int
    components: int
    largest_fraction_phase: float
    spans_any: bool


def tetrahedral_centres(rho: float) -> np.ndarray:
    signs = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]])
    return rho * signs / np.sqrt(3.0)


def generate_vector_field(n: int, padding: int, lc: float, seed: int) -> np.ndarray:
    """Generate on a larger periodic grid, then take a bounded central crop."""
    rng = np.random.default_rng(seed)
    full_n = n + 2 * padding
    crop = slice(padding, padding + n)
    return np.stack([
        generate_gaussian_field(full_n, lc, rng)[crop, crop, crop]
        for _ in range(3)
    ])


def distance_to_lobes(fields: np.ndarray, centres: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    best = np.full(fields.shape[1:], np.inf, dtype=np.float32)
    nearest = np.zeros(fields.shape[1:], dtype=np.int8)
    for label, centre in enumerate(centres):
        d2 = sum((fields[j] - centre[j]) ** 2 for j in range(3))
        update = d2 < best
        best[update] = d2[update]
        nearest[update] = label
    return np.sqrt(best, out=best), nearest


def metrics(phase: np.ndarray) -> tuple[int, int, float, bool]:
    labels, count = ndimage.label(phase, ndimage.generate_binary_structure(3, 1))
    sizes = np.bincount(labels.ravel())
    occupied = int(phase.sum())
    largest = int(sizes[1:].max(initial=0))
    spans = False
    for axis in range(3):
        common = np.intersect1d(
            np.unique(np.take(labels, 0, axis=axis)),
            np.unique(np.take(labels, -1, axis=axis)),
        )
        spans |= bool(np.any(common != 0))
    return (int(euler_number(phase, connectivity=1)), int(count),
            largest / occupied if occupied else 0.0, spans)


def gaussian_volume_mc(centres: np.ndarray, radii: np.ndarray, samples: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    counts = np.zeros(len(radii), dtype=np.int64)
    for start in range(0, samples, 200_000):
        z = rng.normal(size=(min(200_000, samples - start), 3))
        d2 = np.min(np.sum((z[:, None] - centres[None]) ** 2, axis=2), axis=1)
        counts += np.sum(d2[:, None] <= radii[None] ** 2, axis=0)
    return counts / samples


def run(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    centres = tetrahedral_centres(args.rho)
    rc = float(np.linalg.norm(centres[0] - centres[1]) / 2)
    ratios = np.asarray(args.radius_over_rc)
    radii = rc * ratios
    theory = gaussian_volume_mc(centres, radii, args.mc_samples, args.seed + 100_000)
    distances, nearest_all, results = [], [], []

    for realization in range(args.realizations):
        fields = generate_vector_field(args.n, args.padding, args.lc, args.seed + realization)
        distance, nearest = distance_to_lobes(fields, centres)
        distances.append(distance)
        nearest_all.append(nearest)
        for radius, ratio, phi in zip(radii, ratios, theory):
            phase = distance <= radius
            ec, count, largest, spans = metrics(phase)
            results.append(Result(realization + 1, float(radius), float(ratio),
                                  float(phi), float(phase.mean()), ec, count,
                                  largest, spans))
        print(f"completed realization {realization + 1}/{args.realizations}")

    csv_path = args.output_dir / "coalescing_lobes_sweep.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=Result.__annotations__.keys())
        writer.writeheader()
        writer.writerows(asdict(result) for result in results)

    # Winner-coloured central slices. White is outside the union of balls.
    fig, axes = plt.subplots(args.realizations, len(radii), squeeze=False,
                             figsize=(2.7 * len(radii), 2.6 * args.realizations),
                             constrained_layout=True)
    cmap = plt.colormaps["tab10"].with_extremes(bad="white")
    mid = args.n // 2
    for row in range(args.realizations):
        for col, (radius, ratio) in enumerate(zip(radii, ratios)):
            phase = distances[row][:, :, mid] <= radius
            image = np.where(phase.T, nearest_all[row][:, :, mid].T, np.nan)
            ax = axes[row, col]
            ax.imshow(image, origin="lower", cmap=cmap, vmin=0, vmax=9,
                      interpolation="nearest")
            ax.set(xticks=[], yticks=[])
            if row == 0:
                ax.set_title(rf"$r/r_c={ratio:.2f}$" + "\n" +
                             rf"$\phi={theory[col]:.2f}$")
            ax.text(0.03, 0.04, rf"slice $\phi={phase.mean():.2f}$",
                    transform=ax.transAxes, fontsize=7,
                    bbox=dict(facecolor="white", alpha=.8, edgecolor="none"))
        axes[row, 0].set_ylabel(f"realization {row + 1}")
    fig.suptitle(rf"Tetrahedral Gaussian lobes, $N={args.n}^3$, $L_c={args.lc:g}$; "
                 r"hitting-set contact at $r/r_c=1$", fontsize=14)
    fig.savefig(args.output_dir / "coalescing_lobes_slices.png", dpi=220,
                bbox_inches="tight")
    plt.close(fig)

    groups = [[x for x in results if np.isclose(x.radius_over_rc, ratio)] for ratio in ratios]
    def mean(name: str) -> np.ndarray:
        return np.array([np.mean([getattr(x, name) for x in group]) for group in groups])
    def std(name: str) -> np.ndarray:
        return np.array([np.std([getattr(x, name) for x in group], ddof=1) for group in groups])

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8), constrained_layout=True)
    axes[0, 0].plot(ratios, theory, "k-", label="Gaussian-space MC")
    axes[0, 0].errorbar(ratios, mean("measured_phi"), yerr=std("measured_phi"),
                        fmt="o", label="3D RF")
    axes[0, 1].errorbar(ratios, mean("euler"), yerr=std("euler"), fmt="o-")
    axes[0, 1].axhline(0, color="k", lw=1)
    axes[1, 0].plot(ratios, mean("largest_fraction_phase"), "o-")
    axes[1, 1].plot(ratios, mean("spans_any"), "o-")
    titles = ["Volume calibration", "Euler characteristic",
              "Component coalescence", "Finite-window spanning"]
    ylabels = ["volume fraction", "Euler characteristic",
               "largest component / phase", "spanning frequency"]
    for ax, title, ylabel in zip(axes.flat, titles, ylabels):
        ax.axvline(1, color="crimson", ls="--", label=r"value-space contact $r_c$")
        ax.set(xlabel=r"$r/r_c$", ylabel=ylabel, title=title)
        ax.grid(alpha=.3)
    axes[1, 1].set_ylim(-.05, 1.05)
    axes[0, 0].legend(fontsize=8)
    axes[1, 1].legend(fontsize=8)
    fig.savefig(args.output_dir / "coalescing_lobes_diagnostics.png", dpi=220,
                bbox_inches="tight")
    plt.close(fig)

    (args.output_dir / "README.md").write_text(
        f"# Coalescing-lobe pilot\n\n"
        f"rho = {args.rho}; critical radius r_c = {rc:.8f}.\n\n"
        "The dashed line marks geometric contact in Gaussian value space. It is "
        "not imposed as a physical-space percolation threshold.\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--lc", type=float, default=10.)
    parser.add_argument("--padding", type=int, default=20)
    parser.add_argument("--rho", type=float, default=1.5)
    parser.add_argument("--realizations", type=int, default=4)
    parser.add_argument("--radius-over-rc", type=float, nargs="+",
                        default=[.45, .65, .85, .98, 1.02, 1.15, 1.35, 1.60])
    parser.add_argument("--mc-samples", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=9417)
    parser.add_argument("--output-dir", type=Path,
                        default=Path("results/coalescing_lobes"))
    args = parser.parse_args()
    if not 2 <= args.realizations <= 5:
        raise ValueError("realizations must be between 2 and 5")
    if not 0 < args.lc < 15:
        raise ValueError("lc must be strictly below 15 voxels")
    run(args)


if __name__ == "__main__":
    main()
