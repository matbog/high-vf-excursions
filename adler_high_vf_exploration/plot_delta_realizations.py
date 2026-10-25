"""Plot bounded 3D margin-tessellation slices for increasing delta.

The Gaussian fields are synthesized on a larger FFT grid and centrally cropped,
so the reported 100^3 observation window does not identify opposite faces.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from margin_tessellation import gaussian_margin_volume_fraction, generate_gaussian_field


def generate_margin_and_winners(
    q: int,
    observation_n: int,
    padding: int,
    lc: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate q fields while retaining only their two largest values."""
    synthesis_n = observation_n + 2 * padding
    crop = slice(padding, padding + observation_n)
    rng = np.random.default_rng(seed)
    largest = np.full((observation_n,) * 3, -np.inf, dtype=np.float32)
    second = np.full_like(largest, -np.inf)
    winners = np.zeros(largest.shape, dtype=np.int16)

    for label in range(q):
        full = generate_gaussian_field(synthesis_n, lc, rng)
        field = full[crop, crop, crop]
        beats_largest = field > largest
        second = np.where(beats_largest, largest, np.maximum(second, field))
        largest = np.where(beats_largest, field, largest)
        winners[beats_largest] = label

    return largest - second, winners


def create_plots(
    output_dir: Path,
    q: int,
    n: int,
    lc: float,
    realizations: int,
    deltas: list[float],
    seed: int,
    padding: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    theory_phi = [gaussian_margin_volume_fraction(q, delta) for delta in deltas]
    margins: list[np.ndarray] = []
    winners_all: list[np.ndarray] = []
    measured: list[list[float]] = []

    for realization in range(realizations):
        margin, winners = generate_margin_and_winners(
            q, n, padding, lc, seed + realization
        )
        margins.append(margin)
        winners_all.append(winners)
        measured.append([float(np.mean(margin >= delta)) for delta in deltas])
        print(f"generated realization {realization + 1}/{realizations}")

    extent = (0, n, 0, n)
    fig_binary, axes_binary = plt.subplots(
        realizations, len(deltas), figsize=(3.0 * len(deltas), 2.85 * realizations),
        squeeze=False, constrained_layout=True,
    )
    fig_labels, axes_labels = plt.subplots(
        realizations, len(deltas), figsize=(3.0 * len(deltas), 2.85 * realizations),
        squeeze=False, constrained_layout=True,
    )
    label_cmap = plt.colormaps["tab10"].with_extremes(bad="white")
    z_index = n // 2

    for row in range(realizations):
        for column, delta in enumerate(deltas):
            phase = margins[row] >= delta
            phase_slice = phase[:, :, z_index].T
            labels = np.where(
                phase_slice,
                winners_all[row][:, :, z_index].T,
                np.nan,
            )
            axes_binary[row, column].imshow(
                phase_slice, origin="lower", cmap="gray_r", interpolation="nearest",
                extent=extent, vmin=0, vmax=1,
            )
            axes_labels[row, column].imshow(
                labels, origin="lower", cmap=label_cmap, interpolation="nearest",
                extent=extent, vmin=0, vmax=max(9, q - 1),
            )
            for axes in (axes_binary, axes_labels):
                ax = axes[row, column]
                ax.set_xticks([])
                ax.set_yticks([])
                if row == 0:
                    ax.set_title(
                        rf"$\delta={delta:.2f}$" + "\n" +
                        rf"theory $\phi={theory_phi[column]:.3f}$"
                    )
                ax.text(
                    0.03, 0.04, rf"$\hat{{\phi}}={measured[row][column]:.3f}$",
                    transform=ax.transAxes, fontsize=8,
                    bbox={"facecolor": "white", "alpha": 0.78, "edgecolor": "none", "pad": 1.5},
                )
        axes_binary[row, 0].set_ylabel(f"realization {row + 1}")
        axes_labels[row, 0].set_ylabel(f"realization {row + 1}")

    common_title = (
        rf"Central slices of bounded margin fields: $q={q}$, "
        rf"$N={n}^3$, $L_c={lc:g}$ voxels"
    )
    fig_binary.suptitle(common_title + " — binary retained phase", fontsize=14)
    fig_labels.suptitle(common_title + " — winner-labelled phase", fontsize=14)
    fig_binary.savefig(output_dir / "delta_slices_binary.png", dpi=220, bbox_inches="tight")
    fig_labels.savefig(output_dir / "delta_slices_winner_labels.png", dpi=220, bbox_inches="tight")
    plt.close(fig_binary)
    plt.close(fig_labels)

    with (output_dir / "delta_slices_metadata.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["realization", "seed", "q", "n", "lc", "delta", "theory_phi", "measured_phi"])
        for row in range(realizations):
            for column, delta in enumerate(deltas):
                writer.writerow([
                    row + 1, seed + row, q, n, lc, delta,
                    theory_phi[column], measured[row][column],
                ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q", type=int, default=8)
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--lc", type=float, default=10.0)
    parser.add_argument("--realizations", type=int, default=4)
    parser.add_argument("--delta", type=float, nargs="+", default=[0.05, 0.20, 0.40, 0.70, 1.00])
    parser.add_argument("--seed", type=int, default=3107)
    parser.add_argument("--padding", type=int, default=20)
    parser.add_argument("--output-dir", type=Path, default=Path("results/delta_slices_q8"))
    args = parser.parse_args()
    if not 2 <= args.realizations <= 5:
        raise ValueError("realizations must be between 2 and 5")
    if not 0 < args.lc < 15:
        raise ValueError("lc must be strictly between 0 and 15 voxels")
    if args.padding < int(np.ceil(args.lc)):
        raise ValueError("padding should be at least one correlation length")
    create_plots(
        args.output_dir, args.q, args.n, args.lc, args.realizations,
        args.delta, args.seed, args.padding,
    )


if __name__ == "__main__":
    main()
