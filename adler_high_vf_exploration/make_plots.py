"""Create saved PNG figures for the Gaussian margin-tessellation study."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from gaussian_minkowski import margin_m0, margin_m1
from margin_tessellation import (
    delta_for_volume_fraction,
    generate_gaussian_field,
    margin_phase,
)


def plot_gmfs(output: Path) -> None:
    deltas = np.linspace(0.02, 1.6, 100)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for q in (2, 4, 8, 16, 32):
        m0 = np.array([margin_m0(q, delta) for delta in deltas])
        m1 = np.array([margin_m1(q, delta) for delta in deltas])
        axes[0].plot(deltas, m0, label=rf"$q={q}$")
        axes[1].plot(deltas, m1, label=rf"$q={q}$")
    axes[0].set(xlabel=rf"margin $\delta$", ylabel=rf"$\mathcal{{M}}_0^\gamma(D_\delta)$")
    axes[1].set(xlabel=rf"margin $\delta$", ylabel=rf"$\mathcal{{M}}_1^\gamma(D_\delta)$")
    axes[0].set_title("Expected volume fraction")
    axes[1].set_title("Gaussian boundary functional")
    for axis in axes:
        axis.grid(alpha=0.3)
        axis.legend(ncols=2)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def read_sweep(path: Path) -> list[dict[str, float | int | bool]]:
    rows: list[dict[str, float | int | bool]] = []
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            rows.append(
                {
                    "q": int(row["q"]),
                    "phi": float(row["theoretical_phi"]),
                    "largest": float(row["largest_component_fraction_phase"]),
                    "spans": any(row[key] == "True" for key in ("spans_x", "spans_y", "spans_z")),
                }
            )
    return rows


def plot_phase_diagram(csv_path: Path, output: Path) -> None:
    grouped: dict[tuple[int, float], list[dict[str, float | int | bool]]] = defaultdict(list)
    for row in read_sweep(csv_path):
        grouped[(int(row["q"]), float(row["phi"]))].append(row)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for q in sorted({key[0] for key in grouped}):
        keys = sorted(key for key in grouped if key[0] == q)
        x = np.array([key[1] for key in keys])
        p_span = np.array([np.mean([row["spans"] for row in grouped[key]]) for key in keys])
        largest = np.array([np.mean([row["largest"] for row in grouped[key]]) for key in keys])
        axes[0].plot(x, p_span, "o-", label=rf"$q={q}$")
        axes[1].plot(x, largest, "o-", label=rf"$q={q}$")
    axes[0].set(
        xlabel="retained volume fraction",
        ylabel="finite-domain spanning frequency",
        ylim=(-0.05, 1.05),
    )
    axes[1].set(
        xlabel="retained volume fraction",
        ylabel="largest component / retained phase",
        yscale="log",
    )
    axes[0].set_title("Preliminary percolation diagram")
    axes[1].set_title("Component concentration")
    for axis in axes:
        axis.grid(alpha=0.3, which="both")
        axis.legend()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_morphology_slices(
    output: Path,
    n: int = 96,
    lc: float = 8.0,
    target_phi: float = 0.70,
    seed: int = 1729,
) -> None:
    q_values = (4, 8, 16)
    rng = np.random.default_rng(seed)
    fields = np.stack([generate_gaussian_field(n, lc, rng) for _ in range(max(q_values))])
    index = n // 2
    cmap = plt.colormaps["tab20"].with_extremes(bad="white")
    fig, axes = plt.subplots(2, len(q_values), figsize=(12, 7.2), constrained_layout=True)

    for column, q in enumerate(q_values):
        delta = delta_for_volume_fraction(q, target_phi)
        phase, winners = margin_phase(fields[:q], delta)
        binary_slice = phase[:, :, index]
        label_slice = np.where(binary_slice, winners[:, :, index] + 1, np.nan)

        axes[0, column].imshow(binary_slice.T, origin="lower", cmap="gray_r", interpolation="nearest")
        axes[1, column].imshow(
            label_slice.T,
            origin="lower",
            cmap=cmap,
            vmin=1,
            vmax=20,
            interpolation="nearest",
        )
        axes[0, column].set_title(rf"$q={q}$, $\delta={delta:.3f}$")
        axes[1, column].set_title("winner-labelled cores")
        for row in range(2):
            axes[row, column].set_xticks([])
            axes[row, column].set_yticks([])

    axes[0, 0].set_ylabel("binary retained phase")
    axes[1, 0].set_ylabel("separated Gaussian chambers")
    fig.suptitle(
        rf"Gaussian margin tessellations at target $\phi={target_phi:.2f}$ "
        rf"($N={n}$, $L_c={lc:g}$ voxels)"
    )
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", type=Path, default=Path("results/margin_sweep_n64.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("figures"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plot_gmfs(args.output_dir / "gaussian_minkowski_functionals.png")
    plot_phase_diagram(args.results, args.output_dir / "percolation_phase_diagram.png")
    plot_morphology_slices(args.output_dir / "margin_morphology_slices_phi070.png")
    print(f"Saved PNG figures in {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()

