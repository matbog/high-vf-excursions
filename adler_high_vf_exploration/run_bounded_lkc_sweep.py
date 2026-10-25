"""Theoretical and numerical L_i(delta) curves on a bounded 100^3 cube."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bounded_lkcs import label_aware_lkcs, label_aware_mesh_lkcs
from gaussian_minkowski import margin_gmfs
from margin_tessellation import generate_gaussian_field, margin_phase


def expected_lkcs_cube(q: int, delta: float, side: float, lc: float) -> np.ndarray:
    """Roubin System 2.17 with the margin-set Gaussian Minkowski functionals."""
    m0, m1, m2, m3 = margin_gmfs(q, delta)
    spectral = np.sqrt((2.0 / lc**2) / (2.0 * np.pi))
    return np.array(
        [
            side**3 * spectral**3 * m3
            + 3.0 * side**2 * spectral**2 * m2
            + 3.0 * side * spectral * m1
            + m0,
            side**3 * spectral**2 * m2
            + 3.0 * side**2 * spectral * m1
            + 3.0 * side * m0,
            2.0 * side**3 * spectral * m1 + 3.0 * side**2 * m0,
            side**3 * m0,
        ],
        dtype=float,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--side", type=float, default=100.0)
    parser.add_argument("--lc", type=float, default=10.0)
    parser.add_argument("--q", type=int, nargs="+", default=[4, 8, 16])
    parser.add_argument(
        "--delta", type=float, nargs="+",
        default=[0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.55, 0.70, 0.90, 1.15],
    )
    parser.add_argument("--realizations", type=int, default=5)
    parser.add_argument("--estimator", choices=["mesh", "cubical"], default="mesh")
    parser.add_argument("--output-dir", type=Path, default=Path("results/bounded_lkc_sweep"))
    args = parser.parse_args()
    if not args.lc < 15.0:
        raise ValueError("This experiment requires lc < 15")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    dx = args.side / args.n

    rows = []
    max_q = max(args.q)
    for realization in range(args.realizations):
        rng = np.random.default_rng(44000 + realization)
        fields = np.stack(
            [generate_gaussian_field(args.n, args.lc / dx, rng) for _ in range(max_q)]
        )
        for q in args.q:
            for delta in args.delta:
                phase, winners = margin_phase(fields[:q], delta)
                estimator = label_aware_mesh_lkcs if args.estimator == "mesh" else label_aware_lkcs
                measured = estimator(phase, winners, q, dx).as_array()
                theoretical = expected_lkcs_cube(q, delta, args.side, args.lc)
                rows.append((realization, q, delta, *measured, *theoretical))
                print("realization", realization, "q", q, "delta", delta, "L", measured)

    csv_path = args.output_dir / "bounded_lkc_sweep.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["realization", "q", "delta", "L0", "L1", "L2", "L3",
             "expected_L0", "expected_L1", "expected_L2", "expected_L3"]
        )
        writer.writerows(rows)

    fig, axes = plt.subplots(4, len(args.q), figsize=(4.3 * len(args.q), 13), constrained_layout=True)
    labels = [r"$L_0$ Euler characteristic", r"$L_1$", r"$L_2$ half surface", r"$L_3$ volume"]
    for column, q in enumerate(args.q):
        for index in range(4):
            axis = axes[index, column]
            means, sems, theory = [], [], []
            for delta in args.delta:
                subset = [row for row in rows if row[1] == q and row[2] == delta]
                values = np.array([row[3 + index] for row in subset])
                means.append(values.mean())
                sems.append(values.std(ddof=1) / np.sqrt(len(values)))
                theory.append(subset[0][7 + index])
            axis.plot(args.delta, theory, "-", lw=2, label="GKF theory")
            axis.errorbar(args.delta, means, yerr=sems, fmt="o", ms=4, capsize=3,
                          label=f"numerical mean (n={args.realizations})")
            axis.axhline(0.0, color="black", lw=0.7)
            axis.grid(alpha=0.3)
            axis.set_xlabel(r"margin $\delta$")
            axis.set_ylabel(labels[index])
            if index == 0:
                axis.set_title(rf"$q={q}$")
            if index == 0 or (index == 3 and column == len(args.q) - 1):
                axis.legend(fontsize=8)
    fig.suptitle(
        rf"Bounded cube: $N={args.n}^3$, $T={args.side:g}$, "
        rf"$L_c={args.lc:g}$, $L_c/T={args.lc/args.side:.2f}$, estimator={args.estimator}"
    )
    fig.savefig(args.output_dir / "bounded_lkc_delta_sweep.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
