"""Compare value-space GMF estimates with 3D excursion Euler characteristics."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from margin_tessellation import delta_for_volume_fraction, generate_gaussian_field, margin_phase
from periodic_euler import periodic_euler_label_aware
from tube_gmfs import estimate_tube_gmfs, q2_exact_gmfs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=64)
    parser.add_argument("--lc", type=float, default=5.0)
    parser.add_argument("--q", type=int, nargs="+", default=[2, 4, 8, 16])
    parser.add_argument("--phi", type=float, default=0.70)
    parser.add_argument("--realizations", type=int, default=20)
    parser.add_argument("--sobol-power", type=int, default=19)
    parser.add_argument("--output-dir", type=Path, default=Path("results/euler_experiment"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    gmfs: dict[int, np.ndarray] = {}
    tube_estimates = {}
    for q in args.q:
        delta = delta_for_volume_fraction(q, args.phi)
        estimate = estimate_tube_gmfs(q, delta, power=args.sobol_power)
        gmfs[q] = estimate.gmfs
        tube_estimates[q] = estimate
        print("GMFs", q, delta, estimate.gmfs, "tube RMS", estimate.residual_rms)

    rows = []
    max_q = max(args.q)
    for realization in range(args.realizations):
        rng = np.random.default_rng(9000 + realization)
        fields = np.stack(
            [generate_gaussian_field(args.n, args.lc, rng) for _ in range(max_q)], axis=0
        )
        for q in args.q:
            delta = delta_for_volume_fraction(q, args.phi)
            phase, winners = margin_phase(fields[:q], delta)
            chi = periodic_euler_label_aware(phase, winners, q)
            rows.append((realization, q, delta, float(phase.mean()), chi))

    with (args.output_dir / "periodic_euler.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["realization", "q", "delta", "volume_fraction", "periodic_euler"])
        writer.writerows(rows)

    # For C(r)=exp(-r^2/lc^2), derivative variance is lambda_2=2/lc^2.
    # A 3-torus has no boundary LKCs, hence E[chi]/V=(lambda_2/2pi)^(3/2) M3.
    lambda2 = 2.0 / args.lc**2
    prefactor = (lambda2 / (2.0 * np.pi)) ** 1.5
    q_axis = np.array(args.q)
    measured_mean = np.array(
        [np.mean([row[4] for row in rows if row[1] == q]) / args.n**3 for q in args.q]
    )
    measured_sem = np.array(
        [
            np.std([row[4] for row in rows if row[1] == q], ddof=1)
            / np.sqrt(args.realizations)
            / args.n**3
            for q in args.q
        ]
    )
    predicted = np.array([prefactor * gmfs[q][3] for q in args.q])

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
    for q in args.q:
        estimate = tube_estimates[q]
        fit = sum(
            estimate.coefficients[j] * estimate.radii**j
            for j in range(len(estimate.coefficients))
        )
        axes[0].plot(estimate.radii, estimate.probabilities, "o", ms=3, label=rf"$q={q}$ samples")
        axes[0].plot(estimate.radii, fit, "-", lw=1.5)
    axes[0].set(xlabel=rf"tube radius $\rho$", ylabel="Gaussian tube probability")
    axes[0].set_title("Value-space tube expansions")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8, ncols=2)

    axes[1].errorbar(q_axis, measured_mean, yerr=measured_sem, fmt="o", capsize=4, label="3D periodic simulations")
    axes[1].plot(q_axis, predicted, "s--", label=rf"GKF from estimated $\mathcal{{M}}_3^\gamma$")
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set(xlabel=rf"number of Gaussian fields $q$", ylabel=rf"Euler density $E[\chi]/V$")
    axes[1].set_title(rf"Euler characteristic at $\phi={args.phi:.2f}$")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.savefig(args.output_dir / "tube_and_euler_validation.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    # Explicit q=2 diagnostic for the numerical tube differentiation.
    if 2 in args.q:
        exact = q2_exact_gmfs(delta_for_volume_fraction(2, args.phi))
        print("q=2 exact", exact, "estimated", gmfs[2])


if __name__ == "__main__":
    main()

