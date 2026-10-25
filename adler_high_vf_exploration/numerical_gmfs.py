"""Estimate M0..M3 independently from periodic fields in dimensions 0..3.

On a d-dimensional flat torus, all domain LKCs except volume vanish. The GKF
therefore reduces to

  E[chi(A)] / Vol(M) = (lambda2 / 2pi)^(d/2) M_d^gamma(D).

This makes spatial dimension d an independent numerical probe of M_d.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from gaussian_minkowski import margin_m0, margin_m1
from margin_tessellation import delta_for_volume_fraction, margin_phase
from periodic_euler import periodic_euler_nd
from tube_gmfs import q2_exact_gmfs


def generate_gaussian_field_nd(
    shape: tuple[int, ...], lc: float, rng: np.random.Generator
) -> np.ndarray:
    frequencies = [2.0 * np.pi * np.fft.fftfreq(n) for n in shape]
    k2 = sum(component**2 for component in np.meshgrid(*frequencies, indexing="ij"))
    amplitude = np.exp(-lc**2 * k2 / 8.0)
    white = rng.normal(size=shape)
    field = np.fft.ifftn(np.fft.fftn(white) * amplitude).real
    field -= field.mean()
    field /= field.std()
    return field.astype(np.float32)


def estimate_one(
    q: int,
    delta: float,
    dimension: int,
    n: int,
    lc: float,
    realizations: int,
    seed: int,
) -> tuple[float, float, list[int]]:
    chis: list[int] = []
    shape = (n,) * dimension
    for realization in range(realizations):
        rng = np.random.default_rng(seed + 1009 * dimension + realization)
        fields = np.stack([generate_gaussian_field_nd(shape, lc, rng) for _ in range(q)])
        phase, winners = margin_phase(fields, delta)
        chi = sum(periodic_euler_nd(phase & (winners == label)) for label in range(q))
        chis.append(chi)
    density = np.asarray(chis, dtype=float) / n**dimension
    lambda2 = 2.0 / lc**2
    scale = (lambda2 / (2.0 * np.pi)) ** (dimension / 2.0)
    estimates = density / scale
    return float(estimates.mean()), float(estimates.std(ddof=1) / np.sqrt(realizations)), chis


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q", type=int, nargs="+", default=[2, 4, 8, 16])
    parser.add_argument("--phi", type=float, default=0.70)
    parser.add_argument("--lc", type=float, default=8.0)
    parser.add_argument("--n1", type=int, default=4096)
    parser.add_argument("--n2", type=int, default=384)
    parser.add_argument("--n3", type=int, default=64)
    parser.add_argument("--r1", type=int, default=100)
    parser.add_argument("--r2", type=int, default=60)
    parser.add_argument("--r3", type=int, default=30)
    parser.add_argument("--output-dir", type=Path, default=Path("results/numerical_gmfs"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    estimates_by_q: dict[int, np.ndarray] = {}
    errors_by_q: dict[int, np.ndarray] = {}
    for q in args.q:
        delta = delta_for_volume_fraction(q, args.phi)
        estimates = [margin_m0(q, delta)]
        errors = [0.0]
        for dimension, n, realizations in (
            (1, args.n1, args.r1),
            (2, args.n2, args.r2),
            (3, args.n3, args.r3),
        ):
            estimate, sem, chis = estimate_one(
                q, delta, dimension, n, args.lc, realizations, seed=31000 + 100 * q
            )
            estimates.append(estimate)
            errors.append(sem)
            rows.append((q, delta, dimension, n, args.lc, realizations, estimate, sem))
        estimates_by_q[q] = np.array(estimates)
        errors_by_q[q] = np.array(errors)
        print("q", q, "delta", delta, "M0..M3", estimates, "SEM", errors)

    with (args.output_dir / "numerical_gmfs.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["q", "delta", "dimension", "n", "lc", "realizations", "estimate", "sem"])
        writer.writerows(rows)

    fig, axes = plt.subplots(1, len(args.q), figsize=(3.5 * len(args.q), 3.8), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for axis, q in zip(axes, args.q):
        x = np.arange(4)
        axis.errorbar(x, estimates_by_q[q], yerr=errors_by_q[q], fmt="o", capsize=4, label="periodic-field estimate")
        axis.scatter([0, 1], [margin_m0(q, delta_for_volume_fraction(q, args.phi)), margin_m1(q, delta_for_volume_fraction(q, args.phi))], marker="s", label="analytic $M_0,M_1$")
        if q == 2:
            exact = q2_exact_gmfs(delta_for_volume_fraction(2, args.phi))
            axis.plot(x, exact, "x--", label="exact q=2")
        axis.axhline(0.0, color="black", lw=0.7)
        axis.set_xticks(x, [rf"$M_{j}$" for j in x])
        axis.set_title(rf"$q={q}$, $\phi={args.phi:.2f}$")
        axis.grid(alpha=0.3)
        axis.legend(fontsize=7)
    fig.suptitle("Gaussian Minkowski functionals from periodic Euler densities")
    fig.savefig(args.output_dir / "numerical_gmfs.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

