"""Pilot study of an Adler-compatible Gaussian margin tessellation.

Let G=(G_1,...,G_q) be iid smooth Gaussian random fields and define

    A_delta = {x: max_i G_i(x) - secondmax_i G_i(x) >= delta}.

The hitting set in R^q is a union of q separated polyhedral chambers. Therefore
A_delta is a Gaussian-related excursion and remains within the Gaussian
kinematic formula.  The excluded phase is a film around winner-tie surfaces.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr
from skimage.measure import euler_number


@dataclass(frozen=True)
class Metrics:
    q: int
    delta: float
    theoretical_phi: float
    measured_phi: float
    euler: int
    n_components: int
    largest_component_fraction_domain: float
    largest_component_fraction_phase: float
    spans_x: bool
    spans_y: bool
    spans_z: bool


def gaussian_margin_volume_fraction(q: int, delta: float) -> float:
    """P(max(Z)-secondmax(Z) >= delta), Z_i iid N(0,1)."""
    if q < 2:
        raise ValueError("q must be at least 2")
    if delta < 0:
        raise ValueError("delta must be non-negative")

    # Choose the winning coordinate (q choices), condition on its value t.
    normalizer = 1.0 / np.sqrt(2.0 * np.pi)

    def integrand(t: float) -> float:
        return normalizer * np.exp(-0.5 * t * t) * ndtr(t - delta) ** (q - 1)

    value, _ = quad(integrand, -10.0, 10.0, epsabs=2e-11, epsrel=2e-11, limit=200)
    return float(q * value)


def delta_for_volume_fraction(q: int, target_phi: float) -> float:
    if not 0.0 < target_phi < 1.0:
        raise ValueError("target_phi must lie in (0,1)")
    return float(
        brentq(
            lambda delta: gaussian_margin_volume_fraction(q, delta) - target_phi,
            0.0,
            10.0,
        )
    )


def generate_gaussian_field(
    n: int, correlation_length_voxels: float, rng: np.random.Generator
) -> np.ndarray:
    """Periodic unit-variance GRF with squared-exponential spectral filtering."""
    kx = 2.0 * np.pi * np.fft.fftfreq(n)
    ky = 2.0 * np.pi * np.fft.fftfreq(n)
    kz = 2.0 * np.pi * np.fft.rfftfreq(n)
    k2 = sum(a * a for a in np.meshgrid(kx, ky, kz, indexing="ij"))
    # C(r) proportional to exp(-r^2/lc^2), hence S(k) proportional to
    # exp(-lc^2 k^2/4).  sqrt(S) is the amplitude filter.
    amplitude = np.exp(-correlation_length_voxels**2 * k2 / 8.0)
    white = rng.normal(size=(n, n, n))
    field = np.fft.irfftn(
        np.fft.rfftn(white) * amplitude, s=white.shape, axes=(0, 1, 2)
    ).real
    field -= field.mean()
    field /= field.std()
    return field.astype(np.float32)


def margin_phase(fields: np.ndarray, delta: float) -> tuple[np.ndarray, np.ndarray]:
    """Return binary retained phase and winning-field labels."""
    top_two = np.partition(fields, kth=fields.shape[0] - 2, axis=0)[-2:]
    margin = top_two[1] - top_two[0]
    winners = np.argmax(fields, axis=0).astype(np.int16)
    return margin >= delta, winners


def phase_metrics(
    phase: np.ndarray,
    winners: np.ndarray,
    q: int,
    delta: float,
    theoretical_phi: float,
) -> Metrics:
    """Measure connectivity without falsely joining different winner cones.

    For delta>0 the q chambers of the hitting set are separated in value space,
    hence components carrying different winner labels cannot touch in the
    continuum.  A binary voxel labeling can nevertheless bridge across an
    under-resolved matrix film.  We therefore label each winner phase
    separately and combine their component statistics.
    """
    structure = ndimage.generate_binary_structure(rank=3, connectivity=1)
    largest = 0
    count = 0
    total_euler = 0
    span_flags = [False, False, False]
    for winner in range(q):
        labels, label_count = ndimage.label(
            phase & (winners == winner), structure=structure
        )
        count += int(label_count)
        total_euler += int(euler_number(phase & (winners == winner), connectivity=1))
        sizes = np.bincount(labels.ravel())
        largest = max(largest, int(sizes[1:].max(initial=0)))
        for axis in range(3):
            face0 = np.take(labels, 0, axis=axis)
            face1 = np.take(labels, -1, axis=axis)
            common = np.intersect1d(np.unique(face0), np.unique(face1))
            span_flags[axis] |= bool(np.any(common != 0))
    occupied = int(phase.sum())

    return Metrics(
        q=q,
        delta=delta,
        theoretical_phi=theoretical_phi,
        measured_phi=occupied / phase.size,
        euler=total_euler,
        n_components=int(count),
        largest_component_fraction_domain=largest / phase.size,
        largest_component_fraction_phase=(largest / occupied if occupied else 0.0),
        spans_x=span_flags[0],
        spans_y=span_flags[1],
        spans_z=span_flags[2],
    )


def run_sweep(
    n: int,
    lc: float,
    q_values: list[int],
    target_phis: list[float],
    realizations: int,
    seed: int,
) -> list[Metrics]:
    results: list[Metrics] = []
    max_q = max(q_values)
    for realization in range(realizations):
        rng = np.random.default_rng(seed + realization)
        all_fields = np.stack(
            [generate_gaussian_field(n, lc, rng) for _ in range(max_q)], axis=0
        )
        for q in q_values:
            fields = all_fields[:q]
            for target_phi in target_phis:
                delta = delta_for_volume_fraction(q, target_phi)
                theory_phi = gaussian_margin_volume_fraction(q, delta)
                phase, winners = margin_phase(fields, delta)
                results.append(phase_metrics(phase, winners, q, delta, theory_phi))
                print(realization, q, target_phi, results[-1])
    return results


def write_csv(path: Path, results: list[Metrics]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(Metrics.__annotations__))
        writer.writeheader()
        for result in results:
            writer.writerow(result.__dict__)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=64)
    parser.add_argument("--lc", type=float, default=5.0)
    parser.add_argument("--q", type=int, nargs="+", default=[2, 4, 8, 16])
    parser.add_argument(
        "--phi", type=float, nargs="+", default=[0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
    )
    parser.add_argument("--realizations", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--output", type=Path, default=Path("results/margin_sweep.csv"))
    args = parser.parse_args()
    results = run_sweep(args.n, args.lc, args.q, args.phi, args.realizations, args.seed)
    write_csv(args.output, results)


if __name__ == "__main__":
    main()
