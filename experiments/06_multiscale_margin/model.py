"""Independent margin unions: volume calibration, conditional GKF, Gaussian sampling.

This experiment is isolated from the generators used by experiments 01-06.
All LKC arrays are ordered [L0, L1, L2, L3] in Euclidean units.
"""
from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from scipy import fft, ndimage
from skimage.measure import euler_number


def class_probabilities(total_phi: float, shares: Iterable[float]) -> np.ndarray:
    """Raw occupancies giving added fractions total_phi * share, coarse first."""
    shares = np.asarray(list(shares), dtype=float)
    if (not np.isfinite(total_phi) or not 0 < total_phi < 1
            or shares.ndim != 1 or len(shares) == 0
            or not np.all(np.isfinite(shares)) or np.any(shares <= 0)
            or not np.isclose(shares.sum(), 1.0, rtol=0, atol=1e-10)):
        raise ValueError("Require 0 < total_phi < 1 and positive shares summing to 1.")
    previous_share = np.cumsum(shares) - shares
    return total_phi * shares / (1.0 - total_phi * previous_share)


def cube_lkcs(length: float) -> np.ndarray:
    return np.array([1.0, 3 * length, 3 * length**2, length**3])


def intersection_operator(gmfs: Iterable[float], rho: float) -> np.ndarray:
    """GKF applied to an independent random observation domain, then averaged."""
    b0, b1, b2, b3 = np.asarray(list(gmfs)) * rho ** np.arange(4)
    return np.array([
        [b0, b1, b2, b3],
        [0.0, b0, np.pi / 2 * b1, 2 * b2],
        [0.0, 0.0, b0, 2 * b1],
        [0.0, 0.0, 0.0, b0],
    ])


def union_expectations(operators: Iterable[np.ndarray], length: float):
    """Return raw, previous-union intersection, and cumulative-union LKCs."""
    domain = cube_lkcs(length)
    previous = np.zeros(4)
    stages = []
    for operator in operators:
        raw = operator @ domain
        intersection = operator @ previous
        previous = previous + raw - intersection
        stages.append({"raw": raw, "intersection": intersection, "union": previous.copy()})
    return stages


class GaussianSampler:
    """Real-white-noise FFT filter, with deterministic ensemble normalization.

    The full periodic field is Gaussian with variance one. Its central crop is
    used as a bounded observation window. We retain the random DC coefficient;
    no empirical centering or rescaling is performed. Even-grid Nyquist modes
    are removed to give an unambiguous stationary trigonometric interpolant.
    """

    def __init__(self, n: int, length: float, lc: float, padding_length: float):
        if n < 8 or length <= 0 or lc <= 0 or padding_length < 0:
            raise ValueError("Invalid grid, correlation length, or padding.")
        self.n = n
        self.dx = length / n
        pad = int(np.ceil(padding_length / self.dx))
        self.full_n = fft.next_fast_len(n + 2 * pad, real=True)
        start = (self.full_n - n) // 2
        self.crop = (slice(start, start + n),) * 3
        k = 2 * np.pi * fft.fftfreq(self.full_n, d=self.dx)
        amp = np.exp(-lc**2 * k**2 / 8)
        if self.full_n % 2 == 0:
            amp[self.full_n // 2] = 0
        power = amp**2
        self.lambda2 = float(np.sum(k**2 * power) / power.sum())
        self.rho = float(np.sqrt(self.lambda2 / (2 * np.pi)))
        self.fd_variance = float(
            2 * np.sum((1 - np.cos(k * self.dx)) * power) / power.sum() / self.dx**2
        )
        norm = float(np.mean(power) ** 1.5)
        half = amp[:self.full_n // 2 + 1]
        self.filter = (
            amp[:, None, None] * amp[None, :, None] * half[None, None, :] / norm
        ).astype(np.float32)

    def sample(self, rng: np.random.Generator) -> np.ndarray:
        white = rng.standard_normal((self.full_n,) * 3, dtype=np.float32)
        spectrum = fft.rfftn(white, workers=1, overwrite_x=True)
        del white
        spectrum *= self.filter
        full = fft.irfftn(spectrum, s=(self.full_n,) * 3, workers=1, overwrite_x=True)
        return full[self.crop].copy()


def margin_from_fields(fields: Iterable[np.ndarray]):
    """Keep only the two largest values and winner while streaming competitors."""
    top = second = winner = None
    count = 0
    for i, field in enumerate(fields):
        if top is None:
            top = field.copy()
            second = np.full_like(top, -np.inf)
            winner = np.zeros(top.shape, dtype=np.int32)
        else:
            larger = field > top
            second = np.where(larger, top, np.maximum(second, field))
            top = np.maximum(top, field)
            winner[larger] = i
        count += 1
    if count < 2:
        raise ValueError("A margin requires at least two competitors.")
    return top - second, winner


def sample_margin(sampler: GaussianSampler, q: int, seed: int, realization: int,
                  class_index: int):
    """Independent class/component streams; thresholds reuse the same fields."""
    diagnostics = []

    def components():
        for i in range(q):
            rng = np.random.default_rng(np.random.SeedSequence([seed, realization, class_index, i]))
            field = sampler.sample(rng)
            record = {
                "realization": realization + 1, "class_index": class_index,
                "component": i, "mean": float(field.mean(dtype=np.float64)),
                "variance": float(field.var(dtype=np.float64)),
                "lambda2_spectral": sampler.lambda2,
                "forward_difference_variance_expected": sampler.fd_variance,
            }
            for axis, name in enumerate("xyz"):
                diff = np.diff(field, axis=axis) / sampler.dx
                record[f"forward_difference_variance_{name}"] = float(np.mean(diff**2))
            diagnostics.append(record)
            yield field

    margin, winner = margin_from_fields(components())
    return margin, winner, diagnostics


def topology(mask: np.ndarray, connectivity: int, *, return_labels: bool = False):
    """Bounded binary topology. Connectivity is either 6 or 26; no wrapping."""
    if connectivity not in (6, 26):
        raise ValueError("Connectivity must be 6 or 26.")
    rank = 1 if connectivity == 6 else 3
    labels, count = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, rank))
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    occupied = int(sizes.sum())
    spans = []
    for axis in range(3):
        shared = np.intersect1d(np.take(labels, 0, axis=axis),
                                np.take(labels, -1, axis=axis))
        spans.append(bool(np.any(shared > 0)))
    record = {
        f"euler_{connectivity}": int(euler_number(mask, connectivity=rank)),
        f"components_{connectivity}": int(count),
        f"largest_fraction_{connectivity}": float(sizes.max() / occupied) if occupied else 0.0,
        f"spans_{connectivity}": any(spans),
    }
    for axis, value in zip("xyz", spans):
        record[f"spans_{connectivity}_{axis}"] = value
    return (record, labels, sizes) if return_labels else record


def component_table(labels: np.ndarray, sizes: np.ndarray, owner: np.ndarray,
                    dx: float, n_classes: int):
    """Actual union components, with coarse-first class ownership for accounting."""
    boundary = np.zeros(len(sizes), dtype=bool)
    spanning = np.zeros_like(boundary)
    for axis in range(3):
        low = np.unique(np.take(labels, 0, axis=axis))
        high = np.unique(np.take(labels, -1, axis=axis))
        boundary[low] = boundary[high] = True
        spanning[np.intersect1d(low, high)] = True
    class_sizes = [np.bincount(labels[owner == k + 1], minlength=len(sizes))
                   for k in range(n_classes)]
    for label in range(1, len(sizes)):
        volume = float(sizes[label] * dx**3)
        row = {
            "component": label, "voxels": int(sizes[label]), "volume": volume,
            "equivalent_diameter": float((6 * volume / np.pi) ** (1 / 3)),
            "touches_boundary": bool(boundary[label]), "spans": bool(spanning[label]),
            "mixed_classes": sum(counts[label] > 0 for counts in class_sizes) > 1,
        }
        for k, counts in enumerate(class_sizes):
            row[f"class_{k}_volume_fraction"] = float(counts[label] / sizes[label])
        yield row
