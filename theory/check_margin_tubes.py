"""Independent value-space checks for margin_tube_derivation.md.

No fields or spatial grids are generated. Production GMFs are imported only
inside run_checks, for comparison with geometric quadrature. delta=0 in
face_weight/branch_gmfs denotes a right-hand limit, NOT the full-space GMFs.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import gammaincc, log_ndtr, ndtr, roots_legendre

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OMEGA = 3 * np.arccos(1 / 3) - np.pi


@lru_cache(maxsize=None)
def face_weight(q: int, delta: float, k: int) -> float:
    """S_k including S_0=M_0; k>=q is an absent stratum.

    Integrate the entire real line, with the multiplicity inside the
    integrand so the absolute quadrature tolerance applies to the result.
    """
    if q < 2 or delta < 0 or k < 0:
        raise ValueError("Require q>=2, delta>=0, k>=0")
    if k >= q:
        return 0.0
    log_scale = np.log(q * math.comb(q - 1, k) * np.sqrt(k + 1))

    def integrand(t):
        log_value = (log_scale - .5 * ((t + delta)**2 + k * t*t)
                     - .5 * (k + 1) * np.log(2 * np.pi)
                     + (q - 1 - k) * log_ndtr(t))
        return np.exp(log_value)

    value, error = quad(integrand, -np.inf, np.inf,
                        epsabs=3e-13, epsrel=3e-13, limit=250)
    if error > 2e-10 * max(1, abs(value)):
        raise RuntimeError(f"Face integral did not converge: q={q}, k={k}")
    return float(value)


def outward_normals(k: int) -> np.ndarray:
    normals = np.zeros((k + 1, k))
    normals[0, :] = -1 / np.sqrt(2)
    normals[1:, :] = np.eye(k) / np.sqrt(2)
    return normals


def unit_interval(order):
    x, w = roots_legendre(order)
    return (x + 1) / 2, w / 2


def spherical_simplex(vectors, order):
    """Directions and spherical area weights from a positive simplex chart.

    Columns span a pointed cone; k=2 is an arc, k=3 a spherical triangle.
    Returns unit directions, angular weights, and the unnormalized norm.
    """
    k = vectors.shape[1]
    x, w = unit_interval(order)
    if k == 2:
        bary = np.stack((x, 1 - x))
        parameter_weights = w
    elif k == 3:
        u, v = np.meshgrid(x, x, indexing="ij")
        wu, wv = np.meshgrid(w, w, indexing="ij")
        bary = np.stack((u.ravel(), ((1-u)*v).ravel(),
                         ((1-u)*(1-v)).ravel()))
        parameter_weights = (wu * wv * (1-u)).ravel()
    else:
        raise ValueError("Only two- and three-generator cones are needed")
    points = vectors @ bary
    lengths = np.linalg.norm(points, axis=0)
    jacobian = np.sqrt(np.linalg.det(vectors.T @ vectors)) / lengths**k
    return (points / lengths).T, parameter_weights * jacobian, lengths


@lru_cache(maxsize=None)
def normal_angles(k, order=48):
    if k == 1:
        return np.ones(1), np.ones(1)
    _, weights, lengths = spherical_simplex(outward_normals(k), order)
    # Sum of coefficients of a unit normal direction in the normal basis.
    return weights, 1 / lengths


def branch_gmfs(q, delta, order=48):
    """Geometric GMFs, including their positive-margin limit at zero."""
    s = [face_weight(q, delta, k) for k in range(4)]
    a = delta / np.sqrt(2)
    w2, h2 = normal_angles(2, order)
    w3, _ = normal_angles(3, order)
    return np.array([s[0], s[1], a*s[1] + w2.sum()*s[2],
                     (a*a-1)*s[1] + 2*a*np.dot(w2, h2)*s[2]
                     + 2*w3.sum()*s[3]])


def cone_integral(k, radius, a, order=48):
    angular_weights, h = normal_angles(k, order)
    nodes, weights = unit_interval(order)
    rho = radius * nodes
    radial_weights = radius * weights * rho**(k-1)
    # Small matrix: at most order**2 directions by order radial nodes.
    values = np.exp(a*h[:, None]*rho - .5*rho**2) @ radial_weights
    return float(np.dot(angular_weights, values))


def cone_tube(q, delta, radius, order=48):
    """Full tube measure for q=2,3,4, before chamber tubes overlap."""
    if q not in (2, 3, 4) or delta <= 0 or not 0 <= radius < delta/np.sqrt(2):
        raise ValueError("Require q=2,3,4, delta>0, 0<=r<delta/sqrt(2)")
    return face_weight(q, delta, 0) + sum(
        face_weight(q, delta, k) * cone_integral(k, radius, delta/np.sqrt(2), order)
        for k in range(1, q))


def ordered_projection(points, delta):
    """Project descending-sorted rows onto chamber 0, by active constraints."""
    points = np.asarray(points, dtype=float)
    if np.any(np.diff(points, axis=-1) > 1e-12):
        raise ValueError("Rows must be sorted in descending order")
    t = points[..., 0].copy()
    cumulative = np.cumsum(points, axis=-1)
    for j in range(1, points.shape[-1]):
        active = points[..., j] + delta > t
        t = np.where(active, (cumulative[..., j] + j*delta)/(j+1), t)
    projection = np.minimum(points, np.expand_dims(t-delta, -1))
    projection[..., 0] = t
    return projection


@lru_cache(maxsize=None)
def ordered_directions(q, order):
    """One coordinate-order sector in the orthogonal contrast space."""
    rays = np.stack([(np.arange(q) < j).astype(float) - j/q
                     for j in range(1, q)], axis=1)
    rays /= np.linalg.norm(rays, axis=0)
    directions, weights, _ = spherical_simplex(rays, order)
    sphere_area = 2*np.pi if q == 3 else 4*np.pi
    return directions, weights * math.factorial(q)/sphere_area


def distance_tube(q, delta, radius, order=128):
    """Independent radial Gaussian integral using nearest-point distances.

    Uses neither face weights nor normal-cone coefficients. In an ordered
    sector the nearest winner chamber is chamber 0 (swap argument in note).
    Orthogonal contrasts have chi_(q-1) radius and uniform angular law.
    """
    if q not in (2, 3, 4) or delta <= 0 or not 0 <= radius < delta/np.sqrt(2):
        raise ValueError("Require q=2,3,4, delta>0, 0<=r<delta/sqrt(2)")
    if q == 2:
        return float(2*ndtr(radius-delta/np.sqrt(2)))
    directions, weights = ordered_directions(q, order)
    upper = delta / (directions[:, 0] - directions[:, 1])
    if radius == 0:
        boundary = upper
    else:
        lower = np.zeros_like(upper)
        for _ in range(54):
            middle = (lower + upper)/2
            points = directions * middle[:, None]
            projection = ordered_projection(points, delta)
            outside = np.sum((points-projection)**2, axis=1) > radius**2
            lower = np.where(outside, middle, lower)
            upper = np.where(outside, upper, middle)
        boundary = (lower+upper)/2
    return float(np.dot(weights, gammaincc((q-1)/2, boundary**2/2)))


def asymptotic_coefficients(q):
    s1, s2, s3, s4 = [face_weight(q, 0.0, k) for k in range(1, 5)]
    limit = np.array([1, s1, np.pi*s2/3, -s1+2*OMEGA*s3])
    slope = np.array([-s1/np.sqrt(2), -np.sqrt(2/3)*s2,
                      s1/np.sqrt(2)-np.pi*s3/(2*np.sqrt(3)),
                      np.sqrt(6)*s2-4*OMEGA*s4/np.sqrt(5)])
    return limit, slope


def expected_maximum(q):
    return quad(lambda t: q*t*np.exp(-t*t/2)/np.sqrt(2*np.pi)
                * ndtr(t)**(q-1), -np.inf, np.inf,
                epsabs=3e-13, epsrel=3e-13)[0]


def tight_truncated_face_weight(q, delta, k):
    """Precision audit only: production bounds, tighter raw-integral tolerance."""
    phi = lambda t: np.exp(-t*t/2)/np.sqrt(2*np.pi)
    integral = quad(lambda t: phi(t+delta)*phi(t)**k*ndtr(t)**(q-1-k),
                    -10., 10., epsabs=2e-16, epsrel=2e-13, limit=300)[0]
    return float(q*math.comb(q-1, k)*np.sqrt(k+1)*integral)


def csv_write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_checks(output, order):
    from src.adler import margin_gmfs, margin_stratum_measure

    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"Refusing nonempty output directory: {output}")
    output.mkdir(parents=True, exist_ok=True)
    w2, h2 = normal_angles(2)
    w3, _ = normal_angles(3)
    cone_rows = [dict(quantity=name, quadrature=float(value), exact=float(exact),
                      absolute_error=float(abs(value-exact))) for name, value, exact in
                 [("two_normal_angle", w2.sum(), np.pi/3),
                  ("two_normal_first_moment", np.dot(w2, h2), 2/np.sqrt(3)),
                  ("three_normal_solid_angle", w3.sum(), OMEGA)]]
    csv_write(output / "normal_cones.csv", cone_rows)
    gmfs = []
    for q in (2, 3, 4, 8, 16, 32, 64):
        for delta in (.01, .1, .3, .6, 1.2):
            independent = branch_gmfs(q, delta)
            production = margin_gmfs(q, delta)
            for j in range(4):
                gmfs.append(dict(q=q, delta=delta, j=j,
                                 geometric=float(independent[j]),
                                 production=float(production[j]),
                                 absolute_error=float(abs(independent[j]-production[j]))))
    csv_write(output / "gmf_comparison.csv", gmfs)
    precision_rows = []
    for k in (1, 2, 3):
        reference = face_weight(32, .01, k)
        tight = tight_truncated_face_weight(32, .01, k)
        precision_rows.append(dict(q=32, delta=.01, k=k, full_line=reference,
                                   production=margin_stratum_measure(32, .01, k),
                                   tighter_truncated=tight,
                                   tighter_absolute_error=abs(tight-reference)))
    csv_write(output / "quadrature_precision.csv", precision_rows)
    tube_rows = []
    for q in (2, 3, 4):
        for delta in (.3, .6, 1.2):
            for fraction in (0., .15, .35, .65):
                radius = fraction * delta/np.sqrt(2)
                geometric = cone_tube(q, delta, radius)
                coarse = distance_tube(q, delta, radius, order)
                fine = distance_tube(q, delta, radius, 2*order)
                tube_rows.append(dict(q=q, delta=delta, radius=radius,
                                      radius_over_separation_half=fraction,
                                      cone_tube=geometric, distance_coarse=coarse,
                                      distance_fine=fine,
                                      order=order, refined_order=2*order,
                                      refinement_change=abs(fine-coarse),
                                      absolute_error=abs(fine-geometric)))
        print(f"q={q} distance-tube checks complete", flush=True)
    csv_write(output / "distance_tubes.csv", tube_rows)
    remainders = []
    for q in (3, 4):
        delta = .6
        m = branch_gmfs(q, delta)
        for fraction in (.4, .2, .1, .05):
            radius = fraction*delta/np.sqrt(2)
            polynomial = sum(m[j]*radius**j/math.factorial(j) for j in range(4))
            remainder = cone_tube(q, delta, radius)-polynomial
            remainders.append(dict(q=q, delta=delta, radius=radius,
                                   remainder=remainder, remainder_over_r4=remainder/radius**4))
    csv_write(output / "cubic_remainders.csv", remainders)
    limits, high_vf = [], []
    h = 64/(5*np.sqrt(np.pi))
    cube_weights = np.array([1, 3*h, 3*h*h, h**3])
    for q in (2, 3, 4, 8, 16, 32, 64):
        limit, slope = asymptotic_coefficients(q)
        c = limit[1]/np.sqrt(2)
        s2 = face_weight(q, 0., 2)
        limits.append(dict(q=q, L=64, lc=5, vf_loss_slope=c,
                           expected_gaussian_max=expected_maximum(q),
                           M1_right_limit=limit[1], M2_right_limit=limit[2],
                           M3_right_limit=limit[3],
                           L0_right_limit=float(cube_weights@limit),
                           L0_delta_slope=float(cube_weights@slope)))
        for vf in (.9, .95, .99, .999):
            eps = 1-vf
            delta = brentq(lambda d: face_weight(q, d, 0)-vf, 0., 4., xtol=5e-15)
            delta_first = eps/c
            delta_second = delta_first + np.sqrt(2/3)*s2/limit[1]**3*eps**2
            exact = branch_gmfs(q, delta)
            high_vf.append(dict(q=q, target_vf=vf, delta=delta,
                                delta_first_order=delta_first,
                                delta_second_order=delta_second,
                                delta_second_error=abs(delta-delta_second),
                                L0=float(cube_weights@exact),
                                L0_first_order_in_matrix_vf=float(cube_weights@(limit+slope*eps/c)),
                                L3=64**3*vf,
                                gmf_linear_residual=float(np.max(abs(exact-limit-slope*delta)))))
    csv_write(output / "high_vf_limits.csv", limits)
    csv_write(output / "high_vf_comparison.csv", high_vf)
    max_gmf_error = max(row["absolute_error"] for row in gmfs)
    max_tube_error = max(row["absolute_error"] for row in tube_rows)
    max_tube_change = max(row["refinement_change"] for row in tube_rows)
    # Numerical acceptance thresholds are diagnostics, not rigorous error bounds.
    checks = dict(normal_moments=max(r["absolute_error"] for r in cone_rows) < 1e-12,
                  # The raw production integral has epsabs=2e-11 before
                  # combinatorial scaling. The precision audit isolates the
                  # resulting few-billionths error; this is not a sign error.
                  production_agreement=max_gmf_error < 1e-8,
                  tighter_quadrature=max(r["tighter_absolute_error"] for r in precision_rows) < 1e-11,
                  distance_agreement=max_tube_error < 2e-7,
                  distance_refinement=max_tube_change < 2e-6,
                  maximum_identity=bool(max(abs(r["vf_loss_slope"]-r["expected_gaussian_max"])
                                            for r in limits) < 1e-11))
    paths = [Path(__file__), Path(__file__).with_name("test_margin_tubes.py"),
             ROOT / "src/adler.py", ROOT / "src/fields.py"]
    summary = dict(schema_version=1, study="margin value-space tube checks",
                   study_version="2026-11-10", executed_utc=datetime.now(timezone.utc).isoformat(),
                   python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__,
                   angular_order=order, refined_angular_order=2*order,
                   cone_order=48, cube_L=64, cube_lc=5,
                   max_gmf_absolute_error=max_gmf_error,
                   max_tube_absolute_error=max_tube_error,
                   max_tube_refinement_change=max_tube_change,
                   production_absolute_tolerance=1e-8,
                   distance_absolute_tolerance=2e-7, refinement_absolute_tolerance=2e-6,
                   numerical_checks=checks, passed=all(checks.values()),
                   source_sha256={p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in paths},
                   output_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted(output.glob("*.csv"))})
    (output / "summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("max_gmf_absolute_error", "max_tube_absolute_error",
                                             "max_tube_refinement_change", "passed")}, indent=2))
    if not summary["passed"]:
        raise RuntimeError("A numerical check failed; inspect saved results before drawing conclusions")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("results_margin_tubes"))
    parser.add_argument("--order", type=int, default=128, help="Distance angular quadrature order, then doubled")
    args = parser.parse_args()
    if args.order < 16:
        parser.error("--order must be at least 16")
    run_checks(args.output, args.order)


if __name__ == "__main__":
    main()
