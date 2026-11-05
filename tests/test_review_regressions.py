"""Regression controls for the review fixes; no full experiment is launched."""
from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")
from src import adler, fields, measures
from src.fields import GridSpec


def local_module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lobes = local_module("lobe_measurement", "experiments/05_coalescing_lobes/measurement.py")
osc = local_module("oscillatory_controls", "explorations/oscillatory_bochner/run.py")


class TheoryControls(unittest.TestCase):
    grid = GridSpec(n=12, length=20.)

    def same_lkcs(self, a, b, **kwargs):
        for key in ("L0", "L1", "L2", "L3", "phi"):
            np.testing.assert_allclose(a[key], b[key], **kwargs)

    def test_cube_expansion_against_general_flag_formula(self):
        # Independent general-dimension expression catches missing pi/2 flags.
        u, lc, side = .7, 4., self.grid.length
        gmfs = [norm.sf(u), norm.pdf(u), u*norm.pdf(u), (u*u-1)*norm.pdf(u)]
        domain = [1, 3*side, 3*side**2, side**3]
        omega = lambda d: np.pi**(d/2) / gamma(1+d/2)
        result = adler.expected_lkcs_gaussian(self.grid, lc, np.array([u]))
        for i in range(4):
            expected = sum(math.comb(i+j, j)*omega(i+j)/(omega(i)*omega(j))
                           * (1/(lc*np.sqrt(np.pi)))**j * domain[i+j]*gmfs[j]
                           for j in range(4-i))
            self.assertAlmostEqual(float(result[f"L{i}"][0]), expected, places=10)

    def test_gaussian_component_rescaling(self):
        u = np.array([-.5, 0., 1.3])
        self.same_lkcs(adler.expected_lkcs_gaussian(self.grid, 4., 3*u, sigma=3),
                       adler.expected_lkcs_gaussian(self.grid, 4., u))

    def test_student_common_component_scale_cancels(self):
        u = np.array([-.8, 0., 1.5])
        self.same_lkcs(adler.expected_lkcs_t(self.grid, 4., 3, u, sigma=7),
                       adler.expected_lkcs_t(self.grid, 4., 3, u))

    def test_student_zero_threshold_is_numerator_halfspace(self):
        for nu in (3, 5, 10, 1000):
            self.same_lkcs(adler.expected_lkcs_t(self.grid, 4., nu, np.array([0.])),
                           adler.expected_lkcs_gaussian(self.grid, 4., np.array([0.])))

    def test_large_student_degrees_approach_gaussian_without_overflow(self):
        u = np.array([-.5, .3, 1.4])
        self.same_lkcs(adler.expected_lkcs_t(self.grid, 4., 100000, u),
                       adler.expected_lkcs_gaussian(self.grid, 4., u), rtol=1e-4, atol=1e-5)

    def test_student_contract_rejects_unsupported_degrees(self):
        for nu in (0, 1, 2, np.inf, np.nan):
            with self.assertRaises(ValueError):
                adler.expected_lkcs_t(self.grid, 4., nu, np.array([0.]))
        for nu in (3.5, True, 2):
            with self.assertRaises(ValueError):
                fields.generate_student_t_rf(self.grid, 4., nu, np.random.default_rng(1))

    def test_chi_squared_one_equals_two_gaussian_tails(self):
        u = np.array([.2, .7, 2.])
        chi = adler.expected_lkcs_chi2(self.grid, 4., 1, (3*u)**2, sigma=3)
        gaussian = adler.expected_lkcs_gaussian(self.grid, 4., u)
        self.same_lkcs(chi, {k: 2*v for k, v in gaussian.items()})

    def test_margin_two_equals_two_halfspaces(self):
        delta = .7
        margin = adler.expected_lkcs_margin(self.grid, 4., 2, delta)
        gaussian = adler.expected_lkcs_gaussian(self.grid, 4., np.array([delta/np.sqrt(2)]))
        self.same_lkcs(margin, {k: 2*v for k, v in gaussian.items()}, rtol=1e-11)

    def test_margin_amplitude_standardizes_value_margin(self):
        self.same_lkcs(adler.expected_lkcs_margin(self.grid, 4., 8, .9, sigma=3),
                       adler.expected_lkcs_margin(self.grid, 4., 8, .3))

    def test_zero_margin_is_full_domain(self):
        side = self.grid.length
        for q in (2, 4, 16):
            np.testing.assert_array_equal(adler.margin_gmfs(q, 0), [1, 0, 0, 0])
            self.same_lkcs(adler.expected_lkcs_margin(self.grid, 4., q, 0),
                           dict(L0=1, L1=3*side, L2=3*side**2, L3=side**3, phi=1))
            self.assertEqual(adler.margin_m1(q, 0), 0)

    def test_margin_facet_measure_matches_probability_derivative(self):
        h, delta = 1e-5, .4
        for q in (2, 4, 16):
            derivative = (adler.margin_m0(q, delta+h)-adler.margin_m0(q, delta-h))/(2*h)
            self.assertAlmostEqual(adler.margin_m1(q, delta), -np.sqrt(2)*derivative, places=8)

    def test_parameter_contracts(self):
        for q, delta in ((1, .5), (2.5, .5), (True, .5), (4, -.1), (4, np.nan)):
            with self.assertRaises(ValueError):
                adler.margin_gmfs(q, delta)
        with self.assertRaises(ValueError):
            adler.find_percolation_threshold(self.grid, 4., nu=3, k_chi2=1)
        with self.assertRaises(ValueError):
            adler.expected_lkcs_chi2(self.grid, 4., 1, np.array([0.]))
        with self.assertRaises(ValueError):
            adler.expected_lkcs_t(self.grid, 4., 3, np.array([0.]), sigma=0)

    def test_central_ball_derivatives_against_radial_density(self):
        for r in (.5, 1., 2.):
            k = np.sqrt(2/np.pi)*np.exp(-r*r/2)
            expected = [k*r*r, k*(2*r-r**3), k*(2-5*r*r+r**4)]
            np.testing.assert_allclose(adler.ball_lkcs(np.zeros(3), r)[1:], expected,
                                       rtol=2e-6, atol=2e-7)

    def test_ball_scaling_and_contact_guard(self):
        centres = np.array([[-2., 0., 0.], [2., 0., 0.]])
        self.same_lkcs(adler.expected_lkcs_lobes_disjoint(self.grid, 4., centres*3, 1.5, sigma=3),
                       adler.expected_lkcs_lobes_disjoint(self.grid, 4., centres, .5))
        for r in (2., 2.1):
            with self.assertRaises(ValueError):
                adler.expected_lkcs_lobes_disjoint(self.grid, 4., centres, r)


class GeometryControls(unittest.TestCase):
    def test_crofton_full_cube_includes_exact_face_term(self):
        self.assertEqual(measures.measure_surface_area_crofton(np.ones((6,)*3, bool), .5), 27.)

    def test_all_positive_signed_field_still_has_closed_boundary(self):
        full = np.full((5,)*3, 5.)
        self.assertGreater(measures.measure_l1_from_field(full, 1), 0)
        self.assertGreater(measures.measure_surface_area_from_field(full, 1), 0)
        self.assertEqual(measures.measure_surface_area_from_field(-full, 1), 0)

    def test_connectivity_is_an_explicit_choice(self):
        diagonal = np.zeros((3,)*3, bool)
        diagonal[np.arange(3), np.arange(3), np.arange(3)] = True
        self.assertFalse(any(measures.measure_percolation(diagonal, 1).values()))
        self.assertTrue(all(measures.measure_percolation(diagonal, 3).values()))
        self.assertEqual(measures.measure_euler(diagonal, 1), 3)
        self.assertEqual(measures.measure_euler(diagonal, 3), 1)

    def test_overlapping_partition_is_not_union_topology(self):
        grid = GridSpec(n=8, length=8.)
        x, y, z = np.indices((8,)*3, dtype=float)
        distances = np.stack([np.sqrt((x-cx)**2+(y-cy)**2+(z-3.5)**2)
                              for cx, cy in ((1, 1), (1, 6), (6, 1), (6, 6))])
        m = lobes.measure_lobes(distances, 20., grid)
        self.assertEqual(m["union_euler_26"], 1)
        self.assertEqual(m["partition_euler_26"], 4)
        self.assertEqual(m["sum_ball_euler_26"], 4)
        self.assertEqual(m["union_L3"], 8**3)
        self.assertEqual(m["sum_ball_L3"], 4*8**3)
        self.assertGreater(m["partition_L2_crofton"], m["union_L2_crofton"])
        self.assertTrue(m["union_spans_6"])

    def test_duplicate_ball_does_not_double_union(self):
        grid = GridSpec(n=9, length=9.)
        xyz = np.indices((9,)*3, dtype=float)
        d = np.sqrt(((xyz-4)**2).sum(axis=0))
        one = lobes.measure_lobes(d[None], 2.4, grid)
        duplicate = lobes.measure_lobes(np.stack([d, d]), 2.4, grid)
        for key in ("union_L3", "union_L1_field", "union_L2_field", "union_L2_crofton", "union_euler_26"):
            self.assertAlmostEqual(one[key], duplicate[key])
        self.assertEqual(duplicate["sum_ball_L3"], 2*one["union_L3"])

    def test_separated_ball_controls_are_additive(self):
        grid = GridSpec(n=15, length=15.)
        x, y, z = np.indices((15,)*3, dtype=float)
        distances = np.stack([np.sqrt((x-cx)**2+(y-7)**2+(z-7)**2) for cx in (3, 11)])
        m = lobes.measure_lobes(distances, 2., grid)
        self.assertEqual(m["union_euler_26"], 2)
        for union, per_ball in (("union_L3", "sum_ball_L3"), ("union_L2_crofton", "sum_ball_L2_crofton"),
                                ("union_euler_26", "sum_ball_euler_26")):
            self.assertEqual(m[union], m[per_ball])

    def test_empty_union(self):
        m = lobes.measure_lobes(np.full((2, 5, 5, 5), 10.), 1., GridSpec(n=5, length=5))
        for key in ("union_euler_6", "union_euler_26", "union_L3", "union_L1_field",
                    "union_L2_field", "union_L2_crofton", "union_largest_fraction_26"):
            self.assertEqual(m[key], 0)
        self.assertFalse(m["union_spans_26"])


class OscillatoryControls(unittest.TestCase):
    def test_covariance_is_radial_inverse_of_implemented_spectrum(self):
        for ld, lam in ((2., 3.), (8., 8.), (3., 100.)):
            omega = 2*np.pi/lam
            spectrum = lambda k: np.exp(-ld**2*(k-omega)**2/4)+np.exp(-ld**2*(k+omega)**2/4)
            upper = omega + 16/ld
            mass = quad(lambda k: k*k*spectrum(k), 0, upper, epsabs=1e-12)[0]
            for r in (0., 1., 4.):
                expected = quad(lambda k: k*k*spectrum(k)*np.sinc(k*r/np.pi),
                                0, upper, epsabs=1e-12)[0]/mass
                self.assertAlmostEqual(float(osc.oscillatory_covariance(r, ld, lam)), expected, places=9)

    def test_nonzero_complex_angle_has_negative_spectrum(self):
        for theta in (.1, .4, .7):
            k = np.sqrt(2*(3*theta+np.pi)/np.sin(2*theta))
            self.assertLess(osc.re_spectral_density(np.array([k]), 1., theta)[0], 0)
        with self.assertRaises(ValueError):
            osc.re_spectral_density(np.array([0.]), 1., np.pi/4)

    def test_sampler_ensemble_normalization_preserves_random_mean_and_scale(self):
        rng = np.random.default_rng(412)
        samples = [osc.generate_oscillatory_rf(8, 1., 2., 50., rng) for _ in range(250)]
        means = np.array([a.mean() for a in samples])
        second = np.array([np.mean(a*a) for a in samples])
        self.assertLess(abs(means.mean()), .1)
        self.assertLess(abs(second.mean()-1), .08)
        self.assertGreater(means.std(), .1)
        self.assertGreater(second.std(), .1)


class RunnerControls(unittest.TestCase):
    def test_lobe_smoke_outputs_and_overwrite_protection(self):
        with tempfile.TemporaryDirectory(prefix="rf_review_") as tmp:
            output = Path(tmp)/"union"
            command = [sys.executable, "-B", str(ROOT/"experiments/05_coalescing_lobes/run.py"),
                       "--n", "8", "--realizations", "1", "--radii", ".65", "1.15",
                       "--mc-samples", "2000", "--output", str(output)]
            completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True,
                                       timeout=180, env=dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8"))
            self.assertEqual(completed.returncode, 0, completed.stdout+"\n"+completed.stderr)
            with (output/"coalescing_lobes_sweep.csv").open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual(len(rows), 2)
            self.assertTrue(np.isfinite(float(rows[0]["theory_L0"])))
            self.assertTrue(np.isnan(float(rows[1]["theory_L0"])))
            for row in rows:
                self.assertAlmostEqual(float(row["union_L3"])/100**3, float(row["measured_phi"]))
                self.assertGreaterEqual(float(row["sum_ball_L3"]), float(row["union_L3"]))
                self.assertIn("partition_euler_26", row)
            metadata = json.loads((output/"metadata.json").read_text())
            self.assertEqual(metadata["measurement_schema"], "union_v2")
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(len(list(output.glob("*.png"))), 3)
            before = (output/"coalescing_lobes_sweep.csv").read_bytes()
            refused = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=60)
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn("not empty", refused.stderr)
            self.assertEqual((output/"coalescing_lobes_sweep.csv").read_bytes(), before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
