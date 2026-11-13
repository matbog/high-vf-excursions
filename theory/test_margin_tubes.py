"""Analytic and independent geometry controls for the margin theory draft."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_margin_tubes as check
from src.adler import margin_gmfs


class MarginTubeControls(unittest.TestCase):
    def test_q2_halfspaces_and_factorials(self):
        for delta in (.05, .6, 2.):
            a = delta/np.sqrt(2)
            phi = np.exp(-a*a/2)/np.sqrt(2*np.pi)
            expected = [2*ndtr(-a), 2*phi, 2*a*phi, 2*(a*a-1)*phi]
            np.testing.assert_allclose(check.branch_gmfs(2, delta), expected, atol=2e-13)
            for radius in (0., .2*a, .8*a):
                self.assertAlmostEqual(check.cone_tube(2, delta, radius),
                                       2*ndtr(radius-a), places=12)

    def test_normals_orientation_and_gram(self):
        delta = .7
        for k in (1, 2, 3, 4):
            normals = check.outward_normals(k)
            gram = normals.T@normals
            np.testing.assert_allclose(gram, (np.eye(k)+np.ones((k, k)))/2, atol=1e-15)
            self.assertAlmostEqual(np.linalg.det(gram), (k+1)/2**k)
            point = np.full(k+1, -.2)
            point[0] += delta
            np.testing.assert_allclose(point@normals, -delta/np.sqrt(2), atol=1e-15)

    def test_normal_cone_moments(self):
        weights, h = check.normal_angles(2)
        self.assertAlmostEqual(weights.sum(), np.pi/3, places=12)
        self.assertAlmostEqual(weights@h, 2/np.sqrt(3), places=12)
        weights, _ = check.normal_angles(3)
        self.assertAlmostEqual(weights.sum(), check.OMEGA, places=12)

    def test_low_dimensional_closed_strata(self):
        for delta in (0., .1, .7, 2.):
            self.assertAlmostEqual(check.face_weight(3, delta, 1),
                                   3*np.sqrt(2/np.pi)*np.exp(-delta**2/4)*ndtr(-delta/np.sqrt(6)), places=12)
            self.assertAlmostEqual(check.face_weight(3, delta, 2),
                                   3/(2*np.pi)*np.exp(-delta**2/3), places=12)
            self.assertAlmostEqual(check.face_weight(4, delta, 3),
                                   4/(2*np.pi)**1.5*np.exp(-3*delta**2/8), places=12)

    def test_stratum_derivative_recurrence(self):
        delta, step = .4, 1e-4
        for q in (2, 3, 4, 8, 32):
            for k in range(min(q, 5)):
                f = lambda d: check.face_weight(q, d, k)
                derivative = (f(delta-2*step)-8*f(delta-step)+8*f(delta+step)-f(delta+2*step))/(12*step)
                expected = (-k*delta/(k+1)*f(delta)
                            - np.sqrt((k+1)/(k+2))*check.face_weight(q, delta, k+1))
                self.assertAlmostEqual(derivative, expected, delta=2e-9)

    def test_projection_against_constrained_optimizer(self):
        rng = np.random.default_rng(831)
        for q in (3, 4):
            for _ in range(6):
                z = np.sort(rng.normal(size=q))[::-1]
                delta = .7
                projected = check.ordered_projection(z, delta)
                result = minimize(lambda y: .5*np.sum((y-z)**2), z,
                                  jac=lambda y: y-z, method="SLSQP",
                                  constraints=[dict(type="ineq", fun=lambda y: y[0]-y[1:]-delta)],
                                  options=dict(ftol=1e-13, maxiter=100))
                self.assertTrue(result.success)
                np.testing.assert_allclose(projected, result.x, atol=1e-7)
                np.testing.assert_allclose(check.ordered_projection(z+3, delta), projected+3, atol=1e-14)
                self.assertTrue(np.all(projected[0]-projected[1:] >= delta-1e-13))

    def test_ordered_sector_covers_gaussian_directions(self):
        for q in (3, 4):
            directions, weights = check.ordered_directions(q, 64)
            np.testing.assert_allclose(np.linalg.norm(directions, axis=1), 1, atol=1e-14)
            np.testing.assert_allclose(directions.sum(axis=1), 0, atol=1e-14)
            self.assertTrue(np.all(np.diff(directions, axis=1) <= 1e-14))
            self.assertAlmostEqual(weights.sum(), 1, places=12)

    def test_distance_tube_without_face_decomposition(self):
        for q in (3, 4):
            delta = .6
            for fraction in (0., .35, .65):
                radius = fraction*delta/np.sqrt(2)
                self.assertAlmostEqual(check.distance_tube(q, delta, radius, 192),
                                       check.cone_tube(q, delta, radius), delta=3e-7)

    def test_cubic_tube_remainder(self):
        for q in (3, 4):
            delta = .6
            m = check.branch_gmfs(q, delta)
            errors = []
            for radius in (.04, .02, .01):
                polynomial = m@[1, radius, radius**2/2, radius**3/6]
                errors.append(abs(check.cone_tube(q, delta, radius)-polynomial))
            self.assertGreater(errors[0]/errors[1], 14)
            self.assertGreater(errors[1]/errors[2], 14)

    def test_high_vf_expansion_and_maximum_identity(self):
        for q in (2, 3, 4, 8, 32, 64):
            limit, slope = check.asymptotic_coefficients(q)
            np.testing.assert_allclose(check.branch_gmfs(q, 0.), limit, atol=1e-12)
            self.assertAlmostEqual(limit[1]/np.sqrt(2), check.expected_maximum(q), places=11)
            errors = [np.max(abs(check.branch_gmfs(q, delta)-limit-slope*delta))
                      for delta in (.002, .001)]
            self.assertGreater(errors[0]/errors[1], 3.9)
            for delta in (.01, .005):
                epsilon = 1-check.face_weight(q, delta, 0)
                estimate = (np.sqrt(2)/limit[1]*epsilon
                            + np.sqrt(2/3)*check.face_weight(q, 0., 2)/limit[1]**3*epsilon**2)
                self.assertLess(abs(delta-estimate), 2*delta**3)

    def test_zero_margin_is_separate_and_radius_is_restricted(self):
        np.testing.assert_array_equal(margin_gmfs(4, 0.), [1, 0, 0, 0])
        self.assertGreater(check.branch_gmfs(4, 0.)[1], 1)
        for function in (check.cone_tube, check.distance_tube):
            for delta, radius in ((0., 0.), (.6, .6/np.sqrt(2)), (.6, -.1)):
                with self.assertRaises(ValueError):
                    function(4, delta, radius)

    def test_nonempty_results_are_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            marker = Path(folder)/"existing.txt"
            marker.write_text("preserve")
            with self.assertRaises(ValueError):
                check.run_checks(folder, 16)
            self.assertEqual(marker.read_text(), "preserve")


if __name__ == "__main__":
    unittest.main()
