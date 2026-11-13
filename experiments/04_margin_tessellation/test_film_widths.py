"""Analytical, sampling, and direct-field controls for normal film widths."""
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.stats import chi

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_film_widths as study
from theory import margin_film_widths as law
from theory.check_margin_tubes import face_weight


class AffineFields:
    def __init__(self, third=None):
        self.third = third

    def evaluate(self, points):
        x = np.asarray(points).reshape(-1, 3)[:, 0]
        values = [x/2, -x/2]
        if self.third is not None:
            values.append(np.full_like(x, self.third))
        return np.array(values)

    def value_gradient(self, point, label):
        if label == 2:
            return self.third, np.zeros(3)
        sign = 1 if label == 0 else -1
        return sign*np.asarray(point)[0]/2, np.array([sign/2, 0., 0.])


def planar_crossing():
    return dict(x=0., y=0., z=0., nx=1., ny=0., nz=0.,
                winner_i=0, winner_j=1, gradient_norm=1.)


class FilmWidthControls(unittest.TestCase):
    def test_area_size_bias_changes_chi3_to_chi4(self):
        mean = chi.mean(df=3)
        for width in (.3, .6, 1., 2.):
            integral = quad(lambda r: r*chi.pdf(r, df=3)/mean, 1/width, np.inf)[0]
            self.assertAlmostEqual(integral, float(law.width_cdf(width, 1., 1.)), places=11)
        self.assertGreater(abs(chi.sf(2, df=3)-float(law.width_cdf(.5, 1., 1.))), .1)

    def test_closed_cdf_and_inverse(self):
        delta, lc = .2, 5.
        width = np.array([.3, .7, 1.5])
        x = (delta*lc/width)**2/2
        np.testing.assert_allclose(law.width_cdf(width, delta, lc), np.exp(-x)*(1+x), atol=1e-14)
        p = np.array([.01, .1, .5, .9, .99])
        np.testing.assert_allclose(law.width_cdf(law.width_quantile(p, delta, lc), delta, lc), p, atol=1e-14)
        np.testing.assert_array_equal(law.width_cdf([-1., 0., np.inf], delta, lc), [0., 0., 1.])

    def test_moments_scaling_and_volume_consistency(self):
        self.assertAlmostEqual(law.width_moment(1, 1., 1.), np.sqrt(np.pi/8))
        self.assertAlmostEqual(law.width_moment(2, 1., 1.), .5)
        self.assertTrue(np.isinf(law.width_moment(4, 1., 1.)))
        np.testing.assert_allclose(law.width_quantile(.5, .2, 10.),
                                   2*law.width_quantile(.5, .2, 5.))
        for q in (2, 3, 16, 32):
            self.assertAlmostEqual(law.interface_area_density(q, 5.)*law.width_moment(1, .01, 5.),
                                   .01*face_weight(q, 0., 1)/np.sqrt(2), places=12)

    def test_axis_projection_weight_has_no_tangent_singularity(self):
        rng = np.random.default_rng(42)
        normals = rng.normal(size=(100, 3))
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        normals = np.concatenate((normals, np.eye(3)))
        for n in normals:
            weight = study.axis_weight(n)
            self.assertGreaterEqual(weight, 1/np.sqrt(3)-1e-14)
            self.assertLessEqual(weight, 1+1e-14)
            self.assertAlmostEqual(sum(abs(n))*weight, 1.)

    def test_gstools_spectral_values_and_gradients(self):
        fields = study.Fields(3, 5., 4517, modes=128)
        for point in (np.zeros(3), np.array([1.2, -3.1, .8])):
            direct = fields.evaluate(point)[:, 0]
            for label in range(3):
                value, gradient = fields.value_gradient(point, label)
                self.assertAlmostEqual(value, direct[label], places=12)
                step = 1e-4
                finite = np.array([(fields.evaluate(point+step*e, [label])[0, 0]
                                    - fields.evaluate(point-step*e, [label])[0, 0])/(2*step)
                                   for e in np.eye(3)])
                np.testing.assert_allclose(gradient, finite, atol=2e-9)

    def test_affine_tie_location(self):
        fields = AffineFields()
        origins = np.array([[0., .2, .5], [0., .8, .5]])
        t = np.array([-1., -.1, .1, 1.])
        points = np.repeat(origins[:, None, :], len(t), axis=1)
        points[:, :, 0] = t
        values = fields.evaluate(points.reshape(-1, 3)).reshape(2, 2, len(t))
        roots, rejected = study.locate_crossings(fields, origins, 0, t, values)
        self.assertEqual(rejected, 0)
        self.assertEqual(len(roots), 2)
        for root in roots:
            self.assertAlmostEqual(root["x"], 0.)
            self.assertAlmostEqual(root["gradient_norm"], 1.)
            self.assertAlmostEqual(root["area_weight"], 1.)

    def test_affine_normal_width_is_exact(self):
        rows, _ = study.normal_profile(AffineFields(), planar_crossing(), [.1, .2], 5., points=65)
        for row in rows:
            self.assertTrue(row["stable"])
            self.assertTrue(row["pair_controls_boundaries"])
            self.assertAlmostEqual(row["measured_width"], 2*row["delta"], places=11)

    def test_third_runner_changes_width_without_changing_endpoint_labels(self):
        rows, _ = study.normal_profile(AffineFields(-.03), planar_crossing(), [.2], 5., points=65)
        row = rows[0]
        self.assertTrue(row["stable"])
        self.assertTrue(row["pair_preserved"])
        self.assertFalse(row["pair_controls_boundaries"])
        self.assertAlmostEqual(row["measured_width"], .68, places=11)
        self.assertAlmostEqual(row["linearized_width"], .4, places=11)

    def test_censored_profile_is_not_a_measured_width(self):
        rows, _ = study.normal_profile(AffineFields(), planar_crossing(), [10.], 1., points=65)
        self.assertTrue(rows[0]["censored"])
        self.assertFalse(rows[0]["stable"])
        self.assertIsNone(rows[0]["measured_width"])

    def test_refinement_reports_new_crossings(self):
        row = dict(axis=0, line=0, winner_i=0, winner_j=1, coordinate=.3)
        second = dict(row, coordinate=.7)
        result = study.match_crossings([row], [dict(row), second])
        self.assertEqual(result["fine_unmatched"], 1)
        self.assertEqual(result["coarse_unmatched"], 0)

    def test_invalid_inputs(self):
        for delta, lc in ((0., 1.), (1., -1.), (float("nan"), 1.)):
            with self.assertRaises(ValueError):
                law.width_cdf(.1, delta, lc)
        with self.assertRaises(ValueError):
            law.width_quantile(1., .1, 5.)

    def test_existing_output_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            marker = Path(tmp)/"existing.txt"
            marker.write_text("keep")
            with self.assertRaises(ValueError):
                study.run({}, tmp)
            self.assertEqual(marker.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
