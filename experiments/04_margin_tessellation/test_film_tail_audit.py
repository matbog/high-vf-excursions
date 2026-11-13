from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import audit_film_tails as audit
from test_film_widths import AffineFields, planar_crossing


def source_rows(delta, extent=.25):
    return [dict(q=2, realization=0, crossing_id=0, target_vf=.75, margin_factor=1.,
                 delta=delta, extent=extent, area_weight=1., linearized_width=2*delta,
                 measured_width=None, relative_error=None, stable=False, censored=True)]


def measure(fields, delta, extent=.25, **kwargs):
    t = np.linspace(-extent, extent, 129)
    crossing = dict(planar_crossing(), third_clearance=10.)
    profile = audit.Profile(fields, crossing, t, fields.evaluate(np.column_stack((t, np.zeros((129, 2))))), 8193)
    return audit.measure_profile(profile, source_rows(delta, extent), lc=1., **kwargs)


class FilmTailControls(unittest.TestCase):
    def test_extended_affine_profile_resolves_censoring(self):
        rows, log = measure(AffineFields(), .4)
        row = rows[0]
        self.assertTrue(row["stable"])
        self.assertFalse(row["censored"])
        self.assertGreater(row["extent"], .25)
        self.assertAlmostEqual(row["measured_width"], .8, places=10)
        self.assertAlmostEqual(row["curvature_contribution"], 0., places=10)
        self.assertAlmostEqual(row["competition_contribution"], 0., places=10)
        self.assertGreaterEqual(log[-1]["stable_objects"], 2)

    def test_third_runner_is_isolated_from_curvature(self):
        rows, _ = measure(AffineFields(third=-.03), .2)
        row = rows[0]
        self.assertAlmostEqual(row["pair_width"], .4, places=10)
        self.assertAlmostEqual(row["measured_width"], .68, places=10)
        self.assertAlmostEqual(row["curvature_contribution"], 0., places=10)
        self.assertAlmostEqual(row["competition_contribution"], .7, places=10)
        self.assertTrue(row["pair_preserved"])
        self.assertFalse(row["pair_controls_boundaries"])

    def test_curved_pair_has_no_competition_contribution(self):
        class Curved:
            def evaluate(self, points):
                x = np.asarray(points).reshape(-1, 3)[:, 0]
                h = x+x**3
                return np.array([h/2, -h/2])
        rows, _ = measure(Curved(), .4)
        row = rows[0]
        self.assertTrue(row["stable"])
        self.assertLess(row["curvature_contribution"], -.05)
        self.assertAlmostEqual(row["competition_contribution"], 0.)
        self.assertAlmostEqual(row["relative_error"], row["curvature_contribution"]+row["competition_contribution"])

    def test_extent_cap_keeps_censoring_and_bounds(self):
        rows, _ = measure(AffineFields(), 2., extent=.5, max_extent_lc=1.)
        row = rows[0]
        self.assertTrue(row["censored"])
        self.assertFalse(row["stable"])
        self.assertIsNone(row["measured_width"])
        self.assertEqual(row["audit_stop"], "extent_cap")
        self.assertAlmostEqual(row["width_lower"], 2.)
        self.assertTrue(np.isinf(row["width_upper"]))

    def test_unique_point_budget(self):
        fields = AffineFields()
        p = audit.Profile(fields, planar_crossing(), [0.], fields.evaluate([[0., 0., 0.]]), 3)
        p.evaluate([0., 0., .1, .1, .2])
        self.assertEqual(len(p.cache), 3)
        with self.assertRaises(audit.BudgetReached):
            p.evaluate([.3])

    def test_selection_covers_all_failures_without_weight_changes(self):
        rows = []
        for key in range(5):
            row = dict(source_rows(.1)[0], crossing_id=key, stable=True, censored=False, relative_error=0.)
            rows.append(row)
        rows[0]["censored"] = True
        rows[1]["stable"] = False
        rows[2]["relative_error"] = .2
        groups, reasons = audit.select_profiles(rows, controls=1)
        self.assertTrue({0, 1, 2}.issubset(reasons))
        self.assertEqual(list(reasons.values()).count("control"), 1)
        self.assertTrue(all(group[0]["area_weight"] == 1 for group in groups.values()))

    def test_cdf_bounds_include_unresolved_weight(self):
        known = audit.inherited(dict(source_rows(.1)[0], stable=True, censored=False,
                                     measured_width=.2, relative_error=0.))
        unresolved = audit.inherited(source_rows(.1)[0])
        unresolved.update(width_lower=.4, width_upper=float("inf"))
        summaries, cdf = audit.summarize([known, unresolved])
        at_quarter = next(r for r in cdf if r["width_threshold"] == .25)
        self.assertEqual(at_quarter["area_cdf_lower"], .5)
        self.assertEqual(at_quarter["area_cdf_upper"], .5)
        at_half = next(r for r in cdf if r["width_threshold"] == .5)
        self.assertEqual(at_half["area_cdf_lower"], .5)
        self.assertEqual(at_half["area_cdf_upper"], 1.)
        self.assertTrue(np.isinf(summaries[0]["mean_width_upper"]))

    def test_nonempty_outputs_are_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"marker"
            path.write_text("keep")
            with self.assertRaises(ValueError):
                audit.run("missing", "missing", folder)
            self.assertEqual(path.read_text(), "keep")


if __name__ == "__main__":
    unittest.main()
