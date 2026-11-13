from pathlib import Path
import json
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_film_widths as base
import refine_film_widths as refined
from test_film_widths import AffineFields


class RefinedFilmControls(unittest.TestCase):
    def test_stored_spectrum_matches_direct_field_and_gradient(self):
        original = base.Fields(3, 5., 539, 128)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"modes.npz"
            original.save(path)
            stored = refined.StoredFields(path)
            points = np.random.default_rng(6).normal(size=(20, 3))
            np.testing.assert_allclose(stored.evaluate(points), original.evaluate(points), atol=2e-14)
            for i in range(3):
                np.testing.assert_allclose(stored.value_gradient(points[0], i)[1],
                                           original.value_gradient(points[0], i)[1], atol=1e-14)

    def test_hidden_winner_is_split_into_two_valid_ties(self):
        fields = AffineFields(third=.1)
        origins = np.zeros((1, 3))
        t = np.array([-1., 1.])
        points = np.column_stack((t, np.zeros((2, 2))))
        values = fields.evaluate(points).reshape(3, 1, 2)
        old, rejected = base.locate_crossings(fields, origins, 0, t, values)
        self.assertEqual(len(old), 0)
        self.assertEqual(rejected, 1)
        roots, splits = refined.adaptive_crossings(fields, origins, 0, t, values)
        self.assertEqual(splits, 1)
        np.testing.assert_allclose([r["x"] for r in roots], [-.2, .2], atol=1e-12)
        self.assertTrue(all(r["third_clearance"] >= 0 for r in roots))

    def test_pure_pair_agrees_with_original_locator(self):
        fields = AffineFields()
        origins = np.zeros((1, 3))
        t = np.array([-1., 1.])
        values = fields.evaluate(np.column_stack((t, np.zeros((2, 2))))).reshape(2, 1, 2)
        original, _ = base.locate_crossings(fields, origins, 0, t, values)
        current, splits = refined.adaptive_crossings(fields, origins, 0, t, values)
        self.assertEqual(splits, 0)
        self.assertEqual(base.match_crossings(original, current)["matched_count"], 1)

    def test_two_realizations_keep_integer_refinement_factor(self):
        # Regression: a margin-factor loop must not overwrite the line-grid
        # refinement factor before the next independent realization.
        config = dict(qs=[2], realizations=2, lines_per_axis=1, profiles_per_group=1,
                      subdivisions=4, profile_points=17, modes=32, seed=37,
                      lc=2., length=20., vfs=[.75, .85], factors=[1., .5, .25])
        with tempfile.TemporaryDirectory() as folder:
            source, output = Path(folder)/"source", Path(folder)/"refined"
            base.run(config, source)
            refined.run(source, output, 8)
            metadata = json.loads((output/"metadata.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertTrue((output/"profiles_q2_r1.npz").is_file())


if __name__ == "__main__":
    unittest.main()
