"""Geometric, sampling and tiny-scene controls for the local film audit."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_local_films as films
from src import fields
from src.fields import GridSpec


class LocalFilmControls(unittest.TestCase):
    def test_affine_gap_lost_by_margin_interpolation(self):
        endpoints = np.array([[-.5, .5], [.5, -.5]])
        regions, gaps = films.affine_intervals(endpoints, .2)
        np.testing.assert_allclose(gaps, [[.4, .6]])
        self.assertEqual({r[2] for r in regions}, {0, 1})
        margins, _ = films.margin_and_winner(endpoints)
        self.assertEqual(films.sampled_gaps(np.array([0., 1.]), margins-.2), [])
        self.assertAlmostEqual(films.total_width(gaps), .2)

    def test_third_winner_can_appear_between_endpoint_winners(self):
        endpoints = np.array([[0., 0.], [-1.5, .5], [.5, -1.5]])
        regions, gaps = films.affine_intervals(endpoints, .2)
        np.testing.assert_allclose(gaps, [[.15, .35], [.65, .85]])
        middle = next(r for r in regions if r[2] == 0)
        np.testing.assert_allclose(middle[:2], [.35, .65])
        self.assertNotIn(0, np.argmax(endpoints, axis=0))
        # Omitting the interior competitor underestimates the matrix interval.
        _, two_field_gaps = films.affine_intervals(endpoints[1:], .2)
        self.assertLess(films.total_width(two_field_gaps), films.total_width(gaps))

    def test_same_label_endpoints_can_hide_a_curved_field_gap(self):
        t = np.linspace(0., 1., 257)
        margin = (t-.5)**2
        gaps = films.sampled_gaps(t, margin-.1)
        expected = [[.5-np.sqrt(.1), .5+np.sqrt(.1)]]
        np.testing.assert_allclose(gaps, expected, atol=1e-5)
        _, affine = films.affine_intervals(np.array([[.25, .25], [0., 0.]]), .1)
        self.assertEqual(affine, [])
        self.assertGreater(films.symmetric_difference(affine, gaps), .63)

    def test_unresolved_winner_switch_and_interval_disagreement(self):
        self.assertEqual(films.unresolved_switches(np.ones(2), np.array([0, 1]), .2), 1)
        self.assertEqual(films.unresolved_switches(np.array([1., 0., 1.]), np.array([0, 0, 1]), .2), 0)
        self.assertIsNone(films.boundary_shift([], [(0., .1)]))
        self.assertAlmostEqual(films.symmetric_difference([(.2, .4)], [(.3, .5)]), .2)
        self.assertAlmostEqual(films.boundary_shift([(.2, .4)], [(.3, .5)]), .1)
        self.assertEqual(films.complement([(0., .3), (.3, 1.)]), [])

    def test_face_pair_selection_and_eligible_population(self):
        margin = np.ones((4,)*3)
        winners = (np.indices(margin.shape)[0] >= 2).astype(int)
        rows, pool = films.select_segments(margin, winners, .5, 100, "cross_label", np.random.default_rng(9))
        self.assertEqual(pool["pool_total"], 16)
        self.assertEqual(pool["pool_per_axis"], [16, 0, 0])
        self.assertEqual(len(rows), 16)
        self.assertTrue(all(r["axis"] == 0 and r["i"] == 1 for r in rows))
        controls, pool = films.select_segments(margin, winners, .5, 8, "same_label", np.random.default_rng(9))
        self.assertEqual(pool["pool_total"], 3*3*4*4-16)
        self.assertEqual(len({(r["axis"], r["i"], r["j"], r["k"]) for r in controls}), 8)
        self.assertTrue(all(r["winner_start"] == r["winner_end"] for r in controls))

    def test_seed_stream_matches_shared_sampler_and_batching(self):
        grid = GridSpec(4, 4.)
        legacy = fields.generate_margin_fields(grid, 2., 3, np.random.default_rng(321729))
        points = np.indices((4,)*3).reshape(3, -1).T
        sampler = films.OriginalFields(3, 2., 1., 321729, batch=7)
        result = sampler.evaluate(points).astype(np.float32).reshape((3, 4, 4, 4))
        np.testing.assert_array_equal(result, legacy)
        sampler.batch = 100
        np.testing.assert_array_equal(sampler.evaluate(points).astype(np.float32).reshape(legacy.shape), legacy)

    def test_tiny_cli_seed_guard_and_saved_outputs(self):
        with tempfile.TemporaryDirectory(prefix="rf_films_") as directory:
            source, output = Path(directory)/"scene", Path(directory)/"result"
            source.mkdir()
            n, dx, q, lc, seed, modes = 8, .5, 2, 1., 321729, 16
            points = np.indices((n,)*3).reshape(3, -1).T
            sampler = films.OriginalFields(q, lc, dx, seed, mode_no=modes)
            values = sampler.evaluate(points).astype(np.float32)
            margin, winners = films.margin_and_winner(values)
            phase = margin >= films.saved.adler.delta_for_vf(q, .75)
            for filename, name, array in (("margin.vti", "margin", margin),
                                          ("winners.vti", "winner_label", winners),
                                          ("phase.vti", "phase", phase.astype(np.uint8)),
                                          ("field_g1.vti", "G_1", values[0])):
                mesh = pv.ImageData(dimensions=(n,)*3, spacing=(dx,)*3)
                mesh.point_data[name] = array.reshape((n,)*3).flatten(order="F")
                mesh.save(source/filename)
            shaped_margin, shaped_winners, _, _ = films.saved.load_scene(source, q, .75)
            wrong = films.OriginalFields(q, lc, dx, seed+1, mode_no=modes)
            with self.assertRaisesRegex(ValueError, "not reproduced"):
                films.validate_nodes(wrong, points[:32], shaped_margin, shaped_winners, source)
            before = {p.name: films.saved.sha256(p) for p in source.iterdir()}
            command = [sys.executable, "-B", str(HERE/"run_local_films.py"),
                       "--source", str(source), "--output", str(output), "--q", str(q),
                       "--lc", str(lc), "--field-seed", str(seed), "--mode-no", str(modes),
                       "--contacts", "3", "--controls", "3", "--validation-nodes", "16",
                       "--phis", ".75", "--max-points", "65", "--tolerance", ".01"]
            result = subprocess.run(command, capture_output=True, text=True, timeout=120,
                                    env=dict(os.environ, PYTHONIOENCODING="utf-8"))
            self.assertEqual(result.returncode, 0, result.stdout+"\n"+result.stderr)
            metadata = json.loads((output/"metadata.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["segment_count"], 6)
            self.assertEqual(metadata["seed_validation"]["margin_max_abs"], 0.)
            self.assertTrue((output/"profile_examples.png").exists())
            self.assertTrue((output/"profiles.npz").exists())
            self.assertEqual({p.name: films.saved.sha256(p) for p in source.iterdir()}, before)
            repeated = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("Output is not empty", repeated.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
