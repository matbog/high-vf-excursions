"""Small geometric and VTI controls for rethresholding a saved scene."""
from __future__ import annotations

import csv
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
import run_saved_sweep as sweep
from src.fields import GridSpec


def save_array(folder, filename, name, array, origin=(0., 0., 0.)):
    mesh = pv.ImageData(dimensions=array.shape, spacing=(.5,)*3, origin=origin)
    mesh.point_data[name] = array.flatten(order="F")
    mesh.save(folder/filename)


def make_scene(folder):
    folder.mkdir()
    margin = (np.arange(64).reshape((4,)*3)*.03).astype(np.float32)
    winners = (np.indices((4,)*3)[0] % 2).astype(np.int16)
    phase = (margin >= sweep.adler.delta_for_vf(2, .75)).astype(np.uint8)
    for filename, name, array in (("margin.vti", "margin", margin),
                                  ("winners.vti", "winner_label", winners),
                                  ("phase.vti", "phase", phase)):
        save_array(folder, filename, name, array)
    return margin, winners


class SavedSweepControls(unittest.TestCase):
    def test_cross_label_bridge_does_not_make_a_winner_span(self):
        margin = np.zeros((5,)*3)
        margin[:, 2, 2] = 10.
        winners = np.zeros((5,)*3, dtype=int)
        winners[3:] = 1
        rows, labels, components = sweep.measure_threshold(margin, winners, GridSpec(5, 5.), 2, 1., .75)
        for row in rows:
            self.assertTrue(row["union_spans_x"])
            self.assertFalse(row["union_spans_y"])
            self.assertEqual(row["spanning_winners"], 0)
            self.assertEqual(row["per_winner_components"], 2)
            self.assertEqual(row["union_components"], 1)
            self.assertEqual(row["cross_label_contacts"], 1)
            self.assertEqual(row["occupied_neighbor_pairs"], 4)
            self.assertEqual(row["euler_additivity_gap"], -1)
            self.assertEqual(row["largest_separate_component_fraction"], 3/5)
            self.assertEqual(row["separate_spanning_volume_fraction"], 0.)
            selected = [r for r in components if r["connectivity"] == row["connectivity"]]
            self.assertEqual(sum(r["volume"] for r in selected), 5.)

    def test_diagonal_connectivity_and_cross_label_pair_count(self):
        mask = np.zeros((3,)*3, dtype=bool)
        mask[np.arange(3), np.arange(3), np.arange(3)] = True
        six, _ = sweep.component_statistics(mask, 1., 6)
        diagonal, _ = sweep.component_statistics(mask, 1., 26)
        self.assertEqual(six["components"], 3)
        self.assertFalse(six["spans"])
        self.assertEqual(diagonal["components"], 1)
        self.assertTrue(diagonal["spans"])
        winners = np.indices(mask.shape)[0]
        contacts = sweep.cross_label_contacts(mask, winners)
        self.assertEqual(contacts[6]["cross_label_contacts"], 0)
        self.assertEqual(contacts[26]["cross_label_contacts"], 2)

    def test_cavity_euler_boundary_and_volume_weighted_sizes(self):
        mask = np.ones((5,)*3, dtype=bool)
        mask[2, 2, 2] = False
        for c in (6, 26):
            stats, components = sweep.component_statistics(mask, .5, c)
            self.assertEqual(stats["L0"], 2)
            self.assertEqual(stats["components"], 1)
            self.assertEqual(stats["volume"], 124*.5**3)
            self.assertEqual(stats["boundary_volume"], stats["volume"])
            self.assertEqual(stats["spanning_volume"], stats["volume"])
            self.assertTrue(components[0]["touches_boundary"])
            self.assertAlmostEqual(components[0]["equivalent_diameter"],
                                   (6*stats["volume"]/np.pi)**(1/3))
        example = [dict(equivalent_diameter=1., volume=1.),
                   dict(equivalent_diameter=2., volume=8.)]
        self.assertEqual(sweep.size_quantile(example, .5), 2.)

    def test_empty_phase(self):
        rows, labels, components = sweep.measure_threshold(
            np.zeros((4,)*3), np.zeros((4,)*3, dtype=int), GridSpec(4, 4.), 2, 1., .75)
        self.assertFalse(components)
        self.assertEqual(len(labels), 4)
        for row in rows:
            self.assertEqual(row["measured_phi"], 0.)
            self.assertEqual(row["union_L0"], 0)
            self.assertEqual(row["cross_label_fraction"], 0.)
            self.assertEqual(row["diameter_volume_p90"], 0.)

    def test_vti_order_volume_convention_and_baseline_guard(self):
        with tempfile.TemporaryDirectory(prefix="rf_saved_scene_") as directory:
            folder = Path(directory)/"scene"
            expected_margin, expected_winners = make_scene(folder)
            margin, winners, grid, info = sweep.load_scene(folder, 2, .75)
            np.testing.assert_array_equal(margin, expected_margin)
            np.testing.assert_array_equal(winners, expected_winners)
            self.assertEqual(grid.length, 2.)
            self.assertEqual(info["point_coordinate_span"], 1.5)
            self.assertEqual(info["reference_mask_mismatches"], 0)
            save_array(folder, "phase.vti", "phase", np.zeros((4,)*3, dtype=np.uint8))
            with self.assertRaisesRegex(ValueError, "reference phase disagrees"):
                sweep.load_scene(folder, 2, .75)

    def test_mismatched_geometry_and_noninteger_winners_rejected(self):
        with tempfile.TemporaryDirectory(prefix="rf_saved_bad_") as directory:
            folder = Path(directory)/"scene"
            _, winners = make_scene(folder)
            save_array(folder, "winners.vti", "winner_label", winners, origin=(1., 0., 0.))
            with self.assertRaisesRegex(ValueError, "identical grid"):
                sweep.load_scene(folder, 2, .75)
            save_array(folder, "winners.vti", "winner_label", winners.astype(float)+.5)
            with self.assertRaisesRegex(ValueError, "integers"):
                sweep.load_scene(folder, 2, .75)

    def test_small_cli_outputs_provenance_and_overwrite_protection(self):
        with tempfile.TemporaryDirectory(prefix="rf_saved_cli_") as directory:
            source, output = Path(directory)/"scene", Path(directory)/"output"
            make_scene(source)
            before = {p.name: sweep.sha256(p) for p in source.iterdir()}
            command = [sys.executable, "-B", str(HERE/"run_saved_sweep.py"),
                       "--source", str(source), "--output", str(output), "--q", "2",
                       "--lc", "1", "--target-phis", ".6", ".8"]
            result = subprocess.run(command, capture_output=True, text=True, timeout=120,
                                    env=dict(os.environ, PYTHONIOENCODING="utf-8"))
            self.assertEqual(result.returncode, 0, result.stdout+"\n"+result.stderr)
            metadata = json.loads((output/"metadata.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["completed_thresholds"], 2)
            self.assertEqual(metadata["independent_realizations"], 1)
            self.assertEqual(metadata["scene"]["input_sha256"], before)
            self.assertEqual({p.name: sweep.sha256(p) for p in source.iterdir()}, before)
            with (output/"measurements.csv").open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(len(rows), 4)
            self.assertLessEqual(float(rows[0]["measured_phi"]), float(rows[2]["measured_phi"]))
            self.assertGreater(float(rows[0]["delta"]), float(rows[2]["delta"]))
            self.assertTrue((output/"threshold_sweep.png").exists())
            self.assertTrue((output/"components.csv").exists())
            repeated = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("Output is not empty", repeated.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
