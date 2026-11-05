"""Small controls for experiment 05's paired-resolution runner."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1]))
import measurement
import run_resolution as study
from src import excursions, fields
from src.fields import GridSpec


def config(**kwargs):
    result = dict(grids=[8, 16], radii=[.85, .98], realizations=1, seed=9417,
                  length=16., lc=3., rho=1.5, reference_n=8, reference_padding=2,
                  mode_no=32, threads=1, slab=3, origin=4., include_mesh=False)
    result.update(kwargs)
    return result


class ResolutionControls(unittest.TestCase):
    def test_common_samples_match_under_refinement_and_slab_changes(self):
        c = config()
        coarse = study.sample_vector(8, c, 1)
        for n in (12, 16):
            fine = study.sample_vector(n, c, 1)
            error, count = study.common_node_error(coarse, fine)
            self.assertLessEqual(error, 1e-6)
            self.assertEqual(count, (4 if n == 12 else 8)**3)
        whole = study.sample_vector(8, config(slab=8), 1)
        np.testing.assert_array_equal(coarse, whole)

    def test_reference_grid_reproduces_shared_legacy_sampler(self):
        c = config(mode_no=1000)
        reference = fields.generate_vector_field(GridSpec(8, 16.), 3., 3,
                                                 np.random.default_rng(9417), padding=2)
        paired = study.sample_vector(8, c, 1)
        np.testing.assert_allclose(paired, reference, atol=2e-6, rtol=0)
        self.assertNotEqual(study.component_seeds(9417, 1), study.component_seeds(9417, 2))

    def test_distances_match_original_vector_norm(self):
        vector = np.random.default_rng(2).normal(size=(3, 5, 5, 5)).astype(np.float32)
        centres = excursions.tetrahedral_centres(1.5)
        actual = study.ball_distances(vector, centres)
        expected = np.stack([np.linalg.norm(vector-c[:, None, None, None], axis=0) for c in centres])
        np.testing.assert_allclose(actual, expected, atol=1e-14, rtol=0)

    def test_contact_counts_are_undirected_and_connectivity_specific(self):
        phase = np.ones((2, 2, 2), bool)
        labels = np.indices(phase.shape)[0]
        result = study.cross_label_contacts(phase, labels)
        self.assertEqual(result["cross_label_contacts_6"], 4)
        self.assertEqual(result["cross_label_contacts_26"], 16)
        phase[:] = False
        phase[0, 0, 0] = phase[1, 1, 1] = True
        result = study.cross_label_contacts(phase, labels)
        self.assertEqual(result["cross_label_contacts_6"], 0)
        self.assertEqual(result["cross_label_contacts_26"], 1)

    def test_topology_only_retains_voxel_values_without_mesh_calls(self):
        x = np.indices((4,)*3)[0]
        distances = np.stack([np.where(x < 2, 0., 2.), np.where(x >= 2, 0., 2.)])
        with patch.object(measurement, "measure_l1_from_field", side_effect=AssertionError("mesh called")), \
             patch.object(measurement, "measure_surface_area_from_field", side_effect=AssertionError("mesh called")):
            row = study.measure_resolution(distances, 1., GridSpec(4, 4.), include_mesh=False)
        self.assertEqual(row["union_euler_6"], 1)
        self.assertEqual(row["sum_ball_euler_6"], 2)
        self.assertEqual(row["euler_additivity_gap_6"], -1)
        self.assertEqual(row["cross_label_contacts_6"], 16)
        self.assertEqual(row["union_L3"], row["sum_ball_L3"])
        self.assertTrue(np.isnan(row["union_L1_field"]))
        self.assertTrue(np.isnan(row["sum_ball_L2_field"]))

    def test_checkpoint_resume_and_parameter_protection(self):
        c = config(grids=[4, 8], slab=4)
        real_measurement = study.measure_resolution
        calls = 0
        def interrupt_second(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise KeyboardInterrupt("test interruption")
            return real_measurement(*args, **kwargs)
        with tempfile.TemporaryDirectory(prefix="rf_resolution_test_") as directory:
            output = Path(directory)/"run"
            with patch.object(study, "measure_resolution", side_effect=interrupt_second):
                with self.assertRaises(KeyboardInterrupt):
                    study.run(c, output)
            metadata = json.loads((output/"metadata.json").read_text())
            self.assertEqual(metadata["status"], "interrupted")
            self.assertEqual(metadata["completed_cases"], 1)
            with self.assertRaisesRegex(ValueError, "not empty"):
                study.run(c, output)
            with self.assertRaisesRegex(ValueError, "identical parameters"):
                study.run(c | {"mode_no": 64}, output, resume=True)
            with patch.object(study, "measure_resolution", wraps=real_measurement) as measured:
                study.run(c, output, resume=True)
                self.assertEqual(measured.call_count, 3)
            rows = study.read_csv(output/"measurements.csv")
            self.assertEqual(len(rows), 4)
            self.assertEqual(len({(r["N"], r["radius_over_rc"]) for r in rows}), 4)
            metadata = json.loads((output/"metadata.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["expected_cases"], 4)
            self.assertTrue((output/"resolution.png").exists())
            self.assertTrue((output/"summary.md").exists())
            before = (output/"measurements.csv").read_bytes()
            with patch.object(study, "sample_vector", side_effect=AssertionError("unnecessary sampling")):
                study.run(c, output, resume=True)
            self.assertEqual(before, (output/"measurements.csv").read_bytes())

    def test_cli_small_full_mesh_run_and_precontact_contract(self):
        with tempfile.TemporaryDirectory(prefix="rf_resolution_cli_") as directory:
            command = [sys.executable, "-B", str(HERE/"run_resolution.py"),
                       "--grids", "4", "8", "--radii", ".85", ".98", "--mode-no", "16",
                       "--slab", "3", "--output", str(Path(directory)/"full")]
            result = subprocess.run(command, capture_output=True, text=True, timeout=120,
                                    env=dict(os.environ, PYTHONIOENCODING="utf-8"))
            self.assertEqual(result.returncode, 0, result.stdout+"\n"+result.stderr)
            rows = study.read_csv(Path(directory)/"full"/"measurements.csv")
            self.assertEqual(len(rows), 4)
            self.assertTrue(all(np.isfinite(r["union_L1_field"]) for r in rows))
            invalid = subprocess.run(command+["--radii", "1"], capture_output=True, text=True, timeout=30)
            self.assertNotEqual(invalid.returncode, 0)
            self.assertIn("0<r/rc<1", invalid.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
