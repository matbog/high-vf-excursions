"""Small controls for the geometry sweep; no full N=100 study."""
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
from scipy.stats import ncx2

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_geometry as geometry
from src import excursions
from src.fields import GridSpec


def config(**kwargs):
    values = dict(n=8, rhos=[1.5, 2.5], radii=[.85, .98], realizations=2,
                  seed=9417, length=16., lc=3., reference_n=8, reference_padding=2,
                  origin=4., mode_no=16, threads=1, slab=4, theory_only=False)
    values.update(kwargs)
    return values


class GeometrySweepControls(unittest.TestCase):
    def test_cavity_euler_and_component_boundary_flags(self):
        mask = np.ones((5,)*3, dtype=bool)
        mask[2, 2, 2] = False
        for rank in (6, 26):
            stats, rows = geometry.component_statistics(mask, .5, rank)
            self.assertEqual(stats["components"], 1)
            self.assertEqual(stats["euler"], 2)
            self.assertTrue(stats["spans"])
            self.assertTrue(rows[0]["touches_boundary"])
            self.assertEqual(stats["volume"], 124*.5**3)
            self.assertEqual(rows[0]["volume"], stats["volume"])
            self.assertAlmostEqual(rows[0]["equivalent_diameter"], (6*stats["volume"]/np.pi)**(1/3))

    def test_diagonal_has_different_spanning_under_both_conventions(self):
        mask = np.zeros((3,)*3, dtype=bool)
        mask[np.arange(3), np.arange(3), np.arange(3)] = True
        six, rows6 = geometry.component_statistics(mask, 1., 6)
        twenty_six, rows26 = geometry.component_statistics(mask, 1., 26)
        self.assertEqual(six["components"], 3)
        self.assertFalse(six["spans"])
        self.assertEqual(twenty_six["components"], 1)
        self.assertTrue(twenty_six["spans"])
        self.assertEqual(sum(r["volume"] for r in rows6), rows26[0]["volume"])

    def test_distinct_balls_are_not_merged_for_particle_statistics(self):
        x = np.indices((4,)*3)[0]
        distances = np.stack([np.where(x < 2, 0., 2.), np.where(x >= 2, 0., 2.)])
        row, balls, components = geometry.measure_geometry(distances, 1., GridSpec(4, 4.))
        for rank in (6, 26):
            self.assertEqual(row[f"union_components_{rank}"], 1)
            self.assertEqual(row[f"sum_ball_components_{rank}"], 2)
            self.assertEqual(row[f"sum_ball_L0_{rank}"], 2)
            self.assertEqual(row[f"euler_additivity_gap_{rank}"], -1)
            self.assertEqual(row[f"largest_separate_component_fraction_{rank}"], .5)
            selected = [r for r in components if r["connectivity"] == rank]
            self.assertEqual(sum(r["volume"] for r in selected), row["union_L3"])
            for ball in (0, 1):
                detail = next(r for r in balls if r["connectivity"] == rank and r["ball"] == ball)
                self.assertEqual(detail["volume"], 32.)
        self.assertEqual(row["cross_label_contacts_6"], 16)

    def test_empty_masks_and_invalid_overlap(self):
        distances = np.full((4, 4, 4, 4), 10.)
        row, balls, components = geometry.measure_geometry(distances, 1., GridSpec(4, 4.))
        self.assertEqual(row["measured_phi"], 0)
        self.assertEqual(row["largest_separate_component_fraction_6"], 0)
        self.assertFalse(row["any_ball_spans_26"])
        self.assertEqual(len(components), 0)
        self.assertEqual(len(balls), 8)
        with self.assertRaisesRegex(ValueError, "disjoint"):
            geometry.measure_geometry(np.zeros_like(distances), 1., GridSpec(4, 4.))

    def test_matches_resolution_runner_topology_on_same_field(self):
        c = config(realizations=1)
        vector = geometry.resolution.sample_vector(c["n"], c, 1)
        centres = excursions.tetrahedral_centres(1.5)
        distances = geometry.resolution.ball_distances(vector, centres)
        radius = .98*np.linalg.norm(centres[0]-centres[1])/2
        grid = GridSpec(c["n"], c["length"])
        current, _, _ = geometry.measure_geometry(distances, radius, grid)
        previous = geometry.resolution.measure_resolution(distances, radius, grid, include_mesh=False)
        self.assertEqual(current["measured_phi"], previous["measured_phi"])
        for rank in (6, 26):
            for key in ("union_euler", "union_components", "union_spans", "any_ball_spans", "cross_label_contacts"):
                self.assertEqual(current[f"{key}_{rank}"], previous[f"{key}_{rank}"])
            self.assertEqual(current[f"sum_ball_L0_{rank}"], previous[f"sum_ball_euler_{rank}"])

    def test_theory_only_uses_exact_disjoint_ball_probability(self):
        c = config(theory_only=True)
        with tempfile.TemporaryDirectory(prefix="rf_geometry_theory_") as directory:
            output = Path(directory)/"theory"
            with patch.object(geometry.resolution, "sample_vector", side_effect=AssertionError("field generated")):
                geometry.run(c, output)
            rows = geometry.resolution.read_csv(output/"theory.csv")
            self.assertEqual(len(rows), 4)
            for row in rows:
                expected = 4*ncx2.cdf(row["radius"]**2, df=3, nc=row["rho"]**2)
                self.assertAlmostEqual(row["theory_phi"], expected, places=12)
            self.assertEqual(json.loads((output/"metadata.json").read_text())["status"], "complete")
            self.assertFalse((output/"cases").exists())

    def test_field_reuse_accounting_and_complete_resume(self):
        c = config()
        with tempfile.TemporaryDirectory(prefix="rf_geometry_reuse_") as directory:
            output = Path(directory)/"run"
            with patch.object(geometry.resolution, "sample_vector", wraps=geometry.resolution.sample_vector) as sample:
                geometry.run(c, output)
                self.assertEqual(sample.call_count, 2)
            rows = geometry.resolution.read_csv(output/"measurements.csv")
            self.assertEqual(len(rows), 8)
            self.assertEqual(len(list((output/"cases").glob("*.json"))), 8)
            balls = geometry.resolution.read_csv(output/"balls.csv")
            self.assertEqual(len(balls), 64)
            for row in rows:
                for rank in (6, 26):
                    selected = [r for r in balls if geometry.case_key(r) == geometry.case_key(row)
                                and r["connectivity"] == rank]
                    self.assertEqual(len(selected), 4)
                    self.assertAlmostEqual(sum(r["volume"] for r in selected), row["union_L3"])
            with self.assertRaisesRegex(ValueError, "not empty"):
                geometry.run(c, output)
            with self.assertRaisesRegex(ValueError, "matching"):
                geometry.run(c | {"n": 16}, output, resume=True)
            # CSVs are derived views: resume must rebuild them from complete case files.
            (output/"measurements.csv").write_text("stale\n")
            with patch.object(geometry.resolution, "sample_vector", side_effect=AssertionError("field regenerated")):
                geometry.run(c, output, resume=True)
            self.assertEqual(len(geometry.resolution.read_csv(output/"measurements.csv")), 8)
            self.assertEqual(json.loads((output/"metadata.json").read_text())["status"], "complete")

    def test_resume_recovers_when_report_write_fails_after_atomic_case(self):
        c = config(realizations=1)
        original = geometry.write_outputs
        calls = 0
        def fail_after_case(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise OSError("simulated interrupted report")
            return original(*args, **kwargs)
        with tempfile.TemporaryDirectory(prefix="rf_geometry_resume_") as directory:
            output = Path(directory)/"run"
            with patch.object(geometry, "write_outputs", side_effect=fail_after_case):
                with self.assertRaisesRegex(OSError, "simulated"):
                    geometry.run(c, output)
            checkpoints = list((output/"cases").glob("*.json"))
            self.assertEqual(len(checkpoints), 1)
            before = checkpoints[0].read_bytes()
            with patch.object(geometry, "measure_geometry", wraps=geometry.measure_geometry) as measure:
                geometry.run(c, output, resume=True)
                self.assertEqual(measure.call_count, 3)
            self.assertEqual(checkpoints[0].read_bytes(), before)
            self.assertEqual(len(geometry.resolution.read_csv(output/"measurements.csv")), 4)

    def test_cli_and_strict_precontact_contract(self):
        with tempfile.TemporaryDirectory(prefix="rf_geometry_cli_") as directory:
            command = [sys.executable, "-B", str(HERE/"run_geometry.py"), "--n", "4",
                       "--realizations", "1", "--rhos", "1.5", "2.5", "--radii", ".98",
                       "--mode-no", "8", "--output", str(Path(directory)/"run")]
            result = subprocess.run(command, capture_output=True, text=True, timeout=120,
                                    env=dict(os.environ, PYTHONIOENCODING="utf-8"))
            self.assertEqual(result.returncode, 0, result.stdout+"\n"+result.stderr)
            metadata = json.loads((Path(directory)/"run"/"metadata.json").read_text())
            self.assertEqual(metadata["completed_cases"], 2)
            self.assertEqual(metadata["status"], "complete")
            self.assertTrue((Path(directory)/"run"/"geometry_summary.png").exists())
            bad = subprocess.run(command+["--radii", "1.0"], capture_output=True, text=True, timeout=30)
            self.assertNotEqual(bad.returncode, 0)
            self.assertIn("0<r/rc<1", bad.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
