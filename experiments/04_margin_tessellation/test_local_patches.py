"""Analytic geometry and tiny-scene controls for the local 3D patch study."""
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
import run_local_patches as patches
import refine_patch_components as input_refinement


class LocalPatchControls(unittest.TestCase):
    def test_refining_component_inputs_recovers_hidden_cavity(self):
        reference_points = patches.grid_points(16)
        reference = (np.sum((reference_points-.5)**2, axis=1) >= .1).reshape((33,)*3)
        errors = []
        for factor in (1, 2, 4):
            points = patches.grid_points(factor)
            coarse = np.array([np.sum((points-.5)**2, axis=1), np.zeros(len(points))])
            coarse = coarse.reshape((2,)+(2*factor+1,)*3)
            arrays = input_refinement.component_reconstruction(coarse, 16, factor)
            phase = arrays["trilinear_components"][0] >= .1
            errors.append(patches.fraction_integral(phase != reference))
            self.assertEqual(patches.topology(phase, 26)["L0"], 1 if factor == 1 else 2)
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])

    def test_affine_reproduction_grid_nodes_faces_and_outer_boundary(self):
        x, y, z = np.indices((3,)*3).astype(float)
        coarse = np.array([2*x-3*y+z+1, -x+.25*y+2*z])
        points = np.vstack((patches.COARSE_POINTS,
                            np.random.default_rng(31).uniform(0, 2, size=(200, 3)),
                            [[.4, .4, .4], [1., .2, 1.7], [2., 2., 2.]]))
        expected = np.array([2*points[:, 0]-3*points[:, 1]+points[:, 2]+1,
                             -points[:, 0]+.25*points[:, 1]+2*points[:, 2]])
        for method in ("trilinear", "tetra"):
            np.testing.assert_allclose(patches.interpolate(coarse, points, method), expected, atol=2e-15)
            np.testing.assert_array_equal(patches.interpolate(coarse, patches.COARSE_POINTS, method),
                                          coarse.reshape(2, -1))

    def test_trilinear_cross_term_and_tetrahedral_continuity(self):
        x, y, z = np.indices((3,)*3).astype(float)
        coarse = (x*y)[None]
        query = np.array([[.75, .25, .5]])
        self.assertAlmostEqual(patches.interpolate(coarse, query, "trilinear")[0, 0], .1875)
        self.assertAlmostEqual(patches.interpolate(coarse, query, "tetra")[0, 0], .25)
        rng = np.random.default_rng(19)
        arbitrary = rng.normal(size=(4, 3, 3, 3))
        for query in (np.array([1., .35, .6]), np.array([.4, .4, .7])):
            sides = np.array([query-np.array([1e-10, 0, 0]), query+np.array([1e-10, 0, 0])])
            result = patches.interpolate(arbitrary, sides, "tetra")
            np.testing.assert_allclose(result[:, 0], result[:, 1], atol=2e-9)

    def test_interior_winner_absent_from_cell_vertices(self):
        x = np.indices((2,)*3)[0].astype(float)
        coarse = np.array([np.zeros_like(x), -1.5+2*x, .5-2*x])
        self.assertNotIn(0, np.argmax(coarse, axis=0))
        for method in ("trilinear", "tetra"):
            values = patches.interpolate(coarse, np.array([[.5, .5, .5]]), method)
            margin, winner = patches.films.margin_and_winner(values)
            self.assertEqual(int(winner[0]), 0)
            self.assertGreater(float(margin[0]), .2)

    def test_curved_field_cavity_changes_euler_without_changing_face_connections(self):
        coarse_points = patches.COARSE_POINTS
        coarse = np.array([np.sum((coarse_points-.5)**2, axis=1), np.zeros(27)]).reshape(2, 3, 3, 3)
        points = patches.grid_points(8)
        direct = (np.sum((points-.5)**2, axis=1) >= .1).reshape((17,)*3)
        for method in ("trilinear", "tetra"):
            m, w = patches.films.margin_and_winner(patches.interpolate(coarse, points, method))
            reconstructed = (m >= .1).reshape(direct.shape)
            self.assertTrue(reconstructed.all())
            for connectivity in (6, 26):
                a = patches.labelled_topology(direct, np.zeros_like(direct, dtype=int), connectivity)
                b = patches.labelled_topology(reconstructed, w.reshape(direct.shape), connectivity)
                self.assertEqual(a[0]["L0"], 2)
                self.assertEqual(b[0]["L0"], 1)
                self.assertEqual(a[0]["components"], b[0]["components"])
                self.assertEqual(a[0]["face_pairs"], b[0]["face_pairs"])
                self.assertFalse(patches.same_topology(a, b))

    def test_face_signature_and_diagonal_connectivity(self):
        xbar = np.zeros((5,)*3, dtype=bool)
        xbar[:, 2, 2] = True
        ybar = xbar.swapaxes(0, 1)
        x = patches.topology(xbar, 6)
        y = patches.topology(ybar, 6)
        self.assertEqual(x["L0"], y["L0"])
        self.assertEqual(x["components"], y["components"])
        self.assertNotEqual(x["face_pairs"], y["face_pairs"])
        self.assertNotEqual(patches.signature(x), patches.signature(y))
        diagonal = np.zeros((3,)*3, dtype=bool)
        diagonal[np.arange(3), np.arange(3), np.arange(3)] = True
        self.assertEqual(patches.topology(diagonal, 6)["components"], 3)
        self.assertEqual(patches.topology(diagonal, 26)["components"], 1)

    def test_fixed_domain_quadrature_and_nested_reuse(self):
        for n in (3, 5, 9):
            self.assertAlmostEqual(patches.fraction_integral(np.ones((n,)*3)), 1.)
            self.assertAlmostEqual(patches.fraction_integral(np.indices((n,)*3)[0]/(n-1)), .5)

        class Sampler:
            def __init__(self):
                self.count = 0

            def evaluate(self, points):
                self.count += len(points)
                return np.array([points[:, 0]+2, -points[:, 1]])

        sampler = Sampler()
        coarse = sampler.evaluate(patches.COARSE_POINTS)
        m, w = patches.films.margin_and_winner(coarse)
        m, w = m.reshape((3,)*3), w.reshape((3,)*3)
        old = m.copy()
        m, w, added = patches.refine_direct(sampler, np.zeros(3), m, w, 2)
        self.assertEqual(added, 5**3-3**3)
        np.testing.assert_array_equal(m[::2, ::2, ::2], old)
        m, w, added = patches.refine_direct(sampler, np.zeros(3), m, w, 4)
        self.assertEqual(sampler.count, 9**3)
        expected = patches.grid_points(4).sum(axis=1)  # z does not enter this field
        expected -= patches.grid_points(4)[:, 2]
        np.testing.assert_allclose(m.ravel(), expected+2)

    def test_selection_categories_are_disjoint_and_match_coarse_data(self):
        shape = (32,)*3
        xyz = np.indices(shape)
        margin = np.full(shape, .2)
        winners = np.where(xyz[0] < 12, xyz[1] % 3, 0)
        rows = [dict(kind="cross_label", target_phi=.85, i=4, j=4, k=4,
                     direct_gap_fraction=.5, affine_phase_disagreement_fraction=.09, segment_id=0),
                dict(kind="cross_label", target_phi=.85, i=8, j=20, k=20,
                     direct_gap_fraction=.5, affine_phase_disagreement_fraction=.01, segment_id=1)]
        selected, info = patches.select_patches(margin, winners, rows, .1, candidates=8192)
        self.assertEqual(len(selected), 6)
        self.assertGreater(info["junction_candidates"], 0)
        for a, row in enumerate(selected):
            lower = np.array(row["lower"])
            for other in selected[:a]:
                self.assertFalse(np.all(np.abs(lower-np.array(other["lower"])) <= 2))
            sl = tuple(slice(v, v+3) for v in lower)
            if row["kind"] == "junction_candidate":
                self.assertGreaterEqual(len(np.unique(winners[sl])), 3)
            elif row["kind"] == "same_label_near_threshold":
                self.assertEqual(len(np.unique(winners[sl])), 1)
                self.assertTrue(np.all(margin[sl] >= .1))

    def test_tiny_cli_provenance_volumes_and_overwrite_guard(self):
        with tempfile.TemporaryDirectory(prefix="rf_patches_") as directory:
            source, output = Path(directory)/"scene", Path(directory)/"output"
            source.mkdir()
            sampler = patches.films.OriginalFields(3, 1., .5, 321729, mode_no=16)
            points = np.indices((6,)*3).reshape(3, -1).T
            values = sampler.evaluate(points).astype(np.float32)
            margin, winners = patches.films.margin_and_winner(values)
            phase = margin >= patches.saved.adler.delta_for_vf(3, .75)
            for filename, name, array in (("margin.vti", "margin", margin),
                                          ("winners.vti", "winner_label", winners),
                                          ("phase.vti", "phase", phase.astype(np.uint8)),
                                          ("field_g1.vti", "G_1", values[0])):
                mesh = pv.ImageData(dimensions=(6,)*3, spacing=(.5,)*3, origin=(2., 3., 4.))
                mesh.point_data[name] = array.reshape((6,)*3).ravel(order="F")
                mesh.save(source/filename)
            before = {p.name: patches.saved.sha256(p) for p in source.iterdir()}
            fixed = Path(directory)/"patches.json"
            fixed.write_text(json.dumps([dict(patch_id=5, kind="test", lower=[1, 1, 1], cells_per_axis=2)]))
            command = [sys.executable, "-B", str(HERE/"run_local_patches.py"),
                       "--source", str(source), "--output", str(output), "--patches", str(fixed),
                       "--q", "3", "--lc", "1", "--mode-no", "16", "--refinements", "1", "2", "4"]
            result = subprocess.run(command, capture_output=True, text=True, timeout=120,
                                    env=dict(os.environ, PYTHONIOENCODING="utf-8"))
            self.assertEqual(result.returncode, 0, result.stdout+"\n"+result.stderr)
            metadata = json.loads((output/"metadata.json").read_text())
            self.assertEqual(metadata["status"], "complete")
            self.assertEqual(metadata["seed_validation"]["margin_max_abs"], 0.)
            self.assertEqual(metadata["completed_patches"], 1)
            rows = patches.read_csv(output/"measurements.csv")
            self.assertEqual(len(rows), 48)
            for row in rows:
                self.assertGreaterEqual(float(row["local_phi"]), 0.)
                self.assertLessEqual(float(row["local_phi"]), 1.)
                self.assertAlmostEqual(float(row["local_L3_quadrature"]), float(row["local_phi"]))
                if row["method"] == "direct":
                    self.assertEqual(float(row["phase_disagreement_phi"]), 0.)
                    self.assertEqual(row["matches_direct_topology"], "True")
            mesh = pv.read(output/"patch_05.vti")
            np.testing.assert_allclose(mesh.bounds, [2.5, 3.5, 3.5, 4.5, 4.5, 5.5])
            with np.load(output/"patch_05.npz") as arrays:
                np.testing.assert_array_equal(mesh.point_data["direct_margin"].reshape((9,)*3, order="F"),
                                              arrays["direct_margin"])
            self.assertTrue((output/"patch_05_phi85.png").exists())
            refined = Path(directory)/"refined"
            refine_command = [sys.executable, "-B", str(HERE/"refine_patch_components.py"),
                              "--source", str(output), "--output", str(refined), "--input-factors", "2"]
            result = subprocess.run(refine_command, capture_output=True, text=True, timeout=120,
                                    env=dict(os.environ, PYTHONIOENCODING="utf-8"))
            self.assertEqual(result.returncode, 0, result.stdout+"\n"+result.stderr)
            info = json.loads((refined/"metadata.json").read_text())
            self.assertEqual(info["status"], "complete")
            self.assertEqual(len(patches.read_csv(refined/"measurements.csv")), 16)
            self.assertTrue((refined/"input_refinement.png").exists())
            self.assertEqual({p.name: patches.saved.sha256(p) for p in source.iterdir()}, before)
            repeated = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn("Output is not empty", repeated.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
