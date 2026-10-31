"""Analytical controls for the interpolated geometry and GSTools field streams."""
from __future__ import annotations

import tempfile
from pathlib import Path
import unittest

import numpy as np
import pyvista as pv

from gstools_fields import sample_components
from interpolated_geometry import (clip_linear, clip_winner, measure_region,
                                   measure_union, tetrahedral_grid)


def add_fields(mesh, prefix, arrays):
    names = []
    for i, array in enumerate(arrays):
        name = f"{prefix}_{i}"
        mesh.point_data[name] = np.broadcast_to(array, (mesh.n_points,)).copy()
        names.append(name)
    return names


class InterpolatedControls(unittest.TestCase):
    def test_full_domain_geometry(self):
        mesh = tetrahedral_grid(2, 3.0)
        stats, _, _ = measure_region(mesh, 3.0)
        self.assertEqual(mesh.n_cells, 48)
        self.assertEqual(stats["L0"], 1)
        self.assertAlmostEqual(stats["L3"], 27.)
        self.assertEqual(tuple(mesh.bounds), (0., 3., 0., 3., 0., 3.))

    def test_plane_clip_volume_and_boundary(self):
        mesh = tetrahedral_grid(2, 1.0)
        clipped = clip_linear(mesh, mesh.points[:, 0], .37)
        stats, _, _ = measure_region(clipped, 1.0)
        self.assertAlmostEqual(stats["L3"], .63, places=10)
        self.assertEqual(stats["L0"], 1)

    def test_missing_film_between_equal_positive_margin_samples(self):
        mesh = tetrahedral_grid(1, 2.0)
        x = mesh.points[:, 0]-1
        names = add_fields(mesh, "a", [x, -x])
        summary, components, *_ = measure_union(mesh, [names], [1.0], 2.0)
        self.assertEqual(summary["L0"], 2)
        self.assertEqual(summary["components"], 2)
        self.assertAlmostEqual(summary["L3"], 4.0, places=10)
        self.assertAlmostEqual(summary["largest_fraction"], .5)
        # Interpolating the endpoint margin gives a falsely full cube.
        self.assertTrue(np.all(np.abs(2*x) >= 1.0))

    def test_winner_can_exist_without_winning_at_any_original_vertex(self):
        mesh = tetrahedral_grid(1, 1.0)
        x = mesh.points[:, 0]
        names = add_fields(mesh, "a", [0.0, 2*x-1.5, .5-2*x])
        clipped = clip_winner(mesh, names, 0, .2)
        stats, _, _ = measure_region(clipped, 1.0)
        self.assertAlmostEqual(stats["L3"], .3, places=10)
        self.assertEqual(stats["L0"], 1)

    def test_thin_film_is_not_welded_away(self):
        mesh = tetrahedral_grid(1, 1.0)
        x = mesh.points[:, 0]-.5
        names = add_fields(mesh, "a", [x, -x])
        summary, *_ = measure_union(mesh, [names], [2e-11], 1.0)
        self.assertEqual(summary["components"], 2)
        self.assertEqual(summary["L0"], 2)
        self.assertAlmostEqual(summary["L3"], 1-2e-11, places=12)

    def test_cross_class_bridge_and_triple_overlap(self):
        mesh = tetrahedral_grid(2, 1.0)
        x, y, z = mesh.points.T
        groups = [add_fields(mesh, "a", [x-.5, .5-x]),
                  add_fields(mesh, "b", [y, 0.0]), add_fields(mesh, "c", [z, 0.0])]
        pair, components, _, _, added = measure_union(mesh, groups[:2], [.5, .8], 1.0)
        self.assertEqual(pair["L0"], 1)
        self.assertEqual(pair["components"], 1)
        self.assertAlmostEqual(pair["L3"], .6, places=9)
        self.assertTrue(pair["spans"])
        np.testing.assert_allclose(added, [.5, .1], atol=1e-9)
        triple, components, terms, _, added = measure_union(mesh, groups, [.5, .8, .8], 1.0)
        self.assertEqual(triple["L0"], 1)
        self.assertEqual(triple["components"], 1)
        self.assertAlmostEqual(triple["L3"], .68, places=9)
        self.assertTrue(any(t["order"] == 3 for t in terms))
        np.testing.assert_allclose(added, [.5, .1, .08], atol=1e-9)
        reverse, *_ = measure_union(mesh, groups[::-1], [.8, .8, .5], 1.0)
        self.assertEqual(reverse["L0"], triple["L0"])
        self.assertAlmostEqual(reverse["L3"], triple["L3"], places=9)

    def test_tunnel_euler_is_not_component_count(self):
        mesh = tetrahedral_grid(3, 3.0)
        centers = mesh.cell_centers().points
        mask = ~((centers[:, 0] > 1) & (centers[:, 0] < 2)
                 & (centers[:, 1] > 1) & (centers[:, 1] < 2))
        stats, _, _ = measure_region(mesh.extract_cells(mask), 3.0)
        self.assertEqual(stats["L0"], 0)
        self.assertEqual(stats["components"], 1)
        self.assertAlmostEqual(stats["L3"], 24.)

    def test_union_volume_against_independent_point_integration(self):
        rng = np.random.default_rng(138)
        coefficients = rng.normal(size=(3, 3, 4))
        mesh = tetrahedral_grid(1, 1.0)
        groups = [add_fields(mesh, str(k), [mesh.points @ c[:3] + c[3] for c in group])
                  for k, group in enumerate(coefficients)]
        deltas = [.3, .4, .5]
        summary, _, _, _, added = measure_union(mesh, groups, deltas, 1.0)
        points = rng.random((200_000, 3))
        mask = np.zeros(len(points), bool)
        fractions = []
        for group, delta in zip(coefficients, deltas):
            values = np.sort(np.stack([points @ c[:3]+c[3] for c in group]), axis=0)
            phase = values[-1]-values[-2] >= delta
            fractions.append(np.mean(phase & ~mask))
            mask |= phase
        self.assertAlmostEqual(summary["L3"], float(mask.mean()), delta=.004)
        np.testing.assert_allclose(added, fractions, atol=.004, rtol=0)

    def test_sixteen_competitors_with_inactive_regions(self):
        mesh = tetrahedral_grid(1, 1.0)
        x = mesh.points[:, 0]-.5
        names = add_fields(mesh, "a", [x, -x]+[-10.-j for j in range(14)])
        summary, *_ = measure_union(mesh, [names], [.4], 1.0)
        self.assertEqual(summary["L0"], 2)
        self.assertEqual(summary["components"], 2)
        self.assertAlmostEqual(summary["L3"], .6, places=10)

    def test_empty_and_full_margin_union(self):
        mesh = tetrahedral_grid(1, 1.0)
        names = add_fields(mesh, "a", [2., 0.])
        full, *_ = measure_union(mesh, [names], [.5], 1.0)
        empty, *_ = measure_union(mesh, [names], [3.], 1.0)
        self.assertEqual(full["L0"], 1)
        self.assertAlmostEqual(full["L3"], 1.)
        self.assertEqual(empty["L0"], 0)
        self.assertEqual(empty["L3"], 0.)

    def test_exported_volume_matches_measured_volume(self):
        mesh = tetrahedral_grid(1, 1.0)
        names = add_fields(mesh, "a", [mesh.points[:, 0], 0.])
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            out = Path(tmp)
            summary, *_ = measure_union(mesh, [names], [.4], 1.0, export_directory=out)
            exported = pv.read(out / "class_0_winner_0.vtu")
            stats, _, _ = measure_region(exported, 1.)
            self.assertAlmostEqual(stats["L3"], summary["L3"])
            self.assertEqual(stats["L0"], summary["L0"])

    def test_gstools_paired_refinement(self):
        classes = [{"lc": 1.5, "q": 2}]
        coarse, diag = sample_components(2, 2., classes, 123, mode_no=32)
        fine, _ = sample_components(4, 2., classes, 123, mode_no=32)
        for a, b in zip(coarse[0], fine[0]):
            np.testing.assert_allclose(a, b[::2, ::2, ::2], atol=1e-12, rtol=0)
        self.assertNotEqual(diag[0]["seed"], diag[1]["seed"])
        self.assertGreater(abs(coarse[0][0].mean()), 1e-6)


if __name__ == "__main__":
    unittest.main(verbosity=2)
