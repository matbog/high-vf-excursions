"""Small analytical and geometric controls; no experiment outputs are written."""
from __future__ import annotations

import itertools
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src import adler
from src.fields import GridSpec
from model import (GaussianSampler, class_probabilities, component_table,
                   cube_lkcs, intersection_operator, margin_from_fields,
                   topology, union_expectations)


class MultiscaleControls(unittest.TestCase):
    def test_added_volume_calibration(self):
        shares = np.array([.40, .35, .25])
        for phi in (.01, .60, .75, .85, .99):
            raw = class_probabilities(phi, shares)
            free = 1.0
            added = []
            for probability in raw:
                added.append(free * probability)
                free *= 1 - probability
            np.testing.assert_allclose(added, phi * shares, rtol=1e-13)
            self.assertAlmostEqual(1 - free, phi, places=13)
        with self.assertRaises(ValueError):
            class_probabilities(.8, [.4, .3])
        with self.assertRaises(ValueError):
            class_probabilities(1., shares)

    def test_independent_volume_and_priority_accounting(self):
        phi = .8
        shares = [.4, .35, .25]
        raw = class_probabilities(phi, shares)
        uniforms = np.random.default_rng(184).random((250_000, 3))
        remaining = np.ones(len(uniforms), dtype=bool)
        added = []
        for k in range(3):
            mask = uniforms[:, k] < raw[k]
            added.append(np.mean(mask & remaining))
            remaining &= ~mask
        np.testing.assert_allclose(added, phi * np.array(shares), atol=.003, rtol=0)

    def test_single_scale_matches_existing_normalization(self):
        length, lc, q, delta = 40., 5., 8, .4
        actual = intersection_operator(adler.margin_gmfs(q, delta), 1/(lc*np.sqrt(np.pi))) @ cube_lkcs(length)
        expected = adler.expected_lkcs_margin(GridSpec(40, length), lc, q, delta)
        np.testing.assert_allclose(actual, [expected[f"L{j}"] for j in range(4)])

    def test_q2_intersection_against_orthant_tube_coefficients(self):
        # Two independent |H| >= a excursions: four orthants in R^2.
        # Separate facet and right-angle corner integration supplies this control.
        a, lc, length = .6, 4., 32.
        p, Q = norm.pdf(a), norm.sf(a)
        orthant_gmfs = [4*Q**2, 8*p*Q, 8*a*p*Q + 2*np.pi*p**2,
                       8*(a*a-1)*p*Q + 16*a*p**2]
        rho = 1/(lc*np.sqrt(np.pi))
        single = intersection_operator(adler.margin_gmfs(2, np.sqrt(2)*a), rho)
        intersection = single @ single @ cube_lkcs(length)
        expected = intersection_operator(orthant_gmfs, rho) @ cube_lkcs(length)
        np.testing.assert_allclose(intersection, expected, rtol=1e-11, atol=1e-11)

    def test_union_is_order_independent_and_has_exact_volume(self):
        operators = [intersection_operator(adler.margin_gmfs(q, delta), rho)
                     for q, delta, rho in ((2, .4, .1), (4, .6, .2), (8, .9, .3))]
        length = 35.
        reference = union_expectations(operators, length)[-1]["union"]
        for permutation in itertools.permutations(operators):
            np.testing.assert_allclose(union_expectations(permutation, length)[-1]["union"], reference)
        probability = 1 - np.prod([1 - op[0, 0] for op in operators])
        self.assertAlmostEqual(reference[3]/length**3, probability, places=13)

    def test_full_and_empty_phase_limits(self):
        op = intersection_operator(adler.margin_gmfs(4, .5), .2)
        domain = cube_lkcs(20.)
        empty = np.zeros((4, 4))
        full = np.eye(4)
        np.testing.assert_allclose(union_expectations([op, empty], 20.)[-1]["union"], op @ domain)
        np.testing.assert_allclose(union_expectations([op, full], 20.)[-1]["union"], domain)

    def test_sampler_ensemble_variance_and_isotropic_spectrum(self):
        sampler = GaussianSampler(24, 24., 2., 4.)
        n = sampler.full_n
        power = sampler.filter.astype(float)**2
        weights = np.full(n//2+1, 2.)
        weights[0] = 1
        if n % 2 == 0:
            weights[-1] = 1
        weighted = power * weights[None, None, :]
        self.assertAlmostEqual(weighted.sum()/n**3, 1., places=6)
        k = 2*np.pi*np.fft.fftfreq(n, d=sampler.dx)
        kr = 2*np.pi*np.fft.rfftfreq(n, d=sampler.dx)
        variances = [(weighted * wave**2).sum()/weighted.sum()
                     for wave in (k[:, None, None], k[None, :, None], kr[None, None, :])]
        np.testing.assert_allclose(variances, sampler.lambda2, rtol=1e-6)
        draws = [sampler.sample(np.random.default_rng(i)) for i in range(24)]
        self.assertLess(abs(np.mean([np.mean(f*f) for f in draws])-1), .08)
        # Per-realization means must remain random: centering destroys this property.
        self.assertGreater(np.std([f.mean() for f in draws]), .005)

    def test_streamed_margin_matches_order_statistics(self):
        stack = np.random.default_rng(90).normal(size=(7, 9, 8, 6)).astype(np.float32)
        margin, winner = margin_from_fields(iter(stack))
        ordered = np.sort(stack, axis=0)
        np.testing.assert_array_equal(margin, ordered[-1]-ordered[-2])
        np.testing.assert_array_equal(winner, stack.argmax(axis=0))

    def test_components_preserve_cross_class_connections(self):
        mask = np.zeros((7, 7, 7), dtype=bool)
        mask[:, 3, 3] = True
        owner = np.zeros_like(mask, dtype=np.uint8)
        owner[:3, 3, 3] = 1
        owner[3:, 3, 3] = 2
        stats, labels, sizes = topology(mask, 6, return_labels=True)
        self.assertEqual(stats["components_6"], 1)
        self.assertTrue(stats["spans_6"])
        rows = list(component_table(labels, sizes, owner, .5, 2))
        self.assertTrue(rows[0]["mixed_classes"])
        self.assertAlmostEqual(rows[0]["volume"], 7*.5**3)
        empty = topology(np.zeros_like(mask), 26)
        self.assertEqual(empty["components_26"], 0)
        self.assertFalse(empty["spans_26"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
