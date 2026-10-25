import numpy as np
from scipy.special import ndtr

from gaussian_minkowski import margin_gmfs, margin_m0, margin_m1, normal_density


def test_q2_closed_forms() -> None:
    for delta in (0.1, 0.5, 1.0, 2.0):
        assert np.isclose(margin_m0(2, delta), 2.0 * ndtr(-delta / np.sqrt(2.0)))
        assert np.isclose(margin_m1(2, delta), 2.0 * normal_density(delta / np.sqrt(2.0)))


def test_m1_is_margin_derivative() -> None:
    for q in (2, 4, 8):
        delta = 0.4
        step = 1e-4
        derivative = (margin_m0(q, delta + step) - margin_m0(q, delta - step)) / (2 * step)
        assert np.isclose(margin_m1(q, delta), -np.sqrt(2.0) * derivative, rtol=2e-7)


def test_all_q2_gmfs() -> None:
    for delta in (0.2, 0.7, 1.2):
        a = delta / np.sqrt(2.0)
        phi = normal_density(a)
        exact = np.array(
            [2.0 * ndtr(-a), 2.0 * phi, 2.0 * a * phi, 2.0 * (a * a - 1.0) * phi]
        )
        assert np.allclose(margin_gmfs(2, delta), exact, rtol=2e-10, atol=2e-11)
