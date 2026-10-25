import numpy as np

from bounded_lkcs import cubical_lkcs, mesh_l1_l2


def test_single_voxel() -> None:
    phase = np.ones((1, 1, 1), dtype=bool)
    assert np.allclose(cubical_lkcs(phase).as_array(), [1, 3, 3, 1])


def test_solid_rectangular_box() -> None:
    phase = np.ones((2, 3, 4), dtype=bool)
    # Box side lengths 2,3,4: L3=24, L2=ab+ac+bc=26,
    # L1=a+b+c=9 under the intrinsic-volume convention, L0=1.
    assert np.allclose(cubical_lkcs(phase).as_array(), [1, 9, 26, 24])


def test_mesh_sphere_is_close() -> None:
    coordinates = np.indices((48, 48, 48), dtype=float)
    squared_radius = sum((coordinates[axis] - 23.5) ** 2 for axis in range(3))
    phase = squared_radius <= 15.0**2
    l1, l2 = mesh_l1_l2(phase)
    assert np.isclose(l1, 4.0 * 15.0, rtol=0.15)
    assert np.isclose(l2, 2.0 * np.pi * 15.0**2, rtol=0.12)
