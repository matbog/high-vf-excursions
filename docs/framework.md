# Framework: random field excursion sets and the Gaussian Kinematic Formula

## The excursion set morphology model

A random field `f : M -> R` assigns a random real value to every point of a spatial domain
`M` (here, a cube in R^3). The field is characterised by its marginal distribution — what
value a single point is likely to take — and its covariance function `C(r) = Cov(f(x), f(x+r))`,
which controls spatial smoothness and correlation length. In this project the covariance is
always Gaussian: `C(r) = sigma^2 * exp(-r^2 / lc^2)`, with correlation length `lc`.

The **excursion set** at threshold `u` is the binary field

```
A_u = { x in M : f(x) >= u }
```

This partitions the domain into an inclusion phase `A_u` and a matrix phase `M \ A_u`.
Varying `u` sweeps volume fraction continuously: high `u` gives sparse inclusions; low `u`
gives a dense, percolating phase. The two-parameter family `(lc, u)` covers a wide range of
microstructure morphologies within a single generative model.

The **Lipschitz-Killing curvatures** (LKCs) of `A_u` are the standard scalar morphological
descriptors in three dimensions:

| Symbol | Name | Physical meaning |
|--------|------|-----------------|
| `L3(A_u)` | Volume | Total inclusion volume |
| `L2(A_u)` | Half surface area | Specific surface (inclusion-matrix interface area / 2) |
| `L0(A_u)` | Euler characteristic | Topology: #components - #tunnels + #cavities |

`L0` is the most information-rich: positive values indicate a disconnected (island) morphology;
negative values indicate a tunnel- or percolation-dominated topology.

## The Gaussian Kinematic Formula

For a smooth, stationary, isotropic Gaussian-related random field on a 3D cube of side `a`,
the Adler-Taylor formula gives

```math
E[ L_i(A_u) ] = sum_{j=0}^{3-i} C(i,j) * rho^j * L_{i+j}(M) * M_j^gamma(D_u)
```

For `i=0` (Euler characteristic) this expands to

```math
E[ L0(A_u) ] = rho^3 * a^3 * M3 + 3 * rho^2 * a^2 * M2 + 3 * rho * a * M1 + M0
```

**Terms:**

- `rho = sqrt(lambda2 / 2*pi)` where `lambda2 = -C''(0) = 2*sigma^2 / lc^2` is the second
  spectral moment of the covariance. `rho` encodes spatial regularity and scale.
- `a^3`, `3a^2`, `3a`, `1` are the domain LKCs of a cube of side `a` (volume, surface, edge, vertex).
- `M_j^gamma(D_u)` are the **Gaussian Minkowski functionals** of the hitting set
  `D_u = { z in R : z >= u }` — the region of value space that maps to the inclusion phase.

**The central structural insight:**

The spatial geometry (encoded in `rho` and the domain LKCs) is completely separated from the
value-space geometry (encoded in `M_j`). When the marginal distribution changes — Gaussian to
Student-t, or scalar to vector — only `M_j` changes. The spatial side is identical as long as
the covariance structure is unchanged.

**Gaussian M_j (exact closed form):**

```
M0 = P(Z >= u) = (1/2) erfc(u / sqrt(2))          [volume fraction]
M1 = phi(u) = (1/sqrt(2*pi)) exp(-u^2/2)           [boundary density]
M2 = u * phi(u)                                     [mean curvature density]
M3 = (u^2 - 1) * phi(u)                             [Gaussian curvature density]
```

**The Roubin-Colliat percolation criterion:**

The upper zero-crossing of `E[L0(A_u)]` as a function of `u` corresponds to a transition from
positive Euler characteristic (disconnected inclusions dominate) to negative (tunnels and
percolation dominate). The volume fraction at this crossing, `phi_c+`, is the maximum
disconnected inclusion volume fraction consistent with the Adler expectation. For a scalar
Gaussian field: `phi_c+ ≈ 0.16` (threshold `kappa_c+ ≈ 0.60`).

## Using the library

All computations use the `src/` library. Run from the repository root or set `sys.path`:

```python
import numpy as np
import sys, pathlib
sys.path.insert(0, str(pathlib.Path('.').resolve()))  # from repo root

from src.fields import GridSpec, generate_gaussian_rf
from src.measures import measure_all
from src.adler import expected_lkcs_gaussian

# Define a grid: 100^3 voxels, physical side 100.0
grid = GridSpec(n=100, length=100.0)
rng = np.random.default_rng(42)

# Generate one Gaussian RF with correlation length lc=8.0
field = generate_gaussian_rf(grid, lc=8.0, rng=rng)

# Threshold at kappa=1.0 — gives VF ≈ 16%
phase = field >= 1.0
m = measure_all(phase, grid)
print(f"measured: phi={m['vf']:.3f}, euler={m['euler']:.0f}")

# Compare to Adler-Taylor prediction
theory = expected_lkcs_gaussian(grid, lc=8.0, kappas=np.array([1.0]))
print(f"theory:   phi={theory['phi'][0]:.3f}, euler={theory['L0'][0]:.0f}")
```

**`measure_all` returns:**

| Key | Description |
|-----|-------------|
| `vf` | Volume fraction |
| `volume` | Total inclusion volume |
| `euler` | Euler characteristic (cubical complex) |
| `half_surface` | Half the marching-cubes surface area |
| `percolates_x/y/z` | Whether any component spans the domain |

**Threshold sweep (experiment 01 pattern):**

```python
thresholds = np.linspace(-2.5, 2.5, 25)
theory = expected_lkcs_gaussian(grid, lc=8.0, kappas=thresholds)

measured_L0 = []
for kappa in thresholds:
    m = measure_all(field >= kappa, grid)
    measured_L0.append(m['euler'])
# theory['L0'] and measured_L0 match to within single-realisation noise
```

## What the Adler link gives you

1. **Closed-form calibration.** Given a target volume fraction and Euler density, the system
   `(phi_target, chi_target) = (M0(u), E[L0]/a^3)` can be inverted numerically to find `(lc, u)`
   without running any simulations. This is the primary use in materials modelling.

2. **No simulation budget for expectation.** The expected LKCs are available instantly for any
   parameter set. Simulation is only needed to quantify realisation-to-realisation variance or
   to explore models where the analytic formula does not yet hold.

3. **Limits: expectation only.** The formula gives `E[L_i]`. A single realisation will differ.
   Variance decreases with domain size (ergodicity), but for small domains or near topological
   transitions the fluctuations can be large. See [limits.md](limits.md) for a precise account
   of where the formula holds and where it does not.
