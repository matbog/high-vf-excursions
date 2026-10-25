# Exploration — Oscillatory Covariance RFs and the Bochner Boundary

## Status

Exploratory direction. Not yet connected to the main high-VF Adler story.
Planned for future investigation.

## Scientific question

The Feigenbaum bifurcation diagram lives on the real axis of the Mandelbrot set.
Analogously, standard Gaussian RF morphology lives in a real-parameter slice of a
larger covariance parameter space. What happens when the covariance parameter is
extended to the complex plane?

For C(r) = exp(-r²/2Lc²) with complex Lc = R exp(i theta):

    Re[S(k)] = (2pi)^(3/2) R^3 exp(-k²R²cos2theta/2) cos(3theta - k²R²sin2theta/2)

Bochner's theorem: C is a valid covariance iff Re[S(k)] >= 0 for all k.
The boundary of validity in (theta, R) space is the 'Mandelbrot set' of Gaussian RFs.

This script also generates oscillatory covariance fields via a sum-of-Gaussians
spectral density (always valid by construction), which interpolate between standard
Gaussian morphology and quasi-periodic structures.

## Key parameters

| Parameter | Value |
|---|---|
| `SEED` | 7 |
| `N` | 100 (voxels per side) |
| Bochner map: theta range | 0 to pi/4 (120 points) |
| Bochner map: R range | 2 to 20 voxels (100 points) |
| Oscillatory sweep: L_d | 8.0 |
| Oscillatory sweep: lambda | [500, 40, 16, 8, 4] voxels |

## How to run

```bash
mamba run -n rf_explo python run.py
```

Expected runtime: ~2–3 minutes.

## Outputs

| File | Contents |
|---|---|
| `results/bochner_map.png` | Bochner validity in (theta, R) space + example field slices |
| `results/oscillatory_sweep.png` | Morphology from smooth blobs to quasi-periodic structures |

## Potential future directions

- Map the full Bochner boundary for the squared-exponential and Matern covariance families.
- Derive Adler LKC formulas for oscillatory covariance RFs (lambda_2 changes; are higher
  spectral moments needed?).
- Explore whether morphologies near the Bochner boundary exhibit near-fractal topology.
- Connection to Lee-Yang zeros: analytic continuation of the EC formula E[L0](u + iv) in
  the complex threshold plane.
