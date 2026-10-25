# Experiment 04 — Coalescing Gaussian-Space Lobes

## Scientific question

Can a vector Gaussian field with ball-shaped hitting sets in value space produce
high volume fractions with controlled lobe merging, and what are the morphological
transitions relative to the value-space contact event?

## Model

Let G(x) = (G_1(x), G_2(x), G_3(x)) map each voxel to a point in R^3. Place K=4
fixed centres c_a at the vertices of a regular tetrahedron of radius rho. Include
a voxel when its Gaussian-space point is within radius r of any centre:

    A_r = {x : min_a ||G(x) - c_a|| <= r}

The single parameter r/r_c controls morphology continuously, where the contact radius:

    r_c = rho * sqrt(2/3)

is the radius at which adjacent Gaussian-space balls first touch.

## Adler link status: **Partial (open problem after contact)**

Before r/r_c < 1: the K balls are disjoint. M_j(D_r) = sum_a M_j(B(c_a, r)).
Each ball's contribution involves a non-central chi-squared CDF; derivation of
all four M_j is tractable but not yet written up.

After r/r_c > 1: balls overlap. Inclusion-exclusion requires Gaussian-weighted
ball-intersection volumes in R^3. These have not yet been analytically derived.

**Volume fraction** is calibrated by Gaussian-space Monte Carlo throughout.
No analytic Euler prediction is available after contact.

## Key parameters

| Parameter | Value | Effect |
|---|---|---|
| `N` | 100 | Grid resolution |
| `LC` | 10.0 | Correlation length |
| `PADDING` | 20 | Crop padding to reduce periodic BC artefacts |
| `RHO` | 1.5 | Tetrahedron radius in Gaussian-space units |
| `N_REALIZATIONS` | 4 | Number of realisations |
| `RADIUS_OVER_RC` | [0.45 .. 1.60] | r/r_c sweep |
| `SEED` | 9417 | Base seed |
| `MC_SAMPLES` | 2,000,000 | MC samples for volume calibration |

## How to run

```bash
mamba run -n rf_explo python run.py
```

Expected runtime: ~15–25 minutes (MC volume calibration + 4 realisations × 8 radii).

## Outputs

| File | Contents |
|---|---|
| `results/coalescing_lobes_sweep.csv` | r/r_c, measured phi, Euler, largest component, spanning |
| `results/coalescing_lobes_slices.png` | Winner-coloured central slices per realization × radius |
| `results/coalescing_lobes_diagnostics.png` | 4-panel: volume, Euler, coalescence, percolation |

## Key findings (from pilot study, 2026-08-17)

| r/r_c | phi | Euler | Largest/phase | Spanning |
|---:|---:|---:|---:|---:|
| 0.45 | 0.056 | 340 | 0.051 | 0 |
| 0.65 | 0.165 | 288 | 0.055 | 0 |
| 0.85 | 0.355 | 144 | 0.311 | 1 |
| 0.98 | 0.526 | -2870 | 0.999 | 1 |
| 1.02 | 0.584 | -2048 | 1.000 | 1 |
| 1.15 | 0.732 | 220 | 1.000 | 1 |
| 1.35 | 0.852 | 14 | 1.000 | 1 |
| 1.60 | 0.933 | 62 | 1.000 | 1 |

**Main observation:** physical spanning occurred *before* value-space contact (r/r_c < 1).
Near r/r_c = 1, the Euler characteristic showed an extreme negative dip, suggesting a
tunnel-rich transition. The Euler magnitude near contact must be checked for resolution
sensitivity.

## Open questions

1. Derive M_j for a single off-centre ball (non-central chi-squared preimage).
2. Derive M_j for the ball union via Gaussian-weighted inclusion-exclusion after contact.
3. Characterise why physical percolation precedes value-space contact.
4. Perform multi-resolution check near the transition.
