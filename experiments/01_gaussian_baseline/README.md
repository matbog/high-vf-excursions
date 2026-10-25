# Experiment 01 — Gaussian RF Baseline

## Scientific question

For a scalar Gaussian random field with Gaussian covariance C(r) ~ exp(-r²/lc²),
do the Adler-Taylor expected LKCs match numerical measurements across the full
threshold range?

This is the ground truth. Every subsequent experiment is compared against this baseline.

## Adler link status: **Full**

The hitting set D_u = [u, +∞) under the Gaussian marginal has exact closed-form
Minkowski functionals M_j (see theory/student_t_derivation.md §3). The Adler-Taylor
formula gives exact predictions for E[L0], E[L2], E[L3] as functions of threshold κ.
Monte Carlo mean ± std validates the single-realisation curves.

The maximum disconnected (meatball) volume fraction is ~16%, reached at the upper
zero-crossing of E[L0].

## Key parameters

| Parameter | Value | Effect |
|---|---|---|
| `N` | 100 | Grid resolution (100³ voxels) |
| `LENGTH` | 100.0 | Domain side length |
| `LC` | 8.0 | Correlation length (in same units as LENGTH) |
| `SEED` | 42 | Single-realisation seed |
| `SEED_MC` | 10 | Monte Carlo seed (offset per run) |
| `N_RUNS_MC` | 50 | Number of Monte Carlo realisations |
| `KAPPA_VTI` | 1.0 | Threshold for exported ParaView files |

## How to run

```bash
mamba run -n rf_explo python run.py
```

Expected runtime: ~5–10 minutes (dominated by Monte Carlo surface area computation).

## Outputs

| File | Contents |
|---|---|
| `results/threshold_sweep.csv` | kappa, measured L0/L2/L3/phi, theory L0/L2/L3/phi |
| `results/threshold_sweep.png` | 4-panel single realisation vs Adler theory |
| `results/mc_gaussian.csv` | mean ± std over 50 runs vs theory |
| `results/mc_gaussian.png` | 3-panel MC figure |
| `results/random_field.vti` | Scalar Gaussian field for ParaView |
| `results/excursion.vti` | Binary excursion at κ=1.0 for ParaView |
| `results/excursion_surface.vtp` | Surface mesh for ParaView |

## Key findings

- Adler theory matches measurements to within MC sampling noise across the full
  threshold range.
- Volume fraction saturates at ~50% (lower threshold, sponge regime) and ~0% (high
  threshold). The disconnected meatball regime ends at φ ≈ 16%.
- L2 (half surface area) carries a small O(dx) marching-cubes resolution bias;
  this is expected and does not affect the Adler comparison for the purposes of
  this project.
