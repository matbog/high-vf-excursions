# Experiment 03 — Gaussian Margin Tessellation

## Scientific question

Can competing Gaussian random fields produce disconnected morphologies with volume
fractions above 80–90%, while remaining within the Gaussian Kinematic Formula?

## Model

Let G = (G_1, ..., G_q) be q independent smooth Gaussian RFs. Define:

    A_delta = {x : max_i G_i(x) - second_max_i G_i(x) >= delta}

The hitting set D_delta in R^q is a disjoint union of q polyhedral winner chambers.
Since D_delta is a stratified set in Gaussian value space, the morphology is a
Gaussian-related excursion and the GKF applies.

The matrix film (excluded phase) consists of voxels where two or more fields are
nearly tied. Increasing q creates more, smaller winner cells. Increasing delta
thickens the matrix film.

## Adler link status: **Partial**

Closed one-dimensional integral formulas for M_0 through M_3 of D_delta have been
derived from the polyhedral normal-cone structure (see theory/gaussian_minkowski_functionals.md).
They have been numerically validated via periodic Euler density experiments in
dimensions 1, 2, 3 at phi = 0.50, 0.70, 0.80.

**What is missing:** a formal published proof of the stratified-tube decomposition
for this polyhedral hitting set. The formulas are used here with this caveat explicitly
noted in plots and outputs.

## Key parameters

| Parameter | Value | Effect |
|---|---|---|
| `N` | 64 | Grid resolution |
| `LC` | 5.0 | Correlation length (in voxels: LC/dx = 5) |
| `Q_VALUES` | [2, 4, 8, 16] | Number of competing fields |
| `TARGET_PHIS` | [0.30 .. 0.90] | Target volume fractions |
| `N_REALIZATIONS` | 5 | Realisations per (q, phi) |
| `SEED` | 1729 | Base seed |

## How to run

```bash
mamba run -n rf_explo python run.py
```

Expected runtime: ~20–30 minutes.

## Outputs

| File | Contents |
|---|---|
| `results/margin_sweep.csv` | q, phi, delta, measured phi, Euler, spanning per realisation |
| `results/margin_sweep.png` | Spanning frequency and Euler vs phi per q |

## Key findings (from pilot study, 2026-08-17)

- q=4: spanning already at phi=50% (4 winner cells percolate easily).
- q=8: no spanning through phi=60%, mixed at 70%, spanning at 80%.
- q=16: no spanning through phi=90%. At phi=90%, the largest same-winner component
  averages ~1.4% of retained material — morphology is highly disconnected.
- q=32: no spanning through phi=90%, largest component ~0.6%.

**Critical caveat:** as phi → 1, delta → 0 and the matrix film becomes arbitrarily
thin. A publishable study must report film thickness in voxel units and perform a
resolution convergence study. Label-aware connectivity is mathematically correct for
delta > 0, but a physically manufacturable morphology also requires the separating
film to be resolved on the grid.
