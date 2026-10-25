# Random Field Morphology: Exploring the Gaussian Kinematic Formula at High Volume Fractions

The Adler-Taylor Gaussian kinematic formula (GKF) gives closed-form expectations for
the Lipschitz-Killing curvatures (LKCs) of excursion sets of smooth Gaussian-related
random fields. For materials modelling, this means a direct analytical link between
field parameters (correlation length, threshold) and morphological descriptors (volume
fraction, surface area, Euler characteristic). The Roubin-Colliat framework exploits
this link to calibrate microstructure models to experimental data.

This project investigates how far the disconnected inclusion volume fraction can be
pushed while keeping the GKF tractable. A scalar Gaussian field tops out at ~16%
disconnected volume fraction in 3D. We explore a progression of models — Student-t
fields, winner-competition margin tessellations, and Gaussian-space ball unions —
each achieving higher volume fractions at the cost of increasing analytic complexity.

## Experiment progression

| # | Model | Max disconnected VF | Adler link |
|---|---|---|---|
| 01 | Gaussian RF (scalar) | ~16% | **Full** — exact closed-form M_j, MC validated |
| 02 | Student-t_ν fields (ν ≥ 3) | ~35–45% (ν=3) | **Full** — Worsley M_j formulas, MC validated |
| 03 | Margin tessellation (q competitors) | 85–90%+ | **Partial** — M_j formulas derived & numerically validated, formal stratified-tube proof pending |
| 04 | Coalescing Gaussian-space lobes | 93%+ | **Partial** — M_j before ball contact; inclusion-exclusion terms after contact unresolved |

An additional exploration (oscillatory covariance / Bochner boundary) is in
`explorations/oscillatory_bochner/` — a separate direction, not yet connected to the
high-VF story.

## Setup

```bash
mamba env create -f environment.yml
mamba activate rf_explo
```

## Running experiments

Each experiment is self-contained:

```bash
cd experiments/01_gaussian_baseline
mamba run -n rf_explo python run.py

cd experiments/02_student_t
mamba run -n rf_explo python run.py

cd experiments/03_margin_tessellation
mamba run -n rf_explo python run.py

cd experiments/04_coalescing_lobes
mamba run -n rf_explo python run.py
```

All scripts use fixed seeds and write output to their own `results/` subfolder.
See each experiment's README.md for details on parameters, expected runtime, and outputs.

## Repository layout

```
src/                        Shared Python library
  fields.py                   RF generation (Gaussian, Student-t, vector fields)
  excursions.py               Excursion set construction (threshold, margin, lobes)
  measures.py                 LKC/Minkowski measurement (Euler, surface area, percolation)
  adler.py                    Adler formula implementations for all field families
  io.py                       ParaView VTI/VTP export

experiments/
  01_gaussian_baseline/       Gaussian RF — full Adler validation baseline
  02_student_t/               Student-t fields — heavier tails, higher VF ceiling
  03_margin_tessellation/     Winner-competition model — 90%+ disconnected VF
  04_coalescing_lobes/        Gaussian-space ball union — controlled lobe merging

explorations/
  oscillatory_bochner/        Complex covariance / Bochner boundary (future direction)

theory/
  student_t_derivation.md     Derivation of t_nu Minkowski functionals
  gaussian_minkowski_functionals.md   GMF formulas for the margin tessellation
  paper_directions.md         Research roadmap and open problems

references/                   Reference PDFs
adler_high_vf_exploration/    Historical research scripts and logs (archived)
```

## References

- Adler, R.J. & Taylor, J.E. (2007). *Random Fields and Geometry*. Springer.
- Worsley, K.J. (1994). Local maxima and the expected Euler characteristic of excursion
  sets of χ², F and t fields. *Advances in Applied Probability* 26(1), 13–42.
- Roubin, E. & Colliat, J.-B. (2016). Critical probability of percolation over bounded
  region in N-dimensional Euclidean space. *J. Stat. Mech.* P033306.
- Roubin, E. (2013). *Meso-scale FE and morphological modelling of heterogeneous media* (PhD).
