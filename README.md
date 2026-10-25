# high-vf-excursions

Excursion sets of Gaussian-related random fields as microstructure models, with analytical
morphological descriptors through the Gaussian Kinematic Formula (GKF).

A scalar Gaussian random field in 3D tops out at roughly 16% disconnected inclusion volume
fraction — the threshold at which the Euler characteristic of the excursion set crosses zero.
This repository investigates a progression of model extensions (Student-t fields, margin
tessellations, coalescing Gaussian-space lobes) that push the disconnected volume fraction
ceiling to 90%+, and tracks what happens to the analytical Adler-Taylor link at each step.
Each model is parameterised by (correlation length, threshold, model-specific parameters) and
produces a binary microstructure whose expected morphological descriptors — Euler characteristic,
surface area, volume fraction — are predicted by closed or semi-closed formulae derived from the
GKF. The central question at each step is: how much of that analytical control survives?

## Experiment progression

| # | Model | Max disconnected VF | Adler link |
|---|-------|---------------------|------------|
| 01 | Gaussian RF (scalar) | ~16% | Full — exact closed-form M_j, Monte Carlo validated |
| 02 | Student-t_v fields (v >= 3) | ~35–45% (v=3) | Full — Worsley M_j formulae, Monte Carlo validated |
| 03 | Margin tessellation (q competitors) | 85–90%+ | Partial — M_j derived and numerically validated; formal stratified-tube proof pending |
| 04 | Coalescing Gaussian-space lobes | 93%+ | Partial — M_j before ball contact; inclusion-exclusion terms after contact unresolved |

For background on the GKF and library usage see [docs/framework.md](docs/framework.md).  
For a precise account of what is and is not analytically known see [docs/limits.md](docs/limits.md).  
For a concise research narrative see [docs/explorations.md](docs/explorations.md).

## Quick start

```bash
mamba env create -f environment.yml
mamba activate rf_explo

cd experiments/01_gaussian_baseline
mamba run -n rf_explo python run.py
# outputs written to experiments/01_gaussian_baseline/results/
```

Each experiment writes figures and CSV files to its own `results/` subfolder. See the
experiment's `README.md` for parameters, expected runtime, and key findings.

## Repository layout

```
src/
  fields.py          RF generation: Gaussian (FFT), Student-t (Worsley), margin stack, vector field
  excursions.py      Excursion set construction: threshold, margin phase, lobe distances
  measures.py        LKC measurement: Euler characteristic, surface area (marching cubes), percolation
  adler.py           Adler-Taylor formulae: Gaussian, Student-t, margin tessellation M_j
  io.py              ParaView export: VTI (scalar/binary fields), VTP (surface mesh)

experiments/
  01_gaussian_baseline/   Scalar Gaussian RF — full Adler validation baseline
  02_student_t/           Student-t fields — heavier tails, higher VF ceiling
  03_margin_tessellation/ Winner-competition model — disconnected VF to 90%+
  04_coalescing_lobes/    Gaussian-space ball union — controlled lobe merging and contact

explorations/
  oscillatory_bochner/    Oscillatory covariance and Bochner boundary (ongoing, separate direction)

theory/
  student_t_derivation.md            Derivation of t_v Minkowski functionals (Worsley 1994 route)
  gaussian_minkowski_functionals.md  GMF formulae for the margin tessellation hitting set
  paper_directions.md                Research roadmap and open problems

references/            Reference PDFs
adler_high_vf_exploration/   Archived research scripts and logs from the exploratory phase
```

## References

- Adler, R.J. & Taylor, J.E. (2007). *Random Fields and Geometry*. Springer.
- Worsley, K.J. (1994). Local maxima and the expected Euler characteristic of excursion sets
  of chi-squared, F and t fields. *Advances in Applied Probability* 26(1), 13–42.
- Roubin, E. & Colliat, J.-B. (2016). Critical probability of percolation over bounded region
  in N-dimensional Euclidean space. *Journal of Statistical Mechanics* P033306.
- Roubin, E. (2013). *Meso-scale FE and morphological modelling of heterogeneous media*. PhD
  thesis, ENS Cachan.
