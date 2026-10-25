# Research directions and paper planning

## Core scientific question

Can we construct Gaussian-related random morphologies with:
1. Disconnected inclusion volume fractions significantly above the Gaussian ~16% ceiling,
2. While preserving the Adler-Taylor Gaussian kinematic formula (GKF)?

The GKF gives closed-form expectations for all Lipschitz-Killing curvatures (LKCs) of
excursion sets. Preserving it means the morphology parameters (volume fraction, surface
area, Euler characteristic) remain analytically tractable — a critical property for
inverse calibration to experimental materials data.

## The progression

### Track 1: Scalar non-Gaussian fields (fully analytic)

Student-t and chi-squared fields achieve higher VF ceilings (up to ~45% for t_3)
by using heavier-tailed marginals. The Adler formula remains exact; only M_j changes.
See theory/student_t_derivation.md for the derivation.

This is the cleanest extension: fully within Worsley/Taylor theory, no new mathematics
needed beyond computing M_j for known distributions.

### Track 2: Margin tessellation (GKF-compatible, proof outstanding)

The hitting set D_delta = {z : z_(q) - z_(q-1) >= delta} is a disjoint union of q
polyhedral winner chambers in Gaussian value space. Since this is a stratified algebraic
set, the GKF applies via the stratified Gaussian tube formula.

Closed M_j formulas (M_0 through M_3) have been derived from the polyhedral normal-cone
geometry and numerically validated at several volume fractions. A formal proof of the
stratified-tube decomposition needs to be written up for publication.

Volume fraction up to 90%+ is achievable for q >= 16, with entirely disconnected
winner-cell cores separated by a matrix film of controlled thickness.

**For a first paper:** the margin tessellation, its closed M_j formulas, and a
systematic finite-size percolation/RVE study could form a self-contained contribution
to a morphology or probability journal.

### Track 3: Coalescing Gaussian-space lobes (partially open)

The ball-union hitting set D_r = union_a B(c_a, r) is stratified and GKF-admissible.
Before lobe contact (r < r_c), M_j(D_r) is additive over balls. Each ball's M_j
involves a non-central chi-squared CDF — derivable but not yet done.

After contact (r > r_c), Gaussian-weighted ball-intersection volumes enter via
inclusion-exclusion. These are not yet analytically derived. This is the main
open mathematical problem for Track 3.

The morphological interest is the separation between three transitions:
1. Hitting-set contact in Gaussian value space (prescribed at r = r_c).
2. Zero-crossing of E[Euler characteristic] — topological transition.
3. Physical percolation (P_span = 1/2 with periodic wrapping).

The pilot study showed physical spanning before value-space contact — an unexpected
result that deserves a systematic multi-resolution investigation.

### Track 4 (future): Oscillatory covariance and Bochner boundary

Extending the correlation length to the complex plane generates a richer morphological
landscape. The boundary of Bochner validity in complex parameter space is the
'Mandelbrot set' analogue for Gaussian RFs. Morphologies near this boundary may exhibit
near-fractal topology.

See explorations/oscillatory_bochner/ for the current pilot.

---

## Minimum credible numerical campaign (for publication)

1. **Margin tessellation:**
   - Periodic grids with union-find wrapping for percolation detection.
   - At least three voxel resolutions at fixed physical T/Lc.
   - T/Lc sweep large enough to extrapolate finite-size crossings.
   - Hundreds of realisations near each transition, fewer away from it.
   - Confidence intervals for wrapping probability and largest-cluster fraction.
   - Direct validation of all predicted LKCs (L0, L2) against cubical/mesh estimators.
   - Report matrix film thickness in voxel units; establish resolved-film criterion.

2. **Coalescing lobes:**
   - Derive M_j for a single non-central ball (pre-contact).
   - Derive ball-intersection M_j via inclusion-exclusion (post-contact).
   - Multi-resolution check near the physical percolation transition.
   - Separate measurement of the three transition radii (contact, EC=0, P_span=1/2).

---

## Claims that are not yet justified

- "90% VF with disconnected morphology" is demonstrated numerically for the margin
  tessellation but requires a formal proof of the GMF formulas.
- The pilot Euler dip near r/r_c = 1 in the coalescing-lobe experiment is resolution-
  sensitive and has not been systematically checked.
- Novelty has not been established by a comprehensive literature review. The specific
  winner-runner-up margin construction has not been found in the literature, but this
  is only a preliminary check.
- phi_c+ from E[L0] = 0 is not the same as the spatial percolation threshold.
  Their separation has not been measured systematically.

---

## References

- Adler, R.J. & Taylor, J.E. (2007). *Random Fields and Geometry*. Springer.
- Worsley, K.J. (1994). Local maxima and the expected Euler characteristic of excursion
  sets of chi^2, F and t fields. *Advances in Applied Probability* 26(1), 13–42.
- Taylor, J.E. (2006). A Gaussian kinematic formula. *Annals of Probability* 34(1), 122–158.
- Roubin, E. & Colliat, J.-B. (2016). Critical probability of percolation over bounded
  region in N-dimensional Euclidean space. *J. Stat. Mech.* P033306.
- Roubin, E. (2013). *Meso-scale FE and morphological modelling of heterogeneous media*
  (PhD thesis). Chapter 2.
