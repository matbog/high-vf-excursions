# Paper direction: Adler-linked Gaussian margin tessellations

## Proposed central contribution

Introduce a vector-Gaussian morphology model in which (q) correlated random
fields compete locally and a material point is retained only when the largest
field value exceeds the second largest by a margin. The margin creates a matrix
film between random winner cells. Unlike ad hoc watershed post-processing, the
result is exactly the preimage of an explicit stratified hitting set, so its
expected Lipschitz–Killing curvatures remain governed by the Gaussian kinematic
formula.

The preliminary result worth pursuing is not merely “90% disconnected.” It is:

> Increasing the dimension of the latent Gaussian field changes the
> percolation class of an Adler-linked excursion, allowing high retained volume
> fractions while maintaining separated winner-cell cores.

## Potential paper contributions

1. **New morphology family.** Define the margin hitting set and establish its
   stratification and admissibility under the GKF.
2. **Analytic descriptors.** Prove the derived Gaussian Minkowski functionals
   (\mathcal M_0^\gamma,\ldots,\mathcal M_3^\gamma), all expressed through
   one-dimensional integrals and universal polyhedral normal-cone constants.
3. **Percolation phase diagram.** Estimate wrapping probability in
   ((q,\delta,T/L_c)), preferably with periodic topology.
4. **Two RVE notions.** Compare the domain size needed for stable expected LKCs
   with the size needed for stable wrapping probability.
5. **Resolution/physical-film criterion.** Relate the value-space margin to the
   distribution of spatial matrix-film thickness and identify a resolved,
   manufacturable high-volume region.
6. **Benchmarks.** Compare one-sided Gaussian, (\chi_1^2), independent unions,
   and margin tessellations at equal volume fraction and spectral moment.

## Claims that are not yet justified

- The 90% pilot is not a continuum percolation threshold.
- Face-to-face spanning on a finite nonperiodic voxel cube is not wrapping
  probability.
- Positive Euler characteristic is not equivalent to absence of percolation.
- Novelty has not yet been established by a systematic literature review.
- The proposed closed (\mathcal M_2^\gamma) and (\mathcal M_3^\gamma)
  formulas are numerically supported but still require a formal checked proof.

## Minimum credible numerical campaign

- Periodic grids with union-find wrapping detection.
- At least three resolutions per fixed physical parameter set.
- (T/L_c\) sweep large enough to extrapolate finite-size crossings.
- Hundreds of realizations near each transition, fewer away from it.
- Confidence intervals for wrapping probability and largest-cluster fraction.
- Direct validation of all predicted LKCs against cubical-complex estimators.

## Secondary paper option

If the full higher-order GMF derivation becomes too large, a first paper could
focus on the new random tessellation, exact volume/surface expectations, and
finite-size percolation/RVE. A follow-up could develop the complete stratified
Gaussian tube coefficients and inverse morphological calibration.

## Complementary direction: controlled coalescence

The margin tessellation keeps winner chambers separated for every positive
margin. A complementary Adler-compatible family uses expanding balls in
Gaussian value space: `A_r = G^{-1}(D_r)`, with
`D_r = union_a B(c_a,r)`. Symmetrically placed centres give exact hitting-set
contact levels controlled by the single scalar `r`, while retaining the GKF.

The initial pilot uses four tetrahedral centres in `R^3`. Its paper question is:

> How do prescribed mergers of a stratified Gaussian hitting set control the
> Euler-characteristic and percolation transitions of its spatial preimage?

The exact value-space contact radius, expected-EC zero, and finite-domain
`P_span=1/2` crossing are to be compared, not assumed equal. Their offsets
quantify how value-space topology transfers to physical-space connectivity.
