# Theory notes

The analytical objective is to compute expected LKCs from the Gaussian
Minkowski functionals (GMFs) of a hitting set. This controls expectations,
especially $L_0$ and $L_3$, while a scalar threshold selects a morphology
within a fixed model family. It does not determine every geometric feature
or guarantee nonspanning aggregates.

| Hitting set | Status |
|---|---|
| Gaussian half-line | Standard GMFs and GKF |
| Chi-squared exterior | Standard radial GMFs at positive thresholds |
| Student-t cone | Standard ratio-field coefficients; supported formula range $\nu\geq3$ |
| Margin chambers, $\delta>0$ | Detailed tube derivation draft and independent value-space checks; mathematical review remains |
| Disjoint ball union, $r<r_c$ | Additive ball tube coefficients, numerically evaluated |
| Overlapping ball union | Intersection GMFs not implemented |
| Independent multiscale union | Conditional-GKF recursion for intersections and unions, subject to the component formulas and regularity assumptions |

## Notes and evidence

- [Student-t formulas](student_t_derivation.md): corrected flag coefficients,
  cube expansion, common-amplitude invariance, cone coefficients, and EC-zero
  versus spanning. The note gives consistency checks, not a full cone-tube proof.
- [Margin and ball GMFs](gaussian_minkowski_functionals.md): the margin
  stratum formulas and ball-radius derivatives. Margin $\delta=0$ is a
  separate full-space case, not the positive-margin limit.
- [Margin tube derivation](margin_tube_derivation.md): chamber separation,
  normal-cone integration, all four GMFs, a stratum-derivative recurrence,
  and high-VF asymptotics for expected $L_0$ and $L_3$.
- [Independent margin checks](margin_tube_checks.md): exact low-dimensional
  controls, direct-distance tube quadrature, 12 tests, and saved results.
- [Physical film-width theory](margin_film_widths.md): interface-area sampling,
  the inverse-chi-four thin-film law, and consistency with the leading
  Adler matrix volume; [completed spatial pilot](../experiments/04_margin_tessellation/FILM_WIDTHS.md)
  with 16 controls. Gradient statistics agree well. A subsequent
  [bounded tail audit](../experiments/04_margin_tessellation/FILM_TAILS.md)
  with eight new controls identifies additional competitors as the dominant
  contribution to most large width errors in the sample. All full-field
  endpoints are found, but 31 of 1,536 profile/margin cases retain numerical
  stability uncertainty; the thin-film law is not a finite-VF guarantee.
- [Framework](../docs/framework.md): common normalization and assumptions.
- [Numerical limits](../docs/limits.md): evidence and measurement qualifications.
- [Multiscale experiment](../experiments/06_multiscale_margin/README.md):
  the conditional intersection operator for different spatial scales.

The committed margin sweep supports useful agreement in some regimes.
At $N=256,q=16$, target VF 60%, all four ensemble-mean errors are within 2%.
That is not a uniform bound: other configurations have larger errors.
Previous references to periodic Euler tests describe exploratory work;
a reproducible periodic validation suite is not included in the published
experiment. A proof cannot be replaced by either comparison.

Historical coalescing-lobe results have roughly 2–7% pre-contact errors
depending on descriptor and radius. Their post-contact per-ball and
partition measurements are not union LKCs. The corrected full sweep and
subsequent resolution and geometry checks are now complete; see
[experiment 05](../experiments/05_coalescing_lobes/README.md).

## Outstanding mathematical work

1. Independently review the new margin tube derivation and its use of the
   GKF hypotheses. Orientation, cone constants, positive reach, and the
   singular endpoint are now derived explicitly. The new area-weighted
   film-width law also needs review and finite-margin validation near
   junctions and weak gradients.
2. Derive Gaussian-weighted intersection GMFs for overlapping balls.
3. Validate the continuum/digital comparison with controlled refinement,
   repeated fields, and explicit connectivity conventions.
4. Establish physically resolved high-VF aggregates and their grading curve;
   an expected-Euler zero or positive Euler alone does not establish this.

## References

- Adler & Taylor (2007), *Random Fields and Geometry*.
- [Taylor (2006)](https://arxiv.org/abs/math/0602545), Gaussian kinematic formula.
- [Taylor & Adler](https://arxiv.org/abs/math/0612580), LKC GKF and flag coefficients.
- Worsley (1994), *Local maxima and the expected Euler characteristic of excursion sets of chi-squared, F and t fields*.
- [Roubin & Colliat (2016)](https://arxiv.org/abs/1510.08252), bounded-domain EC-zero proxy.
- [Bogdan (2015)](https://theses.hal.science/tel-01346484v1), multiscale cementitious-material modelling.
