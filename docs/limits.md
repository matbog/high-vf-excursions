# Scope and limits of the Adler link

The GKF predicts continuum **expected LKCs** under its Gaussian-component
and geometric regularity assumptions. It does not predict an individual
sample's topology, particle-size distribution, or spanning probability.
The source of spatial smoothness is the underlying Gaussian vector;
nonlinear derived fields need not share its covariance or gradient moments.

## Analytical status

| Family | Available prediction | Main qualification |
|---|---|---|
| Gaussian | All four expected LKCs | Standard smooth-field GKF; numerical estimators remain approximate |
| Chi-squared | All four at positive thresholds | Squared independent components; low-threshold separating films can be unresolved |
| Student-t | All four for the supported $\nu\geq3$ formula API | Ratio of independent components; generator requires integer $\nu$; common amplitude cancels |
| Margin | Implemented GMFs through order 3 for $\delta>0$ | Stratified-tube derivation remains an outline; exact $\delta=0$ is a separate full-space case |
| Disjoint balls | All four for strictly disjoint balls | Additivity stops at contact; intersection terms have not been implemented |
| Independent multiscale union | Conditional intersection-GKF recursion | Requires independent stages and regularity; derived from component GMFs, with their limitations |

The margin higher-order formulas have numerical support, not a completed
publication-grade proof. The committed sweeps do not establish a uniform
“under 2%” accuracy bound. At $N=256,q=16$, target VF 60%, the five-realisation
mean errors are approximately $-1.53\%$, $-0.77\%$, $-1.89\%$, and $-0.09\%$
for $L_0,L_1,L_2,L_3$, respectively. Other settings differ more: at
$q=2$, VF 30%, the mean Euler error is about $-16.4\%$.

## Topology claims

Keep three quantities separate:

1. Expected-Euler zero: a topology proxy derived from an LKC expectation.
2. Contact of hitting-set pieces in Gaussian value space.
3. Spatial spanning, measured on a particular sample and geometry.

Positive Euler characteristic can coexist with spanning. Distinct margin
winner labels are separated in the continuous model for $\delta>0$,
but a single winner label can span. Digital adjacency can bridge an
unresolved film. Neither 6- nor 26-neighbour connectivity is automatically
the correct continuum answer.

The saved $q=16$, target VF 60%, $N=256$ showcase is nonspanning under
6-neighbour connectivity and spans all axes under 26-neighbour connectivity.
It is evidence of a resolution/geometry ambiguity, not proof of disconnected
physical aggregates. Zero observed spans in ten independent fields also
does not imply zero probability: the one-sided 95% binomial upper bound is
about 25.9%. Reusing fields across thresholds does not create independent
realisations.

## Measurement and generation

Voxel volume, digital Euler, Crofton area, and interpolated mesh curvature
are different approximations to continuum descriptors. Thin features can
be missed by all of them. Meshing a sampled minimum of several scalar
fields is not exact clipping of the original component fields.
No universal convergence rate is established by the committed data.

GSTools supplies Gaussian-field realisations with a numerical generator.
Its finite spectral representation and the sampling grid introduce
approximation beyond the ideal covariance model. Do not estimate an
arbitrary nonlinear field's gradient variance and substitute it for the
Gaussian-component metric.

Historical outputs retain their original parameters and estimator schemas.
Correcting code or documentation does not retroactively validate those
outputs. In particular, experiment 05's historical post-contact per-ball
and partition columns are not union LKCs.

## Remaining work

A new resolution/realisation campaign, consistent interpolated geometry,
a full margin tube proof, and intersecting-ball GMFs remain separate tasks.
The present corrections do not claim these are finished. High VF is easy
to prescribe in expectation; high VF with resolved, nonspanning aggregates
and a measured grading curve still needs demonstration.
