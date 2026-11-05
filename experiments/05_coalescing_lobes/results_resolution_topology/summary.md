# Paired pre-contact resolution study

Status: complete. 6 completed cases.
Curves join resolutions of the same realization. GKF curves are ensemble expectations.
One realization cannot validate an expectation; refinement does not remove sampling uncertainty.

| r/rc | N | completed realisations | VF | Euler 6 | Euler 26 | ball sum 26 | cross-label pairs 6 / 26 |
|---|---|---|---|---|---|---|---|
| 0.85 | 50 | 1 | 0.35242 | -213.00 | -1487.00 | 165.00 | 577.0 / 11251.0 |
| 0.85 | 100 | 1 | 0.35219 | 165.00 | -1449.00 | 178.00 | 17.0 / 4329.0 |
| 0.85 | 200 | 1 | 0.35221 | 185.00 | 85.00 | 184.00 | 0.0 / 107.0 |
| 0.98 | 50 | 1 | 0.52212 | -1724.00 | -702.00 | 66.00 | 7806.0 / 63203.0 |
| 0.98 | 100 | 1 | 0.52208 | -2942.00 | -2040.00 | 71.00 | 13606.0 / 135367.0 |
| 0.98 | 200 | 1 | 0.52189 | -4577.00 | -5344.00 | 71.00 | 14786.0 / 200466.0 |

Cross-label contacts are digital contacts forbidden between distinct continuum lobes for r/rc<1.
Raw pair counts need not decrease monotonically: the number of neighbour pairs grows with N.
The CSV also records contact fractions, same-convention additivity gaps, per-ball spanning,
and common-node sampling errors. Neither digital convention certifies continuum connectivity.
Scalar-mesh interpolation and domain padding retain their usual limitations.
