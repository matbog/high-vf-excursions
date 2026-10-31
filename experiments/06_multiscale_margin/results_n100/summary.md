# Multiscale margin pilot results

Grid: 100³; domain side: 64.0 mm.
Realizations: 5; elapsed: 113.8 s.

| Target VF | Mean VF | Mean L0 (6 / 26) | GKF L0 | Spanning (6 / 26) | Mean largest fraction (6) |
| --- | --- | --- | --- | --- | --- |
| 60% | 59.396% | 426.0 / -16222.2 | 3598.1 | 5/5 / 5/5 | 97.129% |
| 70% | 69.558% | -6098.2 / -18946.8 | -558.2 | 5/5 / 5/5 | 99.055% |
| 75% | 74.627% | -10310.6 / -16345.8 | -3028.2 | 5/5 / 5/5 | 99.516% |
| 80% | 79.717% | -14473.6 / -9451.8 | -5754.6 | 5/5 / 5/5 | 99.785% |
| 85% | 84.793% | -16776.0 / 2557.6 | -8729.3 | 5/5 / 5/5 | 99.926% |

This is a finite-resolution experiment, not validation of a concrete mixture.
The sand class is a scale proxy. A grading curve requires calibration of measured grain sizes.
Coarse-first ownership recovers class-volume accounting but does not remove physical contacts.
GKF predictions inherit the margin GMFs' validation status; 6/26 disagreement diagnoses discretization sensitivity.
Spanning is measured on the actual union with bounded face-to-face connectivity.

See metadata.json for the effective configuration and source hashes.
