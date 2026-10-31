# Multiscale margin pilot results

Grid: 128³; domain side: 64.0 mm.
Realizations: 5; elapsed: 216.9 s.

| Target VF | Mean VF | Mean L0 (6 / 26) | GKF L0 | Spanning (6 / 26) | Mean largest fraction (6) |
| --- | --- | --- | --- | --- | --- |
| 60% | 59.791% | 2084.6 / -18342.8 | 3598.1 | 5/5 / 5/5 | 96.545% |
| 70% | 69.842% | -4150.0 / -27780.4 | -558.2 | 5/5 / 5/5 | 98.783% |
| 75% | 74.862% | -9399.4 / -28414.2 | -3028.2 | 5/5 / 5/5 | 99.384% |
| 80% | 79.871% | -16330.8 / -22477.8 | -5754.6 | 5/5 / 5/5 | 99.725% |
| 85% | 84.912% | -23621.4 / -6165.0 | -8729.3 | 5/5 / 5/5 | 99.901% |

This is a finite-resolution experiment, not validation of a concrete mixture.
The sand class is a scale proxy. A grading curve requires calibration of measured grain sizes.
Coarse-first ownership recovers class-volume accounting but does not remove physical contacts.
GKF predictions inherit the margin GMFs' validation status; 6/26 disagreement diagnoses discretization sensitivity.
Spanning is measured on the actual union with bounded face-to-face connectivity.

See metadata.json for the effective configuration and source hashes.
