# Saved q-scene threshold sweep

Status: complete; 5 thresholds, one saved realization.
Thresholds and winner labels are correlated observations, not independent trials.
6 / 26 always denotes digital connectivity; no continuum claim follows from either.

| Target VF | Measured VF | Delta | Spanning winners, 6 / 26 | Union spans, 6 / 26 | Largest separate / phase, 6 / 26 | Largest union / phase, 6 / 26 |
|---|---|---|---|---|---|---|
| 60% | 59.890% | 0.236500 | 0 / 0 | True / True | 0.641% / 0.642% | 76.854% / 99.848% |
| 70% | 69.898% | 0.167189 | 0 / 0 | True / True | 0.618% / 0.630% | 98.576% / 99.992% |
| 75% | 74.921% | 0.135616 | 0 / 0 | True / True | 0.621% / 0.621% | 99.667% / 99.998% |
| 80% | 79.916% | 0.105761 | 0 / 0 | True / True | 0.612% / 0.612% | 99.899% / 99.999% |
| 85% | 84.935% | 0.077423 | 0 / 0 | True / True | 0.605% / 0.605% | 99.978% / 100.000% |

| Target VF | Expected L0 | Per-winner L0, 6 / 26 | Union L0, 6 / 26 | Separate components, 6 / 26 | Cross-label pairs, 6 / 26 |
|---|---|---|---|---|---|
| 60% | 3788.59 | 4406 / 3575 | 1292 / -37635 | 4795 / 3631 | 4833 / 208119 |
| 70% | 3980.00 | 4640 / 3720 | -13925 / -39705 | 5102 / 3819 | 33848 / 633231 |
| 75% | 4059.36 | 4726 / 3827 | -29057 / -22731 | 5241 / 3928 | 74719 / 1018592 |
| 80% | 4129.19 | 4796 / 3870 | -44462 / 5975 | 5347 / 3991 | 146790 / 1555070 |
| 85% | 4190.34 | 4919 / 3946 | -48525 / 41391 | 5496 / 4063 | 262866 / 2281483 |

## Interpretation

The theory columns are ensemble GKF expectations under the declared q and lc;
one saved scene cannot validate these expectations or supply uncertainty intervals.
Winner masks are measured separately to expose cross-label voxel mergers.
Their Euler sums may still have within-label discretization errors.
Cross-label contacts at positive margin indicate unresolved continuum separation.
This rethresholding does not reconstruct the underlying Gaussian components or their films.
Point data use the existing N*dx voxel-volume convention, not the (N-1)*dx VTI point extent.
The size statistics include boundary-cut and spanning objects. components.csv records
their flags; equivalent diameter is volume-based and is not a calibrated sieve diameter.
No spanning is a statement about these masks, not a probability or continuum guarantee.
