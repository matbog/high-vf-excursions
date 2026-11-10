# Refining the component input grids

The direct 33^3 reference grid is unchanged in the default study.
Numbers below are phase-disagreement percentage points of local patch volume.
Remaining reference sampling and topology uncertainties are unchanged.

| Patch | Target VF | Input refinement | Trilinear error | Tetrahedral error |
|---|---|---|---|---|
| 0 | 75% | 1 | 5.58243 | 6.98471 |
| 0 | 75% | 2 | 1.59683 | 1.97830 |
| 0 | 75% | 4 | 0.35553 | 0.45166 |
| 0 | 85% | 1 | 7.73926 | 9.44672 |
| 0 | 85% | 2 | 1.87607 | 2.24228 |
| 0 | 85% | 4 | 0.44708 | 0.53406 |
| 1 | 75% | 1 | 0.19836 | 1.17035 |
| 1 | 75% | 2 | 0.05493 | 0.28687 |
| 1 | 75% | 4 | 0.01831 | 0.07782 |
| 1 | 85% | 1 | 0.16022 | 1.21460 |
| 1 | 85% | 2 | 0.04120 | 0.28534 |
| 1 | 85% | 4 | 0.01984 | 0.07782 |
| 2 | 75% | 1 | 4.67606 | 5.14145 |
| 2 | 75% | 2 | 1.09634 | 1.14365 |
| 2 | 75% | 4 | 0.21439 | 0.22812 |
| 2 | 85% | 1 | 7.44324 | 7.85217 |
| 2 | 85% | 2 | 2.09045 | 2.24915 |
| 2 | 85% | 4 | 0.47913 | 0.51422 |
| 3 | 75% | 1 | 1.52817 | 2.17056 |
| 3 | 75% | 2 | 0.35400 | 0.53711 |
| 3 | 75% | 4 | 0.10681 | 0.16327 |
| 3 | 85% | 1 | 1.95160 | 2.79999 |
| 3 | 85% | 2 | 0.52567 | 0.67062 |
| 3 | 85% | 4 | 0.11749 | 0.16632 |
| 4 | 75% | 1 | 0.27313 | 0.27313 |
| 4 | 75% | 2 | 0.06866 | 0.07782 |
| 4 | 75% | 4 | 0.01678 | 0.01678 |
| 4 | 85% | 1 | 0.01144 | 0.01144 |
| 4 | 85% | 2 | 0.00076 | 0.00076 |
| 4 | 85% | 4 | 0.00076 | 0.00076 |
| 5 | 75% | 1 | 0.18158 | 0.08240 |
| 5 | 75% | 2 | 0.06485 | 0.02823 |
| 5 | 75% | 4 | 0.01221 | 0.01068 |
| 5 | 85% | 1 | 0.00000 | 0.00000 |
| 5 | 85% | 2 | 0.00000 | 0.00000 |
| 5 | 85% | 4 | 0.00000 | 0.00000 |

Input refinement 1 uses the original 3^3 component nodes per patch;
2 and 4 use 5^3 and 9^3 nodes on the same fixed box.
Topology comparisons use Euler, component count and component/face incidence
under both digital adjacencies; agreement is not a continuum certificate.
