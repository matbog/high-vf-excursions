# Local 3D patch reconstruction

Status: complete.
Six purposefully selected patches by default; one realization, paired thresholds.
All topology is digital, on bounded patches. No global LKC or continuum-topology validation.
Volume uses fixed-domain trapezoid quadrature; errors below are percentage points of patch volume.

| Patch | Target VF | Direct local VF | Trilinear components error | Tetra components error | Scalar margin error |
|---|---|---|---|---|---|
| 0 | 75% | 30.7884% | 5.5824 | 6.9847 | 42.6987 |
| 0 | 85% | 62.3909% | 7.7393 | 9.4467 | 38.5323 |
| 1 | 75% | 69.5866% | 0.1984 | 1.1703 | 22.1092 |
| 1 | 85% | 82.3139% | 0.1602 | 1.2146 | 16.8968 |
| 2 | 75% | 24.8127% | 4.6761 | 5.1414 | 13.7634 |
| 2 | 85% | 50.1438% | 7.4432 | 7.8522 | 33.5495 |
| 3 | 75% | 42.6308% | 1.5282 | 2.1706 | 19.6632 |
| 3 | 85% | 62.8990% | 1.9516 | 2.8000 | 26.2878 |
| 4 | 75% | 99.0822% | 0.2731 | 0.2731 | 0.2731 |
| 4 | 85% | 99.9886% | 0.0114 | 0.0114 | 0.0114 |
| 5 | 75% | 99.5674% | 0.1816 | 0.0824 | 0.1816 |
| 5 | 85% | 100.0000% | 0.0000 | 0.0000 | 0.0000 |

Errors measure spatial phase disagreement, not just net volume bias.
The scalar margin baseline uses trilinear-component winners for label bookkeeping only.
Topology and refinement diagnostics are in measurements.csv and winners.csv.
A stable digital signature is not certified continuum topology; small or tangential features can remain.
Patch-face connections are local and do not imply spanning of the original domain.
