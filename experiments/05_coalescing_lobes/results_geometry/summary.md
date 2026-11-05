# Pre-contact tetrahedron geometry sweep

Status: complete; 36 completed 3D cases.
rho is the value-space centre radius; r/rc is the morphology control after rho is fixed.
The four balls share one Gaussian vector: they are not four independent realizations.

| rho | r/rc | repetitions | expected VF | measured VF | expected L0 | ball sum L0, 6 / 26 | any ball spans, 6 / 26 |
|---|---|---|---|---|---|---|---|
| 1.5 | 0.85 | 4 | 0.35507 | 0.35867 | 170.24 | 169.00 / 165.75 | 1.00 / 1.00 |
| 1.5 | 0.95 | 4 | 0.48289 | 0.48767 | 84.92 | 87.00 / 84.00 | 1.00 / 1.00 |
| 1.5 | 0.98 | 4 | 0.52554 | 0.53070 | 57.88 | 61.00 / 61.75 | 1.00 / 1.00 |
| 2 | 0.85 | 4 | 0.42065 | 0.41961 | 172.34 | 168.75 / 163.75 | 1.00 / 1.00 |
| 2 | 0.95 | 4 | 0.58656 | 0.58672 | 89.84 | 80.75 / 81.75 | 1.00 / 1.00 |
| 2 | 0.98 | 4 | 0.64275 | 0.64294 | 61.68 | 63.75 / 57.75 | 1.00 / 1.00 |
| 2.5 | 0.85 | 4 | 0.41225 | 0.40842 | 189.32 | 178.50 / 174.50 | 1.00 / 1.00 |
| 2.5 | 0.95 | 4 | 0.60443 | 0.60127 | 106.75 | 107.50 / 101.00 | 1.00 / 1.00 |
| 2.5 | 0.98 | 4 | 0.67133 | 0.66842 | 76.41 | 80.25 / 76.00 | 1.00 / 1.00 |

Spanning frequencies use independent vector realizations as trials.
Zero spans in four trials is weak evidence, not proof of a disconnected continuum morphology.
Euler is not a component count. Ball sums are continuum-additive before contact,
but their digital estimators can still miss within-ball topology.
Components include boundary-cut and spanning objects; flags permit explicit filtering.
Equivalent diameter is volume-based, not a sieve diameter.
A single N=100 sweep is screening, not a resolution or continuum validation.
aggregate.csv uses standard errors across realizations; SEM is NaN with one realization.
