# Local directional film check

One saved realization; 75%/85% refer to target VF. Segments are sampled
conditionally on occupied face-neighbour endpoints, separately by label relation.
These are directional gap lengths, not normal film thicknesses or global statistics.

Seed validation: 1056 nodes; first-field max error 0, margin max error 0, winner mismatches 0.

| Target VF | Pair type | Segments | Stable references | Gaps: direct / affine / scalar | Direct gap / dx, median | Affine width error / dx, median / max |
|---|---|---|---|---|---|---|
| 75% | cross_label | 100 | 100 | 100 / 100 / 0 | 0.630299 | 0.004207 / 0.019719 |
| 75% | same_label | 100 | 100 | 0 / 0 / 0 | 0.000000 | 0.000000 / 0.000000 |
| 85% | cross_label | 100 | 100 | 100 / 100 / 0 | 0.438984 | 0.004454 / 0.027149 |
| 85% | same_label | 100 | 100 | 0 / 0 / 0 | 0.000000 | 0.000000 / 0.000000 |

The affine-component calculation clips all winner inequalities exactly along
each segment; it is not a linear interpolation of the already-derived margin.
The direct reference evaluates the finite-mode GSTools realization between nodes.
Nested probes and interpolated threshold crossings provide a numerical stability
check, not certified root isolation. Features narrower than the probes may remain.
Reference stability requires two successive equal-gap-count comparisons with
boundary motion below the configured tolerance, at least 65 samples, and no
unresolved winner switches between occupied consecutive probes.
Unstable references remain flagged in measurements.csv and are included in the
summary; they must not be treated as validated reference widths.
The scalar interpolation has no gap when both endpoint margins exceed delta.
All intervals and interpolation disagreements are measured along the segment.
Same-label controls do not establish 3D connectivity. No global VF or L0 is
re-estimated here; the Adler hitting set and underlying field model are unchanged.
