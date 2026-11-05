# Historical experiment 04 outputs

These CSVs and figures predate the review corrections and remain unchanged.
See the [experiment README](../README.md) for the current interpretation.

- The main sweep has five realisations per setting. Its good q=16, VF=60%
  agreement does not establish an under-2% bound across the whole sweep.
- Ceiling CSV `delta_physical` is a margin in Gaussian-component value
  units, and `delta_voxels` is that margin divided by spatial grid spacing.
  Neither is a spatial film-thickness measurement. Historical plots based
  on the latter must not be used to choose a minimum grid resolution.
- Ceiling spanning is per-winner, 6-neighbour spanning. Phase A has ten
  independent fields per q, reused across ten VF thresholds.
- Zero observed spanning is finite-sample evidence, not a guarantee of
  disconnected continuum aggregates.
- In the main sweep, `n_components` and `largest_fraction_phase` are zero
  placeholders, not measured component statistics. Its convergence whiskers
  pool variation across thresholds and realisations, not just sampling noise
  at a fixed threshold.

The corrected ceiling runner writes `results_ceiling_v2/`, with explicit
value-margin and `dx_over_lc` columns, preserving this directory.
