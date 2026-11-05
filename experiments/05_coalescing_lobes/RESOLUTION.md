# Paired resolution study — experiment 05

## Status

The topology-only pilot is complete: six cases at N=50,100,200 and
r/rc=0.85,0.98, using one paired realization. See
[the saved summary](results_resolution_topology/summary.md) and
[measurements](results_resolution_topology/measurements.csv).
The optional run with scalar-mesh L1/L2 and the four-realization extension
have not been completed. Further exploration is deferred.

The completed `results_union_v2/` run contains 32 cases at $N=100$.
It separates measurements of the union, individual balls, and nearest-centre
partitions. It reveals a large pre-contact discrepancy near $r/r_c=1$:

| $r/r_c$ | Expected $L_0$ | Mean ball-sum $L_0$ (26) | Mean union $L_0$ (6) | Mean union $L_0$ (26) |
|---|---:|---:|---:|---:|
| 0.85 | 170.24 | 165.75 | 148.25 | −1325.25 |
| 0.98 | 57.88 | 61.75 | −2973.25 | −2038.25 |

For $r<r_c$, different ball preimages cannot touch in the continuous model.
Digital contacts between their labels therefore indicate unresolved separation.
The discrepancy between the voxel union and the per-ball sums can be studied
on the **same realization** without confusing it with ensemble variability.

## Pilot

`run_resolution.py` uses:

| Parameter | Default |
|---|---|
| Physical domain side / component correlation length | $L=100$, $l_c=10$ |
| Resolutions | $N=50,100,200$ |
| Radius ratios | $r/r_c=0.85,0.98$ |
| Realizations | 1, using the first original stream |
| Tetrahedron circumradius | 1.5 |
| Gaussian components / spectral modes | 3 / 1000 modes per component |
| GSTools threads | 1 |
| Output | `results_resolution/` |

This is **six measurement cases**. The default includes scalar-mesh $L_1,L_2$.
It is a numerical-refinement pilot, not a new ensemble validation.
The saved completed pilot used `--topology-only` instead; the default
mesh-enabled command below describes an optional future computation.

From the repository root:

~~~powershell
mamba run -n rf_explo python -u experiments/05_coalescing_lobes/run_resolution.py
~~~

For a lighter first pass, omit the meshes:

~~~powershell
mamba run -n rf_explo python -u experiments/05_coalescing_lobes/run_resolution.py --topology-only --output experiments/05_coalescing_lobes/results_resolution_topology
~~~

This keeps voxel $L_0,L_3$, Crofton $L_2$, components, spanning, and contact
diagnostics. Mesh columns are explicitly NaN, not zero.
The saved `results_resolution_topology/` folder is already populated;
choose another `--output` to repeat this computation.

To repeat the study over the four original streams:

~~~powershell
mamba run -n rf_explo python -u experiments/05_coalescing_lobes/run_resolution.py --realizations 4 --output experiments/05_coalescing_lobes/results_resolution_4real
~~~

That command runs 24 cases. Review the one-realization pilot and its timings
before extending it. Arrays scale as $N^3$; $N=200$ contains eight times as
many samples as $N=100$. Mesh cost also depends on interface complexity.
Only one resolution's distance stack is retained at a time, along with the
smallest grid's component samples for the pairing check. No VTK volumes or
large field snapshots are exported.

## How pairing works

The sampler uses the same two-stage seed sequence as `src/fields.py` and
keeps the GSTools RandMeth modes fixed as $N$ changes. It evaluates the
components at

$$
x_i=x_0+iL/N,\qquad x_0=P_{\rm ref}L/N_{\rm ref}.
$$

Defaults $P_{\rm ref}=20,N_{\rm ref}=100$ give $x_0=20$ in physical units.
This matches the cropped observation window of the original $N=100$ run.

Keeping padding at 20 *grid cells* separately at every resolution would move
the observation window. The new sampler instead holds this physical origin
fixed. It uses physical $l_c$, no empirical centering or variance rescaling,
and the same sample/voxel convention as the previous experiment.

Common points are checked against the smallest grid; a maximum component
difference above $2\times10^{-5}$ stops the run. At the default resolutions
all smaller-grid sample locations occur on each finer grid.
Sampling in slabs limits temporary memory without changing the realization.
Tests compare both slab sizes and the shared legacy sampler on a small
reference grid.

The nominal GKF remains the smooth Gaussian-covariance expectation.
Finite spectral sampling, digital adjacency, scalar interpolation, and
domain closure still introduce approximation. Refinement alone does not
eliminate realization-to-realization fluctuation or validate every GKF
assumption.

## Diagnostics

`measurements.csv` retains the union, ball and partition columns used in
`union_v2` and adds:

- `sum_ball_euler_6`: the missing same-convention control for union Euler 6.
- `euler_additivity_gap_6/26`: union Euler minus individual-ball Euler sums,
  with the same convention on both sides.
- `cross_label_contacts_6/26`: number of undirected occupied neighbour pairs
  belonging to different balls.
- `occupied_neighbor_pairs_6/26` and `cross_label_fraction_6/26`: denominators
  and normalized fractions for those counts.
- `any_ball_spans_6/26`: whether an individual ball preimage spans.
  Compare with `union_spans_6/26` to distinguish within-ball spanning
  from connections that require multiple labels.
- `L1_additivity_gap` and `L2_crofton_additivity_gap`.
- Resolution, exact pre-contact theoretical VF/LKCs, and measurement time.

The runner accepts **only $0<r/r_c<1$**. Post-contact label crossings would
not have the same diagnostic meaning. At every case it checks that union
voxel volume equals the sum of the disjoint ball volumes.

`field_diagnostics.csv` records component seeds, means, second moments,
common-node errors, and sampling times. `metadata.json` records parameters,
dependency versions, source hashes, completion count, and status.

`resolution.png` shows the paired Euler curves, cross-label contact counts,
VF, and $L_1$ (or Euler additivity gaps in topology-only mode).
`summary.md` reports the completed rows. A partial run is explicitly marked.

## Completed pilot results

The maximum shared-node field discrepancy is exactly zero. N=100 reproduces
the corresponding first realization of `results_union_v2/` for the checked
volume and topology measurements.

| r/rc | N | VF | Union L0, 6 / 26 | Ball-sum L0, 6 / 26 | Cross-label pairs, 6 / 26 |
|---|---:|---:|---:|---:|---:|
| 0.85 | 50 | 0.352424 | -213 / -1487 | 177 / 165 | 577 / 11251 |
| 0.85 | 100 | 0.352188 | 165 / -1449 | 178 / 178 | 17 / 4329 |
| 0.85 | 200 | 0.352209 | 185 / 85 | 185 / 184 | 0 / 107 |
| 0.98 | 50 | 0.522120 | -1724 / -702 | 80 / 66 | 7806 / 63203 |
| 0.98 | 100 | 0.522078 | -2942 / -2040 | 78 / 71 | 13606 / 135367 |
| 0.98 | 200 | 0.521888 | -4577 / -5344 | 76 / 71 | 14786 / 200466 |

VF varies by less than 0.024 percentage points over the three grids at
each radius. At 0.85, face contacts disappear at N=200 and the 6-neighbour
Euler becomes additive. Diagonal contacts remain. At 0.98, the normalized
contact fractions decrease but the union Euler discrepancy does not:
topological convergence has not been demonstrated.

At least one individual ball spans under each convention in all six cases,
including the 0.85, N=200 case with no cross-label face contacts. This
distinguishes within-label spanning from spurious inter-label connections
on these grids; it is not a proof of continuum spanning.

The recorded generation and measurement work totalled about 26 minutes on
the machine used, mostly at N=200. Peak RAM was not logged.

## Interpretation and possible follow-up

Look for the Euler additivity gaps and sensitivity to connectivity to shrink
under refinement, while volume and per-ball measurements stabilize.
Inspect both raw cross-label pair counts and their fractions: counts alone
need not decrease monotonically, because the number of neighbour pairs grows
with $N$. Absence of observed contacts does not prove continuum convergence.

If the union spans while no single ball spans under the same convention,
cross-label digital contacts are necessary for that measured spanning.
An individual ball can also span before value-space contact; this remains
possible in the mathematical model.

Agreement of a single realization with expected $L_0$ is not the criterion.
Use the four-realization extension to assess the mean, after inspecting
resolution effects. Even then, four samples give limited uncertainty estimates.
No particular resolution, including $N=200$, is promised sufficient near
$r/r_c=0.98$.

The subsequent [centre-position sweep](GEOMETRY.md) has also completed.
It improves the achieved VF but still finds within-lobe spanning. See the
[experiment README](README.md#deferred-work) for the deferred research
questions; no larger resolution campaign is currently required or running.

## Checkpoints and tests

After each completed case the CSV, summary, plot, and metadata are updated.
Existing nonempty destinations are protected. To resume the default run:

~~~powershell
mamba run -n rf_explo python -u experiments/05_coalescing_lobes/run_resolution.py --resume
~~~

For a custom run, repeat its original arguments and add `--resume`.
Parameters, sources and dependency versions must match. A partially measured
case is recomputed; completed cases are retained. Use a new output directory
to change resolutions, realization count or mesh policy.

All seven small controls passed before publication, including shared-node
sampling, interruption/resume and a tiny CLI run:

~~~powershell
mamba run -n rf_explo python -B experiments/05_coalescing_lobes/test_resolution.py
~~~

The implementation is local to experiment 05 and imports the shared `src`
utilities. It does not depend on experiment 06 or launch its interpolation
prototype. Historical experiment 05 outputs are left in place.
