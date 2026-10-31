# Interpolated geometry in experiment 07

**Deferred on 31 October 2026 because of computational cost.** The working
experiment remains the voxel workflow. This document and the prototype are kept
for later research; see the final section of [README.md](README.md) and the
[stopped full-run notes](results_pl_full/STOPPED.md). The commands below describe
the retained prototype, not a currently running job.

`run_interpolated.py` is a separate, experimental measurement path for **L0,
L3, connected components and spanning**. It uses GSTools Gaussian components
and constructs margin regions between samples. `run.py` and its saved FFT/voxel
results remain the historical baseline. All new code, configuration and outputs
are inside experiment 07.

## What is interpolated

N denotes **cells per axis**, with N+1 sample nodes including both domain faces.
For N=100 and L=64 mm, the nodes span exactly [0,64] and dx=0.64 mm. This removes
the mismatch between an N-point VTK array spanning only (N-1)dx and an N-voxel
measurement spanning N dx.

Each cube is split into six conforming tetrahedra using the same body diagonal.
Each original Gaussian component is affine inside each tetrahedron. For class k
and winner i, construct the closed solid

```math
C_{ki}=\bigcap_{j\ne i}\{\widetilde G_{ki}-\widetilde G_{kj}\geq\delta_k\}.
```

The tildes denote these component interpolants. Sequential scalar clipping of
the tetrahedra implements every inequality, interpolating the original component
values onto new vertices. All competitors are considered: a winner region can
exist inside a tetrahedron without winning at any original vertex.

We never interpolate a precomputed top-two gap or a binary mask. For example,
G1=x and G2=-x on [-1,1] have margin 2 at both endpoints but zero at the centre.
Interpolating those two margin samples fills the separating film; interpolating
the two components retains it. The automated tests include this case and a
film only 2e-11 of a domain width thick.

## One geometry for volume and topology

Each clipped region is a tetrahedral volume, including its boundary on the
observation cube. No artificial outside padding is used to close the surfaces.
For each region:

- L3 is the sum of tetrahedron volumes.
- L0 is V-E+F-T of the volumetric complex.
- Connectivity is shared mesh vertices, rather than a choice of 6 or 26 voxel
  neighbours.

Every nonempty term must pass audits: no face belongs to more than two volume
cells; every boundary edge belongs to two boundary triangles; and the boundary
Euler characteristic equals twice the volume Euler characteristic. A failed
audit stops the run, records the error in metadata, and saves a failing mesh
when mesh export is enabled. Completed earlier cases remain checkpointed.

Exactly flat tetrahedra emitted by clipping are removed before these audits.
Positive-volume slivers are retained. Vertex cleaning defaults to **zero
tolerance**, merging only identical positions. A positive `merge_tolerance`
can collapse thin regions and should be treated as a separate numerical change.

## Multiscale union and real connections

Positive margins make different winner regions of one class disjoint. Cross-class
regions can overlap. For at most three classes, enumerate all nonempty single,
pair and triple winner-region intersections and apply inclusion-exclusion:

```math
L_j(A)=\sum_{\varnothing\ne S\subseteq\{1,2,3\}}
(-1)^{|S|+1}\sum_{(i_k)_{k\in S}}
L_j\!\left(\bigcap_{k\in S}C_{ki_k}\right),\qquad j\in\{0,3\}.
```

These are measurements of the PL geometry, not GKF predictions. The same
intersection terms give coarse-first added class volumes without double counting.

Raw connected regions are nodes of an intersection graph. A positive-volume
cross-class intersection joins their nodes. Inclusion-exclusion also supplies
each final component's volume. Domain-face contacts propagate through the graph
to measure actual spanning. A coarse winner label does not prevent a connection
through another scale.

This implementation assumes regular thresholds and full-dimensional solids.
Exact tangency-only contacts and lower-dimensional isolated regions are not
handled as separate topological features. Such singular cases require a more
general arrangement implementation; random regular cases are the intended use.

## GSTools and Adler reference

The local sampler matches the new shared backend's covariance and engine:

```python
model = gstools.Gaussian(dim=3, len_scale=lc, rescale=1.0, var=1.0, nugget=0.0)
srf = gstools.SRF(model, generator="RandMeth", mode_no=1000, seed=component_seed)
```

Thus the target covariance is exp[-(r/Lc)^2], with nominal derivative variance
2/Lc^2. There is no per-realization centering or normalization. The existing
shared wrapper uses voxel-oriented coordinates; this experiment supplies its
own boundary-inclusive coordinates to the same GSTools engine. It does not
modify `src/fields.py`.

Seeds depend on master seed, realization, class and component, **not grid size**.
Re-evaluating with the same mode count on a refined grid samples the same GSTools
realization. This is tested at coincident nodes. The historical FFT runs do not
have this pairing property. See [GSTools RandMeth](https://geostat-framework.readthedocs.io/projects/gstools/en/stable/api/gstools.field.generator.RandMeth.html)
for its finite random Fourier-mode construction.

The theoretical columns retain the existing margin GMFs and conditional GKF
with nominal derivative variance. They describe the intended smooth Gaussian
model. They are **not claimed exact for the finite-mode, piecewise-linear
interpolant**. At a point inside a tetrahedron, interpolation changes variance
and derivative statistics. On a coarse mesh it can markedly lower measured VF
even though the continuum threshold calibration is correct. Fixed-grid
triangulation also introduces directional effects. Refinement and mode-count
sensitivity remain necessary, as does independent validation of the margin GMFs.

## Run and resource use

With `rf_explo` active, from the repository root:

```bash
python -B experiments/07_multiscale_margin/test_interpolated.py
python -B experiments/07_multiscale_margin/run_interpolated.py --n 12 --q 4 --realizations 1 --phis .6 .8 --output experiments/07_multiscale_margin/results_pl_n12
python -B experiments/07_multiscale_margin/run_interpolated.py --n 24 --q 4 --realizations 1 --phis .6 .8 --output experiments/07_multiscale_margin/results_pl_n24
```

Those are small computational pilots, not concrete-mixture calibration runs.
The paired N=12 and N=24 commands retain the same seeds and 1000 modes per
component. N=24 still has dx=2.67 mm, larger than the sand-proxy Lc=2 mm.

`config_interpolated.json` retains the requested three-scale parameters and
N=100, q=16, five realizations, and five VF targets. The full configuration is
available through:

```bash
python -B experiments/07_multiscale_margin/run_interpolated.py --output experiments/07_multiscale_margin/results_pl_full
```

This is a correctness-first implementation. It retains component arrays and
clipped meshes, and enumerates up to q1*q2*q3 triple intersections. At N=100 the
base alone has six million tetrahedra; storing all component data and intermediate
meshes can require several GB. The full run has **not** been benchmarked. Start
with the small pilot and increase resolution/competitor count deliberately.
GSTools, VTK and numerical-library thread counts are limited to one in this
runner, so it can be used alongside other work.

Use `--no-meshes` or `--no-fields` to reduce saved files, and `--mode-no` for a
separate spectral-mode sensitivity study. Outputs must remain in a subdirectory
of experiment 07. Existing output folders require `--overwrite`, which only
accepts directories carrying this runner's metadata and replaces all their old
files to avoid stale mesh blocks.

## Outputs and ParaView

- `summary.md` and `comparison.png`: interpolated geometry versus voxel
  measurements from **the same PL component fields** at cube centres, alongside
  the continuum GKF reference.
- `geometry_slices.png`: optional matched physical slices, generated with
  `preview_interpolated.py` from the saved component fields and clipped volumes.
- `measurements.csv`: measured union L0/L3, connected components, largest-volume
  fraction, spanning, and the voxel comparison.
- `class_volumes.csv`: corrected coarse-first volume contributions.
- `intersection_terms.csv`: every nonempty measured term and its sign.
- `union_components.csv`: actual union-component volumes, spanning and boundary
  flags, and equivalent diameters.
- `raw_regions.csv`: raw winner-region measurements.
- `fields_r*.npz`: the original GSTools nodal components, for reinspection.
- `metadata.json` and `field_diagnostics.csv`: parameters, versions, source
  hashes, seeds, field moments, derivative diagnostics and completion status.
- `meshes_r*_vf_*.vtm`: open directly in ParaView; references the clipped
  class/winner `.vtu` volumes. Surface representation shows the interpolated
  geometry, and its outer-domain faces are already closed.

**The VTM blocks overlap across classes.** Display them together as a set union,
or inspect classes separately. Summing their volumes with ParaView's Integrate
Variables would double-count overlaps; the CSVs apply the required intersection
corrections. This version does not export a single conforming Boolean-union mesh
for downstream finite-element analysis. It also does not measure L1/L2; only
their theoretical values enter the conditional-GKF reference. Both are explicit
remaining extensions, separate from the implemented L0/L3 measurement path.

The relevant VTK operation is [scalar clipping](https://vtk.org/doc/nightly/html/classvtkTableBasedClipDataSet.html),
which cuts cells and retains partial volumes. Merely making a rendered surface
appear smooth is not used as a geometry or topology check.

## Completed validation runs

All **12 analytical controls** in `test_interpolated.py` pass. They cover exact
cube and plane volumes, a tunnel with Euler zero, hidden winner regions,
separating films, cross-class bridges, triple-overlap corrections, 16 competitors,
class-order invariance, VTK round trips, paired GSTools sampling, and independent
Monte Carlo integration of a synthetic three-class union.

The paired N=12 and N=24 pilots each use one realization, q=4 per class, 1000
modes/component, and targets 60% and 80%. They took approximately 36 and 199
seconds, respectively, on this machine with one thread. The saved fields agree
exactly at common nodes. Metadata, source hashes, class/component volume sums,
Euler intersection sums and ParaView multiblock readback all pass validation.

| Cells/axis | Target VF | Interpolated VF | Interpolated L0 | Centre-voxel L0 (6 / 26) |
| --- | --- | --- | --- | --- |
| 12 | 60% | 36.794% | 14 | 60 / -111 |
| 24 | 60% | 45.712% | 566 | 97 / -617 |
| 12 | 80% | 61.707% | -765 | -103 / -44 |
| 24 | 80% | 69.337% | -2076 | -583 / -246 |

These results demonstrate that voxel and interpolated topology differ, while
the same-field volume estimates are much closer. They do not establish
convergence to the smooth-field theory: neither grid resolves the sand scale.
Both interpolated pilot unions span the domain at both targets. This is not a
direct comparison with the old q=16 FFT experiments.

Inspect the [N=24 summary](results_pl_n24/summary.md),
[comparison plots](results_pl_n24/comparison.png), and
[matched slices](results_pl_n24/geometry_slices.png). For ParaView, open
[the 60% target geometry](results_pl_n24/meshes_r1_vf_0.6000.vtm) or
[the 80% target geometry](results_pl_n24/meshes_r1_vf_0.8000.vtm).
To regenerate the slice comparison:

```bash
python -B experiments/07_multiscale_margin/preview_interpolated.py experiments/07_multiscale_margin/results_pl_n24
```
