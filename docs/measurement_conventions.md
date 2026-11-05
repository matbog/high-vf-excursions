# Measurement conventions

The GKF refers to a continuous set inside a finite cube. Current experiments
estimate its LKCs from sampled fields. Always record the estimator and
connectivity alongside a result.

## Volume and grid geometry

`GridSpec.dx = length/n`. Voxel VF is the fraction of occupied samples;
$L_3=\mathrm{VF}\,L^3$ treats each sample as a cube of volume $dx^3$.

The VTI exporter stores these samples as point data, with point coordinates
spanning $(N-1)dx$. ParaView interpolation, voxel occupancy, and a padded
marching-cubes mesh therefore do not define identical geometry. Padding
closes a mesh, but its boundary position is not exact clipping at the
theoretical domain faces. The padded VTP export also needs a coordinate
offset to align it exactly with the unpadded VTI. A consistent interpolation
framework is deferred; a smooth rendering alone is not validation.

## Surface, $L_2$

- `measure_surface_area` meshes a binary phase at 0.5 and returns half
  the triangle area.
- `measure_surface_area_from_field` meshes a sampled signed field at the requested
  level and returns half the triangle area.
- `measure_surface_area_crofton` counts internal inside/outside transitions
  along three coordinate directions and adds occupied domain-face area:

$$
\widehat L_2=\frac{(N_x+N_y+N_z)\,dx^2}{3}
+\frac12 A_{\mathrm{occupied\ domain\ faces}}.
$$

The previous documentation's denominator 6 for the internal transition term
was wrong. The separate face term gives $3L^2$ for a fully occupied cube.
Directional averaging relies on isotropy for the random internal interface;
the box boundary is handled explicitly rather than by isotropic averaging. Finite sampling
can miss multiple crossings between adjacent points, so the formula is
not an unbiased guarantee for every grid.

Margin chambers are polyhedral in **Gaussian value space**. Their spatial
preimages need not have planar or axis-aligned interfaces. Missing thin
films is a sampling problem, not evidence that ensemble isotropy fails.

## Curvature, $L_1$

For an oriented closed piecewise-linear surface,

$$
L_1=\frac{1}{2\pi}\sum_e |e|\theta_e,
$$

with signed exterior dihedral angle $\theta_e$ (positive for convex bends).
`measure_l1` uses a binary marching-cubes mesh;
`measure_l1_from_field` uses the sampled scalar field. Passing scalar values
retains some interface-position information lost on binarization, but does
not remove all grid, padding, nonsmooth-field, or topology errors.
No general error bound or convergence rate has been established here.

For a union of balls in value space, experiment 05 meshes
$s(x)=r-\min_a\|G(x)-c_a\|$. This captures the sampled union boundary,
including its crease geometry only to the accuracy of the mesh.
Interpolating sampled $s$ is not equivalent to interpolating $G$ and
then taking distances or performing exact Boolean clipping.

## Euler characteristic and connectivity

`measure_euler` supports 6-neighbour (`connectivity=1`) and 26-neighbour
(`connectivity=3`) digital conventions. They can differ sharply when
features are near the grid scale. Marching cubes is not generally identical
to either digital convention, and choosing whichever agrees with theory
is not a convergence test.

`measure_percolation` reports opposite-face spanning, with 6 neighbours
by default or 26 when explicitly requested. Margin per-winner spanning
tests each label with 6 neighbours and ignores cross-label contacts;
it is a different statistic from spanning of the binary union.

## Actual experiment policies

| Experiment | Euler / topology | Surface and curvature |
|---|---|---|
| 01–03 | Binary phase, default 6 neighbours | Binary marching cubes via `measure_all(compute_l1=True)` |
| 04 | Sum of per-winner 26-neighbour Euler; per-winner spanning uses 6 neighbours | Sampled margin-field $L_1$; Crofton $L_2$ with scalar-mesh comparison |
| 05 corrected runner | Binary-union Euler, components, and spanning under both 6 and 26 neighbours | Union scalar-mesh $L_1,L_2$ and union Crofton $L_2$; ball/partition controls stored separately |
| 06 voxel runner | Binary-union Euler/spanning under both conventions | Numerical $L_1,L_2$ are omitted; their expectations enter the union recursion |

Experiment 05's historical `results/` mixed nearest-centre partition Euler
and Crofton with sums of individual-ball mesh measurements. After contact,
these are not union LKCs. The corrected runner writes `results_union_v2/`.

## Comparing with theory

Report ensemble means, dispersion or standard errors, and discretisation
effects separately. Avoid relative errors where a theoretical LKC is near
zero. Check physical domain size, component covariance, units, boundaries,
and grid resolution before interpreting a discrepancy as a theory failure.
Neither $N\geq128$ nor selecting 26 neighbours guarantees percent-level
agreement or continuum connectivity.
