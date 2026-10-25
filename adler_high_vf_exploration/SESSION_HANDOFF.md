# Session handoff

## Research objective

Develop Gaussian-related random morphologies that retain Adler's Gaussian
kinematic formula while allowing either unusually high disconnected inclusion
volume fractions or controlled component coalescence.

## Model 1: Gaussian margin tessellation

Generate `q` independent smooth Gaussian fields and form

    H(x) = largest_i G_i(x) - second_largest_i G_i(x),
    A_delta = {H(x) >= delta}.

The hitting set is the union of `q` separated polyhedral winner chambers. This
produces persistent matrix films and high disconnected volume fractions, but
different winner families cannot merge for any `delta>0`. Thus this model does
not have gradual cross-family coalescence.

Closed candidate formulas for the four Gaussian Minkowski functionals and the
bounded-cube Adler expectations are recorded in `THEORY_NOTE.md` and
implemented in `gaussian_minkowski.py`. Numerical evidence supports them, but a
formal stratified-tube proof remains necessary. The numerical `L1` estimator is
not yet trustworthy for complex surfaces.

## Model 2: coalescing Gaussian-space lobes

The complementary model introduced in this session is

    R(x) = min_a ||G(x)-c_a||,
    A_r = {R(x) <= r} = G^{-1}(D_r),
    D_r = union_a B(c_a,r).

Here `G(x)` is a vector of Gaussian RF values. In the first experiment it has
three components, so each physical voxel maps to one point in auxiliary
Gaussian space `R^3`. Four fixed targets `c_a` are placed at tetrahedron
vertices. A voxel is inclusion when its Gaussian-space point lies within radius
`r` of at least one target.

For tetrahedral centre norm `rho`,

    pairwise centre distance = rho*sqrt(8/3),
    first ball-contact radius r_c = rho*sqrt(2/3).

Increasing `r` grows the four Gaussian-space balls. Before `r_c`, different
lobe families are disjoint, although the spatial preimage of a single lobe can
already connect or span. At `r_c`, different hitting-set lobes first touch; for
larger `r` they overlap. The ball union is stratified, so the morphology remains
within the Adler/GKF framework.

## Plain-language construction

1. Make a physical `100 x 100 x 100` grid.
2. Generate three smooth correlated Gaussian RFs `G1,G2,G3`.
3. At each voxel collect their values into the auxiliary point
   `G(x)=(G1(x),G2(x),G3(x))`.
4. Place four symmetric target points in this auxiliary space.
5. Compute the distance from `G(x)` to the closest target.
6. Retain the voxel if that distance is at most `r`.

The Gaussian-space balls are selection rules, not physical spheres. Their
preimages are irregular spatial inclusions. The four plot colors identify the
nearest target; all colors represent the same inclusion phase.

## Coalescing-lobe pilot

Script: `coalescing_lobes.py`

Parameters:

- bounded central observation window `100^3`;
- larger FFT synthesis followed by cropping;
- three Gaussian RF components;
- four tetrahedral centres;
- `rho=1.5`, `Lc=10` voxels;
- four realizations;
- `r/r_c = 0.45, 0.65, 0.85, 0.98, 1.02, 1.15, 1.35, 1.60`.

Outputs are under `results/coalescing_lobes/`:

- `coalescing_lobes_slices.png`;
- `coalescing_lobes_diagnostics.png`;
- `coalescing_lobes_sweep.csv`.

Pilot ensemble means:

| r/r_c | volume | Euler | largest/phase | spanning frequency |
|---:|---:|---:|---:|---:|
| 0.45 | 0.056 | 340 | 0.051 | 0 |
| 0.65 | 0.165 | 288 | 0.055 | 0 |
| 0.85 | 0.355 | 144 | 0.311 | 1 |
| 0.98 | 0.526 | -2870 | 0.999 | 1 |
| 1.02 | 0.584 | -2048 | 1.000 | 1 |
| 1.15 | 0.732 | 220 | 1.000 | 1 |
| 1.35 | 0.852 | 14 | 1.000 | 1 |
| 1.60 | 0.933 | 62 | 1.000 | 1 |

Main preliminary observation: physical spanning occurred before Gaussian-space
lobe contact. Near `r/r_c=1`, the voxel Euler characteristic showed an extreme
negative dip, suggesting a tunnel-rich transition, but its magnitude must be
checked for resolution sensitivity near the singular contact radius.

## Potential paper question

Compare three transitions rather than assuming they coincide:

1. prescribed contact/merger of the hitting set in Gaussian value space;
2. zero crossing of the expected Euler characteristic;
3. finite-size physical-space percolation (`P_span=1/2` or periodic wrapping).

A systematic separation between these transitions could quantify how hitting-
set topology transfers to physical-space connectivity.

## Recommended next work

1. Densely sample radii around both the apparent spatial transition and
   `r/r_c=1`.
2. Run at least three voxel resolutions at fixed physical `T/Lc`.
3. Use many more realizations and confidence intervals near transition points.
4. Replace face spanning with periodic wrapping for percolation estimates,
   while retaining bounded domains for direct Adler LKC comparison.
5. Decompose Euler characteristic into Betti numbers where practical.
6. Derive Gaussian Minkowski functionals of off-centre balls before contact,
   then treat ball intersections by Gaussian-weighted inclusion-exclusion.

