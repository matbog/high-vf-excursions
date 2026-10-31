# Measurement conventions

How we estimate the four Lipschitz-Killing curvatures $L_0 \ldots L_3$ of
a voxel realization, which choices carry known bias, and when those
biases matter for comparing against the Adler expectation. All
estimators live in `src/measures.py`.

The GKF predicts expectations in the continuum. A single voxel
realization gives you a discretized set, and every estimator imposes
some structure on how that voxel set is turned into a geometric
quantity. The convention matters: on thin films or near-tangent
surfaces, two defensible estimators can disagree by a factor of two or
more.

## $L_3$ — Volume

Trivial. $L_3 = (\text{\# True voxels}) \cdot dx^3$. Exact on the voxel
grid.

## $L_2$ — Half surface area

Two estimators are implemented. They answer slightly different questions
and have different biases.

**Marching cubes on a binary mask** (`measure_surface_area`): zero-pads
the phase, runs `skimage.measure.marching_cubes` at level 0.5, sums the
triangle areas, divides by two. The classic amplitude bias is roughly
1.088 on an isotropic sphere surface — individual triangles slightly
underestimate curved area and the net surface sits half a voxel offset
from the true level set. At our grid sizes ($N \geq 128$) this bias is
about 1–3% on Gaussian excursions.

**Marching cubes on a signed field** (`measure_surface_area_from_field`):
takes a smooth scalar field whose zero level set is the true surface.
Vertices sit on the continuum level set rather than on voxel edges.
This is what we use whenever the signed field is available (margin
tessellation, coalescing lobes) because it eliminates the sub-voxel
snapping error.

**Crofton with analytic face term** (`measure_surface_area_crofton`):
counts voxel-pair transitions in the three axis directions and
reconstructs

$$
L_2^{\text{interior}} = \frac{(N_x + N_y + N_z)\, dx^2}{6}
$$

by the stereological Cauchy relation (unbiased for isotropic surfaces).
The domain-face contribution (voxels on $\partial M$ that belong to the
phase) is added analytically to match the GKF's $3 a^2 \mathcal M_0$
boundary term.

On isotropic Gaussian excursions, Crofton is unbiased to <0.1% at
$N \geq 128$ on a sphere; marching-cubes-on-binary carries the ~1.09
amplitude bias. On polyhedral surfaces (e.g. the matrix films of the
margin tessellation), marching-cubes-on-the-signed-field is the better
estimator because Crofton's 3-direction isotropy assumption fails along
axis-aligned facets.

Default policy in experiments: Crofton for Gaussian baselines where the
surface is isotropic; signed-field marching cubes for margin and lobe
models.

## $L_1$ — Integrated mean curvature $/\,\pi$

**Dihedral-angle estimator** (`measure_l1`, `measure_l1_from_field`). For
a consistently-oriented PL surface mesh with directed edges,

$$
\int_S H\, dA \;=\; \tfrac{1}{2} \sum_{\text{edges}} \lvert e\rvert\, \theta_e,
$$

where $\theta_e$ is the signed dihedral angle at edge $e$ (positive
convex, negative concave). We extract the mesh from
`skimage.marching_cubes`, pair each half-edge with its reverse, compute
dihedrals from the two face normals, and return the result divided by
$\pi$ to match the convention $L_1 = (\int H\, dA)/\pi$ used in the GKF
formula.

Both variants exist for the same reason as for $L_2$: `measure_l1` runs
marching cubes on the binary phase mask;
`measure_l1_from_field` runs it on the smooth signed field.

**Known bias.** Binarizing before meshing loses sub-voxel interface
positions, which produces a residual amplitude bias of roughly $1.7\times$
at $N = 64$ on isotropic Gaussian surfaces (verified against the GKF
expectation). Passing the signed field directly removes most of this —
the amplitude ratio drops to ~1.00 at $N = 128$. We have not
characterized a convergence rate; treat `measure_l1` as a sanity-check
estimator and `measure_l1_from_field` as the production one when the
signed field exists.

## $L_0$ — Euler characteristic

Two connectivities are supported via the `connectivity` argument to
`measure_euler`, which wraps `skimage.measure.euler_number`:

- `connectivity=1` → 6-connectivity (face-adjacent voxels).
- `connectivity=3` → 26-connectivity (face + edge + vertex adjacency).

These can disagree substantially on morphologies with thin films or
sub-voxel features. In the coalescing-lobes experiment near contact,
6-conn produced large negative Euler values that looked like a
topological transition but were an artifact of thin films being read as
tunnels in 6-conn; 26-conn gave the correct count.

There is no universally "right" convention. Our rule of thumb:

- **Gaussian / Student-$t$ bulk excursions**: either connectivity gives
  near-identical results away from $L_0 = 0$. We use 6-conn
  (`connectivity=1`, the `measure_all` default) because it matches
  `skimage.marching_cubes` by construction.
- **Thin-film morphologies** (margin tessellation at high $q$, lobes near
  contact): use 26-conn for comparison to GKF, since the matrix film is
  often sub-voxel and 6-conn systematically miscounts handles.

The GKF predicts a continuum $\chi(A)$ that is agnostic to voxelization.
Pick whichever estimator matches the morphology's natural scale.

## Percolation

`measure_percolation` reports per-axis spanning: whether any 6-connected
component touches both opposite domain faces along $x$, $y$, or $z$.
`measure_percolation_per_winner` is the margin-tessellation refinement —
spanning of a single winner label, not the union.

Percolation is a *finite-domain* property. It is not predicted by the
GKF; it is a geometric fact about one realization. The relevant
GKF-adjacent quantity is the volume fraction $\varphi_c$ where
$\mathbb{E}[L_0] = 0$, which we call $\varphi_c^+$ throughout. The two
thresholds ($\mathbb{E}[L_0] = 0$ and $P_{\text{span}} = 1/2$) coincide
only in the ergodic limit $L/l_c \to \infty$.

## When this matters

For single-realization comparisons to the Adler expectation on $N = 128$
grids with $L/l_c \geq 8$, every estimator above agrees with theory to
within a few percent. The failure cases where the choice matters:

- Thin matrix films ($\leq 2$ voxel widths) — pick 26-conn Euler and
  signed-field marching cubes.
- Near-tangent isosurfaces (field barely crosses $0$) — Crofton
  underestimates because the Cauchy isotropy assumption fails;
  signed-field MC handles it.
- Reporting mean curvature on a binary field at $N \leq 64$ — don't; use
  `measure_l1_from_field` with the signed field.
