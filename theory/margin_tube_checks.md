# Independent checks of the margin tube derivation

The [derivation draft](margin_tube_derivation.md) is supported by deterministic
calculations in Gaussian value space. The saved run is dated 10 November
2026; this report was completed on 12 November. It generates no spatial
random fields. These checks address tube coefficients and asymptotics,
not digital LKC accuracy, spatial spanning, or aggregate grading.

## Methods and independence

[check_margin_tubes.py](check_margin_tubes.py) uses three calculations:

1. **Geometric coefficients.** Face weights are integrated over the whole
   real line with their multiplicities inside the integrand. Normal-cone
   angles and first moments are computed from the actual normal vectors,
   without inserting the production constants $\pi/3$ or
   $3\arccos(1/3)-\pi$. The production GMFs are imported only for comparison.
2. **Full cone tubes.** Radial and angular quadrature evaluates the exact
   face/normal decomposition for $q=2,3,4$ at finite radii below
   $\delta/\sqrt2$. This includes higher powers of radius, so it can check
   the remainder after the cubic GMF expansion.
3. **Distance tubes.** A separate calculation integrates the Gaussian
   probability of lying within a specified distance of the chambers.
   It uses neither face weights nor cone coefficients. It combines
   nearest-point projection with the radial $\chi_{q-1}$ distribution in
   the orthogonal contrast space.

For angular quadrature, columns of a matrix $V$ generate a cone. Normalizing
$Vw$ for simplex weights $w_j\geq0$, $\sum w_j=1$, maps the simplex onto its
spherical section. The angular Jacobian is
$\sqrt{\det(V^TV)}/\|Vw\|^k$ for $k$ generators, with respect to the first
$k-1$ simplex coordinates. For three generators, the parameterization
$w=(u,(1-u)v,(1-u)(1-v))$ adds a factor $1-u$. Tensor Gauss-Legendre
quadrature then integrates these smooth charts.

The distance calculation uses the sector $z_1\geq\cdots\geq z_q$ and
multiplies its angular integral by $q!$. Its extreme rays in contrast space
are proportional to $\mathbf1_{\{i\leq j\}}-j\mathbf1/q$, $j=1,\ldots,q-1$.
The normalized angular weights sum to one; this is checked separately.

In this sector the nearest chamber is the first one. To see this, take
any point $y$ in chamber $j$ and swap its first and $j$th coordinates.
The result belongs to chamber 1, and the squared distance decreases by
$2(z_1-z_j)(y_j-y_1)\geq0$. Projection onto chamber 1 reduces to choosing
its winning coordinate $t$; every other projected coordinate is
$\min(z_j,t-\delta)$. The active-set condition is

$$
t=\frac{z_1+\sum_{j\in J}(z_j+\delta)}{1+|J|},
\qquad J=\{j>1:z_j+\delta>t\}.
$$

The ordered competitors let the code find this set sequentially. Tests
compare the projection against an independent constrained optimizer.
Along each angular direction, bisection locates the radius where the
distance equals the requested tube radius. Convexity of the chamber makes
the radial sublevel interval start at this boundary and continue to
infinity. Its Gaussian radial probability is the corresponding chi survival
function. At tube radius zero, the entry radius is available directly as
$\delta/(u_1-u_2)$.

Normal-cone quadrature uses order 48. Distance quadrature uses order 128
and repeats at order 256, with 54 radial bisection steps. Differences after
refinement are numerical diagnostics, not certified error bounds.

## Results

All six numerical acceptance checks in
[summary.json](results_margin_tubes/summary.json) pass.

| Check | Scope | Largest absolute discrepancy |
|---|---|---:|
| Cone angles and first moment | Three geometric constants | $7.8\times10^{-16}$ |
| GMFs against production | $q=2,3,4,8,16,32,64$; five margins; four GMFs, 140 comparisons | $4.58\times10^{-9}$ |
| Distance tubes against full cone tubes | $q=2,3,4$; three margins; four radii, 36 comparisons | $1.19\times10^{-9}$ |
| Distance quadrature, order 128 to 256 | Same 36 comparisons | $1.04\times10^{-8}$ |

The GMF margins are 0.01, 0.1, 0.3, 0.6, and 1.2. The distance-tube margins
are 0.3, 0.6, and 1.2, with radii 0, 0.15, 0.35, and 0.65 times
$\delta/\sqrt2$. All finite-radius checks remain strictly below tube contact.

The production integrator applies absolute tolerance $2\times10^{-11}$
to an integral before multiplying by a combinatorial factor. An initially
stricter $10^{-9}$ comparison therefore flagged a few-billionths difference
at $q=32$. Recomputing its three face integrals at $\delta=0.01$ with the
same finite bounds but tighter tolerance agrees with the independent
whole-line calculation within $3.6\times10^{-15}$; see
[quadrature_precision.csv](results_margin_tubes/quadrature_precision.csv).
The saved check uses a $10^{-8}$ absolute production-comparison tolerance
and records the tighter audit separately. This gives no evidence of an
incorrect cone constant or sign, and required no production formula change.

At $\delta=0.6$, the full-tube minus cubic-polynomial residual decreases
approximately as $r^4$ for both $q=3$ and $q=4$. The normalized residuals
are recorded in [cubic_remainders.csv](results_margin_tubes/cubic_remainders.csv).
This checks the GMF factorial convention; the separate distance calculation
is what checks the full face decomposition independently.

## High-volume-fraction results

The positive-margin limiting expected Euler values below are for a cube
with $L=64,l_c=5$, under the GKF assumptions in the derivation. They are
expectations, not observed component counts.

| $q$ | $c_q=\mathbb E\max_i Z_i$ | $\lim_{\delta\downarrow0}\mathbb E L_0$ |
|---:|---:|---:|
| 2 | 0.564190 | -282.2144 |
| 3 | 0.846284 | -345.5938 |
| 4 | 1.029375 | -253.8113 |
| 8 | 1.423600 | 517.8421 |
| 16 | 1.765991 | 2067.9321 |
| 32 | 2.069669 | 4328.5544 |
| 64 | 2.343733 | 7218.6965 |

At exactly zero margin, $L_0=1$ for the full cube. This is a different
endpoint from the positive-margin limit: arbitrarily thin separating
films still exist for every $\delta>0$.

[high_vf_comparison.csv](results_margin_tubes/high_vf_comparison.csv) compares
the asymptotic calibration and expected Euler with direct one-dimensional
quadrature at VF 90%, 95%, 99%, and 99.9%, for all seven $q$ values.
For example, at $q=32$, VF 99%, the true margin is 0.00485169835 and
the quadratic approximation is 0.00485157217. Expected Euler is 4321.2413;
its first-order approximation in matrix VF is 4321.3667. These are
small-margin approximations, with no uniform accuracy claim over all VFs.

## Reproduction and saved evidence

From the repository root, using the existing `rf_explo` environment:

```powershell
mamba run -n rf_explo python -B theory/test_margin_tubes.py
mamba run -n rf_explo python -B -u theory/check_margin_tubes.py --output theory/results_margin_tubes_rerun
```

The runner refuses a nonempty output directory. The saved scientific
calculation took several seconds on the development machine; it needs no
3D grid or large field arrays. Its angular arrays scale with quadrature
order, so increasing `--order` increases memory and time.

The 12 tests cover exact two-half-space tubes, normal orientation and Gram
matrices, cone moments, closed $q=3,4$ strata, the derivative recurrence,
projection against constrained optimization, angular normalization,
independent distance tubes, cubic remainders, high-VF expansions, the
singular endpoint, and output preservation. Passing tests does not replace
review of the mathematical argument.

The saved CSV files contain all numerical comparisons:

- [normal_cones.csv](results_margin_tubes/normal_cones.csv)
- [gmf_comparison.csv](results_margin_tubes/gmf_comparison.csv)
- [quadrature_precision.csv](results_margin_tubes/quadrature_precision.csv)
- [distance_tubes.csv](results_margin_tubes/distance_tubes.csv)
- [cubic_remainders.csv](results_margin_tubes/cubic_remainders.csv)
- [high_vf_limits.csv](results_margin_tubes/high_vf_limits.csv)
- [high_vf_comparison.csv](results_margin_tubes/high_vf_comparison.csv)

The summary records Python/NumPy/SciPy versions, quadrature orders,
acceptance tolerances, and SHA-256 hashes of both scripts, the production
`adler.py` and `fields.py` sources, and all seven CSV files. Production
sources retain their previous bytes, preserving the earlier experiment 04
provenance. The current derivation status is recorded in the theory notes.
