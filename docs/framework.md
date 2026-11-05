# Random-field morphology and the Gaussian Kinematic Formula

Let $G:M\to\mathbb R^q$ have independent, identically distributed Gaussian
components and let $A=G^{-1}(D)$. The hitting set $D$ determines which field
values belong to the inclusion phase. Scalar Gaussian thresholds, chi-squared
exteriors, Student-t cones, margin chambers, and unions of balls are different
choices of $D$.

Our unit-variance components have correlation $C_G(r)=\exp(-r^2/l_c^2)$.
A nonlinear field derived from them need not have this correlation: for
example, the correlation of a sum of squared independent components is
$C_G(r)^2$. The GKF uses the **component-induced metric**.

## Descriptors and conventions

| LKC | Meaning |
|---|---|
| $L_3$ | Inclusion volume |
| $L_2$ | Half the total boundary area, including boundary on the domain faces |
| $L_1$ | Integrated mean curvature divided by $\pi$ |
| $L_0$ | Euler characteristic $\beta_0-\beta_1+\beta_2$ |

Positive Euler characteristic does not imply a nonspanning phase.
$L_0$ equals the number of components only when tunnels and cavities
are absent. Four LKCs do not uniquely determine a morphology.

## GKF and cube expansion

For smooth, independent unit-variance components with an isotropic common
metric, and a suitably regular domain and hitting set,

$$
\mathbb E[L_i(A)]=\sum_{j=0}^{3-i}
\left[\begin{matrix}i+j\\j\end{matrix}\right]
\rho^j L_{i+j}(M)\mathcal M_j^\gamma(D),
\quad
\left[\begin{matrix}n\\k\end{matrix}\right]
=\binom nk\frac{\omega_n}{\omega_k\omega_{n-k}},
\quad
\omega_k=\frac{\pi^{k/2}}{\Gamma(1+k/2)}.
$$

With unit components, $\lambda_2=2/l_c^2$ and
$\rho=\sqrt{\lambda_2/(2\pi)}=1/(l_c\sqrt\pi)$.
For a cube of side $a$, $(L_0,L_1,L_2,L_3)=(1,3a,3a^2,a^3)$:

$$
\begin{aligned}
\mathbb E[L_0]&=\rho^3a^3\mathcal M_3+3\rho^2a^2\mathcal M_2+3\rho a\mathcal M_1+\mathcal M_0,\\
\mathbb E[L_1]&=2\rho^2a^3\mathcal M_2+\tfrac{3\pi}{2}\rho a^2\mathcal M_1+3a\mathcal M_0,\\
\mathbb E[L_2]&=2\rho a^3\mathcal M_1+3a^2\mathcal M_0,\\
\mathbb E[L_3]&=a^3\mathcal M_0.
\end{aligned}
$$

The flag coefficient, for example $[2\;1]=\pi/2$, is essential.
See [Taylor & Adler](https://arxiv.org/abs/math/0612580).

For a standard Gaussian half-line $[u,\infty)$,

$$
(\mathcal M_0,\mathcal M_1,\mathcal M_2,\mathcal M_3)
=(\bar\Phi(u),\varphi(u),u\varphi(u),(u^2-1)\varphi(u)).
$$

For common component SD $\sigma$, either standardize the hitting set and
use the unit-component metric, or consistently pair nonstandard GMFs with
the unstandardized metric. The APIs use these equivalent conventions:

- Gaussian: threshold $u/\sigma$ (implemented with scaled density derivatives).
- Chi-squared: threshold $u/\sigma^2$.
- Margin: margin $\delta/\sigma$.
- Balls: centres $c/\sigma$ and radius $r/\sigma$.
- Student-t: common $\sigma$ cancels exactly in the ratio.

Student-t GMFs require the cone in Gaussian component space, not derivatives
of the scalar t survival function. See [the formula note](../theory/student_t_derivation.md).
Margin formulas apply to separated chambers for $\delta>0$; at exactly
$\delta=0$ the union is all value space and its GMFs are $(1,0,0,0)$.

## EC-zero, contact, and spanning

The positive-threshold root of $\mathbb E[L_0]=0$ is an **EC-zero proxy**.
For the Gaussian bulk term it is $u=1$, VF 15.87%. Cube boundary terms
shift it: $a=100,l_c=8$ gives $u\approx0.779$, VF 21.81%.
These numbers depend on physical $a/l_c$, not directly on grid size $N$.

The historical function name `find_percolation_threshold` is retained,
but it returns an EC-zero, not measured spanning probability.
Value-space contact is a separate event (for example, touching balls).
Finite-domain spanning must be measured with an explicit geometry and
connectivity convention. None of these thresholds is generally guaranteed
to coincide with the others, including in a large-domain limit.

## Minimal use

~~~python
import numpy as np
from src.fields import GridSpec, generate_gaussian_rf
from src.measures import measure_all
from src.adler import expected_lkcs_gaussian

grid = GridSpec(n=100, length=100.0)
field = generate_gaussian_rf(grid, lc=8.0, rng=np.random.default_rng(42))
measured = measure_all(field >= 1.0, grid)
theory = expected_lkcs_gaussian(grid, lc=8.0, kappas=np.array([1.0]))
print(measured["vf"], theory["phi"][0])
print(measured["euler"], theory["L0"][0])
~~~

`measure_all` uses voxel volume, 6-neighbour Euler/spanning, and binary
marching-cubes surface; curvature is enabled with `compute_l1=True`.
It is not the measurement policy of
every experiment. See [measurement conventions](measurement_conventions.md).

## What can be calibrated

Given a fixed model family, a threshold controls its expected VF and the
associated expected LKC curve. Allowing an additional scale parameter may
permit fitting both expected VF and Euler density, but feasibility and
uniqueness require checking. A single scalar does not independently fit
several arbitrary descriptors, grain sizes, and connectivity constraints.

The GKF gives expectations, including finite-domain boundary terms.
It does not itself give variances, convergence rates, or a guarantee that
one realisation is near its expectation. Under suitable mixing assumptions,
spatial averages can concentrate with increasing domain size; the variance
of an extensive quantity need not decrease. Grid refinement addresses
discretisation error, not finite-domain sampling uncertainty.
