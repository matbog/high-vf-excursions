# Candidate: Gaussian margin tessellations

## Complementary model: coalescing Gaussian-space lobes

For a vector Gaussian field `G:M -> R^q` and centres `c_1,...,c_K`, define

    R(x) = min_a ||G(x)-c_a||,
    A_r = {R(x) <= r} = G^{-1}(D_r),
    D_r = union_a B(c_a,r).

The radius `r` is a single continuous morphology control. Balls `a,b` first
touch at `r=||c_a-c_b||/2`, so this hitting set has finite merger levels.
For four tetrahedral centres in `R^3`, all with norm `rho`, their common
pairwise distance is `rho*sqrt(8/3)` and

    r_c = rho*sqrt(2/3).

The ball union is stratified, so the GKF still expresses `E[L_i(A_r)]` through
the Gaussian Minkowski functionals of `D_r`. Before contact, additivity gives

    M_j^gamma(D_r) = sum_a M_j^gamma(B(c_a,r)).

In particular, because `||Z-c_a||^2` is noncentral chi-square,

    M_0^gamma(D_r) = K F_{chi'^2_q(rho^2)}(r^2),  r < r_c.

After contact, Gaussian-weighted ball intersections enter through
inclusion-exclusion. Deriving `M_0,...,M_3` through these overlap regimes is the
main analytic task. Hitting-set contact, expected-EC zero, and spatial spanning
must be measured separately rather than assumed to coincide.

Let (G=(G_1,\ldots,G_q):M\to\mathbb{R}^q) have independent, identically
distributed, smooth, standardized Gaussian components. Define

\[
A_\delta=G^{-1}(D_\delta),\qquad
D_\delta=\{z:z_{(q)}-z_{(q-1)}\geq\delta\},
\]

where (z_{(q)}\) and (z_{(q-1)}\) are the largest and second-largest
coordinates. Equivalently, (D_\delta) is the disjoint union

\[
D_{\delta,i}=\{z:z_i-z_j\geq\delta\quad\forall j\ne i\}.
\]

Each piece is an unbounded polyhedral chamber in Gaussian value space. The
construction is therefore a Gaussian-related field excursion and is admissible
for the Gaussian kinematic formula:

\[
\mathbb E\{\mathcal L_i(A_\delta)\}
=\sum_{j=0}^{N-i}
\begin{bmatrix}i+j\\j\end{bmatrix}
(2\pi)^{-j/2}\mathcal L_{i+j}(M,g)
\mathcal M_j^\gamma(D_\delta).
\]

The zeroth Gaussian Minkowski functional, hence the expected volume fraction,
is already a one-dimensional integral:

\[
\mathcal M_0^\gamma(D_\delta)
=q\int_{-\infty}^{\infty}\varphi(t)\Phi(t-\delta)^{q-1}\,dt.
\]

The first functional is also a one-dimensional integral. Each ordered pair of
winner and runner-up coordinates contributes one facet, giving

\[
\mathcal M_1^\gamma(D_\delta)=
\sqrt{2}\,q(q-1)\int_{-\infty}^{\infty}
\varphi(t)\varphi(t-\delta)\Phi(t-\delta)^{q-2}\,dt
=-\sqrt{2}\,\frac{d}{d\delta}\mathcal M_0^\gamma(D_\delta).
\]

### Closed formulas through order three

The higher functionals require both facet and lower-dimensional edge-stratum
terms. For the present chambers these terms remain tractable because all facet
normals have the same geometry.

Set (a=\delta/\sqrt{2}), and define the total Gaussian measure of the
codimension-(k) strata by

\[
S_k(q,\delta)=q{q-1\choose k}\sqrt{k+1}
\int_{-\infty}^{\infty}
\varphi(t+\delta)\varphi(t)^k
\Phi(t)^{q-1-k}\,dt,
\qquad 1\le k\le\min(3,q-1),
\]

with (S_k=0) when (k>q-1). The factors have direct meanings:

- (q{q-1\choose k}) chooses the winner and the (k) tied runners-up;
- (\sqrt{k+1}) is the induced surface Jacobian;
- (t+\delta) is the winner value and the (k) active runners have value (t);
- the remaining coordinates must be at most (t).

The four functionals are

\[
\boxed{\mathcal M_0^\gamma(D_\delta)=
q\int\varphi(t)\Phi(t-\delta)^{q-1}\,dt,}
\]

\[
\boxed{\mathcal M_1^\gamma(D_\delta)=S_1,}
\]

\[
\boxed{\mathcal M_2^\gamma(D_\delta)=
aS_1+\frac{\pi}{3}S_2,}
\]

and

\[
\boxed{\mathcal M_3^\gamma(D_\delta)=
(a^2-1)S_1+
\frac{4a}{\sqrt3}S_2+
2\left(3\arccos\frac13-\pi\right)S_3.}
\]

The constants arise from the normal cones:

- a facet has one outward normal;
- two active facets form a wedge of angle (\pi/3);
- three active facets form an equilateral spherical triangle with solid angle
  (3\arccos(1/3)-\pi).

For (\delta>0), the (q) chambers are separated. Hence, for a sufficiently
small tube radius, their tubes are disjoint and their coefficients add. The
formulas reduce exactly to the two-half-space result when (q=2):

\[
(\mathcal M_0,\mathcal M_1,\mathcal M_2,\mathcal M_3)
=\left(2\Psi(a),2\varphi(a),2a\varphi(a),
2(a^2-1)\varphi(a)\right).
\]

These formulas have passed algebraic (q=2) tests and numerical periodic-field
checks at several volume fractions. A formal publication proof should still
state the chamber stratification, reach/tube-radius restriction and normal-cone
integration carefully.

### Explicit Adler expectation in three spatial dimensions

Substitution into the GKF gives, for a three-dimensional parameter domain,

\[
\mathbb E\{\chi(A_\delta)\}
=\sum_{j=0}^{3}
\left(\frac{\lambda_2}{2\pi}\right)^{j/2}
\mathcal L_j(M)\mathcal M_j^\gamma(D_\delta).
\]

For a Euclidean cube of side (T), with

\[
(\mathcal L_0,\mathcal L_1,\mathcal L_2,\mathcal L_3)
=(1,3T,3T^2,T^3),
\]

this becomes

\[
\boxed{
\begin{aligned}
\mathbb E\{\chi(A_\delta)\}={}&
\mathcal M_0^\gamma
+3T\left(\frac{\lambda_2}{2\pi}\right)^{1/2}\mathcal M_1^\gamma\\
&+3T^2\left(\frac{\lambda_2}{2\pi}\right)\mathcal M_2^\gamma
+T^3\left(\frac{\lambda_2}{2\pi}\right)^{3/2}\mathcal M_3^\gamma .
\end{aligned}}
\]

Thus the four formulas above can be inserted into the same Adler expectation
structure used by Roubin; only the Gaussian hitting-set functionals change.

## Numerical recovery of all four functionals

On a flat periodic domain (a (d)-torus), all domain LKCs except its volume
vanish. For a stationary isotropic Gaussian vector field, the Euler form of
the GKF becomes

\[
\frac{\mathbb E\{\chi(A_\delta)\}}{|M|}
=\left(\frac{\lambda_2}{2\pi}\right)^{d/2}
\mathcal M_d^\gamma(D_\delta).
\]

Consequently, periodic experiments in dimensions 1, 2 and 3 independently
recover (\mathcal M_1^\gamma,\mathcal M_2^\gamma) and
(\mathcal M_3^\gamma). This avoids unstable numerical differentiation of a
Monte Carlo tube-probability curve.

At target volume fraction (\phi=0.70), using a squared-exponential covariance
with (L_c=8), the first numerical campaign produced:

| q | M0 | M1 | M2 | M3 |
|---:|---:|---:|---:|---:|
| 2 | 0.700 | 0.739 | 0.281 | -0.644 |
| 4 | 0.700 | 1.161 | 1.145 | -0.363 |
| 8 | 0.700 | 1.496 | 2.273 | 1.186 |
| 16 | 0.700 | 1.799 | 3.512 | 4.005 |

For (q=2), all four exact values are available because the hitting set is two
separated half-spaces. They are approximately

\[
(0.700,\ 0.7408,\ 0.2854,\ -0.6308),
\]

which provides a successful control of the periodic numerical estimates. The
remaining differences, especially for (\mathcal M_3), quantify finite-grid
and discretization bias in addition to sampling uncertainty.

The sign of (\mathcal M_3) is especially informative in three dimensions:
on a torus it is proportional to expected Euler density. At (\phi=0.70), the
change from negative for (q=4) to positive for (q\ge8) is consistent with a
transition from handle/percolation-dominated winner domains to
component-dominated winner cores.

The closed formulas above supply the required facet, pairwise-edge and
triple-stratum contributions. Remaining work is formal proof presentation and
systematic convergence validation, rather than discovery of the coefficients.

The central hypothesis is that the rejected near-tie set forms a thin matrix
film along interfaces between the winning labels. Increasing (q) should
produce more, smaller winner cells. A positive margin erodes their interfaces,
potentially yielding a high total fraction of mutually disconnected cell
cores.

## Why this may be new enough to investigate

This is not a monotone remapping, a chi-square field, or an independent union
of ordinary excursions. It is a vector-Gaussian excursion whose hitting-set
geometry directly encodes a random tessellation and an adjustable interphase
thickness. It links:

- analytic expected LKCs through the GKF;
- high-volume matrix/inclusion morphologies;
- percolation and finite-size/RVE analysis;
- a controllable number of latent competing phases (q);
- a geometrically meaningful margin parameter (\delta).

## Required theory before publication

1. Derive or numerically evaluate (\mathcal M_j^\gamma(D_\delta)) for
   (j=1,2,3). Polyhedral edges make the stratified hitting-set version of the
   GKF relevant.
2. Establish whether additivity over the (q) separated chambers applies directly
   for every (\delta>0), and handle the singular limit (\delta\downarrow0).
3. Compare the expected-EC zero with direct wrapping probability.
4. Determine an RVE in both LKC and percolation senses as a function of
   (T/L_c), (q), and (\delta).
5. Compare against Gaussian, (\chi_1^2), and independent-union baselines at
   equal volume fraction and equal local spectral moment.
