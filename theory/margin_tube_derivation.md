# Margin chambers: tube derivation and the high-VF limit

Status, 12 November 2026: this is a detailed derivation draft, with independent
value-space quadrature checks. It derives the coefficients already used in
experiment 04; no coefficient correction is indicated. The argument should
receive independent mathematical review before publication. The spatial GKF
is invoked under its usual hypotheses, rather than proved here.

The main results are an exact finite-radius tube decomposition, the GMFs
through order three, a recurrence for their stratum integrals, and a
high-volume-fraction expansion. The accompanying checks generate no random
fields and do not establish convergence of a digital LKC estimator.
See [the numerical report](margin_tube_checks.md) for results and commands.

## 1. Definition and normalization

Fix a finite integer $q\geq2$ and a value margin $\delta>0$. Let

$$
D_\delta=\bigcup_{i=1}^q D_{\delta,i},\qquad
D_{\delta,i}=\{z\in\mathbb R^q:z_i-z_j\geq\delta\text{ for all }j\ne i\}.
$$

Write $\gamma_q$ for standard Gaussian measure, $\varphi,\Phi$ for the
one-dimensional standard normal density and CDF, and
$T_r(D)=\{z:\operatorname{dist}(z,D)\leq r\}$. Our convention is

$$
\gamma_q(T_r(D))=\mathcal M_0(D)+\mathcal M_1(D)r
+\frac{\mathcal M_2(D)}{2!}r^2+\frac{\mathcal M_3(D)}{3!}r^3+O(r^4).
$$

The factorials matter. This is the Gaussian tube convention underlying the
GKF, not the convention of raw polynomial coefficients.
[Taylor (2006)](https://arxiv.org/abs/math/0602545).

All margins here are standardized Gaussian value differences. If component
standard deviation is $\sigma$, replace the supplied margin by
$\delta/\sigma$. Neither $\delta$ nor a tube radius below is a spatial film
thickness.

## 2. Removing the common Gaussian direction

Let $e=\mathbf1/\sqrt q$ and $H=e^\perp$. Decompose $z=x+te$ orthogonally.
Coordinate differences depend only on $x$, so

$$
D_\delta=\widetilde D_\delta+\mathbb Re,\qquad
\operatorname{dist}(z,D_\delta)=\operatorname{dist}(x,\widetilde D_\delta).
$$

The standard Gaussian splits into independent components on $H$ and
$\mathbb Re$. Integrating the latter gives

$$
\gamma_q(T_r(D_\delta))=\gamma_{q-1}(T_r(\widetilde D_\delta)).
$$

Thus $q=3$ is a planar problem and $q=4$ a three-dimensional one. A raw vector
of differences such as $(Z_1-Z_q,\ldots,Z_{q-1}-Z_q)$ is correlated; treating
it as a standard independent Gaussian would give the wrong geometry.

In $H$, chamber $i$ is the translated zero-margin ordering cone:

$$
\widetilde D_{\delta,i}=v_i+\widetilde D_{0,i},\qquad
v_i=\delta(e_i-\mathbf1/q),\qquad
\|v_i\|=\delta\sqrt{(q-1)/q}.
$$

## 3. Separation, reach, and the admissible tube radius

For $x\in D_{\delta,i}$ and $y\in D_{\delta,j}$, $i\ne j$,

$$
(x-y)\cdot(e_i-e_j)\geq2\delta,
\qquad \|x-y\|\geq\sqrt2\delta.
$$

Equality is attained by taking $x_i=\delta/2$, $x_j=-\delta/2$, reversing
these two coordinates for $y$, and giving all other coordinates the same
value at most $-\delta/2$. A common translation puts both points in $H$.
Consequently,

$$
\operatorname{dist}(D_{\delta,i},D_{\delta,j})=\sqrt2\delta.
$$

Set $a=\delta/\sqrt2$. The chamber tubes are disjoint for $0\leq r<a$.
Each closed convex chamber has a unique nearest-point projection. Two
different chambers cannot both be nearest at a distance less than $a$,
because their projected points would be less than $2a$ apart. The midpoint
of a closest pair has two nearest projections at distance $a$. This proves
that the reach of the union is exactly $a$.

In particular, chamber tube measures add for **every fixed $\delta>0$ and
$0\leq r<\delta/\sqrt2$**. There is no positive radius that works uniformly
as $\delta\downarrow0$.

## 4. Faces and their Gaussian weights

A relative-open face of codimension $k$ is specified by a winner $i$ and a
set $J$ of $k$ active competitors. Its coordinates satisfy

$$
x_i=t+\delta,\qquad x_j=t\ (j\in J),\qquad x_\ell<t\ (\ell\notin J\cup\{i\}).
$$

There are $q\binom{q-1}{k}$ such faces. The tangent vector for $t$ has
$k+1$ entries equal to one, giving surface Jacobian $\sqrt{k+1}$; the
remaining free coordinates are orthogonal to it. Integrating their Gaussian
densities up to $t$ gives the total face weight

$$
S_k(q,\delta)=q\binom{q-1}{k}\sqrt{k+1}
\int_{\mathbb R}\varphi(t+\delta)\varphi(t)^k
\Phi(t)^{q-1-k}\,dt,\qquad 1\leq k\leq q-1.
$$

Set $S_k=0$ for $k\geq q$. We will also use $S_0=\mathcal M_0$:

$$
\mathcal M_0=q\int_{\mathbb R}\varphi(t+\delta)\Phi(t)^{q-1}\,dt
=q\int_{\mathbb R}\varphi(t)\Phi(t-\delta)^{q-1}\,dt.
$$

The outward unit normals on the face are

$$
n_j=\frac{e_j-e_i}{\sqrt2},\quad j\in J,
\qquad x\cdot n_j=-a.
$$

The sign is fixed by the inequality defining the chamber: increasing
$x_j-x_i$ moves outside. With these normals as columns of $N$, their Gram
matrix is

$$
Q_k=N^TN=\frac{I_k+\mathbf1\mathbf1^T}{2},\qquad
\det Q_k=\frac{k+1}{2^k}.
$$

The face's outward normal cone is $\{Ns:s\in\mathbb R_+^k\}$.

## 5. Exact normal-cone tube decomposition

Projection onto a convex polyhedron partitions its exterior, up to
Lebesgue-null boundaries, into points $x+v$ where $x$ belongs to a
relative-open face and $v$ belongs to that face's normal cone. Tangent and
normal spaces are orthogonal. Since the face is flat, this decomposition
has unit Jacobian in orthonormal tangent and normal coordinates.

For $v=Ns$, the Gaussian density ratio is

$$
\frac{\varphi_q(x+v)}{\varphi_q(x)}
=\exp(-x\cdot v-\|v\|^2/2)
=\exp\left(a\sum_j s_j-\tfrac12s^TQ_ks\right).
$$

It is independent of $x$ along the face. Passing from $v$ to $s$ contributes
Jacobian $\sqrt{\det Q_k}$. Define

$$
C_k(r,a)=\sqrt{\det Q_k}
\int_{\substack{s_j\geq0\\s^TQ_ks\leq r^2}}
\exp\left(a\sum_j s_j-\tfrac12s^TQ_ks\right)\,ds.
$$

Summing over the faces and then the separated chambers yields the exact
identity

$$
\boxed{\gamma_q(T_r(D_\delta))=
S_0(q,\delta)+\sum_{k=1}^{q-1}S_k(q,\delta)C_k(r,a),
\qquad 0\leq r<a.}
$$

This also deals with the unbounded faces: each $S_k$ is finite, and the
normal factor is bounded on the compact normal integration region for
fixed $a,r$. Scaling $s=ru$ gives a fixed bounded integration domain and a
smooth exponential integrand. Taylor expansion under the integral is
therefore justified. In particular $C_k=O(r^k)$; faces of codimension four
or more do not contribute through order three. The constants in these
estimates are for fixed finite $q,\delta$.

## 6. Deriving the first three coefficients

For a facet there is one unit normal, so

$$
C_1(r,a)=\int_0^r e^{au-u^2/2}\,du
=r+\frac a2r^2+\frac{a^2-1}{6}r^3+O(r^4).
$$

For two active normals, $n_1\cdot n_2=1/2$. Their cone is a sector of angle
$\pi/3$, with bisector $b=(n_1+n_2)/\sqrt3$. The projection of $x$ onto
their span is $-(2a/\sqrt3)b$. Polar coordinates give

$$
\begin{aligned}
C_2(r,a)
&=\int_{-\pi/6}^{\pi/6}\int_0^r
u\exp\left(\frac{2a}{\sqrt3}u\cos\theta-u^2/2\right)\,du\,d\theta\\
&=\frac\pi6r^2+\frac{2a}{3\sqrt3}r^3+O(r^4).
\end{aligned}
$$

For three active normals, their spherical triangle has side lengths
$\pi/3$. The spherical cosine rule gives interior angle
$\arccos(1/3)$, hence solid angle

$$
\Omega=3\arccos(1/3)-\pi.
$$

The leading term is the volume of the unit-cone sector:

$$
C_3(r,a)=\frac\Omega3r^3+O(r^4).
$$

Collecting coefficients and multiplying by the appropriate factorial gives

$$
\boxed{\begin{aligned}
\mathcal M_0&=S_0,\\
\mathcal M_1&=S_1,\\
\mathcal M_2&=aS_1+\frac\pi3S_2,\\
\mathcal M_3&=(a^2-1)S_1+\frac{4a}{\sqrt3}S_2+2\Omega S_3.
\end{aligned}}
$$

These are precisely the positive-margin formulas in `src/adler.py`.

## 7. Low-dimensional analytic controls

For $q=2$, the contrast coordinate is standard normal and the hitting set
is $(-\infty,-a]\cup[a,\infty)$. Thus, for $r<a$,

$$
\gamma(T_r(D_\delta))=2\Phi(r-a),
$$

giving $(\mathcal M_0,\ldots,\mathcal M_3)
=(2\Phi(-a),2\varphi(a),2a\varphi(a),2(a^2-1)\varphi(a))$.

The highest-codimension weight has the closed form

$$
S_{q-1}(q,\delta)=q(2\pi)^{-(q-1)/2}
\exp\left(-\frac{q-1}{2q}\delta^2\right).
$$

This follows by completing the square in the face integral. In particular,

$$
S_2(3,\delta)=\frac3{2\pi}e^{-\delta^2/3},\qquad
S_3(4,\delta)=\frac4{(2\pi)^{3/2}}e^{-3\delta^2/8}.
$$

For $q=3$, the remaining weight also simplifies:

$$
S_1(3,\delta)=3\sqrt{\frac2\pi}\,e^{-\delta^2/4}
\Phi(-\delta/\sqrt6).
$$

Here the product $\varphi(t+\delta)\varphi(t)$ is proportional to a normal
density with mean $-\delta/2$ and variance $1/2$; integrating $\Phi(t)$
against it gives the last CDF. For $q=3,4$, all boundary strata appear in
the exact tube identity above, allowing a finite-radius check rather than
unstable differentiation of a Monte Carlo probability curve.

## 8. A derivative recurrence for all strata

Let $I_k$ denote the integral in $S_k$, without its multiplicity and
Jacobian. Differentiation in $\delta$ and integration by parts in $t$ give

$$
I_k'=-\int(t+\delta)\varphi(t+\delta)\varphi(t)^k\Phi(t)^{q-1-k}\,dt,
$$

$$
0=-(k+1)\int t\varphi(t+\delta)\varphi(t)^k\Phi(t)^{q-1-k}\,dt
-\delta I_k+(q-1-k)I_{k+1}.
$$

The Gaussian decay kills the boundary terms. Combining these identities
and the binomial factors yields

$$
\boxed{S_k'=-\frac{k}{k+1}\delta S_k
-\sqrt{\frac{k+1}{k+2}}S_{k+1},\qquad 0\leq k\leq q-1.}
$$

The last term vanishes at $k=q-1$. The $k=0$ case gives
$\mathcal M_0'=-S_1/\sqrt2$, so VF decreases strictly with the margin.
This recurrence is also a route to higher derivatives without finite
differences. Its use here extends the face-weight definition beyond the
production API's codimension-three limit.

## 9. High-VF calibration and the singular endpoint

Let $s_k=S_k(q,0^+)$ and $c_q=s_1/\sqrt2$. These are face-weight limits,
not the higher GMFs of $D_0$. For fixed finite $q$ the recurrence gives

$$
\mathcal M_0(q,\delta)=1-c_q\delta
+\frac{s_2}{2\sqrt3}\delta^2+O(\delta^3).
$$

There is a useful probabilistic interpretation:

$$
\boxed{c_q=\mathbb E\max_{1\leq i\leq q} Z_i,
\qquad Z_i\stackrel{\mathrm{iid}}\sim N(0,1).}
$$

Indeed the maximum has density $q\varphi(t)\Phi(t)^{q-1}$, and
integration by parts using $t\varphi(t)=-\varphi'(t)$ gives
$\mathbb E\max Z_i=q(q-1)\int\varphi(t)^2\Phi(t)^{q-2}\,dt=s_1/\sqrt2$.

If $\epsilon=1-\mathrm{VF}$, inversion gives a practical calibration:

$$
\boxed{\delta=\frac{\sqrt2}{s_1}\epsilon
+\sqrt{\frac23}\frac{s_2}{s_1^3}\epsilon^2+O(\epsilon^3).}
$$

The full positive-branch expansion is
$\mathcal M_j=m_j^++b_j\delta+O(\delta^2)$, where

$$
(m_0^+,m_1^+,m_2^+,m_3^+)
=\left(1,s_1,\frac\pi3s_2,-s_1+2\Omega s_3\right),
$$

$$
(b_0,b_1,b_2,b_3)=\left(
-\frac{s_1}{\sqrt2},
-\sqrt{\frac23}s_2,
\frac{s_1}{\sqrt2}-\frac{\pi}{2\sqrt3}s_3,
\sqrt6s_2-\frac{4\Omega}{\sqrt5}s_4
\right).
$$

The $s_4$ term is needed for the slope of $\mathcal M_3$ when $q\geq5$.
It comes from differentiating $S_3$; it does not introduce an order-four
GMF into the three-dimensional GKF.

At **exactly** $\delta=0$, the hitting set is all of $\mathbb R^q$, so
its GMFs are $(1,0,0,0)$. The positive-margin limits above generally differ.
There is no contradiction: as positive margins decrease, their union omits
the tie skeleton, whereas $D_0$ includes it. Gaussian volume ignores that
null set; topology and boundary measures need not. The admissible additive
tube radius also shrinks to zero, so taking the zero-margin limit and taking
tube derivatives at zero are different operations.

These expansions are for fixed $q$. No uniform error estimate as
$q\to\infty$ is claimed.

## 10. Consequences for the spatial Adler link

For the spatial application, assume a compact stratified domain and iid
centered, unit-variance Gaussian components satisfying the smoothness,
nondegeneracy, and regularity requirements of the GKF. The hitting set here
has a finite polyhedral stratification, is locally convex for $\delta>0$,
has positive reach, and has finite Gaussian face integrals. The all-LKC GKF
then connects its GMFs to spatial expected LKCs. The general theorem and
its analytic hypotheses are external inputs to this derivation; a
publication should check them explicitly for its chosen field model.
[Taylor & Adler (2009), Sections 1.4 and 4](https://arxiv.org/abs/math/0612580).

For a cube of side $L$ and covariance $\exp(-r^2/l_c^2)$, the Gaussian
component gradient variance is $2/l_c^2$. Put $\rho=1/(l_c\sqrt\pi)$ and
$h=L\rho$. With the [repository LKC normalization](../docs/framework.md),

$$
\begin{aligned}
\mathbb E L_0&=\mathcal M_0+3h\mathcal M_1+3h^2\mathcal M_2+h^3\mathcal M_3,\\
\mathbb E L_1&=3L\mathcal M_0+\tfrac{3\pi}{2}L^2\rho\mathcal M_1
+2L^3\rho^2\mathcal M_2,\\
\mathbb E L_2&=3L^2\mathcal M_0+2L^3\rho\mathcal M_1,\\
\mathbb E L_3&=L^3\mathcal M_0.
\end{aligned}
$$

In particular, define

$$
K_0=1+3hs_1+\pi h^2s_2+h^3(-s_1+2\Omega s_3),
$$

$$
K_1=-\frac{s_1}{\sqrt2}-\sqrt6hs_2
+3h^2\left(\frac{s_1}{\sqrt2}-\frac{\pi s_3}{2\sqrt3}\right)
+h^3\left(\sqrt6s_2-\frac{4\Omega s_4}{\sqrt5}\right).
$$

Then

$$
\boxed{\mathbb E L_0=K_0+K_1\delta+O(\delta^2)
=K_0+\frac{K_1}{c_q}\epsilon+O(\epsilon^2),
\qquad \mathbb E L_3=L^3(1-\epsilon).}
$$

For $q=32,L=64,l_c=5$, quadrature gives $c_q=2.0696688279$,
$K_0=4328.5543713$, and $K_1=-1487.6001004$. Thus the leading calibration
is $\delta\simeq0.483169088(1-\mathrm{VF})$. At VF 99%, the numerical
margin is 0.00485169835 and the quadratic approximation is 0.00485157217.
At exactly zero margin the excursion is the whole cube and $L_0=1$;
$K_0$ is the limit along positive margins only.

This makes the high-VF expectation curve explicit. It does not turn
Euler characteristic into a particle count or establish nonspanning:
$L_0=\beta_0-\beta_1+\beta_2$ in three dimensions. Nor does a single scalar
independently prescribe a grading curve and all shape descriptors. It
selects a one-parameter family once $q$, covariance, and domain are fixed.

## 11. What remains to establish

The coefficient derivation and endpoint analysis are now written out and
have independent numerical controls. The remaining theoretical review is
to audit this proof and its precise use of the GKF's hypotheses, including
the unbounded hitting set. Separately, continuum/digital convergence,
spanning probabilities, and particle grading need spatial evidence.

The subsequent [physical film-width derivation](margin_film_widths.md)
addresses regular two-winner interfaces using the local relation
$w\simeq2\delta/\|\nabla(G_i-G_j)\|$. Interface-area weighting gives an
inverse-chi-four law for the linearized width and recovers the leading
Adler matrix volume. Its finite-margin validity near junctions and weak
gradients remains distinct from the value-space tube argument.
