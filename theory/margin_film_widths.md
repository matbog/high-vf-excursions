# Physical film widths in Gaussian margin tessellations

This derivation connects the scalar value margin to a physical normal
width. It supplements the [margin tube derivation](margin_tube_derivation.md)
without changing the hitting set or its Adler coefficients. It is a
small-margin result for regular two-winner interfaces, not an exact
finite-margin width law or a guarantee of globally separated particles.
The [spatial pilot](../experiments/04_margin_tessellation/FILM_WIDTHS.md)
tests both the interface sampling law and the finite-margin approximation.

## Model and the meaning of a typical interface

Let $G_1,\ldots,G_q$ be iid, centered, unit-variance, smooth Gaussian fields
in three dimensions, with covariance

$$
C(r)=\exp(-r^2/l_c^2).
$$

The regular zero-margin interface between winners $i$ and $j$ is

$$
\Sigma_{ij}=\{x:G_i(x)=G_j(x)>G_k(x)\text{ for all }k\notin\{i,j\}}.
$$

Count each unordered pair once. Triple ties and critical points are
excluded from this regular surface; under the usual smoothness and
nondegeneracy assumptions they have no contribution to its area measure.
For a bounded observation window $B$, the area-Palm average of a mark $h$
means

$$
\mathbb E_A h=
\frac{\mathbb E\sum_{i<j}\int_{\Sigma_{ij}\cap B}h(x)\,dA}
{\mathbb E\sum_{i<j}\operatorname{Area}(\Sigma_{ij}\cap B)}.
$$

This samples interfaces in proportion to area across the ensemble. It is
not the average obtained by giving every realization equal weight after
normalizing each one's interface area separately. Nor is it a uniform
spatial-point average or a conditional sample of occupied grid edges.

The marked Kac-Rice/coarea formula supplies the gradient-norm weight in
such surface integrals. We use it under the smoothness, density, and
integrability assumptions appropriate to the Gaussian fields here.
[Armentano, Azaïs & León, *On a general Kac-Rice formula for the measure of
a level set*](https://arxiv.org/abs/2304.07424).

## Why the gradient has a chi-four law on the interface

Put $H=G_i-G_j$. Stationarity and the covariance above give

$$
\operatorname{Var}(\partial_aG_i)=2/l_c^2,\qquad
\operatorname{Var}(\partial_aH)=4/l_c^2.
$$

The gradient coordinates are independent. At a fixed point, $\nabla H$
is independent of all the component values, since their cross-covariances
vanish and the joint vector is Gaussian. Thus imposing a tie and the
winner inequalities affects the values, but not this conditional gradient
law. At an ordinary point,

$$
R=\frac{l_c}{2}\|\nabla H\|\sim\chi_3.
$$

An interface-area integral additionally weights this density by
$\|\nabla H\|$, through coarea. Therefore

$$
f_{R,A}(r)=\frac{r f_{\chi_3}(r)}{\mathbb E\chi_3}
=\frac12r^3e^{-r^2/2},\qquad r>0.
$$

This is the **$\chi_4$ density**. The conditioning on the other $q-2$
competitors cancels out of this normalized gradient law. For iid components
sharing this covariance, $q$ changes the amount of interface and the
margin required for a target VF, but not this area-Palm gradient law.

The law applies to the ideal Gaussian ensemble. A finite random-mode
generator approximates that model and can have component-specific,
slightly anisotropic gradient covariances conditional on its modes.

## Local normal width and its limiting distribution

At a regular interface point let $g=\|\nabla H\|>0$ and
$n=\nabla H/g$. Along its fixed normal line,

$$
H(x+sn)=gs+O(s^2).
$$

While these same two fields remain the top competitors, the matrix lies
between $H=-\delta$ and $H=\delta$. Its local linearized full width is

$$
W_{\mathrm{lin}}=\frac{2\delta}{g}
=\frac{\delta l_c}{R},\qquad R\sim\chi_4\text{ under area sampling}.
$$

For a fixed regular point with a positive gap to the third competitor,
the implicit-function theorem gives this leading term as
$\delta\downarrow0$. With sufficient smoothness, the quadratic shifts
of the two endpoints cancel, giving a pointwise $O(\delta^3)$ absolute
width correction. This remainder is not uniform near small gradients
or multiple-winner junctions. We do not claim a uniform rate for ensemble
width moments.

Let $a=\delta l_c$ and $Y=W_{\mathrm{lin}}/a$. For $w>0$,

$$
\boxed{\mathbb P_A(W_{\mathrm{lin}}\leq w)
=\exp\left(-\frac{a^2}{2w^2}\right)
\left(1+\frac{a^2}{2w^2}\right).}
$$

The density and moments are

$$
f_W(w)=\frac{a^4}{2w^5}e^{-a^2/(2w^2)},\qquad
\mathbb E_A W_{\mathrm{lin}}^p
=a^p2^{-p/2}\Gamma((4-p)/2),\quad 0\leq p<4.
$$

In particular,

$$
\boxed{\mathbb E_AW_{\mathrm{lin}}=\sqrt{\pi/8}\,\delta l_c
\simeq0.626657\,\delta l_c,\qquad
\operatorname{Var}_A W_{\mathrm{lin}}
=(1/2-\pi/8)(\delta l_c)^2.}
$$

The quantile at probability $p\in(0,1)$ is
$a/F_{\chi_4}^{-1}(1-p)$. Fourth and higher moments of this *linearized*
law diverge. This is not a statement that actual finite-margin films have
those divergent moments: the linear approximation fails in precisely the
small-gradient regime that creates the long tail.

At a point sampled uniformly from thin matrix volume, the leading
interface sampling measure receives an extra width factor proportional
to $1/R$. This changes $\chi_4$ back to $\chi_3$. Width statistics sampled
by interface area and by matrix volume are consequently different objects.

## Interface area and the Adler consistency identity

Define

$$
I_q=\int_{\mathbb R}\varphi(t)^2\Phi(t)^{q-2}\,dt,
\qquad c_q=q(q-1)I_q=\mathbb E\max_{1\leq i\leq q}Z_i.
$$

The value integral for one winning tie is $I_q$, and
$\mathbb E\|\nabla H\|=(2/l_c)\mathbb E\chi_3
=4\sqrt{2/\pi}/l_c$. Summing over unordered pairs gives the expected
zero-margin interface area per spatial volume, counted once:

$$
\boxed{\mathcal S_q=\frac{2\sqrt{2/\pi}}{l_c}\,c_q.}
$$

Consequently,

$$
\boxed{\mathcal S_q\,\mathbb E_AW_{\mathrm{lin}}=c_q\delta.}
$$

This exactly matches the leading matrix VF derived from the GMFs:
$1-\mathcal M_0(q,\delta)=c_q\delta+O(\delta^2)$.
The geometric interpretation counts the two sides of a film once, via
its full width. Treating both aggregate boundaries as independent films
would introduce a factor-of-two error. Junction overlap and finite-margin
geometry enter beyond this leading volume term.

It also agrees with the bulk GKF surface limit. The aggregate boundary
area per volume is $4\rho\mathcal M_1$ with $\rho=1/(l_c\sqrt\pi)$;
as $\delta\downarrow0$, this tends to $2\mathcal S_q$, since each tie has
two aggregate-facing sides. External observation-window faces are not
part of $\mathcal S_q$.

## Using VF as the single control

For prescribed aggregate fraction $\phi$, keep using the existing
calibration $\mathcal M_0(q,\delta)=\phi$. Once $q$ and $l_c$ are fixed,
this also fixes the predicted width law. With $\epsilon=1-\phi$,

$$
\delta=\epsilon/c_q+O(\epsilon^2),\qquad
\mathbb E_A W_{\mathrm{lin}}
=\sqrt{\pi/8}\,\frac{l_c}{c_q}\epsilon+O(\epsilon^2).
$$

Thus raising $q$ at fixed high VF and $l_c$ makes the leading films
thinner through the margin calibration. Increasing $l_c$ scales widths
up, while also changing the grain morphology and expected LKCs. The
extra prediction introduces no new morphology parameter.

For any physical threshold $h>0$, the CDF gives the expected *area
fraction of interfaces* with linearized normal width below $h$. It is
not the matrix volume fraction in those films, a probability that a
voxel method makes an error, or a guarantee of global Euler accuracy.
Replacing $h$ by the grid spacing gives a resolution diagnostic, not a
universal voxel-count rule. Grid orientation and nearby junctions matter.

There is no positive deterministic minimum width in the Gaussian model:
the gradient magnitude has unbounded support. Quantiles are the meaningful
statistical design quantities here.

## Finite-margin qualifications to test

The full top-two-gap field can depart from the two-field approximation
in several ways. Curvature changes the normal profile. A third field can
become runner-up before either adjacent grain changes label, or can itself
become the first occupied winner on one side. At a weak-gradient point,
the first occupied endpoints may lie beyond the local Taylor region.

The spatial pilot therefore measures the connected matrix interval
containing the original tie along its normal. It records endpoint winners,
endpoint runners, censored profiles, and a probe-refinement diagnostic.
Matching the two endpoint labels is insufficient to certify a pure
two-field film. The measurements are normal-line widths, not minimum
distances between curved surfaces, and are distinct from the directional
grid-edge lengths in the earlier film audit.

The completed [bounded tail audit](../experiments/04_margin_tessellation/FILM_TAILS.md)
compares linear widths, nonlinear pair widths, and full-component widths
along the same normals. Additional competitors contribute more than pair
curvature in 152 of 157 resolved baseline cases with width errors above 10%.
Its extended searches find every full-field endpoint, while 31 of 1,536
profile/margin cases retain numerical stability uncertainty. These finite-
margin limitations do not invalidate the GMF expectation formulas.

Independent mathematical review of the marked-area argument, and broader
spatial tests at small margins and near junctions, remain appropriate
before treating this as a validated finite-VF material model.
