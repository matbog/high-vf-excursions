# Analytical limits of the Adler link across model families

## What "Adler link" means in practice

- **Closed-form prediction.** Given field parameters $(l_c, u)$ and model
  parameters (e.g. $\nu$, $q$, $\delta$), the expected Euler
  characteristic and volume fraction are computable without simulation,
  via one-dimensional integrals or explicit formulae.
- **Direct calibration.** The inverse problem — find parameters that
  produce a target $(\varphi, \chi/V)$ — is a 2D root-finding problem
  over analytic functions, not a simulation loop.
- **Convergence guarantee.** The formula is exact in the ergodic limit
  (large domain relative to $l_c$). On a finite domain it holds in
  expectation, with variance that decreases predictably.

## Full link: Gaussian and Student-$t$ fields

Both families have exact, closed-form $\mathcal M_j$ for all four
functionals.

**Gaussian** (`src/adler.py: expected_lkcs_gaussian`):

$$
\mathcal M_0 = \tfrac{1}{2}\,\mathrm{erfc}(u/\sqrt{2}),\quad
\mathcal M_1 = \varphi(u),\quad
\mathcal M_2 = u\,\varphi(u),\quad
\mathcal M_3 = (u^2 - 1)\,\varphi(u).
$$

These are the exact Gaussian Minkowski functionals of the half-line
$D_u = [u, +\infty)$. Monte Carlo experiments (experiment 01, 10
realisations) confirm that measured $L_0$ tracks $\mathbb{E}[L_0]$ to
within sampling noise across the full threshold range.

**Student-$t$** (`src/adler.py: expected_lkcs_t`, following Worsley
1994, with the cone-geometry correction derived in
`theory/student_t_derivation.md`):

$$
\begin{aligned}
\mathcal M_0 &= \mathbb{P}(T_\nu \geq u), \\[2pt]
\mathcal M_1 &= c_\nu(u) = (2\pi)^{-1/2}\bigl(1 + u^2/\nu\bigr)^{-(\nu-1)/2}, \\[2pt]
\mathcal M_2 &= A_\nu\, u\, c_\nu(u),\quad
    A_\nu = \sqrt{2/\nu}\,\Gamma((\nu+1)/2)/\Gamma(\nu/2), \\[2pt]
\mathcal M_3 &= \bigl(\tfrac{\nu-1}{\nu}\, u^2 - 1\bigr)\, c_\nu(u).
\end{aligned}
$$

These apply to the Worsley construction
$T(x) = Z_0(x) / \sqrt{\chi^2_\nu(x)/\nu}$ where all components share the
same $l_c$. MC experiments (experiment 03) confirm $\mathbb{E}[L_0]$ for
$\nu \in \{3, 5, 10, 30\}$ and Gaussian ($\nu \to \infty$).

**Inverse calibration for both families** is fully operational: given
$(\varphi_{\text{target}}, L_0^{\text{target}})$, solve for $(l_c, u)$ or
$(l_c, \nu, u)$ via `scipy.optimize`.

All four $L_i$ for $t_\nu$ are derived from the cone-geometry
$\mathcal M_j^{(t)}$ (see `theory/student_t_derivation.md` §5) and
validated against MC in `experiments/03_student_t/`.

## Partial link: margin tessellation

The hitting set
$D_\delta = \{ z \in \mathbb{R}^q : \max_i z_i - \mathrm{secondmax}_i\,z_i \geq \delta \}$
is a union of $q$ separated polyhedral chambers. The GKF applies to this
stratified set, but deriving $\mathcal M_j$ requires integrating over
the normal cones of each stratum.

**$\mathcal M_0$ and $\mathcal M_1$** are exact one-dimensional
integrals:

$$
\mathcal M_0 = q \int \varphi(t)\, \Phi(t-\delta)^{q-1}\, dt, \qquad
\mathcal M_1 = \sqrt{2}\, q(q-1) \int \varphi(t)\, \varphi(t-\delta)\, \Phi(t-\delta)^{q-2}\, dt.
$$

These are reliable. Volume fraction predictions match measurements to
better than 0.5% for all tested $(q, \delta)$ combinations.

**$\mathcal M_2$ and $\mathcal M_3$** are derived from codimension-$k$
stratum measures $S_k$:

$$
\mathcal M_2 = a S_1 + \tfrac{\pi}{3} S_2, \qquad
\mathcal M_3 = (a^2 - 1) S_1 + \tfrac{4a}{\sqrt{3}} S_2 + 2(3\arccos(1/3) - \pi) S_3,
$$

where $a = \delta/\sqrt{2}$ and $S_k$ involves a one-dimensional
integral with combinatorial multiplicity and the induced surface
Jacobian $\sqrt{k+1}$.

**Numerical validation status:**
- $q = 2$: exact identities available (hitting set is two separated
  half-spaces). All four $\mathcal M_j$ match the exact values.
- $q \in \{4, 8, 16\}$ at $\varphi \in \{0.50, 0.70, 0.80\}$:
  $\mathcal M_3$ recovered from periodic Euler density experiments in
  dimensions 1, 2, 3. Agreement is consistent with the formulae being
  correct.

**What is missing:** a formal published proof of the stratified-tube
decomposition for the polyhedral hitting set. The geometry is
well-defined and the $q = 2$ control is exact, but citing these formulae
in a paper requires the proof. Until then, the Euler predictions should
be treated as a strong numerical conjecture, not a theorem.

## Partial link: coalescing Gaussian-space lobes

The hitting set is $D_r = \bigcup_a B(c_a, r)$, a union of $K$ balls in
$\mathbb{R}^d$ centred at fixed targets $c_a$. The GKF applies, but the
$\mathcal M_j$ depend on the geometry of the ball union.

**Before contact ($r < r_c$):** The balls are disjoint, so by additivity

$$
\mathcal M_j^\gamma(D_r) = \sum_a \mathcal M_j^\gamma\bigl(B(c_a, r)\bigr).
$$

For a single ball $B(c, r)$ in $\mathbb{R}^d$,
$\mathcal M_0^\gamma(B(c, r)) = \mathbb{P}(\lVert Z - c\rVert \leq r) = F_{\chi'^2}((d, \lVert c\rVert^2), r^2)$.
Higher functionals follow by differentiating the tube probability in the
ball radius; the four $\mathcal M_j$ for a single ball are implemented
in `src/adler.py: ball_lkcs`. Agreement with Monte Carlo is better than
5% on all four LKCs at $N = 100$ for $r/r_c \leq 0.95$ (see
`experiments/05_coalescing_lobes/`).

**After contact ($r \geq r_c$):** Ball intersections
$B(c_a, r) \cap B(c_b, r)$ contribute through inclusion-exclusion. The
required quantities are Gaussian-weighted volumes of spherical lens
regions in $\mathbb{R}^d$, i.e. integrals of the form

$$
\int_{\lVert z - c_a\rVert \leq r,\; \lVert z - c_b\rVert \leq r}
\varphi_d(z)\, dz.
$$

These are nontrivial (no simple closed form in $d \geq 3$) and have not
been derived. There is no analytic Euler prediction in the post-contact
regime.

**Summary of what is available:**

| Quantity | Before contact | After contact |
|----------|---------------|---------------|
| $\mathcal M_0 \ldots \mathcal M_3$ | Closed-form, validated < 5% | Not implemented |
| $\mathbb{E}[L_0 \ldots L_3]$ | Closed-form, validated < 5% | Not available |

## Summary table

| Model | $\mathcal M_0$ (vol frac) | $\mathcal M_1$ | $\mathcal M_2, \mathcal M_3$ | $\mathbb{E}[L_0]$ reliable? |
|-------|---------------|----|---------|-----------------|
| Gaussian | Exact | Exact | Exact | Yes — MC validated |
| Student-$t$ ($\nu \geq 3$) | Exact | Exact | Exact | Yes — MC validated |
| Margin tessellation | Exact | Exact | Derived, proof pending | Strong numerical conjecture |
| Coalescing lobes ($r < r_c$) | Exact | Exact | Exact | Yes — MC validated |
| Coalescing lobes ($r \geq r_c$) | Not implemented | Not implemented | Not implemented | No |

## What would complete the link for each model

**Margin tessellation:** Publish the stratified-tube proof for the
polyhedral hitting set $D_\delta$. The geometry is fully worked out;
what remains is a careful presentation of the tube-radius restriction
and normal-cone integration steps.

**Coalescing lobes, after contact:** Derive or numerically evaluate the
Gaussian-weighted spherical lens volumes in $\mathbb{R}^3$. This is the
primary open analytic problem for this model. The pre-contact case is
already closed (`src/adler.py: ball_lkcs`).
