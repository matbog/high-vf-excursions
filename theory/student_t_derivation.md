# Adler–Taylor LKC Expectations for Student-$t$ Random Fields

## 1. The Adler–Taylor Kinematic Formula

For a smooth, stationary, isotropic Gaussian-related random field $f$ on an $N$-dimensional
compact manifold $M$ (with boundary), the expected Lipschitz–Killing curvatures (LKCs) of
the excursion set $A_u = \{x \in M : f(x) \geq u\}$ are given by

$$
\mathbb{E}\bigl[L_i(A_u)\bigr]
= \sum_{j=0}^{N-i} \binom{i+j}{j}
  \left(\frac{\lambda_2}{2\pi}\right)^{j/2}
  L_{i+j}(M)\; \mathcal{M}_j(D_u),
\tag{1}
$$

where:

- $L_k(M)$ are the LKCs of the domain $M$ (purely geometric constants),
- $\lambda_2 = -C''(0)$ is the second spectral moment of the covariance $C(r)$,
- $D_u = [u, +\infty)$ is the *hitting set* in distribution space,
- $\mathcal{M}_j(D_u)$ are the Minkowski functionals of $D_u$ — they carry all
  information about the marginal distribution.

**What changes between field families:** only $\mathcal{M}_j(D_u)$.
The spatial side $(\lambda_2,\, L_{i+j}(M))$ is identical as long as the correlation
structure is unchanged.

---

## 2. LKCs of a 3D Cube

For a cube of side $a$, the non-zero domain LKCs are

$$
L_3(M) = a^3,\quad
L_2(M) = 3a^2,\quad
L_1(M) = 3a,\quad
L_0(M) = 1.
$$

Substituting into (1) for $N = 3$ and expanding all binomial coefficients:

$$
\mathbb{E}[L_0(A_u)]
= \left(\frac{\lambda_2}{2\pi}\right)^{3/2} a^3\,\mathcal{M}_3
+ 3\left(\frac{\lambda_2}{2\pi}\right) a^2\,\mathcal{M}_2
+ 3\left(\frac{\lambda_2}{2\pi}\right)^{1/2} a\,\mathcal{M}_1
+ \mathcal{M}_0,
\tag{2a}
$$

$$
\mathbb{E}[L_1(A_u)]
= 2\left(\frac{\lambda_2}{2\pi}\right) a^3\,\mathcal{M}_3
+ 3\left(\frac{\lambda_2}{2\pi}\right)^{1/2} a^2\,\mathcal{M}_2
+ 3\,a\,\mathcal{M}_1,
\tag{2b}
$$

$$
\mathbb{E}[L_2(A_u)]
= \frac{3}{2}\left(\frac{\lambda_2}{2\pi}\right)^{1/2} a^3\,\mathcal{M}_3
+ 3\,a^2\,\mathcal{M}_2
+ 3\,a\,\mathcal{M}_1,
\tag{2c — check binomials}
$$

$$
\mathbb{E}[L_3(A_u)] = a^3\,\mathcal{M}_0.
\tag{2d}
$$

> **Note.** $\mathbb{E}[L_3]$ is simply the expected volume, confirming
> $\mathcal{M}_0(D_u) = \mathbb{P}(f \geq u)$.

For a **Gaussian covariance** $C(r) = \sigma^2 \exp(-r^2/l_c^2)$:

$$
\lambda_2 = \frac{2\sigma^2}{l_c^2}.
$$

---

## 3. Minkowski Functionals for the Gaussian Distribution

The hitting set $D_u = [u, +\infty)$ for the standard Gaussian $Z \sim \mathcal{N}(0,\sigma^2)$
has Minkowski functionals (derived from the Steiner tube formula applied to the half-line):

$$
\mathcal{M}_0^{(\gamma)}(u) = \bar\Phi\!\left(\frac{u}{\sigma}\right)
= \frac{1}{2}\operatorname{erfc}\!\left(\frac{u}{\sigma\sqrt{2}}\right),
\tag{3a}
$$

$$
\mathcal{M}_1^{(\gamma)}(u) = \frac{1}{\sigma\sqrt{2\pi}}\exp\!\left(-\frac{u^2}{2\sigma^2}\right),
\tag{3b}
$$

$$
\mathcal{M}_2^{(\gamma)}(u) = \frac{u}{\sigma^2}\,\mathcal{M}_1^{(\gamma)}(u),
\tag{3c}
$$

$$
\mathcal{M}_3^{(\gamma)}(u) = \frac{u^2 - \sigma^2}{\sigma^4}\,\mathcal{M}_1^{(\gamma)}(u).
\tag{3d}
$$

The pattern is $\mathcal{M}_j^{(\gamma)} = \mathcal{M}_1^{(\gamma)}\cdot H_{j-1}(u/\sigma)/\sigma^{j-1}$,
where $H_k$ are (probabilist's) Hermite polynomials.

---

## 4. Student-$t$ Random Fields

### 4.1 Construction

Let $Z_0, Z_1, \ldots, Z_\nu$ be independent stationary Gaussian RFs, all with the same
covariance $C(r) = \sigma^2 \exp(-r^2/l_c^2)$. Define the **Worsley $t_\nu$ field**:

$$
T_\nu(x) = \frac{Z_0(x)}{\sqrt{\chi^2_\nu(x)/\nu}},
\qquad \chi^2_\nu(x) = \sum_{i=1}^{\nu} Z_i(x)^2.
\tag{4}
$$

This is a valid stationary RF with:

- **Marginal distribution** $T_\nu(x) \sim t_\nu$ at every point (numerator and denominator
  are independent because they involve different components $Z_i$),
- **Second spectral moment** $\lambda_2 = 2\sigma^2/l_c^2$ (of the numerator $Z_0$; proved
  in Worsley 1994 that the Adler formula uses this quantity),
- **Spatially varying denominator** — this is the key: a scalar $\chi^2$ construction
  ($T = Z/\sqrt{s/\nu}$ with $s$ a fixed scalar per realisation) would merely rescale
  a Gaussian field and the $t_\nu$ Adler formula would not apply.

The marginal PDF is

$$
p_\nu(u) = \frac{\Gamma\!\left(\frac{\nu+1}{2}\right)}{\sqrt{\nu\pi}\;\Gamma\!\left(\frac{\nu}{2}\right)}
\left(1 + \frac{u^2}{\nu}\right)^{-(\nu+1)/2}.
\tag{5}
$$

For $\nu \geq 2$ the variance is $\sigma^2_\nu = \nu\sigma^2/(\nu-2)$; for $\nu \geq 3$ the
second moment $\lambda_2$ is finite and the Adler formula applies.

### 4.2 Tail behavior

For large $|u|$: $p_\nu(u) \sim C_\nu |u|^{-(\nu+1)}$ — **polynomial**, much heavier than
the Gaussian $\exp(-u^2/2\sigma^2)$.  This shifts the volume fraction still present at
the percolation threshold:

| Field | Tail | EC-zero phi_c+ (finite domain, lc=8, a=100) |
|-------|------|----------------------------------------------|
| Gaussian | $\exp(-u^2/2)$ | $\approx 21.8\%$ |
| $\chi^2(k=1)$ | $u^{-1/2}\exp(-u/2)$ | $\approx 37\text{–}51\%$ (lc-dependent) |
| $t_\nu$, $\nu=10$ | $u^{-11}$ | $\approx 21.7\%$ |
| $t_\nu$, $\nu=3$ | $u^{-4}$ | $\approx 21.1\%$ |

Note: $t_\nu$ ceilings are nearly identical to Gaussian — see Section 7.

---

## 5. Minkowski Functionals for the $t_\nu$ Distribution

**Important:** The M_j are NOT obtained by differentiating the survival function. The
excursion set $\{T_\nu \geq u\}$ is the preimage of a cone in $\mathbb{R}^{\nu+1}$ Gaussian
space, and the GKF requires the Gaussian Minkowski functionals of that cone — not the
scalar $t_\nu$ marginal derivatives. (The PDF-derivative approach gives expressions that
fail the decisive consistency check at $u=0$; see Section 7.)

Define the cone-geometry base function:

$$
c_\nu(u) = \frac{1}{\sqrt{2\pi}}\left(1 + \frac{u^2}{\nu}\right)^{-(\nu-1)/2}.
\tag{6e}
$$

Then:

$$
\mathcal{M}_0^{(t)}(u) = \bar F_\nu(u) = \mathbb{P}(T_\nu \geq u),
\tag{6a}
$$

$$
\mathcal{M}_1^{(t)}(u) = c_\nu(u),
\tag{6b}
$$

$$
\mathcal{M}_2^{(t)}(u) = A_\nu\, u\, c_\nu(u),
\quad A_\nu = \sqrt{\tfrac{2}{\nu}}\,\frac{\Gamma((\nu+1)/2)}{\Gamma(\nu/2)},
\tag{6c}
$$

$$
\mathcal{M}_3^{(t)}(u) = \left(\frac{\nu-1}{\nu}\,u^2 - 1\right) c_\nu(u).
\tag{6d}
$$

### Limiting cases

- $\nu \to \infty$: $c_\nu(u) \to \varphi(u)$, $A_\nu \to 1$, and all
  $\mathcal{M}_j^{(t)} \to \mathcal{M}_j^{(\gamma)}$ — recovers the Gaussian formulas (3).
- $u = 0$: $c_\nu(0) = \varphi(0)$ for **all** $\nu$ — consistent with the identity
  $\{T_\nu \geq 0\} = \{Z_0 \geq 0\}$ (same excursion set for all $\nu$, so all LKCs at
  $\kappa=0$ must agree across field families).
- $\mathcal{M}_3(u) = 0$ at $u_c = \sqrt{\nu/(\nu-1)}$ — the bulk EC-zero threshold.
- $\nu = 1$ (Cauchy): no finite variance, $\lambda_2$ is undefined, Adler formula breaks.
- **Safe range:** $\nu \geq 3$.

---

## 6. Full Expected LKCs for a 3D Cube with a $t_\nu$ Field

Substituting (6) into (2) with $\lambda_2 = 2\sigma^2/l_c^2$
(and writing $\rho = \sqrt{\lambda_2 / 2\pi} = \sigma/(l_c\sqrt{\pi})$):

$$
\mathbb{E}[L_0(A_u)]
= \rho^3 a^3\,\mathcal{M}_3^{(t)}
+ 3\rho^2 a^2\,\mathcal{M}_2^{(t)}
+ 3\rho\, a\,\mathcal{M}_1^{(t)}
+ \mathcal{M}_0^{(t)},
\tag{7a}
$$

$$
\mathbb{E}[L_3(A_u)] = a^3\,\bar F_\nu(u).
\tag{7b}
$$

The volume fraction $\Phi = \bar F_\nu(u)$ is controlled by the threshold $u$ through the
$t_\nu$ survival function, which decays polynomially — hence for the same $u$, a $t_\nu$
field produces significantly more volume than a Gaussian field.

---

## 7. Percolation Criterion

Following Roubin & Colliat (2016), the **percolation thresholds in volume fraction** are
the roots of $\mathbb{E}[L_0(A_u)] = 0$:

$$
\bar F_\nu(\kappa_c^\pm) = \Phi_c^\pm.
\tag{8}
$$

The upper root $\Phi_c^+$ (high threshold, meatball regime — disconnected inclusions)
and lower root $\Phi_c^-$ (sponge regime — connected complement) bound the disconnected
phase. For a Gaussian, $\Phi_c^+ \approx 15.9\%$ at $\beta \to \infty$.

For a $t_\nu$ field, $\mathcal{M}_3(u) = 0$ at $u_c = \sqrt{\nu/(\nu-1)}$, which is larger
than the Gaussian $u_c = 1$ and increases with heavier tails. However, $\Phi_c^+ =
\bar F_\nu(u_c) = \mathbb{P}(T_\nu \geq \sqrt{\nu/(\nu-1)})$ is slightly **smaller** than
the Gaussian value, because the heavier tail raises the threshold faster than the CDF
responds. The two effects nearly cancel:

| $\nu$ | $u_c$ | $\Phi_c^+$ (bulk) | $\Phi_c^+$ (finite domain) |
|--------|--------|-------------------|---------------------------|
| $\infty$ (Gaussian) | 1.000 | 15.87% | 21.81% |
| 10 | 1.054 | 15.83% | 21.68% |
| 5  | 1.118 | 15.72% | 21.48% |
| 3  | 1.225 | 15.40% | 21.09% |

**Conclusion:** $t_\nu$ fields do not meaningfully raise the disconnected VF ceiling.
Marginal tail heaviness is almost irrelevant; the ceiling is controlled by the spatial
correlation structure ($\rho$, $a/l_c$), not the marginal distribution.

$$
\Phi_c^+(t_\nu) \lesssim \Phi_c^+(\text{Gaussian}).
$$

---

## 8. Connection to the Paper Extension

The natural extension of Roubin–Colliat (2016) to non-Gaussian fields follows these steps:

1. **Derive** $\mathcal{M}_j^{(t)}$ (done above; refs: Worsley 1994, Taylor 2006).
2. **Compute** $\mathbb{E}[L_0]$ analytically for $t_\nu$ via (7a).
3. **Find** $\kappa_c^\pm(\nu, \beta)$ numerically (roots of the EC curve).
4. **Map** $\Phi_c^+(\nu, \beta)$ — the achievable disconnected VF as function of tail heaviness.
5. **Compare** against Monte Carlo measurements (test_4.py validates steps 1–4).
6. **Apply** to cement paste / battery electrode morphologies where $\Phi > 16\%$ is required.

---

## References

- Adler, R.J. & Taylor, J.E. (2007). *Random Fields and Geometry*. §11.7 (non-Gaussian EC densities).
- Worsley, K.J. (1994). Local maxima and the expected Euler characteristic of excursion sets of $\chi^2$, $F$ and $t$ fields. *Advances in Applied Probability*, 26(1), 13–42.
- Taylor, J.E. (2006). A Gaussian kinematic formula. *Annals of Probability*, 34(1), 122–158.
- Roubin, E. & Colliat, J.-B. (2016). Critical probability of percolation over bounded region in $N$-dimensional Euclidean space. *Journal of Statistical Mechanics*, P033306.
- Roubin, E. (2013). *Meso-scale FE and morphological modelling of heterogeneous media* (PhD thesis). §2 (LKC framework, $\chi^2$ derivation).
