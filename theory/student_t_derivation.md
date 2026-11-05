# Expected LKCs of Student-t random fields

## 1. GKF convention

Let the Gaussian components be independent, unit variance, stationary,
and isotropic, with correlation $\exp(-r^2/l_c^2)$. For a regular hitting
set $D$, use

$$
\mathbb E[L_i(G^{-1}(D))]
=\sum_{j=0}^{3-i}
\left[\begin{matrix}i+j\\j\end{matrix}\right]
\rho^jL_{i+j}(M)\mathcal M_j^\gamma(D),
\qquad \rho=\frac{1}{l_c\sqrt\pi}.
$$

The bracket is the flag coefficient
$\binom nk\omega_n/(\omega_k\omega_{n-k})$, with
$\omega_k=\pi^{k/2}/\Gamma(1+k/2)$.
It is not an ordinary binomial coefficient.

## 2. Cube expansion

For side $a$, the cube LKCs are $(1,3a,3a^2,a^3)$, and

$$
\begin{aligned}
\mathbb E[L_0]&=\rho^3a^3\mathcal M_3+3\rho^2a^2\mathcal M_2+3\rho a\mathcal M_1+\mathcal M_0,\\
\mathbb E[L_1]&=2\rho^2a^3\mathcal M_2+\tfrac{3\pi}{2}\rho a^2\mathcal M_1+3a\mathcal M_0,\\
\mathbb E[L_2]&=2\rho a^3\mathcal M_1+3a^2\mathcal M_0,\\
\mathbb E[L_3]&=a^3\mathcal M_0.
\end{aligned}
$$

Earlier versions of this note had incorrect GMF indices and coefficients
in the $L_1,L_2$ lines. The equations above agree with the shared implementation.

## 3. Gaussian reference

For a standard half-line threshold $u$,

$$
(\mathcal M_0,\mathcal M_1,\mathcal M_2,\mathcal M_3)
=(\bar\Phi(u),\varphi(u),u\varphi(u),(u^2-1)\varphi(u)).
$$

For component SD $\sigma$, standardize $u$ to $u/\sigma$.
This reference is recovered in the Student-t large-$\nu$ limit.

## 4. Ratio construction and scale

With independent Gaussian fields $Z_0,\ldots,Z_\nu$ of common SD $\sigma$,

$$
T_\nu(x)=\frac{Z_0(x)}
{\sqrt{\sum_{i=1}^\nu Z_i(x)^2/\nu}}.
$$

The common amplitude cancels exactly. The marginal is standard $t_\nu$,
with variance $\nu/(\nu-2)$ only for $\nu>2$; there is no factor $\sigma^2$.
The spatial metric in the GKF belongs to the **standardized Gaussian
components**, not to the nonlinear t field.

In particular, the t field's gradient second moment diverges for
$\nu\leq4$; this is not a reason to insert an infinite spectral factor into
the component-space GKF. The formula API supports finite $\nu\geq3$;
the direct generator requires integer $\nu\geq3$. Lower-degree cases
require separate treatment of singularities and are outside this API.
No blanket claim that all Cauchy-field geometric formulas are impossible
is intended. Use the Gaussian API for $\nu=\infty$.

A single random denominator constant over the whole domain is a different
spatial model and does not justify these same expected-LKC formulas.

## 5. Cone GMFs

The hitting set is the cone
$D_u=\{z_0\geq(u/\sqrt\nu)\|(z_1,\ldots,z_\nu)\|\}$ in Gaussian component
space. Its Gaussian tube coefficients are not obtained by differentiating
the scalar t survival function.

Writing

$$
c_\nu(u)=\frac{1}{\sqrt{2\pi}}
(1+u^2/\nu)^{-(\nu-1)/2},
\qquad
A_\nu=\sqrt{\frac2\nu}
\frac{\Gamma((\nu+1)/2)}{\Gamma(\nu/2)},
$$

the implemented coefficients are

$$
\begin{aligned}
\mathcal M_0&=\bar F_\nu(u),&
\mathcal M_1&=c_\nu(u),\\
\mathcal M_2&=A_\nu u c_\nu(u),&
\mathcal M_3&=\left(\frac{\nu-1}{\nu}u^2-1\right)c_\nu(u).
\end{aligned}
$$

These are the standard t-field EC-density expressions in the present GMF
normalization. This note states the cone formulas and consistency checks;
it is not a complete derivation of every tube integral.
The implementation evaluates gamma ratios logarithmically.

Checks with geometric meaning:

- At $u=0$, the excursion is the numerator's Gaussian half-space, so all
  four expected LKCs agree with the Gaussian formulas.
- As $\nu\to\infty$, $A_\nu\to1$ and $c_\nu\to\varphi$.
- Scaling every Gaussian component by the same positive number leaves
  the ratio and all predictions unchanged.

## 6. Finite-domain expectations

Insert the coefficients above into Section 2 with
$\rho=1/(l_c\sqrt\pi)$, independently of common component amplitude.
In particular $\mathbb E[L_3]=a^3\bar F_\nu(u)$.
This gives an expectation, not a guarantee that one sample or one
voxel estimator agrees with it.

## 7. EC-zero proxy

The bulk positive EC-zero is $u_{\rm EC}=\sqrt{\nu/(\nu-1)}$.
Finite-cube boundary terms shift it. For $a=100,l_c=8$:

| $\nu$ | Bulk threshold | Bulk VF | Finite-cube threshold | Finite-cube VF |
|---|---:|---:|---:|---:|
| Gaussian limit | 1.000 | 15.87% | 0.7787 | 21.81% |
| 10 | 1.054 | 15.83% | 0.8159 | 21.68% |
| 5 | 1.118 | 15.72% | 0.8588 | 21.48% |
| 3 | 1.225 | 15.40% | 0.9280 | 21.09% |

Thus this family barely changes this particular EC-zero proxy at fixed
component correlation and domain. It does not prove that marginal
distributions never matter, or that these VFs bound nonspanning samples.
A spatial percolation threshold requires a separate measurement.

## References

- Adler & Taylor (2007), *Random Fields and Geometry*.
- Worsley (1994), *Local maxima and the expected Euler characteristic of excursion sets of chi-squared, F and t fields*.
- [Taylor (2006), A Gaussian kinematic formula](https://arxiv.org/abs/math/0602545).
- [Taylor & Adler, flag-coefficient GKF](https://arxiv.org/abs/math/0612580).
- [SPM's t-field EC densities](https://github.com/spm/spm/blob/main/spm_ECdensity.m), an independent implementation of the normalization.
