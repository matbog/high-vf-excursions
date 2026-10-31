# Theory notes

Short umbrella for the derivations that live alongside this code. Each
experiment's own `README.md` inlines the formulas it depends on; the files
here are the longer worked-out versions that would otherwise clutter those.

## What this project asks

The Adler–Taylor Gaussian Kinematic Formula (GKF) gives a closed-form
expectation for every Lipschitz–Killing curvature (LKC) of the excursion
set of a Gaussian-related field. We treat it as a *morphology generator*:
pick a hitting set `D`, compute its Gaussian Minkowski functionals
`M_j^γ(D)`, and the LKCs of the random morphology `f^{-1}(D)` follow.

Within that framework we look at four hitting-set families:

| Family | Hitting set | Analytic status of `M_j` |
|--------|-------------|--------------------------|
| **Gaussian half-line** | `[u, ∞)` | fully closed-form |
| **Student-$t_\nu$ cone** | `{z/√(χ²/ν) ≥ u}` in $\mathbb{R}^{\nu+1}$ | closed-form via cone geometry |
| **Margin chambers** | `{z_{(q)} - z_{(q-1)} ≥ δ}` in $\mathbb{R}^q$ | closed-form derived here (proof informal) |
| **Coalescing lobes** | `∪_a B(c_a, r)` in $\mathbb{R}^q$ | closed-form before contact only |

The first two are standard (Worsley 1994, Taylor 2006, Adler–Taylor 2007).
The last two are the small novel pieces of this exploration.

## Files

- **[student_t_derivation.md](student_t_derivation.md)** — Student-$t_\nu$ `M_j`
  from the cone geometry of the ratio field (Section 5 is the one that
  matters; the PDF-derivative shortcut was wrong, see the note there).
  Referenced from `experiments/03_student_t/`.

- **[gaussian_minkowski_functionals.md](gaussian_minkowski_functionals.md)**
  — closed forms for the margin-tessellation chambers (`M_0..M_3`, all
  four with cone-angle constants) and the pre-contact coalescing-lobe
  ball-union additivity. Referenced from `experiments/04_margin_tessellation/`
  and `experiments/05_coalescing_lobes/`.

## What is actually proved vs. numerically supported

- **Gaussian** `M_j`: textbook (Adler–Taylor §4).
- **Student-$t_\nu$** `M_j`: Worsley 1994 for `M_3`; the full set is a
  straightforward cone-tube calculation and matches Monte Carlo for
  `ν ∈ {3, 5, 10, 30}`.
- **Margin tessellation** `M_j`: derived from first principles here; `M_0`
  and `M_1` are elementary; `M_2` and `M_3` use the chamber's codimension-
  `k` stratum measures and known cone-angle constants (π/3 for a wedge,
  3·arccos(1/3) − π for the equilateral spherical triangle). Validated to
  under 2% against periodic Euler-density measurements at `N = 256` across
  `q ∈ {2, 4, 8, 16}`. A formal publication-grade proof of the stratified
  tube decomposition is not written up.
- **Coalescing lobes**, `r < r_c`: ball-additive, `M_j^γ(B(c,r))` derived
  from the non-central $\chi^2$ CDF and its derivatives. Matches Monte
  Carlo to under 5% at `N = 100`.
- **Coalescing lobes**, `r ≥ r_c`: ball-intersection terms via
  inclusion-exclusion are not derived. No analytic Euler prediction in
  this regime.

## Open mathematical items

1. Write up the stratified-tube proof for the margin chamber `D_δ` so the
   `M_2, M_3` formulas can be cited rather than conjectured.
2. Compute `M_j^γ` for two Gaussian-weighted intersecting balls in
   $\mathbb{R}^3$ (post-contact lobes).
3. Separate, for each family, the three thresholds that are often
   conflated: hitting-set contact, `E[L_0] = 0`, and spatial percolation
   `P_span = 1/2` with periodic wrapping.

## References

- Adler, R.J. & Taylor, J.E. (2007). *Random Fields and Geometry*. Springer.
- Worsley, K.J. (1994). Local maxima and the expected Euler characteristic
  of excursion sets of $\chi^2$, $F$ and $t$ fields.
  *Advances in Applied Probability* 26(1), 13–42.
- Taylor, J.E. (2006). A Gaussian kinematic formula.
  *Annals of Probability* 34(1), 122–158.
- Roubin, E. & Colliat, J.-B. (2016). Critical probability of percolation
  over bounded region in $N$-dimensional Euclidean space.
  *J. Stat. Mech.* P033306.
