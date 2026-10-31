# Analytical limits of the Adler link across model families

## What "Adler link" means in practice

- **Closed-form prediction.** Given field parameters `(lc, u)` and model parameters (e.g. `nu`,
  `q`, `delta`), the expected Euler characteristic and volume fraction are computable without
  simulation, via one-dimensional integrals or explicit formulae.
- **Direct calibration.** The inverse problem — find parameters that produce a target
  `(phi, chi/V)` — is a 2D root-finding problem over analytic functions, not a simulation loop.
- **Convergence guarantee.** The formula is exact in the ergodic limit (large domain relative to
  `lc`). On a finite domain it holds in expectation, with variance that decreases predictably.

## Full link: Gaussian and Student-t fields

Both families have exact, closed-form `M_j` for all four functionals.

**Gaussian** (`src/adler.py: expected_lkcs_gaussian`):

```
M0 = (1/2) erfc(u/sqrt(2))
M1 = phi(u)
M2 = u * phi(u)
M3 = (u^2 - 1) * phi(u)
```

These are the exact Gaussian Minkowski functionals of the half-line `D_u = [u, +inf)`.
Monte Carlo experiments (experiment 01, 50 realisations) confirm that measured `L0` tracks
`E[L0]` to within sampling noise across the full threshold range.

**Student-t** (`src/adler.py: expected_lkcs_t`, following Worsley 1994):

```
M0 = P(T_v >= u) = t_v survival function
M1 = p_v(u)                            [t_v PDF]
M2 = p_v(u) * (v+1) u / (v + u^2)
M3 = p_v(u) * [(v+1)(v+3) u^2 / (v+u^2)^2 - (v+1) / (v+u^2)]
```

These apply to the Worsley construction `T(x) = Z_0(x) / sqrt(chi2_field(x)/v)` where all
components share the same `lc`. Monte Carlo experiments (experiment 02) confirm `E[L0]`
for `v = 3, 5, 10, 30` and Gaussian (`v -> inf`).

**Inverse calibration for both families** is fully operational: given `(phi_target, L0_target)`,
solve for `(lc, u)` or `(lc, v, u)` via `scipy.optimize`.

All four `L_i` for `t_v` are derived from the cone-geometry `M_j^{(t)}` (see
`theory/student_t_derivation.md` §5) and validated against MC in `experiments/03_student_t/`.

## Partial link: margin tessellation

The hitting set `D_delta = { z in R^q : max_i z_i - second_max_i z_i >= delta }` is a union
of `q` separated polyhedral chambers. The GKF applies to this stratified set, but deriving
`M_j` requires integrating over the normal cones of each stratum.

**`M0` and `M1`** are exact one-dimensional integrals:

```
M0 = q * integral phi(t) * Phi(t-delta)^(q-1) dt
M1 = sqrt(2) * q*(q-1) * integral phi(t) * phi(t-delta) * Phi(t-delta)^(q-2) dt
```

These are reliable. Volume fraction predictions match measurements to better than 0.5% for
all tested `(q, delta)` combinations.

**`M2` and `M3`** are derived from codimension-`k` stratum measures `S_k`:

```
M2 = a * S1 + (pi/3) * S2
M3 = (a^2 - 1) * S1 + (4a/sqrt(3)) * S2 + 2*(3*arccos(1/3) - pi) * S3
```

where `a = delta/sqrt(2)` and `S_k` involves a one-dimensional integral with combinatorial
multiplicity and the induced surface Jacobian `sqrt(k+1)`.

**Numerical validation status:**
- `q = 2`: exact identities available (hitting set is two separated half-spaces). All four
  `M_j` match the exact values.
- `q = 4, 8, 16` at `phi = 0.50, 0.70, 0.80`: `M3` recovered from periodic Euler density
  experiments in dimensions 1, 2, 3. Agreement is consistent with the formulae being correct.

**What is missing:** a formal published proof of the stratified-tube decomposition for the
polyhedral hitting set. The geometry is well-defined and the `q=2` control is exact, but
citing these formulae in a paper requires the proof. Until then, the Euler predictions should
be treated as a strong numerical conjecture, not a theorem.

## Partial link: coalescing Gaussian-space lobes

The hitting set is `D_r = union_a B(c_a, r)`, a union of `K` balls in `R^d` centred at fixed
targets `c_a`. The GKF applies, but the `M_j` depend on the geometry of the ball union.

**Before contact (`r < r_c`):** The balls are disjoint, so by additivity

```
M_j^gamma(D_r) = sum_a M_j^gamma(B(c_a, r))
```

For a single ball `B(c, r)` in `R^d`, `M_0^gamma(B(c,r)) = P(||Z - c|| <= r) =
F_{nc chi2}((d, ||c||^2), r^2)`. Higher functionals follow by differentiating the tube
probability in the ball radius; the four `M_j` for a single ball are implemented in
`src/adler.py: ball_lkcs`. Agreement with Monte Carlo is better than 5% on all four LKCs
at `N = 100` for `r/r_c <= 0.95` (see `experiments/05_coalescing_lobes/`).

**After contact (`r >= r_c`):** Ball intersections `B(c_a, r) ∩ B(c_b, r)` contribute through
inclusion-exclusion. The required quantities are Gaussian-weighted volumes of spherical lens
regions in `R^d`, i.e. integrals of the form

```
integral_{||z-c_a||<=r, ||z-c_b||<=r} phi_d(z) dz
```

These are nontrivial (no simple closed form in `d >= 3`) and have not been derived. There is
no analytic Euler prediction in the post-contact regime.

**Summary of what is available:**

| Quantity | Before contact | After contact |
|----------|---------------|---------------|
| `M0..M3` | Closed-form, validated < 5% | Not implemented |
| `E[L0..L3]` | Closed-form, validated < 5% | Not available |

## Summary table

| Model | M0 (vol frac) | M1 | M2, M3 | E[L0] reliable? |
|-------|---------------|----|---------|-----------------|
| Gaussian | Exact | Exact | Exact | Yes — MC validated |
| Student-t (v >= 3) | Exact | Exact | Exact | Yes — MC validated |
| Margin tessellation | Exact | Exact | Derived, proof pending | Strong numerical conjecture |
| Coalescing lobes (r < r_c) | Exact | Exact | Exact | Yes — MC validated |
| Coalescing lobes (r >= r_c) | Not implemented | Not implemented | Not implemented | No |

## What would complete the link for each model

**Margin tessellation:** Publish the stratified-tube proof for the polyhedral hitting set
`D_delta`. The geometry is fully worked out; what remains is a careful presentation of the
tube-radius restriction and normal-cone integration steps.

**Coalescing lobes, after contact:** Derive or numerically evaluate the Gaussian-weighted
spherical lens volumes in `R^3`. This is the primary open analytic problem for this model.
The pre-contact case is already closed (`src/adler.py: ball_lkcs`).
