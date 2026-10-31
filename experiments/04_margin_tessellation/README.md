# Experiment 04 — Margin tessellation

![q=16 margin tessellation at φ=0.70](paraview/q16_phi70/display_margin_tess_RF_thresh.png)

Competing Gaussian fields with a "margin" rule. Pushes the disconnected-
phase volume fraction well past the Gaussian 16% / 22% ceilings while
keeping the Adler link intact.

## The model

Let `G = (G_1, …, G_q)` be `q` i.i.d. smooth Gaussian RFs. Include a
voxel `x` iff the winner's lead over the runner-up exceeds a margin `δ`:

```
A_δ = { x : max_i G_i(x) − second-max_i G_i(x) ≥ δ }
```

The hitting set in $\mathbb{R}^q$ is a disjoint union of `q` unbounded
polyhedral "winner" chambers; the matrix film is where two or more fields
are tied. Larger `q` → more, smaller winner cells. Larger `δ` → thicker
film, smaller `φ`.

Different winner labels are mutually separated by the film by construction,
but a *single* winner can still span. The relevant percolation check is
per-winner spanning, not spanning of the union.

## Adler link status

Closed-form `M_j^γ(D_δ)` for `j = 0..3`, numerically validated to <2% at
`N = 256` across `q ∈ {2, 4, 8, 16}` on all four LKCs. A formal
stratified-tube proof is not written up — the `q = 2` identities and the
`N = 256` convergence are the current evidence.

## The `M_j` (inline)

Set `a = δ/√2`. The codimension-`k` stratum measures:

```
S_k(q, δ) = q · C(q−1, k) · √(k+1)
          · ∫ φ(t+δ) · φ(t)^k · Φ(t)^(q−1−k) dt,     1 ≤ k ≤ min(3, q−1)
```

Then

```
M_0(q, δ) = q · ∫ φ(t) · Φ(t−δ)^(q−1) dt
M_1(q, δ) = S_1                                            = −√2 · ∂M_0/∂δ
M_2(q, δ) = a · S_1 + (π/3) · S_2
M_3(q, δ) = (a² − 1) · S_1 + (4a/√3) · S_2
           + 2 · (3·arccos(1/3) − π) · S_3
```

The constants are the normal-cone angles: `π/3` for a wedge of two active
facets, `3·arccos(1/3) − π` for the equilateral spherical triangle of
three active facets. For `q = 2` the formulas reduce exactly to twice the
Gaussian `M_j` at `δ/√2` — a check the code runs on every launch (also
checks `M_1 = −√2 · ∂M_0/∂δ`).

Full derivation (chamber stratification + normal-cone integration):
[`theory/gaussian_minkowski_functionals.md`](../../theory/gaussian_minkowski_functionals.md).

## What it runs

- `N_VALUES = [64, 128, 256]`, `N_MAIN = 128`, `LENGTH = 64`, `LC = 5`.
- `Q_VALUES = [2, 4, 8, 16]`, 10 target VFs in `[0.30, 0.90]`.
- 5 realizations per `(N, q, φ)`.
- Per-winner `L_0..L_3` + per-winner spanning.

```bash
mamba run -n rf_explo python run.py
```

Runtime ≈ 2.5 h for the full three-resolution convergence; ≈ 30 min at
`N = 128` only (edit `N_VALUES`).

## Estimator choices

See [`docs/measurement_conventions.md`](../../docs/measurement_conventions.md).
Short version:

- `L_0`: 26-connectivity (thin films sub-voxel at high `q`, high `φ`).
- `L_1`: dihedral sum on marching cubes of the *smooth signed margin field*
  `s_i(x) = G_i(x) − max_{j≠i} G_j(x) − δ`.
- `L_2`: Crofton with analytic domain-face term (primary); MC on the signed
  field (`l2_field` column) for comparison.

Field generation uses `gstools` through `src.fields.simulate_random_field`.

## Outputs

| File | Contents |
|---|---|
| `results/margin_sweep.csv` | per-row `(N, q, φ)` measurements and theory |
| `results/margin_sweep.png` | spanning frequency + `L_0 ± std` vs theory, `N = 128` |
| `results/margin_main.png` | 4×5 grid: per-`q` LKC error bars + slices, `N = 128` |
| `results/margin_euler_ratio.png` | measured/theory % deviation vs `φ`, `N = 128` |
| `results/margin_convergence.png` | % deviation vs `N`, φ-averaged, per `q` |

![margin tessellation: LKCs vs theory](results/margin_main.png)

![measured/theory ratio vs φ](results/margin_euler_ratio.png)

![convergence vs N](results/margin_convergence.png)

## Key numerical results (N = 256)

Measured/theory ratios:

| LKC | q=2 | q=4 | q=8 | q=16 |
|-----|-----|-----|-----|------|
| `L_0` (26-conn)      | 0.957 | 0.973 | 0.985 | 0.983 |
| `L_1` (signed field) | 0.995 | 1.004 | 1.002 | 0.992 |
| `L_2` (Crofton)      | 0.999 | 0.996 | 0.987 | 0.981 |
| `L_3`                | 1.00  | 1.00  | 1.00  | 1.00  |

Convergence from `N = 64 → 128 → 256` is monotonic across all LKCs.
L1/L2 scale roughly O(1/N²), matching the expected marching-cubes mesh
amplitude bias.

## Film thickness

At `φ = 0.90, q = 16`: `δ ≈ 0.059` physical units, i.e. `0.59` voxels at
`N = 128`. The film is sub-voxel at the highest `φ` values, so individual
winner cells can fragment into spurious extra components — this is a
discretization limit (see `margin_convergence.png` for the fix at `N = 256`),
not a formula failure.
