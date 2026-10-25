# Experiment 02 — Student-t Random Fields

## Scientific question

Can heavier-tailed marginal distributions push the disconnected volume fraction
ceiling above the Gaussian ~16%, while keeping the Adler formula exact?

## Adler link status: **Full**

The Worsley t_nu field T(x) = Z_0(x) / sqrt(chi2_field(x) / nu) has a spatially
varying denominator. This is the construction for which Worsley (1994) proved that
the Adler formula applies with t_nu Minkowski functionals M_j (see
theory/student_t_derivation.md §5). The second spectral moment lambda_2 = 2/lc²
is the same as for the underlying Gaussian Z_0.

What changes: M_j have polynomial tails instead of Gaussian tails. This shifts the
volume fraction at the Adler percolation threshold (upper zero-crossing of E[L0])
to higher values.

**The disconnected VF ceiling rises from ~16% (Gaussian) to ~40–45% for nu=3.**

## Key parameters

| Parameter | Value | Effect |
|---|---|---|
| `N` | 80 | Grid resolution |
| `LC` | 8.0 | Correlation length |
| `NU_VALUES` | [None, 30, 10, 5, 3] | None = Gaussian; lower nu = heavier tails |
| `SEED` | 42 | Single-realisation seed |
| `N_RUNS_MC` | 20 | MC runs per case |

Valid range for nu: >= 3. At nu=2, variance is infinite; at nu=1 (Cauchy), lambda_2
is undefined and the Adler formula breaks.

## How to run

```bash
mamba run -n rf_explo python run.py
```

Expected runtime: ~15–20 minutes.

## Outputs

| File | Contents |
|---|---|
| `results/student_t_sweep.csv/.png` | Threshold sweep for all nu, measured vs Adler |
| `results/mc_student_t.csv/.png` | MC mean ± std for Gaussian, t_5, t_3 |
| `results/percolation_diagram.csv/.png` | phi_c+ vs nu: the VF ceiling curve |

## Key findings

- Adler formula with t_nu M_j matches MC measurements for nu >= 3.
- phi_c+ (max disconnected VF from EC criterion) increases from 16% (Gaussian, nu=inf)
  to ~40–45% at nu=3.
- The percolation diagram shows a smooth monotone curve phi_c+(nu), confirming that
  heavier tails systematically raise the disconnected VF ceiling.
- Caveat: phi_c+ from the EC criterion (E[L0]=0) is a theoretical estimate; the
  actual spatial percolation transition may differ. A systematic comparison between
  E[L0]=0 and finite-size wrapping probability has not been performed.
