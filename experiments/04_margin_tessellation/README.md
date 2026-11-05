# Experiment 04 — Margin tessellation

Let $G_1,\ldots,G_q$ be independent unit-variance Gaussian fields with
correlation $\exp(-r^2/l_c^2)$. Write $G_{(q)}$ and $G_{(q-1)}$ for the
largest and second-largest values. Include a point when

$$
G_{(q)}(x)-G_{(q-1)}(x)\geq\delta.
$$

The scalar margin $\delta$ controls expected VF. For $\delta>0$, different
winner labels are separated in the continuous model. This does not prevent
a single label from spanning, and does not guarantee that a sampled binary
phase resolves the separating matrix film.

![Margin tessellation](paraview/q16_phi70/display_margin_tess_RF_thresh.png)

## GMFs and their scope

The expected VF is

$$
\mathcal M_0=q\int_{\mathbb R}\varphi(t)\Phi(t-\delta)^{q-1}\,dt.
$$

Define $a=\delta/\sqrt2$ and, for $1\leq k\leq\min(3,q-1)$,

$$
S_k=q\binom{q-1}{k}\sqrt{k+1}
\int_{\mathbb R}\varphi(t+\delta)\varphi(t)^k
\Phi(t)^{q-1-k}\,dt.
$$

Set absent strata to zero. For $\delta>0$, the implemented coefficients are

$$
\begin{aligned}
\mathcal M_1&=S_1,\\
\mathcal M_2&=aS_1+\tfrac\pi3 S_2,\\
\mathcal M_3&=(a^2-1)S_1+\tfrac{4a}{\sqrt3}S_2
+2(3\arccos(1/3)-\pi)S_3.
\end{aligned}
$$

Insert these in the [cube GKF](../../docs/framework.md). These are
one-dimensional integral formulas with an outlined stratified-tube
derivation; the complete proof remains outstanding. At exactly
$\delta=0$, the union is all value space, so the correct GMFs are
$(1,0,0,0)$ and the expected LKCs are those of the domain. Higher GMFs
need not approach these values as $\delta\downarrow0$.

## What the committed sweep shows

`results/margin_sweep.csv` contains 480 rows:
three resolutions (64, 128, 256), four values of $q$ (2, 4, 8, 16),
eight target VFs (30, 40, 50, 60, 70, 80, 85, 90%), and five realisations.
The domain has $L=64,l_c=5$.

At $N=256,q=16$, target VF 60%, ensemble-mean relative errors are:

| Descriptor | Error |
|---|---:|
| Per-winner Euler, 26 neighbours | −1.53% |
| Sampled margin-field $L_1$ | −0.77% |
| Crofton $L_2$ | −1.89% |
| Volume | −0.09% |

This setting supports good mean agreement; it does not establish a
global under-2% bound. For example, the mean Euler error at $N=256,q=2$,
target VF 30%, is about −16.4%. Percent errors are unstable near theoretical
zeros. The committed results do not establish a universal convergence rate.

The runner sums per-winner 26-neighbour Euler values. Its spanning diagnostic
uses per-winner 6-neighbour components. Neither is identical to applying
both conventions to the binary union. See
[measurement conventions](../../docs/measurement_conventions.md).

## Saved scenes and connectivity

The saved `paraview/q16_phi60_n256/phase.vti` has actual VF 60.77%.
Rechecking that same binary phase gives:

| Neighbours | Components | Largest fraction of inclusion volume | Spanning |
|---|---:|---:|---|
| 6 | 3,206 | 2.04% | None |
| 26 | 1,838 | 44.86% | All axes |

At target VF 75%, the saved scene's actual VF is 75.47% and it spans all
axes under both conventions. Largest-component fractions are about
47.95% (6 neighbours) and 99.62% (26 neighbours).
The historical images and summary files retain their original 6-neighbour
statistics. They are not proof of continuum-disconnected grains.

## Ceiling study: finite evidence, not a VF guarantee

The historical Phase A CSV has 600 rows: six $q$ values, ten target VFs,
and ten independent fields per $q$. A field is reused across VF thresholds;
100 rows per $q$ are not 100 independent trials at a fixed threshold.

- $q=2,4,8$: all ten fields span at every tested VF from 70% to 98%.
- $q=16$: spanning counts rise from 0/10 at 70%, to 1/10 at 75%,
  2/10 at 80–92%, 3/10 at 94%, and 4/10 at 96–98%.
- $q=32,64$: no spanning observed in these ten fields at the tested thresholds.

These are **per-winner, 6-neighbour** results. Zero spans in ten trials has
a one-sided 95% binomial upper bound of about 25.9% for spanning probability.
It does not establish a safe operating ceiling.

Phase B has 90 rows: two $q$ values, three resolutions, three VFs
(85, 90, 95%), and five fields per setting. For $q=16$, counts were 0/5
at $N=64$ and 2/5 at each higher resolution; $q=32$ had 0/5.
Seeds differ across resolutions, so this is not a paired refinement study.

There is no demonstrated EC-zero at VF 85–90% for $q=16$ in this setup:
the implemented expected Euler remains positive through that range and
near VF 100% on the positive-margin branch.

## Margin units and resolution

$\delta$ has Gaussian-component **value units**. It is not a physical
thickness, and $\delta/dx$ is not a number of film voxels.
Near a regular two-field tie, with competing fields $G_i,G_j$ and no
other winner, the local full film thickness is approximately

$$
w(x)\simeq\frac{2\delta}{\|\nabla(G_i-G_j)(x)\|}.
$$

The gradient varies in space, so neither a single thickness nor a minimum
grid resolution follows from $\delta$ alone. Example unit-component margins
are about 0.05897, 0.05044, and 0.04462 at VF 90% for $q=16,32,64$.

Historical ceiling CSV columns `delta_physical` and `delta_voxels` are
misnamed: they store $\delta$ and $\delta/dx$, respectively.
Historical film-resolution figures derived from those columns are not valid
spatial-thickness estimates. They are preserved as historical artifacts.
The corrected ceiling runner writes `delta` and `dx_over_lc` and replaces
those plots with value-margin and measured-VF diagnostics.

## Running

From the repository root:

~~~bash
python experiments/04_margin_tessellation/run.py
python experiments/04_margin_tessellation/run_ceiling.py --help
~~~

The main sweep writes its historical `results/` location; preserve or move
those outputs before a deliberate rerun. The corrected ceiling runner
defaults to `results_ceiling_v2/` and refuses a nonempty destination.
Use `--output` for another destination. No new ceiling or convergence
campaign was run as part of the review corrections.

The scene exporters still report their explicitly labelled 6-neighbour
statistics. A full resolution study under both conventions, and ultimately
a consistent continuum geometry, remain necessary for physical aggregate
connectivity claims.
