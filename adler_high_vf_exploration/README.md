# Adler-linked high-volume disconnected excursions

Research sandbox for Gaussian-related random-field constructions that aim to
increase the volume fraction of disconnected components while preserving the
Gaussian kinematic formula (GKF), called the Adler expectation link in the
source thesis.

Hard constraints:

1. The morphology is an excursion/preimage of a smooth vector-valued Gaussian
   field under an explicitly defined hitting set.
2. The Gaussian Minkowski functionals of that hitting set must be analytic or
   systematically computable, so expected LKCs remain available through GKF.
3. Claims about disconnection are checked with component and wrapping
   statistics, not inferred from expected Euler characteristic alone.
4. Comparisons are made at equal volume fraction and adequate scale ratio.

Python commands in this folder must run through the `rf_explo` mamba
environment.

Generate the saved PNG figures with:

```powershell
mamba run -n rf_explo python make_plots.py
```

The figures are written to `figures/`.

Estimate all four Gaussian Minkowski functionals from periodic Euler-density
experiments with:

```powershell
mamba run -n rf_explo python numerical_gmfs.py
```

Results are written to `results/numerical_gmfs/`.

Generate bounded-cube (L_i(\delta)) comparisons on a (100^3) grid with
(L_c=10<15) using:

```powershell
mamba run -n rf_explo python run_bounded_lkc_sweep.py --n 100 --side 100 --lc 10 --q 4 8 16 --estimator mesh
```

The script rejects correlation lengths greater than or equal to 15.
