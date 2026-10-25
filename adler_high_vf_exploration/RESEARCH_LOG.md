# Research log

## 2026-08-17 — Scope

Objective: identify a defensible, potentially publishable extension of the
Roubin–Colliat excursion-set framework for high-volume, non-percolating
inclusion morphologies, without losing the Gaussian kinematic formula.

Initial candidate families:

- scalar Gaussian fields with multi-interval hitting sets;
- chi-square / radial vector-Gaussian fields;
- unions of orthants or half-spaces in Gaussian value space;
- Gaussian vector fields with competing-label hitting sets;
- altered spectra as a secondary percolation-control parameter.

Important negative control: strictly monotone transformations of a scalar
Gaussian field do not change excursions at fixed volume fraction.

## 2026-08-17 — Gaussian margin tessellation pilot

Introduced the vector-Gaussian hitting set

\[
D_\delta=\{z:\max_i z_i-\operatorname{secondmax}_i z_i\geq\delta\}.
\]

It is the union of (q) mutually separated polyhedral chambers. Its preimage
retains the cores of winner cells and removes a matrix film around tie
interfaces. The exact volume fraction is a one-dimensional Gaussian integral;
the Gaussian surface functional is another one-dimensional integral.

Five (64^3), (L_c=5)-voxel realizations gave:

- (q=4): spanning already at 50% in all five realizations;
- (q=8): no spanning through 60%, mixed at 70%, spanning at 80%;
- (q=16): no spanning through 90%;
- (q=32): no spanning through 90%.

At (q=16,\phi=0.90), the largest same-winner component averaged about 1.4%
of retained material. At (q=32,\phi=0.90), it averaged about 0.6%.

A three-realization (96^3,L_c=8) check also found no spanning for (q=16)
through 85%. This supports, but does not establish, resolution robustness.

Critical caveat: as \(\phi\to1\), \(\delta\to0\) and the matrix film becomes
arbitrarily thin. A publishable study must report film thickness in physical
or voxel units and perform a convergence study. Label-aware connectivity is
mathematically appropriate for \(\delta>0\), but a manufacturable morphology
also requires the separating film to be resolved.

Literature searches found standard GKF/tube-formula work and Gaussian-field
tessellations, but no direct match for using a winner–runner-up margin hitting
set to create high-volume disconnected material with analytic expected LKCs.
This is only a preliminary novelty check, not a definitive claim.

## 2026-08-17 — Four Gaussian Minkowski functionals

Directly fitting third derivatives of Monte Carlo tube probabilities was found
to be numerically unstable and failed the exact (q=2) control. The reliable
alternative uses periodic fields in spatial dimensions 1, 2 and 3, for which
Euler density isolates (M_1,M_2,M_3), respectively.

At target fraction 0.70, the periodic estimates for (q=2) agreed closely with
all four exact half-space-union values. Estimates for (q=4,8,16) show a change
of (M_3), hence periodic Euler density, from negative at (q=4) to positive at
(q=8,16). This independently supports the proposed connectivity transition.

The polyhedral normal-cone calculation then yielded closed one-dimensional
integral formulas for all four GMFs. They use total codimension-(k) stratum
measures (S_k), a pairwise wedge angle (\pi/3), and the triple-normal solid
angle (3\arccos(1/3)-\pi). The formulas reduce exactly to the (q=2)
half-space identities and agree with independent periodic simulations at target
fractions 0.50, 0.70 and 0.80, within the expected finite-grid bias.

## 2026-08-17 — Bounded 100^3 LKC sweep

Ran a nonperiodic (100\times100\times100) cube with side (T=100),
(L_c=10<15), (q=4,8,16), ten margin values and five cubical / two
mesh-based realizations. The theoretical curves use the complete boundary
terms of Roubin System 2.17.

(L_3) and (L_0) closely follow theory. A marching-cubes surface estimator
also gives reasonable (L_2) agreement, with residual voxel-resolution bias.
The VTK integrated-mean-curvature estimator for (L_1) remains systematically
high and must not yet be treated as a validation failure or success; a proper
Crofton or signed-dihedral estimator and a resolution study are required.
