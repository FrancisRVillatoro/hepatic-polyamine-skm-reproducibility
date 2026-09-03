# Validation record for v1.0.1

Release v1.0.1 preserves the v1.0.0 scientific core and adds the post-v1.0.0
audits used through manuscript v10.2. No source BioModels kinetic equation is
changed.

## Frozen v1.0.0 closure workflow

The archived v1.0.0 validation included:

- exact reproduction of the 2006 Monte Carlo stable fractions;
- six directional-ensemble rows agreeing with canonical numerical quantities to below `1e-16`;
- regeneration of 279 states over nine continuation paths;
- 125/125 stable sampled equilibria on the 5x5x5 control cube;
- 81/81 stable sampled equilibria on the refined `r_H=1` face;
- pseudo-arclength branch columns agreeing with archived values to `4.44e-16` and fold quantities to `1.42e-14`;
- exact reproduction of the sampled `tau_Y` sweep;
- repetition of all eight dynamic up/down sweeps, with largest change `4.63e-9` in the switching amplitude `a`.

Adaptive BDF internal step counts are not scientific invariants.

## Substrate-input cross-check

The repository substrate-input core was cross-checked against canonical CSVs at
21 values of `rho`. Maximum absolute discrepancies were `1.14e-13` micromolar in
SAM, `4.97e-14` micromolar in intracellular methionine, `8.92e-10` in the
methionine-specific logarithmic gain, `7.04e-10` in the common-amino-acid gain,
and `8.95e-14 h^-1` in the spectral abscissa.

## MAT-capacity / throughput controls

At the proliferative endpoint:

| MAT policy | SAM (uM) | G_Met | G_AA | alpha_sp (h^-1) |
| --- | ---: | ---: | ---: | ---: |
| physiological reference | 65.061382 | 1.779231 | 1.342716 | -0.015876 |
| published | 41.339452 | 0.738656 | 0.533651 | -0.012826 |
| constant nominal total 480 | 115.746153 | 2.020572 | 1.504599 | -0.011412 |
| baseline-flux matched MATII | 62.685528 | 1.073305 | 0.793140 | -0.011835 |

The physiological basal MAT flux is `205.60042387875146`; matching MATII to that
throughput gives `Vmax_MATII = 261.9019248346911`.

## Simple-pole substrate-saturation audit

The asymptotic calculation gives:

- `mu_inf = 2.151802539828513` from the nonlinear tail fit;
- tail-fit coefficient `C = 131.9096636580272`;
- tail-fit exponent `p = 1.0003121159227146`;
- reduced-balance derivative `B0'(mu_inf) = 82.01063066735516`;
- effective **`B1 = +10763.730351293078`**;
- balance coefficient `K = B1/B0' = 131.24798899489073` micromolar.

The positive sign of `B1` is the convention consistent with the written reduced
balance `B0(mu) + B1/SAM + ... = 0` and the released calculation.

The finite high-SAM nonsingularity diagnostic at `SAM=1e10` has
`sigma_min = 1.6423e-6` and scaled `kappa_2 = 5.2399e10`; this supports, but does
not analytically prove, nonsingularity of the limiting `epsilon=0` subsystem.
The root-specific audit in `mu_inf_precision_audit.py` gives a union of
cross-method estimates spanning `1.20e-8`, motivating the displayed
`mu_inf = 2.15180254` with numerical method spread of order `1e-8`.

## Compatibility-space spectral projection

`compatibility_spectral_projection.py` removes the three exact conservation
directions geometrically by constructing an orthonormal basis `Q` of the
compatibility tangent space and diagonalizing `Q.T @ J @ Q`.

On the extreme physiological high-SAM tail, representative projected physical
rightmost real parts are approximately:

- `-1.33e-8 h^-1` near `SAM=8.97e5 uM`;
- `-1.14e-11 h^-1` near `SAM=3.08e7 uM`;
- `-1.43e-14 h^-1` near `SAM=1.06e9 uM`.

This resolves the three tiny positive values in the legacy magnitude-filtered
column as conservation-mode identification artifacts rather than physical
instabilities.

## Reduced-cusp nondegeneracy

The archived cusp location is
`rho0=0.08180259323371566`, `Y=0.8813154700722712`,
`a=2.007110174251687`. Reconstructing the scalar equilibrium equation with the
same PCHIP metabolic branch gives `F=F_Y=F_YY=0` to numerical precision, with
`F_YYY approximately -4.804` and the `(a,rho0)` unfolding determinant
approximately `-0.812`. Both nonzero values support a nondegenerate codimension-two
cusp of the reduced scalar equilibrium problem.

## Native physiological meal-cycle control

The native BIOMD0000000674 schedule is retained:

- fasting multiplier `0.25`;
- breakfast `1.75`, 07:00--10:00;
- lunch `1.75`, 12:00--15:00;
- dinner `3.25`, 18:00--21:00.

Although dinner exceeds the sustained common-input boundary `eta_inf=2.44486`
for three hours, direct integration remains bounded and approaches a 24-hour
cycle. After 26 warm-up cycles, the tighter repeat has full-state start/end
relative mismatch about `9.5e-8`; SAM ranges from `14.476245` to `301.925780 uM`
and returns to `36.944782 uM` at midnight. Tightening BDF tolerances from
`2e-8/2e-10` to `2e-10/2e-12` changes the SAM maximum by less than `3e-7 uM`.
This establishes that the sustained-input boundary is not a threshold that can
be directly applied to the finite dinner pulse.

## Static release integrity

`python code/verify_release.py` verifies:

- no manuscript/submission `.tex`, `.pdf`, `.bib`, `.doc`, or `.docx` files;
- presence of the three source SBML files and v1.0.1 audit scripts/data;
- `VERSION == 1.0.1`;
- all entries in `SHA256SUMS_v1.0.1.txt`.
