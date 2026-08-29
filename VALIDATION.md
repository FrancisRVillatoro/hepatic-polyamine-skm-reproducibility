# Validation record for v1.0.0

## Frozen closure workflow

The scientific components of the pre-release workflow were rerun from a clean
extraction after the calculations were frozen. The archived validation included:

- exact reproduction of the 2006 Monte Carlo stable fractions;
- six directional-ensemble rows agreeing with canonical numerical quantities to below `1e-16`;
- regeneration of 279 states over nine continuation paths;
- 125/125 stable sampled equilibria on the 5x5x5 control cube;
- 81/81 stable sampled equilibria on the refined `r_H=1` face;
- pseudo-arclength branch columns agreeing with archived values to `4.44e-16` and fold quantities to `1.42e-14`;
- exact reproduction of the sampled `tau_Y` sweep;
- repetition of all eight dynamic up/down sweeps, with largest change in a physical switching value `4.63e-9` in the feedback amplitude `a`.

Adaptive BDF internal step counts are not scientific invariants.

## Substrate-input cross-check

The repository implementation of the substrate-input core was cross-checked
against the frozen canonical CSVs at 21 values of `rho`. Maximum absolute
discrepancies were `1.14e-13` micromolar in SAM, `4.97e-14` micromolar in
intracellular methionine, `8.92e-10` in the methionine-specific logarithmic
gain, `7.04e-10` in the common-amino-acid gain, and `8.95e-14 h^-1` in the
spectral abscissa.

## Final v9 MAT-control regeneration

The v9-specific workflows were regenerated from the release code. At the
proliferative endpoint the results are:

| MAT policy | SAM (uM) | G_Met | G_AA | alpha_sp (h^-1) |
| --- | ---: | ---: | ---: | ---: |
| physiological reference | 65.061382 | 1.779231 | 1.342716 | -0.015876 |
| published | 41.339452 | 0.738656 | 0.533651 | -0.012826 |
| constant nominal total 480 | 115.746153 | 2.020572 | 1.504599 | -0.011412 |
| baseline-flux matched MATII | 62.685528 | 1.073305 | 0.793140 | -0.011835 |

The physiological basal MAT flux is `205.60042387875146`; matching MATII to
that throughput gives `Vmax_MATII = 261.9019248346911`. Relative to the
physiological reference, the published gains fall by 58.5%/60.3%, the
constant-total gains rise by 13.6%/12.1%, and the flux-matched gains fall by
39.7%/40.9% (G_Met/G_AA).

The three rho=1 methionine PALC branches were regenerated. No sign change of
`d log(mu)/ds` occurred through the common range to `mu=12`; extensions reached
`mu=18.72` (published), `15.65` (constant total), and `17.88` (flux matched),
with negative spectral abscissa throughout the reported points.

## Simple-pole substrate-saturation audit

The final asymptotic audit regenerated:

- `mu_inf = 2.151802539828513`;
- tail-fit coefficient `C = 131.9096636580272`;
- tail-fit exponent `p = 1.0003121159227146`;
- reduced-balance derivative `B0'(mu_inf) = 82.01063066735516`;
- effective `-B1 = 10763.730351293078`;
- balance coefficient `K = 131.24798899489073` micromolar.

At the limiting balance, input and output fluxes both equal approximately
`252.0012071`. Input is 91.52% MATIII and 8.48% MATI; output is 71.43% DNMT,
23.02% GNMT, and 5.55% SAMDC. The simple-pole statement is conditional on the
local nonsingularity of the remaining constrained non-SAM subsystem and concerns
the connected equilibrium branch only.

## SBML identity and compatibility-class audit

The physiological and proliferative SBML files have identical exact conserved
pool totals: cytosolic folate `13.4`, mitochondrial folate `40.2`, and
CoA+acetyl-CoA `199.5`. The automated SBML diff classifies the substantive
autonomous differences, after freezing meal bookkeeping at basal input, as the
MAT reparameterization, inverse-SAM ODC/SAMDC feedback, and H2O2 change.
