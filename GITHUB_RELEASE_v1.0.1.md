# v1.0.1 — post-v1.0.0 audit consolidation

This reproducibility release preserves the v1.0.0 source-model core and adds the
numerical audits incorporated through manuscript v10.2.

## Added

- compatibility-space spectral projection near the extreme high-SAM tail;
- reduced-cusp nondegeneracy / unfolding check;
- root-specific numerical-spread audit for `mu_inf`;
- direct native physiological meal-cycle integration;
- one-command post-v1.0.0 audit runner.

## Corrected

- validation text now uses the correct simple-pole convention
  `B1 = +10763.730351...` and `K = B1/B0'`.

## Scope

No BioModels kinetic equation or canonical v1.0.0 source input is changed. The
new files document and reproduce post-v1.0.0 audit calculations. The repository
continues to exclude the article and Supplementary Material.
