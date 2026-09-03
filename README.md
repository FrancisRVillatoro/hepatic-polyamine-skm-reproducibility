# Hepatic polyamine SKM — reproducibility code and data

This repository is the code-and-data reproducibility package for a computational
study of hepatic polyamine metabolism. **It intentionally does not contain the
article manuscript, Supplementary Material, or journal-submission files.**

## Versions and archived releases

- **v1.0.0** is the first frozen archive (Zenodo version DOI `10.5281/zenodo.22162043`).
- **v1.0.1** consolidates the subsequent mathematical/numerical audits used in manuscript v10.2 (Zenodo version DOI `10.5281/zenodo.22274449`).
- Concept DOI for all versions: **10.5281/zenodo.22162042**.

The v1.0.1 additions are reproducibility checks, not a change of the source BioModels kinetics: compatibility-space spectral projection near the extreme high-SAM tail; reduced-cusp nondegeneracy; root-specific numerical spread for the simple-pole boundary; and direct integration of the native physiological meal schedule.

## Contents

- `code/` — executable Python workflows and the three source BioModels SBML files.
- `data/` — canonical numerical outputs and compact audit summaries.
- `requirements.txt` / `requirements-lock.txt` — dependency specifications.
- `SHA256SUMS.txt` — v1.0.0 baseline checksums.
- `SHA256SUMS_v1.0.1.txt` — checksums of canonical v1.0.1 audit data.
- `V1_0_1_ADDITIONS_MANIFEST.tsv` — v1.0.1 delta inventory.

BioModels inputs: `BIOMD0000000190`, `BIOMD0000000674`, and `BIOMD0000000450`.

## Core workflows

```bash
python code/reproduce_closures.py
python code/reproduce_substrate_inputs.py
python code/reproduce_v9_controls.py
```

## v1.0.1 audit workflows

Fast post-v1.0.0 audits:

```bash
python code/reproduce_v10_2_audits.py --skip-meal
```

Complete audit set, including the native 26-day meal-cycle convergence calculation:

```bash
python code/reproduce_v10_2_audits.py
```

Individual workflows:

```bash
python code/compatibility_spectral_projection.py
python code/cusp_nondegeneracy_check.py
python code/mu_inf_precision_audit.py
python code/native_meal_cycle_control.py
```

The native-meal calculation retains the original `floor`/`piecewise` schedule in BIOMD0000000674 and splits integration at every schedule discontinuity. The root-spread calculation uses frozen v1.0.0 asymptotic tables and reports a cross-method numerical spread, not a statistical confidence interval.

## Canonical v1.0.1 audit outputs

- `data/compatibility_spectral_projection_tail.csv`
- `data/compatibility_spectral_projection_summary.json`
- `data/cusp_nondegeneracy_check.json`
- `data/mu_inf_precision_estimates_v10_2.csv`
- `data/mu_inf_precision_summary_v10_2.json`
- `data/native_meal_cycle_day26_keypoints.csv`
- `data/native_meal_cycle_summary_v10_2.json`

`native_meal_cycle_control.py` regenerates the denser day-26 trajectory when rerun. The pre-existing v1.0.0 canonical outputs remain unchanged.

## Validation

See `VALIDATION.md`. Static integrity:

```bash
python code/verify_release.py
```

## AI-assisted development

See `AI_ASSISTANCE.md`. Scientific outputs remain subject to human verification and responsibility.

## License and third-party inputs

Original repository code and metadata are MIT licensed. The three BioModels SBML files are distributed by BioModels under CC0 1.0; see `THIRD_PARTY_NOTICES.md`.

## Funding

Development of this research was supported by the Universidad de Málaga through its II Plan Propio de Investigación, Transferencia y Divulgación Científica (PPRO).

## Citation

For v1.0.0, use Zenodo DOI `10.5281/zenodo.22162043`. For v1.0.1, use Zenodo DOI `10.5281/zenodo.22274449`. The concept DOI `10.5281/zenodo.22162042` always resolves to the latest archived repository version.
