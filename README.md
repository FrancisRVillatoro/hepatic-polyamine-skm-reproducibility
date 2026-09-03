# Hepatic polyamine SKM — reproducibility code and data

This repository is the code-and-data reproducibility package for a computational
study of hepatic polyamine metabolism. **It intentionally does not contain the
article manuscript, Supplementary Material, or journal-submission files.**

## Versions and archived releases

- **v1.0.0** is the first frozen archive (Zenodo version DOI
  `10.5281/zenodo.22162043`).
- **v1.0.1** consolidates the subsequent mathematical/numerical audits used in
  manuscript v10.2. Its GitHub tag should be archived by Zenodo as a new version.
- Concept DOI for all versions: **10.5281/zenodo.22162042**.

The v1.0.1 additions are reproducibility checks, not a change of the source
BioModels kinetics. They include:

1. compatibility-space spectral projection near the extreme high-SAM tail;
2. nondegeneracy checks for the reduced regulatory cusp;
3. root-specific numerical-spread audit for the simple-pole input boundary;
4. direct integration of the native physiological meal schedule.

## Contents

- `code/` — executable Python workflows and the three source BioModels SBML files.
- `data/` — canonical CSV/JSON numerical outputs used for reported quantities.
- `requirements.txt` — supported dependency ranges.
- `requirements-lock.txt` — exact package versions used for the frozen workflow.
- `environment_versions.txt` — recorded software environment.
- `SHA256SUMS.txt` — v1.0.0 baseline checksums.
- `SHA256SUMS_v1.0.1.txt` — checksums of v1.0.1 additions/updated metadata.
- `V1_0_1_ADDITIONS_MANIFEST.tsv` — v1.0.1 delta inventory.

The BioModels inputs are:

- `BIOMD0000000190` — mammalian polyamine model;
- `BIOMD0000000674` — physiological integrated hepatic model;
- `BIOMD0000000450` — proliferative integrated hepatic model.

## Core reproducibility workflows

Create a clean Python environment and install the locked dependencies when
possible:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows PowerShell: .venv\\Scripts\\Activate.ps1
pip install -r requirements-lock.txt
```

Fast deterministic closure checks:

```bash
python code/reproduce_closures.py
```

Original full closure campaign:

```bash
python code/reproduce_closures.py --full
```

Substrate-input workflow:

```bash
python code/reproduce_substrate_inputs.py --smoke
python code/reproduce_substrate_inputs.py
```

v9 MAT-policy and asymptotic controls:

```bash
python code/reproduce_v9_controls.py
python code/reproduce_v9_controls.py --full
```

## v1.0.1 audit workflows

Run the fast post-v1.0.0 audits while skipping the 26-day meal integration:

```bash
python code/reproduce_v10_2_audits.py --skip-meal
```

Run the complete v1.0.1 audit set, including the native meal schedule:

```bash
python code/reproduce_v10_2_audits.py
```

The individual workflows are:

```bash
python code/compatibility_spectral_projection.py
python code/cusp_nondegeneracy_check.py
python code/mu_inf_precision_audit.py
python code/native_meal_cycle_control.py
```

The native-meal calculation retains the original `floor`/`piecewise` schedule in
BIOMD0000000674 and splits numerical integration at every schedule discontinuity.
The root-spread calculation uses the frozen v1.0.0 asymptotic tables and reports a
cross-method numerical spread; it is not a statistical confidence interval.

## Canonical v1.0.1 audit outputs

- `data/compatibility_spectral_projection_tail.csv`
- `data/compatibility_spectral_projection_summary.json`
- `data/cusp_nondegeneracy_check.json`
- `data/mu_inf_precision_estimates_v10_2.csv`
- `data/mu_inf_precision_summary_v10_2.json`
- `data/native_meal_cycle_day26_tight.csv`
- `data/native_meal_cycle_summary_v10_2.json`

The pre-existing v1.0.0 canonical outputs remain unchanged.

## Validation

See `VALIDATION.md` for the frozen v1.0.0 closure results and the v1.0.1 audit
checks. Run

```bash
python code/verify_release.py
```

to verify the v1.0.1 release-delta checksums and required files.

## AI-assisted development

See `AI_ASSISTANCE.md` for disclosure of supportive use of OpenAI ChatGPT and
Anthropic Claude. Scientific outputs remain subject to human verification and
responsibility.

## License and third-party inputs

Original repository code and metadata are MIT licensed. The three BioModels SBML
files are distributed by BioModels under CC0 1.0; see `THIRD_PARTY_NOTICES.md`.

## Funding

Development of this research was supported by the Universidad de Málaga through
its II Plan Propio de Investigación, Transferencia y Divulgación Científica
(PPRO).

## Citation

For v1.0.0, use Zenodo DOI `10.5281/zenodo.22162043`.
For v1.0.1, use the new Zenodo **version DOI** assigned after the GitHub tag is
archived. The concept DOI `10.5281/zenodo.22162042` always resolves to the latest
archived repository version.
