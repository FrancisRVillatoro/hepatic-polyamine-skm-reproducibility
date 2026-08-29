# Hepatic polyamine SKM — reproducibility code and data

This repository is the code-and-data reproducibility package for a computational
study of hepatic polyamine metabolism. **It intentionally does not contain the
article manuscript, Supplementary Material, or journal-submission files.**

## Contents

- `code/` — executable Python workflows and the three source BioModels SBML files.
- `data/` — canonical CSV outputs used for reported numerical quantities.
- `requirements.txt` — supported dependency ranges.
- `requirements-lock.txt` — exact package versions used for the frozen workflow.
- `environment_versions.txt` — recorded software environment.
- `SHA256SUMS.txt` — cryptographic checksums for the release contents.

The BioModels inputs are:

- `BIOMD0000000190` — mammalian polyamine model;
- `BIOMD0000000674` — physiological integrated hepatic model;
- `BIOMD0000000450` — proliferative integrated hepatic model.

## Main reproducibility workflows

Create a clean Python environment and install the locked dependencies when
possible:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows PowerShell: .venv\\Scripts\\Activate.ps1
pip install -r requirements-lock.txt
```

Run the fast deterministic closure checks:

```bash
python code/reproduce_closures.py
```

Run the original full closure campaign (Monte Carlo, path/cube calculations,
47-state continuation, timescale sweep, controls, and dynamic hysteresis):

```bash
python code/reproduce_closures.py --full
```

The substrate-input extension is kept separate because the high-SAM
pseudo-arclength branches are comparatively expensive. A smoke test is:

```bash
python code/reproduce_substrate_inputs.py --smoke
```

The full substrate-input workflow is:

```bash
python code/reproduce_substrate_inputs.py
```

It recomputes basal logarithmic gains and the methionine-only/common-amino-acid
PALC branches into `data/recomputed_substrate/`, leaving canonical CSV files
untouched.


## v9 control calculations

The final pre-submission audit added three MAT-policy controls and an analytical
check of the physiological high-SAM boundary. Regenerate the v9-specific fast
calculations with:

```bash
python code/reproduce_v9_controls.py
```

This recomputes the published, constant-nominal-capacity, and baseline-flux-matched
MAT gain sweeps, the automated SBML diff/conserved-pool audit, and the simple-pole
SAM-balance reduction. To regenerate the three rho=1 methionine PALC branches as
well, run:

```bash
python code/reproduce_v9_controls.py --full
```

Use `--figures` to regenerate the four v9 control figures into `figures_v9/`.
The release retains canonical numerical CSVs in `data/`; these files are the
source for the manuscript's reported v9 control values.

## Canonical results represented in `data/`

The archive contains, among other outputs:

- 2006 SBML audit and augmented SKM elasticities/robustness;
- physiological-to-proliferative equilibrium continuation;
- exact positive-direction decomposition of 19 signed reversible rates;
- anchored structural-robustness ensembles and cancellation diagnostics;
- factorial MAT / inverse-SAM-feedback / H2O2 controls;
- asynchronous continuation paths, a 5x5x5 control cube, and the refined H2O2 face;
- constant-total-nominal-MAT-capacity, baseline-flux-matched MATII, and fixed-H2O2 controls;
- reduced regulatory-gain sensitivity and cusp calculations;
- full 47-state pseudo-arclength branch, folds, and timescale sweep;
- finite-rate hysteresis data;
- methionine-only and common-amino-acid substrate-input PALC data;
- logarithmic SAM input gains across the proliferative coordinate;
- v9 MAT-policy gain/PALC controls, automated SBML identity/pool audit, and the simple-pole substrate-saturation calculation.

The legacy CSV column names `a_SN_on` and `a_SN_off` in some archived data mean
the lower- and upper-parameter folds, respectively. The final analysis uses the
unambiguous notation `a_-` and `a_+`; canonical CSV headers were not changed so
that frozen numerical files remain byte-stable.

## Validation

The frozen workflow was rerun component by component from a clean extraction.
See `VALIDATION.md` for the archived closure tolerances, substrate-input cross-check, and final v9-control regeneration. Canonical CSV files are retained as numerical
artifacts; reproduction scripts write new files rather than silently replacing
them.

## AI-assisted development

See `AI_ASSISTANCE.md` for the disclosure of supportive use of OpenAI ChatGPT
5.6 and Anthropic Claude Opus 5. All scientific outputs remain subject to human
verification and responsibility.

## License and third-party inputs

Original repository code and metadata are MIT licensed. The three BioModels
SBML files are distributed by BioModels under CC0 1.0; see
`THIRD_PARTY_NOTICES.md`.

## Funding

Development of this research was supported by the Universidad de Málaga through its II Plan Propio de Investigación, Transferencia y Divulgación Científica (PPRO).

## Citation

For the journal submission, cite the immutable **Zenodo v1.0.0 DOI** corresponding
to this exact code-and-data release. After the first GitHub/Zenodo release is
created, place the version DOI in the article's Data and Code Availability
section. The concept DOI may additionally be displayed in this README to point
to the latest repository version.
