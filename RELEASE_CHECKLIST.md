# GitHub → Zenodo release checklist for v1.0.1

The release object contains code, SBML inputs, numerical data, and reproducibility
metadata only. It contains no article or Supplementary Material.

## 1. Verify the release-ready commit

From a clean checkout:

```bash
python code/verify_release.py
python code/reproduce_v10_2_audits.py --skip-meal
```

Optionally rerun the complete native meal-cycle audit:

```bash
python code/reproduce_v10_2_audits.py
```

Record the exact commit:

```bash
git rev-parse HEAD
```

## 2. Create GitHub release v1.0.1

Tag the exact audited commit:

```bash
git tag -a v1.0.1 -m "v1.0.1: post-v1.0.0 audit consolidation"
git push origin v1.0.1
```

Create a GitHub Release for tag `v1.0.1` using
`GITHUB_RELEASE_v1.0.1.md` as release notes.

## 3. Verify the Zenodo archive

After Zenodo imports the GitHub release, verify:

- version `1.0.1`;
- author/title metadata;
- open access and MIT license for original repository content;
- all three BioModels SBML inputs;
- `V1_0_1_ADDITIONS_MANIFEST.tsv` and `SHA256SUMS_v1.0.1.txt`;
- no manuscript or Supplementary Material files.

Record the new **version DOI**. The concept DOI remains
`10.5281/zenodo.22162042`.

## 4. Post-release metadata update

After Zenodo assigns the v1.0.1 DOI:

1. add the version DOI to `CITATION.cff` and the README on `main`;
2. update the article Data and Code Availability section from v1.0.0 to the
   exact v1.0.1 commit and Zenodo version DOI;
3. do not alter the tagged v1.0.1 scientific release object.
