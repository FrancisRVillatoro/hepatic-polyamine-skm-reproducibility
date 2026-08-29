# GitHub → Zenodo release checklist

Suggested GitHub repository name:

`hepatic-polyamine-skm-reproducibility`

The release object contains only code, SBML inputs, numerical data, and
reproducibility metadata. It contains no article or supplementary manuscript.

## 1. Create and push the GitHub repository

From the directory containing this file:

```bash
git init
git branch -M main
git add .
git commit -m "Reproducibility code and data for hepatic polyamine SKM"
git remote add origin git@github.com:<USER>/hepatic-polyamine-skm-reproducibility.git
git push -u origin main
```

Record the immutable commit:

```bash
git rev-parse HEAD
```

## 2. Connect the repository to Zenodo

In Zenodo, connect the GitHub account and enable archiving for this repository.
Do this before creating the GitHub release.

## 3. Create GitHub release v1.0.0

Tag the exact audited commit:

```bash
git tag -a v1.0.0 -m "v1.0.0: code and data supporting the JTB submission"
git push origin v1.0.0
```

Then create a GitHub Release for tag `v1.0.0`. Suggested release title:

`v1.0.0 — hepatic polyamine SKM reproducibility code and data`

Suggested release notes are in `GITHUB_RELEASE_v1.0.0.md`.

## 4. Verify the Zenodo archive

After Zenodo imports the GitHub release, verify:

- title and author metadata;
- version `1.0.0`;
- open access;
- MIT license for original repository content;
- uploaded ZIP contains no manuscript `.tex`/`.pdf` files;
- the three BioModels SBML inputs are present and their CC0 status is documented;
- `SHA256SUMS.txt` is present.

Record both identifiers:

- **version DOI** — cite this exact v1.0.0 DOI in the article;
- **concept DOI** — optional README badge/link for all versions.

## 5. Final manuscript replacement

Replace the three placeholders in the article's Data and Code Availability text
with:

1. public GitHub URL;
2. immutable v1.0.0 commit SHA;
3. Zenodo **version DOI**.

Do not upload the article or Supplementary Material to this reproducibility
repository merely to make those links resolvable.

## 6. Optional post-release README update

After Zenodo assigns the DOI, the default `main` branch README may be updated to
show the DOI. The archived v1.0.0 object remains immutable. If the scientific
code/data themselves change, create a new versioned release rather than silently
altering v1.0.0.
