# Publishing shannonpore to bioconda

Goal: `conda install -c bioconda shannonpore` installs the CLI **and**
modkit + pysam in one step (something pip cannot do, since modkit is a
Rust binary distributed on bioconda as `ont-modkit`).

## Package naming (resolved)

The package previously installed a top-level Python module named `src`,
which bioconda review would have rejected (`import src` clashes with
any other package that made the same mistake). The package directory is
now `shannonpore/` and the entry point is `shannonpore.cli:main` — no
further action needed.

## Release checklist (per version)

1. Merge outstanding feature PRs into `main`.
2. Bump `version` in `pyproject.toml` and `shannonpore/__init__.py`;
   move the `## Unreleased` CHANGELOG section under the new version
   heading.
3. Tag and release:
   ```bash
   git tag -a v0.2.0 -m "shannonpore 0.2.0"
   git push ebenstein v0.2.0
   gh release create v0.2.0 --repo ebensteinLab/ShannonPore \
       --title "shannonpore 0.2.0" --notes-file <(changelog section)
   ```
4. Compute the tarball checksum and fill it into the recipe:
   ```bash
   curl -sL https://github.com/ebensteinLab/ShannonPore/archive/refs/tags/v0.2.0.tar.gz | sha256sum
   ```

## Bioconda submission (first time)

1. Fork https://github.com/bioconda/bioconda-recipes and branch
   (e.g. `add-shannonpore`).
2. Copy `conda-recipe/meta.yaml` → `recipes/shannonpore/meta.yaml`
   (sha256 filled in).
3. Lint + build locally (optional but saves review round-trips):
   ```bash
   conda create -n bioconda-build -c conda-forge -c bioconda bioconda-utils
   bioconda-utils lint recipes config.yml --packages shannonpore
   bioconda-utils build recipes config.yml --packages shannonpore
   ```
4. Open the PR against bioconda-recipes `master`; respond to the
   automated lint bot and the human reviewer. Once merged, the package
   appears on the bioconda channel within a few hours.

## Subsequent version updates

The BiocondaBot usually opens the version-bump PR automatically when a
new GitHub release appears (it re-computes the sha256). If it doesn't,
edit `recipes/shannonpore/meta.yaml` (version + sha256, reset
`build.number` to 0) in a small PR.

## User-facing install (after acceptance)

```bash
conda install -c conda-forge -c bioconda shannonpore
# or
mamba install -c conda-forge -c bioconda shannonpore
```
