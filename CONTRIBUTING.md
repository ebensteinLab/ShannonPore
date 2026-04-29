# Contributing to nanoentropy

Thanks for your interest in contributing to **nanoentropy**! This document
describes how to set up a development environment, the standards we follow,
and how to submit changes.

## Getting started

1. **Clone the repository**

   ```bash
   git clone https://github.com/uribertocchitau/nanoentropy.git
   cd nanoentropy
   ```

2. **Install** the conda + Python + R environment:

   ```bash
   bash install.sh
   ```

3. **Verify the install:**

   ```bash
   nanoentropy doctor      # checks runtime dependencies
   nanoentropy selftest    # runs a small end-to-end smoke pipeline
   ```

If both commands report a clean status, you're ready to develop.

## Development workflow

We follow a standard fork-and-pull-request workflow:

1. **Fork** the repository on GitHub.
2. **Create a branch** off `main` using one of these prefixes:
   - `feat/<short-description>` — new feature
   - `fix/<short-description>` — bug fix
   - `docs/<short-description>` — documentation only
   - `refactor/<short-description>` — internal refactor, no behavior change
   - `test/<short-description>` — tests only
3. **Write tests first** (TDD) — see [Testing requirements](#testing-requirements).
4. **Implement** until tests pass.
5. **Run the full check suite locally** before pushing.
6. **Open a pull request** against `main` describing the change and
   linking any related issue.

## Code style

- **Formatter:** [`black`](https://github.com/psf/black) (line length 100,
  configured in `pyproject.toml`).
- **Linter:** [`ruff`](https://github.com/astral-sh/ruff) (rules `E`, `F`,
  `W`, `I`, `B`, `UP`, `SIM`, configured in `pyproject.toml`).
- **Type hints** are required on all new public functions and methods.
- **Docstrings** are required on all public APIs (modules, classes,
  functions). Use a consistent style (Google or NumPy) within a file.
- **File size:** keep individual source files under **800 lines**. Split by
  feature/domain when files grow large.
- **Function size:** target under 50 lines per function; extract helpers
  rather than nesting more than 4 levels deep.
- **Immutability:** prefer returning new objects over mutating inputs.

Run the formatter and linter before committing:

```bash
black src tests
ruff check src tests
```

## Testing requirements

Every new feature **must** include tests. Bug fixes **must** include a
regression test that fails before the fix and passes after.

- **Coverage targets:** 80%+ on `src/pipelines/`, `src/plots/`, `src/io/`.
- **Test runner:**

  ```bash
  pytest -v -p no:anyio
  ```

- **Markers** (declared in `pyproject.toml`):
  - `smoke` — fast smoke tests
  - `unit` — unit tests
  - `integration` — integration tests with fixtures
  - `e2e` — end-to-end Streamlit tests
  - `slow` — long-running tests

- **Coverage report:**

  ```bash
  pytest --cov=src --cov-report=term-missing
  ```

## Commit messages

We use [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>: <short description>

<optional body>
```

Types we use: `feat`, `fix`, `docs`, `test`, `refactor`, `chore`, `perf`,
`ci`.

Examples:

- `feat: add ternary entropy mode to pipeline`
- `fix: handle missing CpG context in modkit output`
- `docs: clarify install steps for macOS`
- `test: add regression test for segmentation edge case`

## Pull requests

A PR is ready to merge when **all** of the following hold:

- [ ] All tests pass: `pytest -v -p no:anyio`
- [ ] `nanoentropy doctor` runs clean
- [ ] `nanoentropy selftest` runs clean
- [ ] `ruff check src tests` is clean
- [ ] `black --check src tests` is clean
- [ ] No merge conflicts with `main`
- [ ] PR description explains the **why** (motivation), not just the
      **what** (diff)
- [ ] New / changed behavior is reflected in `CHANGELOG.md` under
      `[Unreleased]`

For larger changes, please open a **Discussion** or an **Issue** first to
agree on approach.

## Issues and questions

- **Bug reports** and **feature requests:** open a
  [GitHub Issue](https://github.com/uribertocchitau/nanoentropy/issues).
  Include `nanoentropy doctor` output, OS / Python / R versions, a minimal
  reproducer, and the expected vs. actual behavior.
- **Questions, ideas, "how do I…":** use
  [GitHub Discussions](https://github.com/uribertocchitau/nanoentropy/discussions)
  rather than the issue tracker.

## Code review

We use an internal **adversarial review pipeline** to catch regressions and
hallucinated behavior in PRs that touch core pipelines or plotting code.
You can run it locally over a branch:

```bash
bash clawteam/orchestrate_review.sh v4-rev 3 5
```

See [`clawteam/README.md`](clawteam/README.md) for what the five-agent
pipeline does (code QA, hallucination QA, brutal review, defender, fixer)
and how to interpret its output.

For routine PRs, a single human reviewer plus passing CI is sufficient.

## License

By contributing, you agree that your contributions will be licensed under
the project's [MIT License](LICENSE).
