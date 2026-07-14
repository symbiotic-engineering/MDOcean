# MDOcean Copilot Agent Onboarding

## Repository quick orientation
- MDOcean is a MATLAB-first engineering optimization codebase with reproducible pipelines managed through **Calkit + DVC**, plus LaTeX paper builds in `pubs/`.
- Core areas:
  - `mdocean/`: simulation, optimization, analysis, plotting code.
  - `tests/`: MATLAB unit/integration-style validation tests (`run_tests.m`).
  - `pubs/`: LaTeX manuscripts (`applied-ocean-research-model`, `renewable-energy-mdo`, dissertation, etc.).
  - `docs/`: Sphinx documentation.
  - `calkit.yaml`: canonical pipeline stage definitions.

## Default way to run and validate work
- Prefer **targeted Calkit stages** over full-pipeline runs.
- Use exactly one stage with `-s` for focused checks unless explicitly asked otherwise.
- Use `calkit xr ...` for commands that should run in a Calkit-managed environment but do not have a dedicated stage.

Common targeted checks:
- AOR paper build: `calkit run -s build-AOR-paper`
- RE paper build: `calkit run -s build-RE-paper`
- AOR LaTeX lint: `calkit run -s paperlint-aor`
- RE LaTeX lint: `calkit run -s paperlint-re`
- AOR tests: `calkit run -s test-aor`
- RE tests: `calkit run -s test-re`
- Duplicate JSON key check for paper numbers: `calkit run -s check-dup-json-keys`
- Docs build (from `docs/`): `make html`

## Editing and contribution conventions
- MATLAB docstrings should use:
  ```matlab
  function [out1,out2] = func_name(in1,in2)
  % Function func_name
  %
  % :param in1: Input parameter 1
  % :param in2: Input parameter 2
  % :returns: Output parameter 1
  % :returns: Output parameter 2
  ```
- Update `CHANGELOG.md` under `## Unreleased` for completed tasks:
  ```md
  ### Added (or Changed, Removed, Fixed)
  - Category: Short description.
  ```
  Categories: Model, Validation, Speedup, CI, Paper, Figures, Tables, Pipeline, Analysis, Dev, Tests, Docs, Submodules.
- Merge-ready CI expects both:
  - `mdocean/__version__.py`
  - `CHANGELOG.md`
  to be updated in PRs intended for release/merge.

## Git and merge behavior
- If asked to merge, fetch/pull first and produce a true merge commit (two parents), not a squash-like normal commit.
- This repo uses a custom merge driver for `dvc.lock`:
  - `git config merge.dvclock.driver "python dev/dvc_lock_merge.py %O %A %B"`
- If instructed to use "ours/theirs" for pipeline-managed conflicts, this refers to the custom merge driver prompt keys (`o`/`t`) unless otherwise specified.
- When resolving changelog merge conflicts, keep both entries and place **ours** first.

## Known recurring errors and workarounds
- `calkit pull` transient cloud errors:
  - Errors like `failed to pull data from the cloud - <n> files failed to download` or `unexpected error` can be transient.
  - **Workaround:** re-run the same `calkit pull` command until it succeeds.
- `calkit pull` checkout warning:
  - `failed to pull data from the cloud - Checkout failed for following targets`
  - **Workaround:** this is documented as safe to ignore in maintainer guidance.
- Submodule checkout/auth errors in CI (`No url found` for nested submodule):
  - Caused by stale nested submodule metadata under `.git/modules`.
  - **Workaround pattern:** deinit submodules, remove stale module cache, checkout with `submodules: false`, then run explicit `git submodule update --init --recursive`.
- During automated `git pull`, `not our ref` can occur with submodule pointers:
  - **Workaround pattern in CI:** continue flow; only create fallback branch for other pull failures.

## Efficiency tips for cloud agents
- Start with `calkit.yaml`, `tests/run_tests.m`, and relevant paper folder before changing code.
- Avoid running `calkit run` without stage filters unless explicitly requested.
- Keep changes surgical and scoped to the touched subsystem (model, pipeline, papers, docs).