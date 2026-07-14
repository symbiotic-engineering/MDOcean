- This repository performs a multidisciplinary engineering optimization study in matlab, and generates figures and compiles papers in latex.
- When testing code, run the code using calkit to ensure that the correct environments are applied.
  For example, to test latex code in the pubs/applied-ocean-research-model directory, run the following command:
  ```
  calkit run -s build-aor-paper
  ```
  Only run one pipeline stage (`-s` flag) rather than the entire pipeline, as the full pipeline can take a long time to run.
  If you are testing code that does not have a corresponding calkit stage, but would run in a calkit environment, use the `calkit xr` command.
- Write matlab docstrings in the following format:
  ```
  function [out1,out2] = func_name(in1,in2)
  % Function func_name
  %
  % :param in1: Input parameter 1
  % :param in2: Input parameter 2
  % :returns: Output parameter 1
  % :returns: Output parameter 2
  ```
- When you complete a task, add a line to the `CHANGELOG.md` under the `## Unreleased` section, or update the existing entry. Use the following format:
  ```
  ### Added (or changed, removed, fixed)
  - Category: Short description of the change.
  ```
  Available change categories are Model, Validation, Speedup, CI, Paper, Figures, Tables, Pipeline, Analysis, Dev, Tests, Docs, Submodules.
- If I ask you to merge, you should pull/fetch first to ensure the most updated version.
  Ensure it is a proper merge commit with two parents, not a regular commit.
  If I say to use ours or theirs to resolve merge conflicts, this refers to the keyboard input
  of `o` or `t` required by the custom merge driver for pipeline-managed files, 
  unless otherwise specified. When merging the changelog, you should combine the entries,
  putting the entry from ours at the top.