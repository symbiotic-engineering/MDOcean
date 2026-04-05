# Changelog
## Unreleased

## [v1.1.11](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.11) - 2026-4-5
### Changed
- CI: pin pubs conda environment to specific package versions (python=3.14.3, libxml2=2.15.2, pip=26.0.1)
- CI: fix Dependabot conda manifest directory to `/pubs/` and add pip updater for `/docs`

## [v1.1.10](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.10) - 2026-4-4
### Added
- CI: Dependabot configuration for conda, GitHub Actions, and git submodules

## [v1.1.9](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.9) - 2026-4-4
### Changed
- CI: upgrade GitHub Actions to Node.js 24 compatible versions (checkout v6, upload-artifact v6, setup-matlab v3, run-command v3, codecov-action v5, action-junit-report v6, action-download-artifact v20, setup-python v6)

## [v1.1.8](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.8) - 2026-4-1
### Changed
- Readme: links now go to calkit draft publications rather than old google docs drafts

## [v1.1.7](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.7) - 2026-4-1
### Fixed
- Pipeline: remove extraneous calkit outputs that caused two postpro stages to fail
- Pipeline: add xmllint to pubs conda env to avoid issue with mermaid reformatting stage
- CI: avoid unnecessary double dvc pull

## [v1.1.6](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.6) - 2026-3-30
### Added
- Pipeline: python script to update calkit.yaml automatically with matlab auto deps
- CI: enforce calkit.yaml to be up to date with matlab auto deps

## [v1.1.5](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.5) - 2026-3-30
### Changed
- Plot: added frequency domain to dynamics runtime bar plot
### Fixed
- Analysis: add forgotten add_wecsim_path in runtime analysis class

## [v1.1.4](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.4) - 2026-3-29
### Changed
- Dev: extracted WEC-Sim path setup into standalone `add_wecsim_path` function, called only when needed rather than during general path initialization

## [v1.1.3](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.3) - 2026-3-29
### Added
- Model: drag computed via strip theory integral. Note that this increases runtime by over 40%.
- Pipeline: stages to precompute drag lookup table and analysis for drag integral plots
- Validation: WecSim saves fundamental of drag force and phases via FFT
- CI: check to avoid accidental submodule downgrades 
### Changed
- Validation: WecSim uses Morison drag rather than quadratic drag
- Pipeline: dvc.lock merge driver no longer requires unchanged dvc.yaml
- Plots: runtime bar chart numbers rounded to nearest ms
- Plots: end figures saved as .fig in addition to .pdf
- Paper: moved meem appendix back to this repo from OpenFLASH
### Fixed
- Model: sign error in phase of d'alembert force
- Pipeline: track text outputs in git for a few stages forgotten last PR

## [v1.1.2](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.2) - 2026-3-16
### Added
- Paper: mermaid diagram for modeling choice taxonomy
- Paper: explanation of power multiplication and drag integral implementation
- Paper: stage for building UMERC grid paper
### Changed
- Paper: introduction focuses on modeling not optimization
- Paper: updates to structures, design load cases, econ, fixed point iteration, design variable coupling sections
- Paper: removed damping vs reactive appendix
- Paper: moved force saturation commentary from introduction to discussion
- Pipeline: proper storage of text file outputs in git, and dvc images stored individually instead of by folder
- Paper: updated graphical abstract
### Fixed
- Paper: MEEM appendix plots no longer placeholders with improved aesthetics
- Paper: fixed citation/reference warnings

## [v1.1.1](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.1) - 2026-3-13
### Added
- Model: ability to consider float+spar as one body in storm case
### Fixed
- Model: Recommended change in control to ensure stability properly extended to singlebody case 

## [v1.1.0](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.0) - 2026-3-5
### Added
- Model: checks for open loop and closed loop instabilities and, when the system is unstable but stabilizable, calculates the required change to the controller to maintain stability
- Model: saturate coupling hydro coefficients to enforce positive definite radiation matrices
- Model: new slamming diameter constraints for float and spar
- Analysis: `contourx` function that avoids warnings for constant contours
- Dev: small script to facilitate retuning heave force and average power
- Dev: scraper that checks out various commits and runs MEEM to facilitate hydro debugging
- Dev: custom merge strategy for `dvc.lock`
### Changed
- Model: major change to slamming model - now considers both minimum and maximum slamming amplitudes and appropriately differentiates between large and small waves
- Model: `power_scale_multibody`, `m_scale`, and `F_heave_mult` values increased to align validation
- Model: refactor of functions within `get_response_drag`. The new organization will make the comparison between single-loop and nested-loop solver easier, and reduces the passing of many parameters between functions by using a `control_evaluation_fcn` anonymous function.
- Model: the dynamics code in a few places starts to differentiate between `Z_l` and `Z_p` (load on electrical vs mechanical side), but most code still neglects electrical dynamics
- Pipeline: separate `move-results` stage to allow running stages for only one paper not both
- CI: Calkit workflow only runs stages needed for AOR, not RE, paper
- CI: Calkit workflow upgrades to latest calkit version
- Analysis: removed pareto from design space exploration for time savings in AOR development
- Analysis: figure saving in postpro uses `exportgraphics` instead of `print` to avoid cutting off title/etc in large figures
- Analysis: figure saving in intermediate results uses `savefig` instead of storing in struct, with saved position if needed when figure is larger than the screen, to avoid figure size issues 
- Optimization/analysis: renaming of slamming constraints for clarity
### Fixed
- Model: error in spar dynamics where wrong draft was used for wamit hydro coefficient interpolation. Now uses T_s 29 m instead of 35 m.
- Analysis: `is_feasible` not respecting ignored constraints
- Analysis/pipeline: Wecsim stages appropriately split between simulation and figures and save all 100+ figure outputs
- Analysis: figure aesthetics for Runtime, Slamming, and Wecsim stages
- Analysis: fix `power_matrix_compare` output sizes for `report=true` in Wecsim validation plots
- Analysis: corrected edge case logic in figure validity for when object does not have a 'Type' property
- CI: Calkit workflow issue where a failed dvc pull would be falsely misinterpreted as a git merge conflict
- CI: Avoid dvc timeouts by lengthening time and retrying dvc push in workflow
- CI: Add `numeric-results.tex` as part of the `calkit save`
- Readme: Zenodo badge now always points to newest version instead of old version

## [v1.0.4](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.0.4) - 2026-3-1
### Added
- Pipeline: new stage for saving end results as a json and formatting it for use in latex
- Paper: updates numbers to read from this json
### Changed
- Analysis: runtimes for each analysis class are saved
### Fixed
- Pipeline: xdsm stage uses conda environment rather than venv to get around venv activate error

## [v1.0.3](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.0.3) - 2026-2-27
### Added
- Paper: new section on QP for optimal control, PTO multiport dynamic equations
- Pipeline: XDSM diagram is generated by a python pipeline step rather than manually
### Changed
- Paper: notation consistency and matrix equation of motion in dynamics
- Paper: rewrite of drag section with integral formulation
- Paper: rewrite of slamming section recognizing both min and max amplitudes
- Analysis: updated slamming plot code to be consistent with new math
### Fixed
- Validation: issue in add_mdocean_path where WecSim h5 read file was removed from path (cherry picked from revert-dynamics-solver branch)
- Paper: some undefined citations
- Analysis: issue in module runtime comparison where times from profiler were not accurate. Now uses timeit to scale instead.

## [v1.0.2](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.0.2) - 2026-1-7
### Added
- Added config for readthedocs
### Changed
- Move usage instructions from readme to own docs page
- Minor docs changes
- Release notes use info from CHANGELOG.md
### Fixed
- Problem in deployment of docs to gh-pages branch

## [v1.0.1](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.0.1) - 2026-1-6

### Added
- Two new CI workflows related to versioning and releasing
- __version__.py to store version information
### Changed
- Trigger behavior of existing CI workflows
### Fixed
- Incorrect variable naming in damping plate structures
- Calkit stage dependencies so postpro depends on analysis
- Reduce likelihood of ld.so error by loading patch earlier

## [v1.0.0](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.0.0) - 2025-12-19 - Journal
Major changes to structures, dynamics, optimization formulation, caching, testing, and reproducibility made over the last year in preparation for upcoming journal paper submissions.

### Model
#### Added
- Initial environmental and CEM modules
- Stiffeners considered in structures
#### Changed
- Force saturation uses iteration rather than explicit quadratic formula
- Significant damping plate and float structural model updates
- Ability to perform iteration with fsolve solver
#### Fixed
- Multibody dynamics correction
- Reactive control now runs rather than failing to converge

### Optimization and Analysis 
#### Added
- Local parameter sensitivity
- Amplitude constraints
#### Changed
- Structural thicknesses added to design variables
- Geometric design variables are dimensional rather than ratios
- Objectives are output as vector J rather than individual variables
- Perform post processing without repeating analysis using the caching in GenericAnalysis class
- Pareto front axes are power and cost rather than LCOE and coefficient of variation

### Tests, Validation, and Supporting Code
#### Added
- Dynamics validation in WEC-Sim with reports for single and multi body
- Publications are part of the codebase (pubs folder)
- Testcase features: generated files get deleted at end, option to not run slow tests
- Full reproducibility with calkit pipeline
- Sphinx documentation
#### Changed
- Test case regressions show on GitHub actions GUI
- Speed up CI: parallelization and get rid of repeats

## [v2.1](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v2.1) - 2024-10-26 - Zenodo fix
### Fixed
- Citation.cff file valid
### Changed
- Upgraded to matlab actions tests v2

## [v2.0](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v2.0) - 2024-05-23 - first stable version using MEEM hydrodynamics
### Added
- Use MEEM instead of analytical Froude Krylov approximations to get hydro coeffs
- Unit tests, continuous integration, and code coverage
- Maximum capture width check
### Changed
- Updated folder structure
### Fixed
- Opex scaled incorrectly with number of WECs
- Powertrain wasn't included in available volume calculation
- Incorrect float waterplane area

## [v1.4](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.4) - 2022-08-16 - used for conference presentation
### Fixed
- Incorrect buckling end condition
- Incorrect implementation of force saturation multiplier
- Incorrect scaling of cost with number of WECs
- Mass of damping plate support tubes not accounted for
- Powertrain and transmission loss not accounted for
- Wave period and wave height were swapped
- Incorrect nominal controller design variables and objective values
- Nominal simulation did not use power saturation
- Incorrect scale factor when summing power over JPD
- Froude Krylov force coefficient too high
- Float maximum displacement constraint not being enforced

## [v1.3](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.3) - 2022-06-13 - used for the conference video
### Added
- Location sensitivity
- Check for which constraints are active
### Fixed
- Depth of submergence
- Parameter sweep slope calculation failed for NaNs

## [v1.2](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2) - 2022-05-17 - used for conference proceedings paper
### Added
- Comparison of power matrix for 4 designs
- New design variables: h_f_ratio, T_f_ratio, T_s_ratio
- Constrain that float cannot exceed spar height
### Fixed
- Incorrect value for young's modulus of steel

## [v1.0](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.0) - 2022-02-15 - used for conference draft submission
