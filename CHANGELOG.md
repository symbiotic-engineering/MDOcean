# Changelog
## Unreleased
### Changed
- Dev: manage external MATLAB dependencies (OpenFLASH, WEC-Sim, SAFE) with the [mip](https://mip.sh) package manager instead of git submodules: dependencies are declared in `mip.yaml`, pinned in `mip.lock`, and installed into a project-local `./.mip` environment by `setup_mip` / `mip project sync`

## [v1.3.3](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.3.3) - 2026-7-14
### Changed
- Dev: expand `.github/copilot-instructions.md` with cloud-agent onboarding guidance, targeted Calkit commands, and known error/workaround notes.

## [v1.3.2](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.3.2) - 2026-7-14
### Added
- Dev: Create copilot agent instructions
- Dev: Customize copilot setup to include calkit and matlab

## [v1.3.1](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.3.1) - 2026-7-13
### Added
- Pipeline: feature in dvc_lock_merge.py to checkout git-tracked pipeline outputs consistent with dvc.lock merge

## [v1.3.0](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.3.0) - 2026-7-11
Papers at the state of journal submission, plus a few post-submission edits.
### Added
- Paper: custom table column types P, M, L
- Paper: plain language abstract in dissertation
- Figures: stacked bar charts for Comparison outputs
- Figures: create real graphical abstract for RE
- Figures: create `mod_freq_domain_ctrl_combined` mermaid flowchart
- Paper: macro to switch between one and two columns
- Paper: declaration of interests and declaration of generative AI
- AOR paper: describe spar coefficients Olaya fit
- AOR paper: use stars instead of low/med/high for method comparison table
- RE paper: add back design of experiments section
### Changed
- Figures: `DesignSpaceExploration` marks infeasibility with x/o
- Figures: add color to mermaid control flowcharts
- Figures/tables: format tweaks to `LocationSensitivity`, `FitOlaya`, `PtoSweep` fig/table outputs
- Figures: fix aspect ratio of `constraint_active_plot.m` 
- Analysis: `SweepGeoms` real timing as output, increase sweep fidelity from 3 to 5 points per dimension
- Pipeline: Tikz figures are standalone and part of `ReadNonMatlabFigs` instead of externalized
- Pipeline: remove pareto search dependence on constraint active figure
- Pipeline: update overleaf sync paths
- Paper: use `cas-dc` document class
- AOR paper: add overlines to hydrostatic stability line segments
- AOR paper: bring back definition of storm hydro coeffs
- AOR paper: change modified frequency domain to QLPS and describe dimension of QLPS problem in `module_details.tex`
- Dissertation: appendices occur after each chapter instead of all at the end
### Fixed
- Tables: fix extra math `$` that occurred in `table2latex` entries with subscripts
- Paper: use pipeline outputs instead of manual figures for power matrix multiplication
- Paper: appropriate `citet/citep`, capitalization of `Cref`
- Pipeline: removed unnecessary dependencies in optimization folder from `make-calkit-stages` stage
### Removed
- Pipeline: removed overleaf sync for beamer slides
- Pipeline: removed `sim_runtime.pdf` figure
- Paper: deleted `numbers.tex` in AOR
- Paper: commented out multistart bar chart, runtime bar chart, epsilon constraint seed figures and some discussion about sensitivity to drag coefficient in RE

## [v1.2.17](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.17) - 2026-7-8
### Changed
- License: updated Copyright to include UMich and all contributors

## [v1.2.16](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.16) - 2026-7-8
### Added
- Dev: create latex linter to check common mistakes
### Changed
- Paper: Maha's word choice changes from Overleaf - note that these are not yet fully reviewed but are merged anyway to maintain accurate history
- Pipeline: comment out svg2png calkit stage
- Deps: bump `OpenFLASH` submodule by 2 commits
- Pipeline: move `matlab_figs` plot functions into folders, regenerate calkit.yaml, and force commit dependent stages
- Paper: update UMERC paper to use shared citation references and cref for appendices
- Dev: Change `.mlx` to `.m` for file changed scraper
### Fixed
- Pipeline/paper: fixing citations, refs, json2latex, numbers, and addressing linter-identified issues in Maha's wording changes
- Figures: Improve readability of Mermaid control flowchart diagrams, PTO sweep (remove hashing), and `GradientOptim`'s `delta_x`

## [v1.2.15](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.15) - 2026-7-7
### Added
- Docs: instructions for using calkit with overleaf
- Dev: interactive script to help execute a calkit overleaf sync

## [v1.2.14](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.14) - 2026-7-7
### Added
- Dev: script to check convexity with cvxpy
- Dev: script to return which github copilot models were used from json log, for AI acknowledgement statement

## [v1.2.13](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.13) - 2026-7-7
### Fixed
- CI: fixed codecov upload by using OIDC instead of token

## [v1.2.12](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.12) - 2026-6-26
- Paper: portion of Maha's Overleaf changes that move things to appendix and introduce if/else structure, without changing wording.
- Submodule: update OpenFLASH to reflect JFM submission changes (`jfm-thesis` branch with `main` merged into it)

## [v1.2.11](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.11) - 2026-6-19
### Changed
- Tests: comment out assert success command in `run_tests.m`
### Fixed
- Tests: allow stages with tables to pass tests, just have to rerun their postpro

## [v1.2.10](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.10) - 2026-6-17
### Added
- Pipeline: OpenFLASH subproject
- Docs: add defense video to readme and sphinx site
- Presentation: upload pptx and pdf of dissertation slides
- Plots: `raster_text()` matlab utility required for displaying text with certain fonts
- Plots: `unicode()` utility for handling supplementary plane unicode codepoints
- Dev: `bibtex-conda` environment for manual latex citation checking
- Dissertation: create sections for contributions, list of research outputs, conclusion
- Dissertation: add quote to beginning of each chapter
- Calkit status: version sent to committee 6/9
### Changed
- Paper: Use subimport for Overleaf-friendly path inclusion
- CI: Initialize submodules recursively
- CI: Turn off quarto stage (since no longer used)
- CI: stop using retry script for calkit pull
- Docs: update repo status badge from in progress to active
- Docs: better details on licensing of external code
- Dissertation and paper: updates for readability and clarity
- Dissertation: use submitted JFM paper but with US spelling
- Tests: delete unused AI-generated tests
- Pipeline: use new ck DVC remote
- Submodule: move `sea-lab-utils` to `analysis` folder in `OpenFLASH`
- Paper: create line breaks for each sentence
### Fixed
- Analysis: comment out broken AI-generated `compute_pareto_shape_heuristics()` in `ParetoFigFunc`
- Test: fix wecsim filenames being inconsistent after previous shortening for struct fieldnames
- Test: further shorten wecsim filenames (multibody)
- Test: fix string vs char handling in `save_fig_with_diagnostic.m`
- Test: avoid rerunning postprocessing in tests (but this introduces test errors in stages that have tables) 
- Plots: fix filename in `ReadNonMatlabFigs` that caused raster instead of vector images to be used in paper
- Plot: fix multistart parallel axis plot visibility, including latex font
- Tables: fix first column formatting in matlab latex table generation
- Tables: introduce, and then fix, an issue with string escape parsing in `table2latex()`
- Tables: fix multi-character underscores in `Comparison`, `DesignVars`, and `GradientOptimFigFunc` tables
- Paper: fixing citations, including making bibcop compatible with multiple .bib files
- Paper: fix references, including prefacing labels with paper name to avoid collisions
- Citations: fix JFM citations, title of JFM paper in self-citation, and RE bibcop
- CI: ignore `not our ref` error in calkit save branching logic
- CI: turn on write permission for checks

## [v1.2.9](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.9) - 2026-5-7
### Added
- Paper: dissertation writeup
### Changed
- Submodules: update OpenFLASH submodule to include jfm paper in thesis
### Fixed
- Fix error in circle constraint center/radius for positve power
- Location sensitivity wrong class name in analysis outputs

## [v1.2.8](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.8) - 2026-5-5
### Added
- Pipeline: create `test-aor` and `test-re` stages
### Changed
- CI: combine calkit and test workflow, triggers on push to main
- Inputs: break out `make_drag_integral()` into its own function
- Pipeline: more string/char filename sanity checks in `save_fig_with_diagnostic()`
### Fixed
- Analysis: make `LocationSensitivity` run by adding forgotten analysis outputs, fixing figure concatenation,  
- Analysis: make `ParameterSensitivity` run by removing material column from linear constraint A-matrix
- Analysis: fix relative path issue in `GenericAnalysis`
- Tests: fixed `Wecsim` tests by ensuring filenames are within character limit for struct fields
- Tests: fixed relative path of `save_folder` in `all_figures()`

## [v1.2.7](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.7) - 2026-5-5
### Added
- Paper: discussion on scaling and convexity of power and LCOE with respect to dynamic constraint limits

## [v1.2.6](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.6) - 2026-5-4
## Added
- Add report wecsim validation figure to AOR paper appendix
## Changed
- Phase error validation plots now wrap from -pi to pi for better readability
- Wecsim validation plots use log scale x-axis for wecsim geometry, no drag, multibody
- Wecsim validation plots label colorbar axis with % symbol
- Wecsim visualization notebook only shows first page of pdf to reduce filesize
- Wecsim validation cases only turn amplitude saturation on for multibody cases
- Better aesthetics on signed log plot: tick label format, level selection, show contours
## Fixed
- CI: prevent dvc pull issues by retrying up to 10x

## [v1.2.5](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.5) - 2026-5-3
## Fixed
- setup: `add_mdocean_path.m` ensures local parallel pool is used to avoid interference between runners
- CI: use runner-specific temp directory and add cleanup function to avoid orphaned processes 

## [v1.2.4](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.4) - 2026-5-3
### Fixed
- Pipeline: add RunSingleFigFunc and Runtime as deps of `re-results-to-latex` so `nominalLCOE` and `simRuntime` format resolves
### Changed
- CI: use latest version of calkit again

## [v1.2.3](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.3) - 2026-4-29
## Added
- Slides: create first draft defense slides in both quarto and beamer
- Pipeline: stages for rendering quarto slides in html and pdf using quarto and playwright chromium docker containers respectively

## [v1.2.2](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.2) - 2026-4-29
### Changed
- Speedup: `find_nominal_inputs` uses unsaturated simulation to calculate `F_max_nom` instead of optimization
- Validation: suppress warnings in Meem validation
- Analysis: parallelize and remove print statements in `PtoSweep` analysis

## [v1.2.1](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.1) - 2026-4-25
### Added
- Paper: dedicated economics appendix, including beta table moved from main text
- Validation: average percent error rows for CAPEX, OPEX, and LCOE added to validation LaTeX table
- Wecsim: `wecsimAvgPowerErrorBestCase` and `wecsimAvgPowerErrorWorstCase` exported to `end.json`
- Paper: create flowcharts and explanations for modified-frequency-domain method
### Changed
- Validation: include units in the validation table
- Pipeline: Mermaid outputs switched from SVG to PNG and calkit iteration is used for multiple diagrams
- Paper: moved most of PTO section to appendix
- Paper: removed various less-important figures, combined related figures to subfigures, and updated AOR figure mapping
- Paper: various wording fixes in benchmarking.tex, discussion.tex, module-details.tex
- Plots: removed titles from `power_matrix_compare` contour and tiled-layout plots
### Fixed
- Analysis: `Parameters` class normalization now handles function handles
- Analysis: add forgotten intermediate outputs to `MEEM` class
- CI: more comprehensive fix to stale submodules: deinit all before clearing, checkout submodules false, and update submodules separately

## [v1.2.0](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.2.0) - 2026-4-25
Major new functionality: constrained optimal control via brute force search, solver, and analytical QP
### Added
- Model: solver-based controller can now obey amplitude and power constraints, not just force limits, via revamped `control_errors_from_sat_results()`
- Model: ability to determine control via nested brute force sweep instead of simultaneously with the drag solver
- Model: ability to determine control via quadratically-constrained quadradratic program (QCQP) implemented via analytical circle intersections
- Model: parameter `control_solve_type` to change betweeen solver, brute force, and analytical optimal control formulations
- Plots: four new plots in `RunSingleFigFunc` class for optimal control sweep (multiplier and reflection coefficient space) and QCQP (circles and circles with brute-force overlay in reflection coefficient space)
- Validation: create test case for verifying/benchmarking `multibody_response.m`
- Model: `control_type='none'` option for storm case
### Changed
- Speedup: change dynamics argument output order to avoid unnecessary recomputation of phase
- Speedup: reuse shared helper terms within and across controllers in `multibody_response.m`
- Speedup: reuse `D_sys` determinant from `multibody_impedance.m` in `multibody_response.m`
- Speedup: compute angles via `atan2()` rather than `angle()`, and use `complex()` or identities in place of some complex operations
- Speedup: `response_and_ctrl_err_from_ctrl_guess()` uses guessed phase to avoid computing resultant phase
- Analysis: runtime analysis bar chart includes combinations of various constraints on/off
- Analysis: change `GenericAnalysis` to handle class
- Analysis/Speedup: `GenericAnalysis` caches computation of dependencies between analysis and postpro stages
- Model: retune report power and heave force multiplier
- Plots: `improvePlot.m` now detects polar axes, not just Cartesian
- File structure: split control functions out into new `constrained_opt_ctrl.m` file
### Fixed:
- Model: output the stabilized PTO coefficients from `control_evaluation_fcn()` to avoid inconsistencies in downstream calculations
- Model: prevent negative PTO damping

## [v1.1.29](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.29) - 2026-4-24
### Changed
- Validation: plots for damping plate hydro fit now overlay case 2, 4, and WAMIT data in log space and with improved aesthetics

## [v1.1.28](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.28) - 2026-4-20
### Added
- Validation: plots showing acceleration harmonics, position THD, and drag force waveform to validate describing function

## [v1.1.27](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.27) - 2026-4-20
### Fixed
- Model: add pi to multibody phases to make wecsim phase validation match
- Pipeline: add forgotten meem intermediate results figure outputs
- Validation: show all phases between -pi and pi

## [v1.1.26](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.26) - 2026-4-20
### Added
- Simulation: create `run_and_catch_warnings.m` to capture output and count warnings
### Changed
- Model: include identifiers on warnings in `get_response_drag.m`, `make_drag_integral.m`, `run_MEEM.m`, and `simulation.m`
- Analysis: SweepGeoms analysis uses warning suppression and counting

## [v1.1.25](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.25) - 2026-4-17
### Fixed
- Validation: fix PTO force extraction from WEC-Sim so it includes stiffness (previously only included damping component)

## [v1.1.24](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.24) - 2026-4-17
### Fixed
- Model: remove sign discrepancy between MEEM and WAMIT excitation phase 

## [v1.1.23](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.23) - 2026-4-17
### Changed
- Analysis: add jupyter notebooks showing figures in each stage

## [v1.1.22](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.22) - 2026-4-15
### Changed
- Analysis: cache `var_bounds()` struct in `var_bounds.mat` to avoid recomputation

## [v1.1.21](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.21) - 2026-4-15
### Fixed
- Validation: extract one wave period of WEC-Sim output starting at a wave-period boundary so that FFT phase is relative to the wave (`runRM3Parallel.m`)

## [v1.1.20](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.20) - 2026-4-13
### Fixed
- Model: phase sign error in `second_order_transfer_fcn` in `get_response_drag.m`
- Validation: copy-paste error in `runRM3Parallel.m` where `float_pos` was used to compute spar and relative phase instead of `spar_pos` and `rel_pos`
- Validation: mislabeled singlebody/multibody rows in WEC-Sim error table (`post_process_fcn.m`)

## [v1.1.19](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.19) - 2026-4-13
### Added
- Analysis: new `FitOlaya` class that loads Olaya et al. digitized damping plate data, computes derived signals, and generates fit and exploratory plots for hydrodynamic coefficient fitting

## [v1.1.18](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.18) - 2026-4-13
### Changed
- CI: pin calkit to v0.37.3 to avoid corrupted mat files in 0.37.4

## [v1.1.17](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.17) - 2026-4-13
### Added
- Analysis: new `PtoSweep` class that performs a 2-D sweep of F_max and P_max for the nominal RM3 geometry, producing contourf plots of average power, design cost, and LCOE with a hatched infeasibility region and markers for the maximum-power and minimum-LCOE operating points

## [v1.1.16](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.16) - 2026-4-10
### Added
- Analysis: new class for sweeping geometries and plotting radiation eff, surface area, CWR, etc
- Model: cache m_k_h to reduce compute time since it doesn't change across sims
### Changed
- Model: turn on stabilization within optimal control solver
- Model: increase max kappa of drag LUT from 8 to 120
### Fixed
- Model: handle edge case geometries that cause finite precision errors
- Model: fix indexing error in make_drag_integral_LUT that led to out-of-bounds queries being kept NaN instead of numerically computed

## [v1.1.15](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.15) - 2026-4-6
### Changed
- CI: bump paths-filter to v4 in merge-ready.yml workflow

## [v1.1.14](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.14) - 2026-4-5
### Added
- Pipeline: stage that shows figures in jupyter notebook to facilitate code review
 
## [v1.1.13](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.13) - 2026-4-5
### Changed
- CI: bump `astral-sh/setup-uv` from v5 to v7

## [v1.1.12](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v1.1.12) - 2026-4-5
### Fixed
- Model: wave velocity depth evaluated at full draft `T` instead of half-draft `T/2` in `get_dynamic_coeffs.m`
- Paper: drag coefficient formula corrected to use velocity amplitude `|\hat{\dot{\xi}}|` instead of displacement `|\hat{\xi}}|`
- Paper: spar column ζ formula fraction corrected from `A·σ_Y/F_crit` to `σ_Y/(F_crit/A)`
- Paper: float waterplane area corrected to annular formula `π/4·(D_f²−D_s²)` instead of solid disk `π/4·D_f²`
- Paper: structural cost price `p_s` redefined as $/m³ (was $/kg) with table values updated accordingly
- CI: remove stale `sea-lab-utils` nested-submodule gitdir cache before checkout to fix `git submodule foreach --recursive` auth-removal failure on self-hosted runner

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
