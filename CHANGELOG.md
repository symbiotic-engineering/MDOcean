
## Unreleased

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
