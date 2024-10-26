
## Unreleased

## [v2.1](https://github.com/symbiotic-engineering/MDOcean/releases/tag/v2.1) - 2024-10-26 - Zenodo fix
### Fixed
- Citation.cff works in Zenodo

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
