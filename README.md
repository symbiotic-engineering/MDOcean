[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Unit Tests](https://github.com/symbiotic-engineering/MDOcean/actions/workflows/tests.yml/badge.svg)](https://github.com/symbiotic-engineering/MDOcean/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/symbiotic-engineering/MDOcean/graph/badge.svg?token=PQNFQ72IC8)](https://codecov.io/gh/symbiotic-engineering/MDOcean)
[![GitHub](https://img.shields.io/github/license/symbiotic-engineering/MDOcean)](https://github.com/symbiotic-engineering/MDOcean/blob/main/LICENSE.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13997244.svg)](https://doi.org/10.5281/zenodo.13997244)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=symbiotic-engineering/MDOcean)
[![View on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/180694-mdocean)

# MDOcean
This is an open source codebase that uses Multidisciplinary Design Optimization (**MDO**) to optimize an **ocean** wave energy converter (WEC). 

More specifically, it is both a multidisciplinary model and a design optimization framework for two-body point absorber WECs.

MDOcean uses the SQP and pattern search algorithms to find the geometry, PTO, and structural design which optimizes the levelized cost of energy, capital cost, and average power of the Reference Model 3 (RM3) WEC, using a fast simplified frequency domain WEC model.

## Features

### Model features:
- semi-analytical hydrodynamic model, using the matched eigenfunction expansion method to model linear potential flow
- semi-analytical frequency-domain dynamic model, using describing functions to model drag and saturation nonlinearities
- analytical structural model, using tabulated solutions and the equivalent-thickness method to model stiffened plate ultimate and fatigue limits
- algebraic cost model, with costs scaling with PTO force and power and structural material volume
- algebraic geometry model, calculating common areas and volumes
- (in development) algebraic eco-cost model, with eco-costs scaling with structural material volume, hull area, and maintenance frequency
- (in development) integration with grid model, capturing energy market prices and grid-wide emissions

For 210 sea states, the model takes 39 ms to run, which is around a 5 order of magnitude improvement compared to the equivalent ~1 hour set of parallel WEC-Sim MCR simulations.

### Optimization features:
-  SQP (gradient-based) single-objective optimization
-  pattern search (gradient-free) and epsilon constraint SQP (gradient-based) multi-objective optimization
-  multi-start to see how starting point affects the optimal result
-  derivative-based local sensitivity analysis to see how parameter values affect the optimal result with very little additional computation
-  re-optimization-based local sensitivity analysis to see how parameter values affect the optimal result with higher accuracy

Running the single-objective optimization with typical parameters reduces LCOE by 57% compared to the standard RM3 design.

### Benchmarking and Validation features:
- Validation of power production and device amplitude against WEC-Sim
- Validation of hydrodynamic coefficents against WAMIT (for RM3 geometry) and existing matched eigenfunction expansion method results (for a benchmark geometry from Chau and Yeung 2012)
- Validation of structural model against FEA results for RM3 geometry from Reference Model report (Neary et. al 2012)
- Validation of economic and geometric outputs against results for RM3 geometry from Reference Model report (Neary et. al 2012)
- Validation that power production does not violate the theoretical radiation limit
- Unit tests, code coverage, and automatic report generation to monitor status of above checks on demand
- Continuous integration with GitHub Actions to monitor status of above checks on every push to GitHub

## Academic Context

The project is part of research in the [Symbiotic Engineering Analysis (SEA) Lab](https://sea.mae.cornell.edu/).

Model journal paper citation (in prep): R. McCabe, M. Dietrich, and M. N. Haji, “Development, Validation, and Benchmarking of a Multidisciplinary Semi-Analytical Model for Wave Energy Converters,” in preparation, 2025. [Link to draft paper manuscript](https://drive.google.com/file/d/1egObpiE9zCCdZb5DkxvAwoFdwA6j_a-9/view?usp=sharing).

Optimization journal paper citation (in prep): R. McCabe, M. Dietrich, and M. N. Haji, “Leveraging Multidisciplinary Design Optimization to Advance Wave Energy Converter Viability,” in preparation, 2026. [Link to draft paper manuscript](https://drive.google.com/file/d/1hL8iIyhYBUG73HbeGz_ji5bXIvY4WKt0/view?usp=sharing).

Conference paper citation: R. McCabe, O. Murphy, and M. N. Haji, “Multidisciplinary Optimization 
to Reduce Cost and Power Variation of a Wave Energy Converter,” 
*International Design Engineering Technical Conferences & Computers and 
Information in Engineering Conference*, St. Louis, MO, August 14-17, 2022.
[https://doi.org/10.1115/DETC2022-90227](https://doi.org/10.1115/DETC2022-90227).

A video recording of the conference presentation is available [here](https://www.youtube.com/watch?v=LjpfXvujUGY).

## Code Documentation

Documentation for the function API for this code is in progress at [this Sphinx site](https://symbiotic-engineering.github.io/MDOcean/).


## Installation

Clone the repository via Git. Use the `--recursive` flag to include submodules.
```
git clone --recursive https://github.com/symbiotic-engineering/MDOcean.git
```

If you are unfamiliar with Git, click "Code > Download ZIP" to get a .zip file, or try the "Open in MATLAB Online" button above to use the MATLAB Online IDE instead.

## File Structure and Usage
- `tests`: continuous integration tests for validation as well as generating a report to reproduce all figures and results. Running `run_tests.m` is the easiest way to reproduce all results. The `run_slow_tests` flag in `test.m` can be toggled to enable or disable the tests which take longer than a few minutes to run.
- `mdocean/inputs`: numerical inputs needed to run the optimiztion, simulation, and validation, including wave data, parameters, design variable bounds, and validation values.
- `mdocean/simulation`: the simulation that takes design variables and parameters as inputs and returns objective and constraint values as outputs, and its validation.
The script `run_single.m` is a good starting point if you want to run the simulation without optimizing.
- `mdocean/optimization`: scripts and functions to perform single objective and multi-objective optimization and sensitivities. Start with the script `gradient_optim.m`
if you want to run single objective optimization for each of the two objectives. For multi-objective, run `pareto_search.m` followed by `pareto_curve_heuristics.m`.
- `mdocean/plots`: helper functions to visualize outputs. Start with the script `all_figures.m` if you want to generate specific figures from the paper.
- `mdocean/analysis`: classes for running analyses of various sorts, typically wrappers around functions in the `optimization` and `plots` folders. Uses a special cached workflow defined in abstract class `GenericAnalysis.m` so that changes to cheap post-processing scripts (typically plots) do not require rerunning the expensive analysis (typically optimization).
- `dev`: miscellaneous scripts not core to the codebase that were used to inform the development of the simulation.

If you are running individual scipts/functions, you will need to `cd` to the `mdocean` folder and add all subfolders here to the matlab path. This is done automatically if you are running everything at once via `run_tests.m`.

## Reproducibility

This project uses [Calkit](https://docs.calkit.org/) pipelines to ensure all figures and publications are fully reproducible.
To reproduce, [download calkit](https://docs.calkit.org/installation/) and then enter the command `calkit run`. 
You can also view the project artifacts on [calkit.io](https://calkit.io/symbiotic-engineering/mdocean). Larger files like images and pdfs are stored in DVC rather than Git.

## Software Authors

- Rebecca McCabe, rgm222@cornell.edu (Project lead and point of contact, 2021-present) [@rebeccamccabe](https://github.com/rebeccamccabe)
- Madison Dietrich, mjd429@cornell.edu (Project contributor, 2023-present) [@MadisonDietrich](https://github.com/MadisonDietrich)
- Anthony Long, anl58@cornell.edu (Project contributor, 2025) [@anthnlong](https://github.com/anthnlong)
- Olivia Murphy, om66@cornell.edu (Project contributor, 2021-22) [@ommurphy](https://github.com/ommurphy)
- Iris Ren, zr92@cornell.edu (Project contributor, 2024) [@irin0012](https://github.com/irin0012)
- Maha Haji, maha@cornell.edu (Advisor, 2021-present) [@maha-haji](https://github.com/maha-haji)

## License

This project is released open-source under the MIT License. The validation folder contains code taken from NREL's WEC-Sim. 
The Apache 2.0 license for this open source WEC-Sim code is included.

## Dependencies

The following packages are used in this code:
| Package | Required? |
| ---    | ---      |
| MATLAB | Required for simulation |
| Statistics and Machine Learning Toolbox | Required for simulation |
| [OpenFLASH](https://github.com/symbiotic-engineering/OpenFLASH) | Required for simulation |
| Optimization Toolbox | Required for optimization |
| Global Optimization Toolbox | Required for optimization |
| Symbolic Math Toolbox | Optional for simulation code generation |
| Parallel Computing Toolbox | Optional for speedup |
| MATLAB Report Generator | Optional for WEC-Sim validation |
| Simulink | Optional for WEC-Sim validation | 
| Simscape | Optional for WEC-Sim validation | 
| Simscape Multibody | Optional for WEC-Sim validation |
| [WEC-Sim](https://github.com/WEC-Sim/WEC-Sim/) | Optional for WEC-Sim validation |
| [SAFE](https://github.com/rebeccamccabe/SAFE-matlab) | Optional for sensitivity analysis |

The external codes (OpenFLASH, WEC-Sim, and SAFE) are included as submodules for ease of use.
The code has been tested on MATLAB R2022a (Windows) and R2024b (Linux), and likely works on other versions and operating systems.

## Contributing

Suggestions, questions, bug reports, and contributions are welcome. Open an issue or pull request. To discuss the possibility of broader collaborations, please email rgm222@cornell.edu.

## Funding Acknowledgement

This material is based upon work supported by the 
National Science Foundation Graduate Research Fellowship under 
Grant No. DGE–2139899, and the Cornell Engineering Fellowship.
Any opinion, findings, and conclusions or recommendations 
expressed in this material are those of the authors(s) and do not 
necessarily reflect the views of the National Science Foundation.
