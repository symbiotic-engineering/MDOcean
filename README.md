[![Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Unit Tests](https://github.com/symbiotic-engineering/MDOcean/actions/workflows/tests.yml/badge.svg)](https://github.com/symbiotic-engineering/MDOcean/actions/workflows/tests.yml)
[![Code coverage with codecov](https://codecov.io/gh/symbiotic-engineering/MDOcean/graph/badge.svg?token=PQNFQ72IC8)](https://codecov.io/gh/symbiotic-engineering/MDOcean)
[![GitHub](https://img.shields.io/github/license/symbiotic-engineering/MDOcean)](https://github.com/symbiotic-engineering/MDOcean/blob/main/LICENSE.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13997243.svg)](https://doi.org/10.5281/zenodo.13997243)
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

Dissertation (in prep): R. McCabe, "Leveraging Semi-Analytical Modeling, Multidisciplinary Design Optimization, 
        and System Value Metrics to Advance Wave Energy Converter Viability," 2026. [Link to draft dissertation manuscript](https://calkit.io/symbiotic-engineering/mdocean/publications?ref=thesis-update&path=pubs%2Fdissertation%2FsampleThesis.pdf)

A video recording of the dissertation defense is available below. [Link to slides](https://calkit.io/symbiotic-engineering/mdocean/publications?ref=thesis-update&path=pubs%2Fdefense%2Fppt%2FMcCabe_Dissertation_Defense.pdf).
[<img width="618" height="318" alt="image" src="https://github.com/user-attachments/assets/30ce6df6-a7f8-4ed6-bd16-48fd492e7036" />](https://vod.video.cornell.edu/media/t/1_d7sd33c2)

Model journal paper citation (in prep): R. McCabe, M. Dietrich, and M. N. Haji, “Development, Validation, and Benchmarking of a Multidisciplinary Semi-Analytical Model for Wave Energy Converters,” in preparation, 2026. [Link to draft paper manuscript](https://calkit.io/symbiotic-engineering/mdocean/publications#pubs/applied-ocean-research-model/main.pdf).

Optimization journal paper citation (in prep): R. McCabe, M. Dietrich, and M. N. Haji, “Leveraging Multidisciplinary Design Optimization to Advance Wave Energy Converter Viability,” in preparation, 2026. [Link to draft paper manuscript](https://calkit.io/symbiotic-engineering/mdocean/publications#pubs/renewable-energy-mdo/mdocean.pdf).

Conference paper citation: R. McCabe, O. Murphy, and M. N. Haji, “Multidisciplinary Optimization 
to Reduce Cost and Power Variation of a Wave Energy Converter,” 
*International Design Engineering Technical Conferences & Computers and 
Information in Engineering Conference*, St. Louis, MO, August 14-17, 2022.
[https://doi.org/10.1115/DETC2022-90227](https://doi.org/10.1115/DETC2022-90227).

A video recording of the conference presentation is available below.

[<img width="528" height="297" alt="image" src="https://github.com/user-attachments/assets/81934806-9870-44bf-96a6-da8eb279a491" />](https://www.youtube.com/watch?v=LjpfXvujUGY)


## Code Documentation

Documentation for the function API for this code is in progress at [this Sphinx site](https://symbiotic-engineering.github.io/MDOcean/).


## Installation

Clone the repository via Git.
```
git clone https://github.com/symbiotic-engineering/MDOcean.git
```

Then, in MATLAB, run the one-time setup from the repository root to install the
external MATLAB dependencies (OpenFLASH, WEC-Sim, and SAFE):
```
>> setup_mip
```
This uses the [mip](https://mip.sh) package manager to build a project-local
environment (`./.mip`) containing the exact package versions pinned in
`mip.lock`, so every machine gets an identical dependency set. The declared
dependencies live in `mip.yaml`; see the [mip projects
documentation](https://mip.sh/docs/projects) for details. After setup, run
MATLAB commands with the dependencies loaded via

```
>> mip project run '<script, function, or expression>'
```

which is also how the calkit pipeline stages run.

If you are unfamiliar with Git and want to download the code another way, click "Code > Download ZIP" to get a .zip file, or try the "Open in MATLAB Online" button above to use the MATLAB Online IDE instead.

## How to Use

Please see the documentation [here](https://symbiotic-engineering.github.io/MDOcean/usage.html) for instructions on how to run the simulation, optimization, or entire paper-reproduction pipeline.
The package [`calkit`](https://docs.calkit.org/) is used for single-button reproducibility.
Once the installation steps are complete, the command `calkit run --pull` will reproduce all artifacts (figures, tables, papers).
Additionally, [this page](https://symbiotic-engineering.github.io/MDOcean/modification.html) has guidance on how to modify the code to be used for other problems.


## Software Authors

- Rebecca McCabe, rgm222@cornell.edu (Project lead and point of contact, 2021-present) [@rebeccamccabe](https://github.com/rebeccamccabe)
- Madison Dietrich, mjd429@cornell.edu (Project contributor, 2023-present) [@MadisonDietrich](https://github.com/MadisonDietrich)
- Anthony Long, anl58@cornell.edu (Project contributor, 2025) [@anthnlong](https://github.com/anthnlong)
- Olivia Murphy, om66@cornell.edu (Project contributor, 2021-22) [@ommurphy](https://github.com/ommurphy)
- Iris Ren, zr92@cornell.edu (Project contributor, 2024) [@irin0012](https://github.com/irin0012)
- Maha Haji, maha@cornell.edu (Advisor, 2021-present) [@maha-haji](https://github.com/maha-haji)

## License

This project is released open-source under the MIT License. 

As exceptions, the following files have been adapted from other code and have different licenses, reproduced in the relevant folders:
- `mdocean/inputs/validation/WECSim/`: Apache 2.0
- `mdocean/plots/util/MatlabSVG/`: BSD-2-Clause
- `mdocean/plots/util/RBGCube_a.svg`: CC BY-SA 3.0
- `mdocean/plots/util/bluewhitered.m`: BSD-2-Clause
- `mdocean/plots/util/hatchfill2.m` and `hatchfill2_demo.m`: BSD-2-Clause
- `mdocean/plots/util/table2latex.m`: BSD-2-Clause
- `dev/latex/Paper-Linter/`: MIT License, different copyright holder
- `pubs/shared/els-cas/`: LaTeX Project Public License v1.3c
- `pubs/applied-ocean-research-model/stars-cc-by-sa.tex`: CC BY-SA 4.0

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
| [SAFE](https://github.com/SAFEtoolbox/SAFE-matlab) | Optional for sensitivity analysis |

The external codes (OpenFLASH, WEC-Sim, and SAFE) are declared in `mip.yaml`, pinned in `mip.lock`,
and installed by the [mip](https://mip.sh) package manager (run `setup_mip` once; see Installation).
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
