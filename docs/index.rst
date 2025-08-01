.. MDOcean documentation master file, created by
   sphinx-quickstart on Mon Jul 28 14:20:24 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MDOcean documentation
=====================

This is an open source codebase that uses Multidisciplinary Design Optimization (**MDO**) to optimize an **ocean** wave energy converter (WEC). 

More specifically, it uses the SQP and pattern search algorithms to find the geometry, PTO, and structural design which optimizes the levelized cost of energy, capital cost, and average power of the Reference Model 3 (RM3) WEC, using a fast simplified frequency domain WEC model.

For 210 sea states, the model takes 39 ms to run, which is around a 5 order of magnitude improvement compared to the equivalent ~1 hour set of parallel WEC-Sim MCR simulations.

**Context**

The project is part of research in the [Symbiotic Engineering Analysis (SEA) Lab](https://sea.mae.cornell.edu/).

Journal paper citation (in prep): R. McCabe, M. Dietrich, and M. N. Haji, “Leveraging Multidisciplinary Design Optimization and Semi-Analytical Modeling to Advance Wave Energy Converter Viability,” in preparation, 2025.

Conference paper citation: R. McCabe, O. Murphy, and M. N. Haji, “Multidisciplinary Optimization 
to Reduce Cost and Power Variation of a Wave Energy Converter,” 
*International Design Engineering Technical Conferences & Computers and 
Information in Engineering Conference*, St. Louis, MO, August 14-17, 2022.
[https://doi.org/10.1115/DETC2022-90227](https://doi.org/10.1115/DETC2022-90227).

A video recording of the conference presentation is available [here](https://www.youtube.com/watch?v=LjpfXvujUGY).

**Software Authors**
- Rebecca McCabe, rgm222@cornell.edu (Project lead and point of contact, 2021-present) @rebeccamccabe
- Madison Dietrich, mjd429@cornell.edu (Project contributor, 2023-present) @MadisonDietrich
- Olivia Murphy, om66@cornell.edu (Project contributor, 2021-22) @ommurphy
- Iris Ren, zr92@cornell.edu (Project contributor, 2024)
- Maha Haji, maha@cornell.edu (Advisor, 2021-present) @maha-haji

**License**

This project is released open-source under the MIT License. The validation folder contains code taken from NREL's WEC-Sim. 
The Apache 2.0 license for this open source WEC-Sim code is included.

**File Structure and Usage**

Clone the repository via Git. If you are unfamiliar with Git, click "Code > Download ZIP" to get a .zip file, or try the "Open in MATLAB Online" button above to use the MATLAB Online IDE instead.

- `tests`: continuous integration tests for validation as well as generating a report to reproduce all figures and results. Running `run_tests.m` is the easiest way to reproduce all results. The `run_slow_tests` flag in `test.m` can be toggled to enable or disable the tests which take longer than a few minutes to run.
- `mdocean/inputs`: numerical inputs needed to run the optimiztion, simulation, and validation, including wave data, parameters, design variable bounds, and validation values.
- `mdocean/simulation`: the simulation that takes design variables and parameters as inputs and returns objective and constraint values as outputs, and its validation.
The script `run_single.m` is a good starting point if you want to run the simulation without optimizing.
- `mdocean/optimization`: scripts and functions to perform single objective and multi-objective optimization and sensitivities. Start with the script `gradient_optim.m`
if you want to run single objective optimization for each of the two objectives. For multi-objective, run `pareto_search.m` followed by `pareto_curve_heuristics.m`.
- `mdocean/plots`: helper functions to visualize outputs. Start with the script `all_figures.m` if you want to generate specific figures from the paper.
- `dev`: miscellaneous scripts not core to the codebase that were used to inform the development of the simulation.

If you are running individual scipts/functions, you will need to `cd` to the `mdocean` folder and add all subfolders here to the matlab path. This is done automatically if you are running everything at once via `run_tests.m`.

**Dependencies**

The following packages are used in this code:
| Package | Required? |
| ---    | ---      |
| MATLAB | Required for simulation |
| Statistics and Machine Learning Toolbox | Required for simulation |
| Optimization Toolbox | Required for optimization |
| Global Optimization Toolbox | Required for optimization |
| Symbolic Math Toolbox | Optional for simulation code generation |
| Parallel Computing Toolbox | Optional for speedup |
| MATLAB Report Generator | Optional for WEC-Sim validation |
| Simulink | Optional for WEC-Sim validation | 
| Simscape | Optional for WEC-Sim validation | 
| Simscape Multibody | Optional for WEC-Sim validation |
| [WEC-Sim](https://github.com/WEC-Sim/WEC-Sim/) | Optional for WEC-Sim validation |

The code has been tested on R2022a (Windows) and R2024b (Linux), and likely works on other versions and operating systems.

**Contributing**

Suggestions, questions, bug reports, and contributions are welcome. Open an issue or pull request. To discuss the possibility of broader collaborations, please email rgm222@cornell.edu.

**Funding Acknowledgement**

This material is based upon work supported by the 
National Science Foundation Graduate Research Fellowship under 
Grant No. DGE–2139899, and the Cornell Engineering Fellowship.
Any opinion, findings, and conclusions or recommendations 
expressed in this material are those of the authors(s) and do not 
necessarily reflect the views of the National Science Foundation.


.. toctree::
   api

.. note::
   This project is under active development.
