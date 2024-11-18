[![Unit Tests](https://github.com/symbiotic-engineering/MDOcean/actions/workflows/tests.yml/badge.svg)](https://github.com/symbiotic-engineering/MDOcean/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/symbiotic-engineering/MDOcean/graph/badge.svg?token=PQNFQ72IC8)](https://codecov.io/gh/symbiotic-engineering/MDOcean)
![GitHub](https://img.shields.io/github/license/symbiotic-engineering/MDOcean)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13997244.svg)](https://doi.org/10.5281/zenodo.13997244)

# MDOcean
This is an open source codebase that uses Multidisciplinary Design Optimization (**MDO**) to optimize an **ocean** wave energy converter (WEC). 

More specifically, it uses the SQP and pattern search algorithms to find the geometry and controller design which minimizes the energy cost and power variation 
of the Reference Model 3 (RM3) WEC, using a fast simplified frequency domain WEC model.

**Context**

The project is part of research in the [Symbiotic Engineering Analysis (SEA) Lab](https://sea.mae.cornell.edu/) 
and has been accepted for publication/presentation in the 2022 ASME IDETC-CIE.
At this conference, the work was presented at the DAC-6 session and is publication number 90227.
A recording of the conference presentation is available [here](https://www.youtube.com/watch?v=LjpfXvujUGY).
The project began as an effort in Cornell course [MAE 5350](https://classes.cornell.edu/browse/roster/FA21/class/MAE/5350).
Known areas for improvement are listed as GitHub issues. If you find any additional errors, please let us know.

Citation: R. McCabe, O. Murphy, and M. N. Haji, “Multidisciplinary Optimization 
to Reduce Cost and Power Variation of a Wave Energy Converter,” 
*International Design Engineering Technical Conferences & Computers and 
Information in Engineering Conference*, St. Louis, MO, August 14-17, 2022.
https://doi.org/10.1115/DETC2022-90227.

**Authors**
- Rebecca McCabe, rgm222@cornell.edu (Project lead and point of contact) @rebeccamccabe
- Olivia Murphy, om66@cornell.edu (Project contributor) @ommurphy
- Maha Haji, maha@cornell.edu (Advisor) @maha-haji

**License**

This project is released open-source under the MIT License. The validation folder contains code taken from NREL's WEC-Sim. 
The Apache 2.0 license for this open source WEC-Sim code is included.

**File Structure**

- `inputs`: numerical inputs needed to run the optimiztion, simulation, and validation, including wave data, parameters, design variable bounds, and validation values.
- `simulation`: the simulation that takes design variables and parameters as inputs and returns objective and constraint values as outputs, and its validation.
The script `run_single.m` is a good starting point if you want to run the simulation without optimizing.
- `optimization`: scripts and functions to perform single objective and multi-objective optimization and sensitivities. Start with the script `gradient_optim.m`
if you want to run single objective optimization for each of the two objectives.
- `plots`: helper functions to visualize outputs. Start with the script `all_figures.m` if you want to try out the entire pipeline by running all relevant 
optimizations to generate every figure in the paper.
- `dev`: miscellaneous scripts not core to the codebase that were used to inform the development of the simulation.

**Dependencies**

The following packages are used in this code:
- MATLAB
- Optimization Toolbox
- Global Optimization Toolbox
- Statistics and Machine Learning Toolbox
- Symbolic Math Toolbox

All are required except the symbolic math toolbox, which is used only for code generation and exploratory scripting, not core functionality.

**Funding Acknowledgement**

This material is based upon work supported by the 
National Science Foundation Graduate Research Fellowship under 
Grant No. DGE–2139899, and the Cornell Engineering Fellowship.
Any opinion, findings, and conclusions or recommendations 
expressed in this material are those of the authors(s) and do not 
necessarily reflect the views of the National Science Foundation.
