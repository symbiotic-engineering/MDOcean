# MDOcean
This is an open source codebase that uses Multidisciplinary Design Optimization (**MDO**) to optimize an **ocean** wave energy converter (WEC). 
More specifically, it uses the SQP and pattern search algorithms to find the geometry and controller design which minimizes the energy cost and power variation 
of the Reference Model 3 (RM3) WEC, using a fast simplified frequency domain WEC model.

**Context**

The project is part of research in the [Symbiotic Engineering Analysis (SEA) Lab](https://sea.mae.cornell.edu/) and has been accepted for publication/presentation 
in the [2022 ASME IDETC-CIE](https://event.asme.org/IDETC-CIE).
At this conference, the work was presented at the [DAC-6](https://www.designautomationconference.org/dac-6) session and is publication number 90227.
The project began as an effort in Cornell course [MAE 5350](https://classes.cornell.edu/browse/roster/FA21/class/MAE/5350).

**Authors**
- Rebecca McCabe, rgm222@cornell.edu (Project lead and point of contact) @rebeccamccabe
- Olivia Murphy, om66@cornell.edu (Project contributor) @ommurphy
- Maha Haji, maha@cornell.edu (Advisor) @maha-haji

**Disclaimer**

The version of the simulation used in the conference proceedings paper and in the conference video was later found to have a number of errors. 
These errors have since been corrected and the current code is correct to the best of the authors' knowledge, within the limitations of the stated assumptions. 
Known areas for improvement are listed as GitHub issues. If you find any additional errors, please let us know.

**License**

This project is released open-source under the MIT License. The validation folder contains code taken from NREL's BEMIO module, which is part of WEC-Sim. 
The license for this open source WEC-Sim code is included.

**File Structure**

- `inputs`: numerical inputs needed to run the optimiztion, simulation, and validation, including wave data, parameters, design variable bounds, and validation values.
- `optimization`: scripts and functions to perform single objective and multi-objective optimization and sensitivities. Start with the script `gradient_optim.m`
if you want to run single objective optimization for each of the two objectives.
- `plots`: helper functions to visualize outputs. Start with the script `all_figures.m` if you want to try out the entire pipeline by running all relevant 
optimizations to generate every figure in the paper.
- `simulation`: the simulation that takes design variables and parameters as inputs and returns objective and constraint values as outputs, and its validation.
The script `run_single.m` is a good starting point if you want to run the simulation without optimizing.
- `dev`: miscellaneous scripts not core to the codebase that were used to inform the development of the simulation.
