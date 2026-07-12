Modification
============

In addition to running MDOcean as-is, we welcome and encourage users to create a fork and modify the code to suit their needs.

- To change **parameter values** including simulation/optimization settings:

  - For a one-off change, simply call ``p.param_name = new_param_value`` before running the simulation/optimization.
  - To change the default everywhere, edit the values in ``mdocean/inputs/parameters.m``. 

  If the parameter you are changing is the JPD, additional JPD data should be uploaded to ``mdocean/inputs/wave_conditions/``.
- To add a new **reference design**, which is a set of parameters and design variable values, define a new "mode" in 
  ``mdocean/inputs/parameters.m`` and ``mdocean/inputs/var_bounds.m`` following the structure provided. 
  The current reference designs are different versions of RM3: 'wecsim' and 'report'. 
  You may also want to change ``mdocean/optimization/find_nominal_inputs.m`` 
  (which currently sets the force limit just high enough to avoid saturation, since a force limit is not specified for the nominal RM3) 
  or disable it from being called in ``mdocean/optimization/var_bounds.m``.
  If you have validation data available for your reference design, put this in ``mdocean/inputs/validation_inputs.m``.
- To adjust the **design variable bounds** (upper and lower), modify the values in ``mdocean/inputs/var_bounds.m``.
- To change the **problem formulation**, such as to add/remove or change the definition of a design variable (x), constraint (g), or objective (J), 
  modify the names in ``mdocean/inputs/var_bounds.m`` and the code in ``mdocean/simulation/simulation.m``.
- To modify the **model** used for simulation, adjust the code in ``mdocean/simulation/modules``. 
  If you need to add inputs or outputs to the modules, this will also require adjusting ``mdocean/simulation/simulation.m``.

  - If you want to **override the hydrodynamics** model, for example to use hydro data you have previously obtained from a BEM solver,
    set ``p.use_MEEM = false`` and ``p.hydro = hydro_struct`` where ``hydro_struct`` is a struct with the hydro coeffs vs frequency 
    in the same format as WecSim expects. Note that this override only makes sense to use for simulation, not optimization, because the 
    overridden hydro data will not change with device geometry. The override could be used for optimization if the only design variables 
    are PTO or structural rather than bulk geometry.

  Note that the MEEM and dynamics modules use code that is auto-generated from MATLAB symbolic. 

  - The MEEM symbolic code is generated on-demand if the generated filename is not found. If you change any of the symbolic code in
    ``mdocean/simulation/modules/MEEM/run_MEEM.m`` (not recommended unless you really know what you're doing), you should delete the files in 
    ``mdocean/simulation/modules/MEEM/generated`` so they auto-regenerate upon running the simulation. Regeneration may take several minutes. 
  - On the other hand, the dynamics generation is currently done manually by running ``dev/dynamics/multibody/two_body.mlx``.
    For example, if you wanted to change the hydrodynamic degrees of freedom for the dynamics model, you would modify and rerun this mlx file.
    But conveniently, the majority of dynamics model changes can be done by modifying the regular (not auto-generated) files in the
    ``mdocean/simulation/modules/dynamics/`` folder without touching any of the auto-generated files.

- To modify the **optimization or sensitivity analysis methods**, for example to try out a genetic algorithm or implement Sobol analysis, 
  adjust the code in ``mdocean/optimization``. This code should be kept general enough to handle changes to the simulation and problem formulation,
  so avoid hardcoding.

If you encounter problems with any of these modifications, or if you want to modify the code in a way that is not covered by these instructions,
please submit a GitHub issue and the developers are happy to help.
We are also eager to see how you are using MDOcean even if you don't need assistance.