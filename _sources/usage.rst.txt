Usage
=====


- ``tests``: continuous integration tests for validation as well as generating a report to reproduce all figures and results. Running ``run_tests.m`` is the easiest way to reproduce all results. The ``run_slow_tests`` flag in ``test.m`` can be toggled to enable or disable the tests which take longer than a few minutes to run.
- ``mdocean/inputs``: numerical inputs needed to run the optimiztion, simulation, and validation, including wave data, parameters, design variable bounds, and validation values.
- ``mdocean/simulation``: the simulation that takes design variables and parameters as inputs and returns objective and constraint values as outputs, and its validation.
The script ``run_single.m`` is a good starting point if you want to run the simulation without optimizing.
- ``mdocean/optimization``: scripts and functions to perform single objective and multi-objective optimization and sensitivities. Start with the script ``gradient_optim.m``
if you want to run single objective optimization for each of the two objectives. For multi-objective, run ``pareto_search.m`` followed by ``pareto_curve_heuristics.m``.
- ``mdocean/plots``: helper functions to visualize outputs. Start with the script ``all_figures.m`` if you want to generate specific figures from the paper.
- ``mdocean/analysis``: classes for running analyses of various sorts, typically wrappers around functions in the ``optimization`` and ``plots`` folders. Uses a special cached workflow defined in abstract class ``GenericAnalysis.m`` so that changes to cheap post-processing scripts (typically plots) do not require rerunning the expensive analysis (typically optimization).
- ``dev``: miscellaneous scripts not core to the codebase that were used to inform the development of the simulation.

If you are running individual scipts/functions, you will need to ``cd`` to the ``mdocean`` folder and add all subfolders here to the matlab path. This is done automatically if you are running everything at once via ``run_tests.m``.

Reproducibility
---------------

This project uses ``calkit <https://docs.calkit.org/>``_ pipelines to ensure all figures and publications are fully reproducible.
To reproduce, `download calkit <https://docs.calkit.org/installation/>`_ and then enter the command ``calkit run``. 
You can also view the project artifacts on `calkit.io <https://calkit.io/symbiotic-engineering/mdocean>`_. Larger files like images and pdfs are stored in DVC rather than Git.