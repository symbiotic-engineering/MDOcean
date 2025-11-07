from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)

x.add_system("opt", OPT, r"\text{Optimizer}")
x.add_system("geom", FUNC, r"\text{Geometry}")
x.add_system("hydro", FUNC, r"\text{Hydrodynamics}")
x.add_system("solver", SOLVER, r"\text{Iteration}")
x.add_system("dynam", FUNC, (r"\text{Dynamics}", r"\text{and Control}"))
x.add_system("struct", FUNC, r"\text{Structures}")
x.add_system("econ", FUNC, r"\text{Economics}")
x.add_system("F", FUNC, r"\text{Objective}")
x.add_system("G", FUNC, r"\text{Constraints}")

x.connect("opt", "geom", (r"\text{Dimensions,}",r"\text{Thicknesses}"))
x.connect("opt", "hydro", r"\text{Dimensions}")
x.connect("opt", "dynam", (r"\text{Generator}",r"\text{ratings}"))
x.connect("opt", "struct", (r"\text{Dimensions,}", r"\text{thicknesses}"))
x.connect("opt", "econ", (r"\text{Generator}",r"\text{ratings}"))

x.connect("solver","dynam", (r"\text{Dynamic}",r"\text{response guess}"))
x.connect("dynam","solver", (r"\text{Dynamic}",r"\text{response residual}"))

x.connect("geom", "dynam", r"\text{Mass}")
x.connect("geom", "econ", r"\text{Material volume}")
x.connect("hydro", "dynam", (r"\text{Hydrodynamic}",r"\text{coefficients}"))
x.connect("dynam", "struct", r"\text{Loads}")
x.connect("dynam", "econ", r"\text{Power}")

x.connect("geom", "G", (r"\text{Stability and}",r"\text{hydrostatic constraints}"))
x.connect("dynam", "G", r"\text{Amplitude constraints}")
x.connect("struct", "G", r"\text{Structural constraints}")

x.connect("dynam","F", r"\text{Power}")
x.connect("econ","F", (r"\text{LCOE,}",r"\text{Capital cost}"))

x.connect("F", "opt", "J")
x.connect("G", "opt", "g")

x.add_output("opt", "x^*", side=LEFT)
x.add_output("F", "J^*", side=LEFT)
x.add_output("G", "g^*", side=LEFT)
x.write("mdocean/plots/non_matlab_figs/xdsm")
