from pyxdsm.XDSM import XDSM, OPT, SOLVER, FUNC, LEFT

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)

x.add_system("opt", OPT, r"\text{Optimizer}")
x.add_system("geom", FUNC, r"\text{Geometry}")
x.add_system("hydro", FUNC, r"\text{Hydrodynamics}")
x.add_system("solver", SOLVER, r"\text{Iteration}")
x.add_system("dynam", FUNC, (r"\text{Dynamics}", r"\text{and Control}"))
x.add_system("struct", FUNC, r"\text{Structures}")
x.add_system("cost", FUNC, r"\text{Cost}")
x.add_system("grid", FUNC, r"\text{Grid}")
x.add_system("env", FUNC, r"\text{Environment}")
x.add_system("F", FUNC, r"\text{Objective}")
x.add_system("G", FUNC, r"\text{Constraints}")

x.connect("opt", "geom", (r"\text{Dimensions,}",r"\text{Thicknesses}"))
x.connect("opt", "hydro", r"\text{Dimensions}")
x.connect("opt", "dynam", (r"\text{Generator}",r"\text{ratings}"))
x.connect("opt", "struct", (r"\text{Dimensions,}", r"\text{thicknesses}"))
x.connect("opt", "cost", (r"\text{Generator}",r"\text{ratings}"))

x.connect("solver","dynam", (r"\text{Dynamic}",r"\text{response guess}"))
x.connect("dynam","solver", (r"\text{Dynamic}",r"\text{response residual}"))

x.connect("geom", "dynam", r"\text{Mass}")
x.connect("geom", "cost", r"\text{Material volume}")
x.connect("geom", "env", r"\text{Material volume}")
x.connect("hydro", "dynam", (r"\text{Hydrodynamic}",r"\text{coefficients}"))
x.connect("dynam", "struct", r"\text{Loads}")
x.connect("dynam", "cost", r"\text{Power}")
x.connect("dynam", "env", r"\text{Power}")
x.connect("dynam", "grid", r"\zeta,\omega_n")
x.connect("cost", "grid", r"\text{LCOE}")
x.connect("grid", "env", r"\text{Grid Emissions}")

x.connect("geom", "G", (r"\text{Stability and}",r"\text{hydrostatic constraints}"))
x.connect("dynam", "G", r"\text{Amplitude constraints}")
x.connect("struct", "G", r"\text{Structural constraints}")

x.connect("grid","F", r"\text{Net Grid Value}")
x.connect("env","F", r"\text{Net Eco-Value}")

x.connect("F", "opt", "J")
x.connect("G", "opt", "g")

x.add_output("opt", "x^*", side=LEFT)
x.add_output("F", "J^*", side=LEFT)
x.add_output("G", "g^*", side=LEFT)
x.write("mdocean/plots/non_matlab_figs/xdsm_grid")
