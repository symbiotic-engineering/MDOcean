clear;close all;clc
p = parameters();
b = var_bounds(p);


X = [b.X_noms; 1];

[LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM, P_elec, D_d, ~, g] = simulation(X,p)

FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
[feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),g(17),g(18),p)

num_outputs = 9;
runtime = timeit(@()simulation(X,p),num_outputs);

plot_power_matrix(X,p)

figure
power_PDF(X,p)
