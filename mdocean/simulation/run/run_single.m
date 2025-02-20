clear;close all;clc
p = parameters();
b = var_bounds();


X = [b.X_noms; 1];

[LCOE, P_var, ~, g] = simulation(X,p)

[feasible,~,failed] = is_feasible(g,X,p,b)

num_outputs = 4;
runtime = timeit(@()simulation(X,p),num_outputs);

plot_power_matrix(X,p)

figure
power_PDF(X,p)
