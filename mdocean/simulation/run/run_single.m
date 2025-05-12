clear;close all;clc
p = parameters();
b = var_bounds();


X = [b.X_noms; 1];

[LCOE, P_var, ~, g, val] = simulation(X,p)

[feasible,~,failed] = is_feasible(g,X,p,b)

plot_power_matrix(X,p,b,'')

figure
power_PDF(X,p)

visualize_geometry(X,p)
