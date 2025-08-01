clear;close all;clc
p = parameters();
b = var_bounds();


X = [b.X_noms; 1];

[J, ~, g, val] = simulation(X,p)

LCOE = J(1);
P_var = J(2);

[feasible,~,failed] = is_feasible(g,X,p,b)

plot_power_matrix(X,p,b,'')

figure
power_PDF(X,p)

visualize_geometry(X,p)
