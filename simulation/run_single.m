clear;close all;clc
p = parameters();
b = var_bounds(p);

X = [b.D_sft_nom, b.D_i_ratio_nom, b.D_or_nom, b.M_nom, b.N_WEC_nom, b.D_int_nom, b.w_n_nom];
[LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM, P_elec] = simulation(X,p)

plot_power_matrix(X,p)