clear;close all;clc
p = parameters();
b = var_bounds(p);

file={'Humboldt_California_Wave Resource _SAM CSV.csv','PacWave-North_Oregon_Wave-Resource.csv','PacWave-South_Oregon_Wave-Resource.csv','WETS_Hawaii_Wave-Resource.csv'};
n=length(file);
for i=1:n
    
  jpd=readmatrix(file{i}, 'Range', 'A3');
  p.JPD = jpd(2:end,2:end);
  p.Hs = jpd(2:end,1);
  p.T = jpd(1,2:end);
  
X = [b.D_f_nom, b.D_s_ratio_nom, b.h_f_ratio_nom, b.T_s_ratio_nom, ...
    b.F_max_nom, b.D_int_nom, b.w_n_nom, b.M_nom];


[LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM, P_elec, D_d, ~, g] = simulation(X,p)

FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
[feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),p)

num_outputs = 9;
runtime = timeit(@()simulation(X,p),num_outputs);

plot_power_matrix(X,p)

power_PDF(X,p)

end 

