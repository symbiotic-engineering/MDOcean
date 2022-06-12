clear;close all;clc
p = parameters();
b = var_bounds(p);

files = {'Humboldt_California_Wave Resource _SAM CSV.csv',...
    'PacWave-North_Oregon_Wave-Resource.csv',...
    'PacWave-South_Oregon_Wave-Resource.csv','WETS_Hawaii_Wave-Resource.csv'};

for i=1:length(files)
    
    jpd=readmatrix(files{i}, 'Range', 'A3');
    p.JPD = jpd(2:end,2:end);
    p.Hs = jpd(2:end,1);
    p.T = jpd(1,2:end);
  
    X = [b.X_noms; 1];

    [LCOE, P_var, B, FOS1Y, FOS2Y, FOS3Y, FOS_buckling, GM, P_elec, D_d, ~, g] = simulation(X,p)
    
    FOS = min([FOS1Y FOS2Y FOS3Y FOS_buckling]);
    [feasible,failed] = is_feasible(B,FOS,GM,P_elec,D_d,g(16),g(17),g(18),p)
    
    plot_power_matrix(X,p)
    figure(2)
    power_PDF(X,p)
    hold on

end

figure(2)
legend('Humboldt, CA','PacWave North, OR', 'PacWave South, OR','WETS, Hawaii')

