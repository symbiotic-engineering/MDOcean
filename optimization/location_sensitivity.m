clear;%close all;clc
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
  
    X = b.X_start_struct;

    which_obj = 1; % only optimize LCOE
    [X_opts, obj_opts, flags]  = gradient_optim(X,p,b,which_obj);
    
    plot_power_matrix(X_opts,p)
    figure(2)
    power_PDF(X_opts,p)
    hold on
end
figure(2)
legend('Humboldt, CA','PacWave North, OR', 'PacWave South, OR','WETS, Hawaii')
