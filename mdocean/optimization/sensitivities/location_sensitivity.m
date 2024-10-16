function [] = location_sensitivity()

p = parameters();
b = var_bounds();

files = {'Humboldt_California_Wave Resource _SAM CSV.csv',...
    'PacWave-North_Oregon_Wave-Resource.csv',...
    'PacWave-South_Oregon_Wave-Resource.csv','WETS_Hawaii_Wave-Resource.csv'};

X_opts = zeros(8,length(files));
obj_opts = zeros(1,length(files));
flags = zeros(1,length(files));

for i=1:length(files)
    
    jpd=readmatrix(files{i}, 'Range', 'A3');
    p.JPD = jpd(2:end,2:end);
    p.Hs = jpd(2:end,1);
    p.T = jpd(1,2:end);
  
    X = b.X_start_struct;
    
    which_obj = 1; % only optimize LCOE
    [X_opts(:,i), obj_opts(i), flags(i)]  = gradient_optim(X,p,b,which_obj);
    
    plot_power_matrix(X_opts(:,i),p)
    figure(2)
    power_PDF(X_opts(:,i),p)
    hold on
end
figure(2)
locs = {'Humboldt, CA','PacWave North, OR', 'PacWave South, OR','WETS, Hawaii'};
legend(locs)
array2table([X_opts;obj_opts;flags],'VariableNames',locs,'RowNames',[b.var_names,{'LCOE','flag'}])
assert(all(flags>0))

%% try Hawaii with California design
X_cali = X_opts(:,1);

LCOE_hawaii_with_cali_design = simulation(X_cali,p)
pct_diff = (LCOE_hawaii_with_cali_design - obj_opts(4)) / obj_opts(4)

end