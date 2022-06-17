% Generate all figures used in the report

% figure 1 - RM3 image
% Created in powerpoint

% table 1 - design variables table
p = parameters();
b = var_bounds(p);
DV_table = array2table([b.X_mins b.X_noms b.X_maxs], 'VariableNames',{'Mins','Noms','Maxs'})

% table 2 - constraints table

% table 3 - parameters table

% figure 2 - N2 diagram
% Created in powerpoint

% figure 3 - saturation time signal
% figure 4 - saturation alpha
sin_saturation_demo()

% figure 5 - JPD multiplication

% table 4 - validation

% table 5 - optimal DVs for 4 designs
compare()

% table 6 - optimal DVs for 4 locations
location_sensitivity()

% figure 6 - pareto front
% figure 7 - design heuristics
pareto_bruteforce()

% figure 8 - parameter sensitivities
param_sweep()

% figure 9 - overlaid geometry

% figure 10 - probability CDF