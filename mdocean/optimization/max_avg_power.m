function [X_opt_2,val_2,flag_2] = max_avg_power(p,b)

% hack that makes cost constant, so minimizing LCOE is actually maximizing average power
p_mod = p;
p_mod.cost_perN_mult = 0;
p_mod.cost_perW_mult = 0;
p_mod.cost_perkg_mult = [0 0 0]; 


x0 = b.X_start_struct;

% run LCOE minimization (effectively power maximization due to hack above)
[X_opt,~,flag] = gradient_optim(x0,p_mod,b,1); 

% plug back into unmodified simulation with regular params to get average power
[~, ~, ~, ~, val] = simulation(X_opt,p);

% since some DVs don't affect power, they are undetermined and we need
% another optimization to minimize cost while maintaining the same power
p_mod_2 = p;
p_mod_2.avg_power_min = val.power_avg;
x0_2 = cell2struct(num2cell(X_opt(1:end-1)),b.var_names(1:end-1)',1);
[X_opt_2,~,flag_2] = gradient_optim(x0_2,p_mod_2,b,2); 

% plug back into unmodified simulation with regular params to get val
[~, ~, ~, ~, val_2] = simulation(X_opt_2,p);

end