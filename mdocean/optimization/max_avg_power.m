function [X_opt,val,flag] = max_avg_power(p,b)

% hack that makes cost constant, so minimizing LCOE is actually maximizing average power
p_mod = p;
p_mod.cost_perN_mult = 0;
p_mod.cost_perW_mult = 0;
p_mod.cost_perkg_mult = [0 0 0]; 


x0 = b.X_start_struct;

% run LCOE minimization (effectively power maximization due to hack above)
[X_opt,~,flag] = gradient_optim(x0,p_mod,b,1); 

% plug back into unmodified simulation with regular params to get val
[~, ~, ~, ~, val] = simulation(X_opt,p);

end