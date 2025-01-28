function [X_opt,val] = max_avg_power(p,b)

% hack that makes cost constant, so minimizing LCOE is actually maximizing average power
p.cost_perN_mult = 0;
p.cost_perW_mult = 0;
p.cost_perkg_mult = [0 0 0]; 


x0 = b.X_start_struct;

% run LCOE minimization (effectively power maximization due to hack above)
X_opt = gradient_optim(x0,p,b,1); 

% plug back into simulation to get unsaturated power
[~, ~, ~, ~, val] = simulation(X_opt,p);

end