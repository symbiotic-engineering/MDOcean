function [X_opt,val,flag] = max_avg_power(p,b)

b.X_start_struct.P_max = 2000;
b.X_start_struct.F_max = 10;

[X_opt,val,flag,num_occurrences_max_pwr] = max_pwr_then_min_cost(p, b)


% if power has been saturated, meaning the cost-minimization found a better design than the power maximization
iters = 1;
max_iters = 4;
while num_occurrences_max_pwr > 1 && iters < max_iters

    warning('Power maximization failed on try %d of %d, trying again. Average power is %.2f kW.',iters,max_iters, val.power_avg/1e3)

    new_P_max = max(val.P_mech(:)) * p.eff_pto / 1e3;
    X_opt(strcmp(b.var_names,'P_max')) = new_P_max;

    fprintf('Setting new power limit to %d kW.',new_P_max)

    b.X_start_struct = cell2struct(num2cell(X_opt),b.var_names',1);
    [X_opt,val,flag,num_occurrences_max_pwr] = max_pwr_then_min_cost(p, b)

    iters = iters+1;

end

if num_occurrences_max_pwr>1
    error('Power maximization failed after the maximum number of iterations (%d).',max_iters)
end

end

function [X_opt_2,val_2,flag_2,num_occurrences_max_pwr] = max_pwr_then_min_cost(p, b)
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
    [~, ~, P_matrix, ~, val_2] = simulation(X_opt_2,p);

    num_occurrences_max_pwr = sum(P_matrix(:)==max(P_matrix,[],'all'));
end
