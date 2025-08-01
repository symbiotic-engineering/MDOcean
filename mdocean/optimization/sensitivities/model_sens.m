clc
p = parameters();
b = var_bounds();
x0_input = b.X_start_struct;
%[X_opt, ~, ~, ~, ~, ~, ~, val] = gradient_optim(x0_input,p,b,1);
X_opt_mod = X_opt; 

allow_constr_viol = false; % true will run simulation, false will run optimization

if allow_constr_viol
    % force saturation off
    % find this number manually by adding breakpoint
    % in dynamics and doing max(mag_U_unsat,[],'all')/1e6 
    F_unsat = 11.815;
    X_opt_mod(6) = F_unsat;
    val_mod = print_saturation_effect(X_opt_mod,p,b,val);
    Pavg_factor_force_sat = min(val.P_sat_ratio,[],'all') * 100
    Fmax_factor_force_sat = X_opt(6)/F_unsat * 100
    
    % power saturation off
    P_unsat = max(val.P_mech * p.eff_pto,[],'all');
    X_opt_mod = X_opt;
    X_opt_mod(7) = P_unsat;
    print_saturation_effect(X_opt_mod,p,b,val);
    Pavg_factor_power_sat = val.power_max / P_unsat * 100
    
    % both force and power saturation off
    P_unsat_both = max(val_mod.P_mech * p.eff_pto,[],'all');
    X_opt_mod(6) = F_unsat;
    X_opt_mod(7) = P_unsat_both;
    print_saturation_effect(X_opt_mod,p,b,val);
    Pavg_factor_both_sat = val.power_max / P_unsat_both * 100
else
    x0_new = cell2struct(num2cell(X_opt(1:end-1)),b.var_names(1:end-1)',1);

    % force and power saturation both on (default)
    [X_opt_2,val_2] = opt_and_print(p,b,val,x0_new);
    x0_new = cell2struct(num2cell(X_opt_2(1:end-1)),b.var_names(1:end-1)',1);
    opt_and_print(p,b,val_2,x0_new);

    % force saturation off
    p_mod = p;
    p_mod.use_force_sat = false;
    opt_and_print(p_mod,b,val_2,x0_new)

    % power saturation off
    p_mod = p;
    p_mod.use_power_sat = false;
    opt_and_print(p_mod,b,val_2,x0_new)

    % force and power saturation off
    p_mod.use_force_sat = false;
    opt_and_print(p_mod,b,val_2,x0_new)

end

function val_mod = sim_and_print(X_opt_mod,p,b,val)
    [~,~,g,val_mod] = simulation(X_opt_mod,p); 
    constraint_violated = b.constraint_names(g<0)
    print_saturation_effect(val_mod,val)
end
function [X_opt_mod,val_mod] = opt_and_print(p_mod,b,val,x0_input)
    disp(['Force saturation ' num2str(p_mod.use_force_sat) ', power saturation ' num2str(p_mod.use_power_sat)])
    [X_opt_mod, ~, flag, ~, ~, ~, ~, val_mod] = gradient_optim(x0_input,p_mod,b,1);
    X_opt_mod
    flag
    print_saturation_effect(val_mod,val)
end


function print_saturation_effect(val_mod,val)

    power_diff = val_mod.power_avg - val.power_avg
    pct_power_diff = power_diff / val.power_avg * 100
    pct_PTO_cost = (val_mod.capex_PTO - val.capex_PTO) / val.capex_PTO * 100
    pct_struct_cost = (val_mod.capex_struct - val.capex_struct) / val.capex_struct * 100
    pct_LCOE = (val_mod.LCOE - val.LCOE) / val.LCOE * 100

    disp('=============================')
end