function [net_econ_value,net_eco_value,eco_cost_per_wec,LCOE_econ,capex_design] = net_value(F_heave_mat,P_matrix_elec,...
    N_WEC,FCR,eff_array,JPD,marginal_price,eco_cost_steel,m_m, ...
    eco_cost_fiberglass,A_hull,eco_cost_distance,distance_from_shore, ...
    marginal_carbon,M,cost_perkg_mult,cost_perN_mult,cost_perW_mult,F_max,P_max)

%scale force
factor = 1./sum(F_heave_mat,'all','omitnan');
force_matrix = F_heave_mat.*factor;

%assemble the inputs for environmental value
years = 2010 + (0:20);
P_weighted = P_matrix_elec .* JPD / 100 * eff_array; % taken from dynamics.m and LCOE_from_capex_design_power.m
hr_per_yr = 8766;
AEP_matrix = P_weighted * N_WEC * hr_per_yr / 1000; % W to kWh per year, all wecs
avoided_co2_tons = marginal_carbon .* AEP_matrix / 1000; %kgCO2 to tonsCO2
avoided_co2_tons = sum(avoided_co2_tons,'all','omitnan');
scc_case = '3pct';

%calculate capex_design_dep for econ cost
[capex_design,~, ~] = design_cost_model(m_m, M, cost_perkg_mult, N_WEC, ...
                                      cost_perN_mult, cost_perW_mult, F_max, P_max);

%call the cost and value wrappers
[~,LCOE_econ] = econ_cost(force_matrix,P_matrix_elec,capex_design,N_WEC,...
    FCR,eff_array,JPD);
[~, LVOE_econ] = econ_value(force_matrix, P_matrix_elec, FCR, eff_array,...
    JPD, N_WEC, marginal_price);
[~,LCOE_env,eco_cost_per_wec] = env_cost(force_matrix,P_matrix_elec,eco_cost_steel,...
    m_m,eco_cost_fiberglass,A_hull,eco_cost_distance,...
    distance_from_shore,JPD,N_WEC,FCR,eff_array);
[~, LVOE_env] = env_value(force_matrix, P_matrix_elec, FCR, eff_array, years,...
    avoided_co2_tons, scc_case, JPD, N_WEC);
   
net_econ_value = (LVOE_econ - LCOE_econ) * sum(AEP_matrix,'all','omitnan');
net_eco_value = (LVOE_env - LCOE_env) * sum(AEP_matrix,'all','omitnan');
end


function [capex_design_dep, ...
          pto, devicestructure] = design_cost_model(mass_material, M, cost_perkg_mult, N_WEC, ...
                                  cost_perN_mult, cost_perW_mult, F_max, P_max)

    % from RM3 CBS 1.4 and 1.5, with curve fits done in dev/design_cost_scaling.m
    
    % structural cost per wec
    alpha_struct = 0.481;
    cost_per_kg = ( 1.64e6 + 1.31e6 * N_WEC^(-alpha_struct) ) / 687000 * cost_perkg_mult(M);
    devicestructure = cost_per_kg * mass_material;
    
    % PTO cost per wec
    alpha_pto = 0.206;
    pto_const =    92593 + 1051  *N_WEC^(-alpha_pto);
    pto_power = ( 0.4454 + 0.9099*N_WEC^(-alpha_pto) ) * P_max * cost_perW_mult;
    pto_force = ( 0.0086 + 0.0118*N_WEC^(-alpha_pto) ) * F_max * cost_perN_mult;
    pto = pto_const + pto_power + pto_force;
    
    % sum design-dependent cost per wec
    capex_design_dep = devicestructure + pto; 

end