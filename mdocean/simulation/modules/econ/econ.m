
function [LCOE, capex_design_dep, ...
            capex, opex, pto, devicestructure] = econ(mass_material, M, cost_perkg_mult, N_WEC, P_elec, FCR, ...
                                cost_perN_mult, cost_perW_mult, F_max, P_max, efficiency)
    
    [capex_design_dep, ...
     pto, devicestructure] = design_cost_model(mass_material, M, cost_perkg_mult, N_WEC, ...
                                      cost_perN_mult, cost_perW_mult, F_max, P_max);
    
    [LCOE,capex,opex] = LCOE_from_capex_design_power(capex_design_dep, N_WEC, P_elec, FCR, efficiency);

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
    pto_power = ( 0.4454 + 0.9099*N_WEC^(-alpha_pto) ) * P_max * cost_perN_mult;
    pto_force = ( 0.0648 + 0.0893*N_WEC^(-alpha_pto) ) * F_max * cost_perW_mult;
    pto = pto_const + pto_power + pto_force;
    
    % sum design-dependent cost per wec
    capex_design_dep = devicestructure + pto; 

end