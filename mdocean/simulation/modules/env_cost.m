function [LCOE_matrix,LCOE] = env_cost(force_matrix,P_elec)

% calculate capex and opex
capex = 0;

%fix opex: $/yr of energy revenue on the grid
alpha_opex = 0.5567;
opex_per_wec = 1.193e6 * N_WEC^-alpha_opex;
opex = opex_per_wec * N_WEC;

%calculate the levelized cost matrix
levelized_cost_matrix = matrix_levelization(capex, opex, FCR, force_matrix);
[LCOE_matrix,LCOE] = LXOE_func(levelized_cost_matrix,P_elec,efficiency);

function [eco_cost_per_wec, eco_cost_total] = env_cost(steel_mass, distance_from_shore, fiberglass_area, ...
                                                      eco_cost_steel, eco_cost_fiberglass, eco_cost_distance, ...
                                                      num_wecs)

eco_cost_per_wec = eco_cost_steel * steel_mass + eco_cost_fiberglass * fiberglass_area + eco_cost_distance * distance_from_shore;
eco_cost_total   = eco_cost_per_wec * num_wecs;

end