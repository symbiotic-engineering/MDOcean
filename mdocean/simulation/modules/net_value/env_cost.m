function [LCOE_matrix,LCOE] = env_cost(force_matrix,P_elec,eco_cost_steel,...
    steel_mass,eco_cost_fiberglass,fiberglass_area,eco_cost_distance,...
    distance_from_shore,JPD,N_WEC,FCR,efficiency)

% calculate capex and opex
capex_per_wec = eco_cost_steel * steel_mass + eco_cost_fiberglass * fiberglass_area;
opex_per_wec = eco_cost_distance * distance_from_shore;
capex = capex_per_wec * N_WEC;
opex = opex_per_wec * N_WEC;

%opex should be per year
maintenances_per_year = 2; %from RM3 report page 170
opex = opex * maintenances_per_year;

%calculate the levelized cost matrix
levelized_cost_matrix = matrix_levelization(capex, opex, FCR, force_matrix);
[LCOE_matrix,LCOE] = LXOE_func(levelized_cost_matrix,P_elec,efficiency,JPD,N_WEC);