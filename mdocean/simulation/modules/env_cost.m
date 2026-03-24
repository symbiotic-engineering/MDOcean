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