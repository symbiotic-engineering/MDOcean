function [LCOE_matrix,LCOE] = econ_cost(force_matrix,P_elec,capex_design_dep,N_WEC,FCR,efficiency,JPD)

% total capex
alpha_non_design = 0.741;
capex_non_design_dep = 12.68e6 * N_WEC^(-alpha_non_design) + 1.24e6;
capex_per_wec = capex_design_dep + capex_non_design_dep;
capex = capex_per_wec * N_WEC;

% opex = operation, postinstall, replacement, consumables, and insurance
alpha_opex = 0.5567;
opex_per_wec = 1.193e6 * N_WEC^-alpha_opex;
opex = opex_per_wec * N_WEC;

%calculate the levelized cost matrix
levelized_cost_matrix = matrix_levelization(capex, opex, FCR, force_matrix);
[LCOE_matrix,LCOE] = LXOE_func(levelized_cost_matrix,P_elec,efficiency,JPD,N_WEC);