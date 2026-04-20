function [LVOE_matrix, LVOE_env] = env_value(force_matrix, P_elec, FCR, efficiency, years, avoided_co2_tons, scc_case, JPD, N_WEC)

if nargin < 7 || isempty(years) || isempty(avoided_co2_tons) || isempty(scc_case)
    error("env_value(force_matrix,P_elec,FCR,efficiency,years,avoided_co2_tons,scc_case) required.");
end

[~, ~, avoided_dollars_per_year] = carbon_avoided_value(years, avoided_co2_tons, scc_case);

capex = 0;
opex  = sum(avoided_dollars_per_year);

levelized_value_matrix = matrix_levelization(capex, opex, FCR, force_matrix);
[LVOE_matrix, LVOE_env] = LXOE_func(levelized_value_matrix, P_elec, efficiency,JPD,N_WEC);

end