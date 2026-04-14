function [LVOE_matrix, LVOE] = econ_value(force_matrix, P_elec, FCR, efficiency, JPD, N_WEC, price_contour)

capex = 0;
hr_per_yr = 8766;

P_weighted = P_elec .* (JPD/100) .* efficiency;
AEP_matrix = P_weighted * N_WEC * hr_per_yr / 1000;  

revenue_matrix = price_contour .* AEP_matrix;       
opex = sum(revenue_matrix, 'all', 'omitnan');        

levelized_value_matrix = matrix_levelization(capex, opex, FCR, force_matrix);
[LVOE_matrix, LVOE] = LXOE_func(levelized_value_matrix, P_elec, efficiency, JPD, N_WEC);

end