function levelized_X_matrix = matrix_levelization(capex, opex, FCR, force_matrix)

levelized_x_matrix_unscaled = FCR * capex * force_matrix + opex; % $/year
matrix_scale_factor = levelized_cost / mysum(levelized_cost_matrix_unscaled);
levelized_X_matrix = levelized_x_matrix_unscaled * matrix_scale_factor;