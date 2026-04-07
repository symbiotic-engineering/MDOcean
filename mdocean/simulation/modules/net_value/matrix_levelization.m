function levelized_X_matrix = matrix_levelization(capex, opex, FCR, force_matrix)

levelized_x = FCR * capex + opex; %$/year
levelized_x_matrix_unscaled = FCR * capex * force_matrix + opex; % $/year
matrix_scale_factor = levelized_x / sum(levelized_x_matrix_unscaled,'all','omitnan');
levelized_X_matrix = levelized_x_matrix_unscaled * matrix_scale_factor;