function sigma = get_plate_stress(moment_per_length, y_max, h_eq)
    M = moment_per_length;

    % see notebook p21 12/15/24
    sigma = 12 * M .* y_max ./ h_eq.^3;

end
