function [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs(r, k, draft)
    
    % Froude Krylov force coefficient (diffraction is neglected)

    tuning_factor =  4.5; % tune to more closely match WAMIT results which include diffraction

    r_k_term = r^2 - (k.^2 * r^4)/8 + (k.^4 * r^6)/192 - (k.^6 * r^8)/9216 ...
                + (k.^8 * r^10)/737280 - (k.^10 * r^12)/88473600;

    r_k_term = abs(r_k_term); % get rid of any negatives that result at high frequencies

    gamma_over_rho_g = pi * exp(-k * draft * tuning_factor) .* r_k_term;
    
    % Added mass
    A_over_rho = 1/2 * 4/3 * pi * r^3 * 0.63;

    % Radiation damping
    B_over_rho_w = k/2 .* gamma_over_rho_g.^2; % Haskind relationship, using the deep water group velocity

end

