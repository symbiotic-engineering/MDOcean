function [A_f_over_rho, A_s_over_rho, A_c_over_rho, ...
    B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
    gamma_f_over_rho_g, gamma_s_over_rho_g, ...
    gamma_phase_f, gamma_phase_s] = get_hydro_coeffs(r, k, draft, hydro)

    % Froude Krylov force coefficient (diffraction is neglected)

    % tuning_factor =  4.5; % tune to more closely match WAMIT results which include diffraction
    %
    % r_k_term = r^2 - (k.^2 * r^4)/8 + (k.^4 * r^6)/192 - (k.^6 * r^8)/9216 ...
    %             + (k.^8 * r^10)/737280 - (k.^10 * r^12)/88473600;
    %
    % r_k_term = abs(r_k_term); % get rid of any negatives that result at high frequencies
    %
    % gamma_over_rho_g = pi * exp(-k * draft * tuning_factor) .* r_k_term;
    %
    % % Added mass
    % A_over_rho = 1/2 * 4/3 * pi * r^3 * 0.63;

    k_wamit = hydro.w.^2 / 9.8;

    % added mass
    A_f_wamit = hydro.A(3,3,:);
    A_s_wamit = hydro.A(9,9,:);
    A_c_wamit = hydro.A(3,9,:);

    % excitation
    gamma_f_wamit = hydro.ex_ma(3,1,:);
    gamma_phase_f_wamit = -hydro.ex_ph(3,1,:);
    gamma_s_wamit = hydro.ex_ma(9,1,:);
    gamma_phase_s_wamit = -hydro.ex_ph(9,1,:);

    % radiation damping
    B_f_wamit = hydro.B(3,3,:);
    B_s_wamit = hydro.B(9,9,:);
    B_c_wamit = hydro.B(3,9,:);

    % interpolation
    A_f_over_rho       = interp1(k_wamit,A_f_wamit(:),k);
    A_s_over_rho       = interp1(k_wamit,A_s_wamit(:),k);
    A_c_over_rho       = interp1(k_wamit,A_c_wamit(:),k);
    B_f_over_rho_w     = interp1(k_wamit,B_f_wamit(:),k);
    B_s_over_rho_w     = interp1(k_wamit,B_s_wamit(:),k);
    B_c_over_rho_w     = interp1(k_wamit,B_c_wamit(:),k);
    gamma_f_over_rho_g = interp1(k_wamit,gamma_f_wamit(:),k);
    gamma_s_over_rho_g = interp1(k_wamit,gamma_s_wamit(:),k);
    gamma_phase_f      = interp1(k_wamit,gamma_phase_f_wamit(:),k);
    gamma_phase_s      = interp1(k_wamit,gamma_phase_s_wamit(:),k);

    %B_over_rho_w = k/2 .* gamma_over_rho_g.^2; % Haskind relationship, using the deep water group velocity

end
