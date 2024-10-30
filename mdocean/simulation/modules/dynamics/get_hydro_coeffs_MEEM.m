function [A_f_over_rho, A_s_over_rho, A_c_over_rho, ...
        B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
        gamma_f_over_rho_g, gamma_s_over_rho_g, ...
        gamma_phase_f, gamma_phase_s] = get_hydro_coeffs_MEEM(a2, m0, d2, a1, d1, a3, h, harmonics, spar_excitation_coeffs)

heaving_IC = false;
heaving_OC = true;
auto_BCs = false;
spatial_res = 0;
show_A = false;
plot_phi = false;

N_num = harmonics;
M_num = harmonics;
K_num = harmonics;

% run MEEM only for the first row of m0, which should capture all unique freqs
[mu_nondim, lambda_nondim] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, ...
                                               N_num, M_num, K_num, ...
                                               a1/h, a2/h, d1/h, ...
                                               d2/h, 1, m0(1,:)*h, ...
                                               spatial_res, show_A, plot_phi);

% expand matrix for each Hs
num_Hs = size(m0,1);
mu_nondim = repmat(mu_nondim,[num_Hs 1]);
lambda_nondim = repmat(lambda_nondim, [num_Hs 1]);

normalize = pi * a2^3;
A_f_over_rho   = mu_nondim     * normalize;
B_f_over_rho_w = lambda_nondim * normalize;

% Radiation damping
gamma_f_over_rho_g = sqrt(2 * B_f_over_rho_w ./ m0); % Haskind relationship, using the deep water group velocity

% since MEEM currently doesn't take the damping plate into account,
% get the spar values with approximations that include damping plate
[A_s_over_rho,gamma_s_over_rho_g,B_s_over_rho_w] = spar_dynamics(a1/a3, 2*a3, d1, ...
                                                    spar_excitation_coeffs, m0);

% set the remaining coeffs to zero since not sure how to approximate them
A_c_over_rho = 0;
B_c_over_rho_w = 0;
gamma_phase_f = zeros(size(mu_nondim)); % I would maybe be able to get this from Haskind of MEEM results?
gamma_phase_s = 0;

end