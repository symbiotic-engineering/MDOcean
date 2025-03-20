function [A_f_over_rho, A_s_over_rho, A_c_over_rho, ...
        B_f_over_rho_w, B_s_over_rho_w, B_c_over_rho_w, ...
        gamma_f_over_rho_g, gamma_s_over_rho_g, ...
        gamma_phase_f, gamma_phase_s] = get_hydro_coeffs_MEEM(a2, m0, d2, a1, d1, a3, h, g, w, harmonics, spar_excitation_coeffs)

heaving_IC = false;
heaving_OC = true;
auto_BCs = false;
spatial_res = 0;
show_A = false;
plot_phi = false;

N_num = harmonics;
M_num = harmonics;
K_num = harmonics;

% run MEEM only unique freqs
m0_meem = unique(m0(~isnan(m0)));
if ~isrow(m0_meem)
    m0_meem = m0_meem.'; % must be row vector
end
[mu_nondim, lambda_nondim, gamma_phase_f] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, ...
                                               N_num, M_num, K_num, ...
                                               a1/h, a2/h, d1/h, ...
                                               d2/h, 1, m0_meem*h, ...
                                               spatial_res, show_A, plot_phi);

% expand matrix for each Hs - this assumes that m0 is identical (except nans) in each column
num_Hs = size(m0,1);
mu_nondim     = repmat(mu_nondim,     [num_Hs 1]);
lambda_nondim = repmat(lambda_nondim, [num_Hs 1]);
gamma_phase_f = repmat(gamma_phase_f, [num_Hs 1]);

normalize = pi * a2^3;
A_f_over_rho   = mu_nondim     * normalize;
B_f_over_rho_w = lambda_nondim * normalize;

% Excitation
[~, mult] = group_velocity(w, m0, g, h); % finite depth multiplier to group velocity
gamma_f_over_rho_g = sqrt(2 * mult .* B_f_over_rho_w ./ m0); % Haskind relationship, using the finite depth group velocity (see notebook p129 2/26/25)

% since MEEM currently doesn't take the damping plate into account,
% get the spar values with approximations that include damping plate
[A_s_over_rho,gamma_s_over_rho_g,B_s_over_rho_w] = spar_dynamics(a1/a3, 2*a3, d1, ...
                                                    spar_excitation_coeffs, m0, mult);

% set the remaining coeffs to zero since not sure how to approximate them
A_c_over_rho = 0;
B_c_over_rho_w = 0;
gamma_phase_s = 0;

end