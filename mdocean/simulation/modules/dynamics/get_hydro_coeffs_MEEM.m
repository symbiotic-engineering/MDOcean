function [A_over_rho, B_over_rho_w, gamma_over_rho_g] = get_hydro_coeffs_MEEM(a2, m0, d2, a1, d1, h, harmonics)

heaving_IC = false;
heaving_OC = true;
auto_BCs = false;
spatial_res = 0;
show_A = false;
plot_phi = false;

N_num = harmonics;
M_num = harmonics;
K_num = harmonics;

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
A_over_rho   = mu_nondim     * normalize;
B_over_rho_w = lambda_nondim * normalize;

% Radiation damping
gamma_over_rho_g = sqrt(2 * B_over_rho_w ./ m0); % Haskind relationship, using the deep water group velocity

end