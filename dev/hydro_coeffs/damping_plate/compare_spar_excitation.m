% compare getting radiation damping and excitation force for the spar
% using froude krylov and using MEEM and using WAMIT

w_max = 1.5;
w = linspace(0.01, w_max, 100);
a1 = 3;
a2 = 10;
d2 = 2;
d1 = 35;
g = 9.8;
h = 3000;

k = w.^2 / g;

% Froude Krylov
gamma_over_rho_g = pi * exp(-k * d1) * (a1 / 2)^2;
B_over_rho_w = k / 2 .* gamma_over_rho_g.^2;

% Froude Krylov take 2
gamma_over_rho_g_2 = 2 * pi * (a1 / 2) * exp(-k * d1) .* besselj(1, k * a1 / 2) ./ k;
B_over_rho_w_2 = k / 2 .* gamma_over_rho_g_2.^2;

% MEEM
heaving_IC = true;
heaving_OC = false;
auto_BCs = false;
harmonics = 10;
N_num = harmonics;
M_num = harmonics;
K_num = harmonics;
spatial_res = 0;
show_A = false;
plot_phi = false;
[~, lambda_nondim] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, ...
                              N_num, M_num, K_num, ...
                              a1 / h, a2 / h, d1 / h, ...
                              d2 / h, h / h, k * h, ...
                              spatial_res, show_A, plot_phi);

% [~, lambda_nondim] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, ...
%                                                N_num, M_num, K_num, ...
%                                                a1, a2, d1, ...
%                                                d2, h, k, ...
%                                                spatial_res, show_A, plot_phi);

normalize = 1; % pi * a1^3;
B_over_rho_w_MEEM = lambda_nondim * normalize;
gamma_over_rho_g_MEEM = sqrt(2 * B_over_rho_w_MEEM ./ k);

% WAMIT
% borrowing code from hydro_coeff_err.m
hydro = struct();
hydro = readWAMIT(hydro, 'rm3.out', []); % function from WECSim
w_WAMIT = hydro.w;
w_WAMIT = w_WAMIT(hydro.w < w_max);
B = hydro.B(9, 9, hydro.w < w_max);
B = B(:);
gamma = hydro.ex_ma(9, 1, hydro.w < w_max);
gamma = gamma(:);

% plot
figure;
plot(w, gamma_over_rho_g, w, gamma_over_rho_g_2, '--', w, gamma_over_rho_g_MEEM, w_WAMIT, gamma);
legend('Froude Krylov', 'Froude Krlyov Take 2', 'MEEM', 'WAMIT');
% title('Spar Radiation Damping B_{33}/\rho\omega')
title('Spar Excitation \gamma_{3}/\rho g');
xlabel('\omega (rad/s)');
improvePlot;

figure;
plot(w_WAMIT, B ./ gamma.^2 * 10000, w_WAMIT, w_WAMIT.^2 / (2 * 9.8) * 10000, '--');
legend('B/\gamma^2*\rho g^2/\omega', '\omega^2/2g');

%% WAMIT scaling to estimate other geometries
D_d = 30;
d1_sweep = [15 25 35 45];
D_d_sweep = [10 20 30 40];
D_s_sweep = [3 6 9];

k_WAMIT = w_WAMIT'.^2 / g;
k_D_d_WAMIT = k_WAMIT * D_d;
gamma_over_ekz = gamma ./ exp(-d1 * k_WAMIT);
figure;
plot(k_D_d_WAMIT, gamma);
xlabel('k D_d');
ylabel('\gamma');
title('Nondimensionalized WAMIT Excitation');
