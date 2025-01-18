% clear all
% close all

%% settings 
auto_BCs = false;

num_harmonics = 10;
N_num = num_harmonics;
M_num = num_harmonics;
K_num = num_harmonics;

heaving_OC = true;
heaving_IC = false;

plot_phi = false;
show_A = false;

%% set numerical values
a1_num = .5;
a2_num = 1;
d1_num = .5;
d2_num = .25;
h_num = 1.05;
m0_num = 1;
spatial_res = 30;

%% run validation of potential for single frequency
tic
[mu_nondim, lambda_nondim]= run_MEEM(heaving_IC, heaving_OC, auto_BCs, N_num, M_num, K_num, ...
                       a1_num, a2_num, d1_num, d2_num, h_num, m0_num, spatial_res, show_A, plot_phi);
sim_time = toc

mu_nondim*1025*h_num^3

% figure(7)
% plot_potential_validation()
% 
% %% run validation of hydro coeffs for various frequencies
% m0_nums = [linspace(0,0.1,10), linspace(0.1,6,100)];
% plot_phi = false;
% [mu_nondim, lambda_nondim] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, N_num, M_num, K_num, ...
%                        a1_num, a2_num, d1_num, d2_num, h_num, m0_nums, spatial_res, show_A, plot_phi);
% 
% figure
% plot(m0_nums, mu_nondim, m0_nums, lambda_nondim)
% ylabel('Nondimensional Hydro Coeff')
% xlabel('Wavenumber m0')
% legend('Added Mass','Damping')
% grid on
% hold on
% plot_hydro_coeff_validation()
% 
% %% validation functions
% function plot_potential_validation()
%     a1_potential = readmatrix("inputs/validation/MEEM_validation/potential_a1.csv");
%     a2_potential = readmatrix("inputs/validation/MEEM_validation/potential_a2.csv");
%     plot(a1_potential(:,1),a1_potential(:,2),'m-*','DisplayName','Yeung 2012 at a_1')
%     plot(a2_potential(:,1),a2_potential(:,2),'b-*','DisplayName','Yeung 2012 at a_2')
% end
% 
% function plot_hydro_coeff_validation()
%     mu_nondim = readmatrix("inputs/validation/MEEM_validation/added_mass.csv");
%     lambda_nondim = readmatrix("inputs/validation/MEEM_validation/damping.csv");
%     excitation_nondim = readmatrix("inputs/validation/MEEM_validation/excitation.csv");
%     excitation_phase_nondim = readmatrix("inputs/validation/MEEM_validation/excitation_phase.csv");
% 
%     plot(mu_nondim(:,1),     mu_nondim(:,2),    'c--','DisplayName','Added Mass Yeung 2012')
%     plot(lambda_nondim(:,1), lambda_nondim(:,2),'m--','DisplayName','Damping Yeung 2012')
%     improvePlot
% end