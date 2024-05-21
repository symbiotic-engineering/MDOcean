% todo:
% - get actual hydro coeffs (check)
% - loop to get coeffs as a function of frequency (check)
% - write symbolic to file so it runs the solve fast (check)
% - check that the same dimensions scaled by some constant give the same
% nondim hydro coeffs (checked, doesn't match, must resolve)
% - do convergence study on hydro coeff vs N (check)
% - compare hydro coeffs against Yeung (check) and WAMIT (unofficially
% checked, A is ~2x too big and B is reasonable)
% - get excitation force
% - do convergence study on matching quality vs N
% - integrate into mdocean (check)

% extensions for later:
% - time profiling and speedup (check)
% - h = 1 gives divide by zero error, investigate
% - update auto BCs to work with different N,M,K
% - find the difference between auto BCs and regular BCs
% - try the more convenient i2 eigenfunctions mentioned in the 2012 paper
% - condense the A-matrix to be 2M x 2M from p7 of notebook
% - get surge coeffs
% - do low-frequency approximations for cummins equation
% - take the limit in deep water
% - get coupling hydro coeffs for when moving at diff velocities
% - extend to many regions

clear all
close all

%% settings 
auto_BCs = false;

num_harmonics = 3;
N_num = num_harmonics;
M_num = num_harmonics;
K_num = num_harmonics;

heaving_OC = true;
heaving_IC = true;

plot_phi = false;
show_A = false;

%% set numerical values
a1_num = .5;
a2_num = 1;
d1_num = .5;
d2_num = .25;
h_num = 1.001;
m0_nums = linspace(0.1,5,100);
spatial_res = 30;

%% run
tic
[mu_nondim, lambda_nondim] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, N_num, M_num, K_num, ...
                       a1_num, a2_num, d1_num, d2_num, h_num, m0_nums, spatial_res, show_A, plot_phi);

toc

%% plot hydro coeffs vs wavenumber
figure
plot(m0_nums,mu_nondim, m0_nums,lambda_nondim)
xlabel('Wavenumber m_0')
ylabel('Nondimensional Hydro Coeff')
legend('Added Mass','Damping')
grid on
