
clear all
close all

%% settings 
auto_BCs = false;

heaving_OC = true;
heaving_IC = false;

plot_phi = false;
show_A = false;

geometry = 'RM3'; % RM3 or default

%% set numerical values
if strcmp(geometry,'default')
    a1_num = .5;
    a2_num = 1;
    d1_num = .5;
    d2_num = .25;
    h_num = 1.001;
elseif strcmp(geometry,'RM3')
    a1_num = 3/100;
    a2_num = 10/100;
    d1_num = 35/100;
    d2_num = 2/100;
    h_num = 1;
end
m0_nums = linspace(0.1,5,100);
spatial_res = 30;

num_harmonics = [3 5 10 20 30];
run_multiple_harmonics(num_harmonics, heaving_IC, heaving_OC, auto_BCs, ...
                       a1_num, a2_num, d1_num, d2_num, h_num, m0_nums, spatial_res, show_A, plot_phi);

function run_multiple_harmonics(num_harmonics, heaving_IC, heaving_OC, auto_BCs, ...
                       a1_num, a2_num, d1_num, d2_num, h_num, m0_num, spatial_res, show_A, plot_phi)

    numels = [numel(a1_num), numel(a2_num), numel(d1_num), numel(d2_num), ...
                         numel(h_num), numel(m0_num)];
    mu_nondim = zeros(length(num_harmonics), max(numels));
    lambda_nondim = zeros(length(num_harmonics), max(numels));
    
    for i=1:length(num_harmonics)
        N_num = num_harmonics(i);
        M_num = num_harmonics(i);
        K_num = num_harmonics(i);
        [mu_nondim(i,:), lambda_nondim(i,:)] = run_MEEM(heaving_IC, heaving_OC, auto_BCs, N_num, M_num, K_num, ...
                           a1_num, a2_num, d1_num, d2_num, h_num, m0_num, spatial_res, show_A, plot_phi);
    end

    %% plot hydro coeffs vs wavenumber
    if i==1
        figure
    end
    plot(m0_num,mu_nondim, m0_num,lambda_nondim,'--')
    xlabel('Wavenumber m_0')
    ylabel('Nondimensional Hydro Coeff')
    legend('Added Mass','Damping')
    grid on
    hold on
    legend(cellstr(num2str([num_harmonics(:); num_harmonics(:)])))
    improvePlot

    figure
    plot(num_harmonics, mu_nondim(:,1), '*-', num_harmonics, lambda_nondim(:,1),'*-')
    xlabel('Number of Harmonics')
    ylabel('Nondimensional Hydro Coefficient')
    legend('Added Mass','Damping')
    improvePlot
end
