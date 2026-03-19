
function [fig_vs_omega,fig_vs_NMK] = convergence_study()

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

num_harmonics = [5 10 20 30];
[fig_vs_omega,fig_vs_NMK] = run_multiple_harmonics(num_harmonics, heaving_IC, heaving_OC, auto_BCs, ...
                       a1_num, a2_num, d1_num, d2_num, h_num, m0_nums, spatial_res, show_A, plot_phi);

function [fig1,fig2] = run_multiple_harmonics(num_harmonics, heaving_IC, heaving_OC, auto_BCs, ...
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
    fig1 = figure;
    rgb = get(groot,"FactoryAxesColorOrder");
    markers = {'-','--',':','.-'};
    plot(NaN,NaN,'Color',rgb(1,:),'DisplayName','Added Mass')
    hold on
    plot(NaN,NaN,'Color',rgb(2,:),'DisplayName','Damping')
    for i=1:length(num_harmonics)
        h1 = plot(m0_num, mu_nondim(i,:),   markers{i},'Color',rgb(1,:));
        h2 = plot(m0_num,lambda_nondim(i,:),markers{i},'Color',rgb(2,:));
        h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        plot(NaN,NaN,[markers{i} 'k'],'DisplayName',num2str(num_harmonics(i))) % dummy for legend
    end
    xlabel('Nondimensional Wavenumber $m_0h$','Interpreter','latex')
    ylabel('Nondimensional Radiation Coefficient','Interpreter','latex')
    legend('Interpreter','latex')
    grid on
    improvePlot

    fig2 = figure;
    plot(num_harmonics, mu_nondim(:,1), '*-', num_harmonics, lambda_nondim(:,1),'*-')
    xlabel('Number of Harmonics')
    ylabel('Nondimensional Hydro Coefficient')
    legend('Added Mass','Damping')
    improvePlot
end

end
