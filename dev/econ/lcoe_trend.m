% Script used to determine which cost categories scale with N_WEC in econ module
% Shows line of best fit, actual cost, and linear approximation used in sim

N_WEC = [1 10 50 100];
capex_per = [280 98 68 65]; % from figure 5-34 on p176 of RM3 report
opex_per = [165 47 19 13];  % from figure 5-35 on p177 of RM3 report
capex = capex_per .* N_WEC;
opex = opex_per .* N_WEC;

% cost = cost_fixed + cost_mult * N_WEC;

% capex: least squares fit
close all
figure
ft_capex = fit(N_WEC',capex','poly1');
plot(ft_capex,N_WEC,capex)

% capex: MDOcean simulation
cost_fixed = 4.55+0.99+0.53+0.36+5.91+1.59;
cost_mult = 0.53+2.94+0.62;
capex_sim = (cost_fixed + N_WEC*cost_mult) * capex(1)/(cost_fixed+cost_mult);
hold on
plot(N_WEC,capex_sim)
title('Capex')
xlabel('N WEC')

% opex: least squares fit
figure
ft_opex = fit(N_WEC',opex','poly1');
plot(ft_opex,N_WEC,opex)

% opex: MDOcean simulation
cost_fixed = 710+142+227;
cost_mult = 27+54+8;
opex_sim = (cost_fixed + N_WEC*cost_mult) * opex(1)/(cost_fixed+cost_mult);
hold on
plot(N_WEC,opex_sim)
title('Opex')
xlabel('N WEC')