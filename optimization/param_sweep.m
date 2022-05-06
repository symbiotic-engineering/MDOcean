% Brute force parameter sensitivity sweep (reoptimize for each param value)
%% Setup
clear;clc;close all

var_names = {'Hs','Hs_{struct}','T','T_{struct}','\sigma_y','\rho_m','E',...
            'cost_m','t_{sft}','t_{sf}','t_{sfb}','t_r', 't_{vc}','B_{min}',...
            'FOS_{min}','FCR','N_{WEC}'};   % list of parameters to sweep
vars = regexprep(var_names,'[{}\\]','');    % remove the curly braces and slashes

ratios = .8 : .1 : 1.2;
p = parameters();
b = var_bounds(p);

% use the optimal x as x0 to speed up the sweeps
x0 = struct('D_sft',b.D_sft_nom,'D_i_ratio',b.D_i_ratio_nom,'D_or_ratio',...
        b.D_or_ratio_nom,'F_max',b.F_max_nom,'D_int',b.D_int_nom,'w_n',b.w_n_nom);
x0_vec = gradient_optim(x0,p,b);
x0 = struct('D_sft',x0_vec(1),'D_i_ratio',x0_vec(2),'D_or_ratio',x0_vec(3),...
    'M',x0_vec(4),'F_max',x0_vec(5),'D_int',x0_vec(6),'w_n',x0_vec(7));
   
LCOE  = zeros(length(vars),length(ratios));
P_var = zeros(length(vars),length(ratios));

%% Run optimization
for i=1:length(vars)
    p = parameters();
    var_nom = p.(vars{i});
    for j=1:length(ratios)
        p.(vars{i}) = ratios(j) * var_nom;
        [~, obj_opt, flag] = gradient_optim(x0,p,b);
        if flag >= 1
            LCOE(i,j) = obj_opt(1);
            P_var(i,j) = obj_opt(2);
        else
            [LCOE(i,j), P_var(i,j)] = deal(NaN);
        end
    end   
end
%% Plot each sensitivity
col_nom = find(ratios==1);
LCOE_nom = LCOE(1,col_nom);
Pvar_nom = P_var(1,col_nom);

figure
plot(ratios,LCOE/LCOE_nom)
xlabel('Parameter Ratio from nominal')
ylabel('LCOE ratio from nominal')
legend(var_names)
improvePlot
grid on

figure
plot(ratios,P_var/Pvar_nom)
xlabel('Parameter ratio from nominal')
ylabel('Power Variation ratio from nominal')
legend(var_names)
improvePlot
grid on
%% Plot overall slope as tornado chart

slope_LCOE = (LCOE(:,end) - LCOE(:,1))/LCOE_nom;
slope_Pvar = (P_var(:,end) - P_var(:,1))/Pvar_nom;

% separate charts for each objective
figure
subplot 121
barh(categorical(var_names),slope_LCOE)
title('LCOE')
subplot 122
barh(categorical(var_names),slope_Pvar)
title('c_v')
sgtitle('Normalized Sensitivities')
improvePlot

% both objectives on the same chart
figure
barh(categorical(var_names),[slope_LCOE slope_Pvar])
legend('LCOE','P_{var}')
title('Sensitivities')
